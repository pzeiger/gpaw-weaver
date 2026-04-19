import inspect
import json
from pathlib import Path

from gpaw.calculator import GPAW
from gpaw.new.ase_interface import GPAW as _NewGPAW

from .log import extract_scf_convergence

_NEW_GPAW_PARAMS = set(inspect.signature(_NewGPAW).parameters)

_DEFAULT_GPW_DIR = Path('gpw_files')
_DEFAULT_GPW_LOGS = Path('gpw_logs')
_DEFAULT_DB = Path('calculations.db')

try:
    from ase.db.core import reserved_keys as _ase_reserved_keys
    _ASE_RESERVED = set(_ase_reserved_keys)
except ImportError:
    _ASE_RESERVED = {
        'magmoms', 'magmom', 'charges', 'energy', 'free_energy', 'forces',
        'stress', 'dipole', 'numbers', 'positions', 'cell', 'pbc',
        'masses', 'tags', 'momenta', 'constraints', 'calculator',
        'calculator_parameters', 'initial_magmoms', 'initial_charges',
    }


def _resolve_db(db):
    import ase.db
    if db is None:
        return ase.db.connect(str(_DEFAULT_DB))
    if isinstance(db, (str, Path)):
        p = Path(db)
        if not p.suffix:
            p = p.with_suffix('.db')
        return ase.db.connect(str(p))
    return db


def _safe_db_key(key):
    return f'gpaw_{key}' if key in _ASE_RESERVED else key


def _serialize_calc_params(calc_params):
    out = {}
    for key, val in calc_params.items():
        safe_key = _safe_db_key(key)
        if type(val) in (bool, int, float, str):
            out[safe_key] = val
        else:
            try:
                out[safe_key] = json.dumps(val, sort_keys=True)
            except Exception:
                out[safe_key] = json.dumps(val.todict(), sort_keys=True)
    return out


def run_and_store_gpaw_calculation(atoms_initial, calc_params, system,
                                   db=None,
                                   save_gpw=False,
                                   save_gpw_mode='calculation',
                                   legacy_gpaw=True,
                                   gpw_dir=_DEFAULT_GPW_DIR,
                                   gpw_logs=_DEFAULT_GPW_LOGS):
    """Run a GPAW calculation and store results in the ASE database.

    Log and GPW files are named after the database row IDs of the stored
    entries so they are always uniquely identifiable:

    * Log  → ``<gpw_logs>/<converged_id>.txt``
    * GPW  → ``<gpw_dir>/<converged_id>.gpw``

    The initial DB entry is written *before* the calculation so its ID can
    be used for the log file during the run; the log is then renamed to the
    converged ID once that entry is written.

    Parameters
    ----------
    atoms_initial : ase.Atoms
        Starting geometry.
    calc_params : dict
        All GPAW parameters for this calculation. Use ``make_pw_params`` or
        ``make_fd_params`` to build this, then add system-specific settings
        (``nbands``, ``charge``, ``setups``, ``poissonsolver``, …) directly.
        Every key is serialised and stored in the database.
    system : str
        Human-readable label stored alongside the calculation.
    db : ase.db.core.Database or str or Path or None
        Database to write results into.  Accepts an already-connected ASE
        database object, a file path (str or Path) to connect to, or
        ``None`` to use the default ``calculations.db`` in the working
        directory.
    save_gpw : bool
        Whether to write a ``.gpw`` restart file (default False).
    save_gpw_mode : str
        Passed as the ``mode`` argument to ``calc.write()``.
        ``'calculation'`` (default) omits wavefunctions; ``'all'`` saves them.
    legacy_gpaw : bool
        Use the old GPAW implementation (``True``, default) or the new
        refactored one (``False``). Stored in the DB for use by
        ``load_gpaw_calculation``.
    gpw_dir : Path
        Directory for ``.gpw`` restart files (default ``gpw_files/``).
    gpw_logs : Path
        Directory for GPAW log files (default ``gpw_logs/``).

    Returns
    -------
    atoms : ase.Atoms
        Converged atoms with attached calculator.
    initial_id : int
        Database row ID of the initial structure entry.
    converged_id : int
        Database row ID of the converged structure entry.
    """
    db = _resolve_db(db)
    db_params = _serialize_calc_params(calc_params)
    db_params['legacy_gpaw'] = legacy_gpaw

    # Write initial structure first so we have a DB ID for naming the log file.
    initial_id = db.write(atoms_initial, system=system, **db_params)

    log_path = Path(gpw_logs) / f'{initial_id}.txt'
    log_path.parent.mkdir(parents=True, exist_ok=True)

    atoms = atoms_initial.copy()
    if legacy_gpaw:
        gpaw_params = {k: v for k, v in calc_params.items()
                       if k in GPAW.default_parameters}
        magmoms = calc_params.get('magmoms')
        if magmoms is not None:
            atoms.set_initial_magnetic_moments(magmoms)
        calc = GPAW(**gpaw_params, txt=str(log_path))
    else:
        gpaw_params = {k: v for k, v in calc_params.items()
                       if k in _NEW_GPAW_PARAMS}
        calc = _NewGPAW(**gpaw_params, txt=str(log_path))
    atoms.calc = calc
    atoms.get_potential_energy()
    atoms.get_forces()
    if calc.get_number_of_spins() > 1:
        atoms.get_magnetic_moments()

    convergence_data = extract_scf_convergence(log_path)

    data = {
        'system': system,
        'initial_id': initial_id,
        'data': {
            'scf_iter': [d['iter'] for d in convergence_data],
            'scf_energies': [d['energy'] for d in convergence_data],
            'scf_log10_eigst': [d['log10_eigst'] for d in convergence_data],
            'scf_log10_dens': [d['log10_dens'] for d in convergence_data],
            'scf_magmoms': [d['magmom'] for d in convergence_data],
        },
        **db_params,
    }

    converged_id = db.write(atoms, **data)

    # Rename log from initial_id → converged_id now that both are known.
    log_path.rename(Path(gpw_logs) / f'{converged_id}.txt')

    if save_gpw:
        gpw_file = Path(gpw_dir) / f'{converged_id}.gpw'
        gpw_file.parent.mkdir(parents=True, exist_ok=True)
        calc.write(str(gpw_file), mode=save_gpw_mode)
        db.update(converged_id, gpw_file=str(gpw_file))
        atoms.info.setdefault('key_value_pairs', {})['gpw_file'] = str(gpw_file)

    db.update(initial_id, converged_id=converged_id)

    return atoms, initial_id, converged_id


def load_gpaw_calculation(selector, db=None, system=None, legacy_gpaw=None,
                          gpw_logs=_DEFAULT_GPW_LOGS):
    """Load a previously stored calculation from the ASE database.

    Parameters
    ----------
    selector : dict
        calc_params dict (or subset) used to identify the calculation.
        Serialised the same way as in ``run_and_store_gpaw_calculation``.
    db : ase.db.core.Database or str or Path or None
        Database to search.  Accepts an already-connected ASE database
        object, a file path (str or Path) to connect to, or ``None`` to
        use the default ``calculations.db`` in the working directory.
    system : str, optional
        System label to narrow the search.
    legacy_gpaw : bool or None
        Filter by old (``True``) or new (``False``) GPAW implementation.
        When ``None`` (default) the value stored in the DB is used, defaulting
        to ``True`` for rows written before this field was introduced.
    gpw_logs : Path
        Directory where log files are stored (default ``gpw_logs/``).

    Returns
    -------
    atoms_converged : ase.Atoms
        Converged atoms with additional info attached.
    calc : GPAW
        Calculator loaded from the GPW file.
    """
    db = _resolve_db(db)
    initial_id = None
    converged_id = None
    atoms_converged = None

    extra = _serialize_calc_params(selector)
    if system is not None:
        extra['system'] = system
    if legacy_gpaw is not None:
        extra['legacy_gpaw'] = legacy_gpaw
    rows = db.select(**extra)

    for row in rows:
        if 'initial_id' in row.key_value_pairs:
            assert initial_id is None
            initial_id = row.initial_id
            atoms_converged = row.toatoms(add_additional_information=True)
        if 'converged_id' in row.key_value_pairs:
            assert converged_id is None
            converged_id = row.converged_id

    if atoms_converged is None:
        raise LookupError(
            f'No calculation found in DB for system={system!r}, '
            f'selector={selector}'
        )

    atoms_converged.info['key_value_pairs']['db_id'] = converged_id

    try:
        gpw_file = atoms_converged.info['key_value_pairs']['gpw_file']
    except KeyError:
        raise FileNotFoundError(
            f'Calculation (DB id={converged_id}) was stored without a GPW '
            'file. Re-run with save_gpw=True to save the calculator.'
        )

    if not Path(gpw_file).exists():
        raise FileNotFoundError(
            f'GPW file recorded in DB no longer exists on disk: {gpw_file}'
        )

    kv = atoms_converged.info['key_value_pairs']
    use_legacy = kv.get('legacy_gpaw', True)

    log_path = Path(gpw_logs) / f'{Path(gpw_file).stem}.txt'
    if use_legacy:
        calc = GPAW(str(gpw_file), txt=str(log_path))
    else:
        calc = _NewGPAW(str(gpw_file), txt=str(log_path))

    return atoms_converged, calc
