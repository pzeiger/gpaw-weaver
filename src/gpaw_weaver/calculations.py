import hashlib
import inspect
import json
from pathlib import Path

import numpy as np

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


def _calculation_hash(atoms, calc_params=None):
    """Return a hex digest that uniquely identifies a calculation.

    Hashes atomic numbers, positions (rounded to 6 decimal places),
    cell vectors, periodic boundary conditions, magnetic moments, and
    the full serialised calc_params so that different structures, magnetic
    configurations, and calculation settings (XC functional, k-points, …)
    all produce distinct values.

    If calc_params contains a 'magmoms' key it takes precedence over
    atoms.get_initial_magnetic_moments() for the magnetic moment contribution.
    """
    if calc_params is not None and calc_params.get('magmoms') is not None:
        magmoms = np.asarray(calc_params['magmoms'], dtype=float)
    else:
        magmoms = atoms.get_initial_magnetic_moments()
    h = hashlib.sha256()
    h.update(atoms.numbers.tobytes())
    h.update(np.round(atoms.positions, decimals=6).tobytes())
    h.update(np.round(atoms.cell[:], decimals=6).tobytes())
    h.update(atoms.pbc.tobytes())
    h.update(np.round(magmoms, decimals=6).tobytes())
    if calc_params is not None:
        serialized = _serialize_calc_params(calc_params)
        for key in sorted(serialized):
            h.update(key.encode())
            h.update(str(serialized[key]).encode())
    return h.hexdigest()


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


def run_and_store_gpaw_calculation(atoms_initial, calc_params,
                                   db=None,
                                   label=None,
                                   save_gpw=False,
                                   save_gpw_mode='calculation',
                                   legacy_gpaw=True,
                                   gpw_dir=_DEFAULT_GPW_DIR,
                                   gpw_logs=_DEFAULT_GPW_LOGS):
    """Run a GPAW calculation and store results in the ASE database.

    Log and GPW files are named after the calculation hash so they are
    uniquely identifiable without requiring a database row ID:

    * Log  → ``<gpw_logs>/<hash>.txt``
    * GPW  → ``<gpw_dir>/<hash>.gpw``

    Database entries are written only after the calculation converges
    successfully, so a failed run leaves no partial entries in the database.

    Parameters
    ----------
    atoms_initial : ase.Atoms
        Starting geometry.
    calc_params : dict
        All GPAW parameters for this calculation. Use ``make_pw_params`` or
        ``make_fd_params`` to build this, then add system-specific settings
        (``nbands``, ``charge``, ``setups``, ``poissonsolver``, …) directly.
        Every key is serialised and stored in the database.
    db : ase.db.core.Database or str or Path or None
        Database to write results into.  Accepts an already-connected ASE
        database object, a file path (str or Path) to connect to, or
        ``None`` to use the default ``calculations.db`` in the working
        directory.
    label : str or None
        Optional human-readable label stored alongside the calculation.
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
    if label is not None:
        db_params['label'] = label
    calc_hash = _calculation_hash(atoms_initial, calc_params)
    db_params['atoms_hash'] = calc_hash

    log_path = Path(gpw_logs) / f'{calc_hash}.txt'
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

    # Write DB entries only after successful convergence.
    initial_id = db.write(atoms_initial, **db_params)

    data = {
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

    if save_gpw:
        gpw_file = Path(gpw_dir) / f'{calc_hash}.gpw'
        gpw_file.parent.mkdir(parents=True, exist_ok=True)
        calc.write(str(gpw_file), mode=save_gpw_mode)
        db.update(converged_id, gpw_file=str(gpw_file))
        atoms.info.setdefault('key_value_pairs', {})['gpw_file'] = str(gpw_file)

    db.update(initial_id, converged_id=converged_id)

    return atoms, initial_id, converged_id


def query_gpaw_calculations(atoms_initial, db=None, calc_params=None):
    """Return all DB rows whose atoms structure and calc_params match the supplied values.

    Parameters present in *calc_params* must match exactly; parameters stored in
    the database but absent from *calc_params* are ignored (partial match).
    The atoms structure (atomic numbers, positions, cell, pbc, and magnetic
    moments) must also match exactly.

    Parameters
    ----------
    atoms_initial : ase.Atoms
        Reference structure to match against.
    db : ase.db.core.Database or str or Path or None
        Database to search.
    calc_params : dict or None
        Subset of calculation parameters to filter by.  Only keys present here
        are required to match.  When ``None`` only the atomic structure is used
        as a filter.

    Returns
    -------
    list of ase.db.row.AtomsRow
        All rows whose atoms structure and supplied calc_params match.
    """
    db = _resolve_db(db)

    extra = _serialize_calc_params(calc_params) if calc_params is not None else {}
    rows = list(db.select(**extra))

    ref_numbers = atoms_initial.numbers
    ref_positions = np.round(atoms_initial.positions, decimals=6)
    ref_cell = np.round(atoms_initial.cell[:], decimals=6)
    ref_pbc = atoms_initial.pbc
    if calc_params is not None and calc_params.get('magmoms') is not None:
        ref_magmoms = np.round(np.asarray(calc_params['magmoms'], dtype=float), decimals=6)
    else:
        ref_magmoms = np.round(atoms_initial.get_initial_magnetic_moments(), decimals=6)

    matched = []
    for row in rows:
        row_atoms = row.toatoms()
        if not np.array_equal(row_atoms.numbers, ref_numbers):
            continue
        if not np.array_equal(np.round(row_atoms.positions, decimals=6), ref_positions):
            continue
        if not np.array_equal(np.round(row_atoms.cell[:], decimals=6), ref_cell):
            continue
        if not np.array_equal(row_atoms.pbc, ref_pbc):
            continue
        row_magmoms = np.round(row_atoms.get_initial_magnetic_moments(), decimals=6)
        if not np.array_equal(row_magmoms, ref_magmoms):
            continue
        matched.append(row)

    return matched


def load_gpaw_calculation(atoms_initial, calc_params,
                          db=None, legacy_gpaw=None,
                          gpw_logs=_DEFAULT_GPW_LOGS,
                          txt='gpaw_log.txt'):
    """Load a previously stored calculation from the ASE database.

    Parameters
    ----------
    atoms_initial : ase.Atoms
        The initial structure passed to ``run_and_store_gpaw_calculation``.
        Its hash is used to identify the matching database entry, so
        different phases of the same element are distinguished correctly.
    calc_params : dict or None
        The same ``calc_params`` dict used in
        ``run_and_store_gpaw_calculation``.  All keys (XC functional,
        k-points, magmoms, …) are included in the hash, so passing the
        correct dict is required to uniquely identify the calculation when
        multiple runs on the same structure exist.  When ``None`` (default)
        only the atomic structure and ``atoms_initial.get_initial_magnetic_moments()``
        contribute to the hash.
    db : ase.db.core.Database or str or Path or None
        Database to search.  Accepts an already-connected ASE database
        object, a file path (str or Path) to connect to, or ``None`` to
        use the default ``calculations.db`` in the working directory.
    legacy_gpaw : bool or None
        Filter by old (``True``) or new (``False``) GPAW implementation.
        When ``None`` (default) the value stored in the DB is used, defaulting
        to ``True`` for rows written before this field was introduced.
    gpw_logs : Path
        Directory where log files are stored (default ``gpw_logs/``).
    txt : str
        Logfile name for the loaded calculator (default ``gpaw_log.txt``).

    Returns
    -------
    atoms_converged : ase.Atoms
        Converged atoms with additional info attached.
    calc : GPAW
        Calculator loaded from the GPW file.
    """
    db = _resolve_db(db)

    atoms_hash = _calculation_hash(atoms_initial, calc_params)
    extra = {'atoms_hash': atoms_hash}
    if legacy_gpaw is not None:
        extra['legacy_gpaw'] = legacy_gpaw
    rows = list(db.select(**extra))

    converged_rows = [r for r in rows if 'initial_id' in r.key_value_pairs]
    initial_rows = [r for r in rows if 'converged_id' in r.key_value_pairs]

    if len(converged_rows) > 1:
        raise ValueError(
            f'{len(converged_rows)} converged rows match '
            f'atoms_hash={atoms_hash!r}. '
            'Narrow the search with legacy_gpaw or a label stored on the entries.'
        )
    if not converged_rows:
        raise LookupError(
            f'No calculation found in DB for atoms_hash={atoms_hash!r}'
        )

    converged_row = converged_rows[0]
    atoms_converged = converged_row.toatoms(add_additional_information=True)
    converged_id = initial_rows[0].converged_id if initial_rows else None

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

    if use_legacy:
        calc = GPAW(str(gpw_file), txt=txt)
    else:
        calc = _NewGPAW(str(gpw_file), txt=txt)

    return atoms_converged, calc


def delete_gpaw_calculation(atoms_initial=None, calc_params=None,
                            db=None, legacy_gpaw=None, row_id=None):
    """Delete a previously stored calculation from the ASE database.

    Deletes all database rows (initial and converged) that match the
    atoms structure and calc_params, or deletes a calculation by its
    row ID.

    Parameters
    ----------
    atoms_initial : ase.Atoms or None
        The initial structure passed to ``run_and_store_gpaw_calculation``.
        Its hash is used to identify the matching database entries.
        Required unless ``row_id`` is provided.
    calc_params : dict or None
        The same ``calc_params`` dict used in
        ``run_and_store_gpaw_calculation``. All keys are included in the hash.
        When ``None`` only the atomic structure contributes to the hash.
    db : ase.db.core.Database or str or Path or None
        Database to search. Accepts an already-connected ASE database
        object, a file path (str or Path) to connect to, or ``None`` to
        use the default ``calculations.db`` in the working directory.
    legacy_gpaw : bool or None
        Filter by old (``True``) or new (``False``) GPAW implementation.
        When ``None`` (default) matches both. Only used when ``row_id`` is not
        provided.
    row_id : int or None
        Delete the row with this database ID. If the row belongs to a
        stored calculation pair, the linked row is also deleted.

    Returns
    -------
    int
        Number of rows deleted.
    """
    db = _resolve_db(db)

    if row_id is not None:
        row = db.get(row_id)
        if row is None:
            raise LookupError(f'No database row found with id={row_id!r}')
        rows = [row]
    else:
        if atoms_initial is None:
            raise ValueError('atoms_initial is required when row_id is not provided')
        atoms_hash = _calculation_hash(atoms_initial, calc_params)
        extra = {'atoms_hash': atoms_hash}
        if legacy_gpaw is not None:
            extra['legacy_gpaw'] = legacy_gpaw
        rows = list(db.select(**extra))

    deleted_ids = set()
    for row in rows:
        deleted_ids.add(row.id)
        if 'converged_id' in row.key_value_pairs:
            deleted_ids.add(row.key_value_pairs['converged_id'])
        if 'initial_id' in row.key_value_pairs:
            deleted_ids.add(row.key_value_pairs['initial_id'])

    deleted_count = 0
    for row_id in deleted_ids:
        db.delete([row_id])
        deleted_count += 1

    return deleted_count


def list_gpaw_calculations(atoms_initial=None, calc_params=None, db=None,
                           columns=None):
    """List stored calculations in the database with specified columns.

    Returns a list of dictionaries, each containing the database ID and
    the requested columns from the calculation metadata.

    Parameters
    ----------
    atoms_initial : ase.Atoms or None
        Reference structure to filter by. If None, list all calculations.
    calc_params : dict or None
        Subset of calculation parameters to filter by. Only keys present here
        are required to match. If None, no parameter filtering.
    db : ase.db.core.Database or str or Path or None
        Database to search. Accepts an already-connected ASE database
        object, a file path (str or Path) to connect to, or None to
        use the default ``calculations.db`` in the working directory.
    columns : list of str or None
        List of column names to include in the output. If None, includes
        all available key-value pairs.

    Returns
    -------
    list of dict
        Each dict contains 'id' (database row ID) and the requested columns.
    """
    db = _resolve_db(db)

    if atoms_initial is None and calc_params is None:
        # List all rows
        rows = list(db.select())
    else:
        # Use query_gpaw_calculations to filter
        rows = query_gpaw_calculations(atoms_initial, db=db, calc_params=calc_params)

    result = []
    for row in rows:
        entry = {'id': row.id}
        kv = row.key_value_pairs
        if columns is None:
            entry.update(kv)
        else:
            for col in columns:
                if col in kv:
                    entry[col] = kv[col]
        result.append(entry)

    return result
