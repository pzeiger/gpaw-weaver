"""Integration tests for run_and_store_gpaw_calculation / load_gpaw_calculation.

Each test is parametrized over legacy_gpaw=True/False so every scenario runs
against both GPAW implementations.  The ASE database and file paths are
isolated in pytest's tmp_path.

Scenarios covered (× both legacy_gpaw values):
  - non-magnetic
  - collinear spin polarization
  - non-collinear magnetism
  - save_gpw=True  (GPW file written, DB updated)
  - load roundtrip (run then load_gpaw_calculation)
"""
import json
from unittest.mock import patch

import pytest
from ase.calculators.calculator import Calculator, all_changes

from gpaw_weaver.calculations import (
    delete_gpaw_calculation,
    list_gpaw_calculations,
    load_gpaw_calculation,
    query_gpaw_calculations,
    run_and_store_gpaw_calculation,
)
from helpers import make_fake_gpaw_class, make_log


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

LEGACY = pytest.mark.parametrize("legacy_gpaw", [True, False], ids=["legacy", "new"])


def _run(fe_atom, pw_params, db, work_dirs, *, n_spins, magmom_str, legacy_gpaw,
         save_gpw=False):
    """Patch GPAW and call run_and_store_gpaw_calculation."""
    gpw_dir, gpw_logs = work_dirs
    FakeGPAW = make_fake_gpaw_class(
        n_spins=n_spins,
        log_content=make_log(n_iters=3, magmom=magmom_str),
    )
    with patch("gpaw_weaver.calculations.GPAW", FakeGPAW), \
         patch("gpaw_weaver.calculations._NewGPAW", FakeGPAW):
        return run_and_store_gpaw_calculation(
            fe_atom, pw_params, db=db,
            legacy_gpaw=legacy_gpaw, save_gpw=save_gpw,
            gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )


def _assert_common(atoms, initial_id, converged_id, db, gpw_logs, legacy_gpaw):
    """Assertions that hold for every scenario."""
    assert initial_id != converged_id
    assert initial_id >= 1 and converged_id >= 1

    init_row = db.get(id=initial_id)
    assert init_row.key_value_pairs["converged_id"] == converged_id

    conv_row = db.get(id=converged_id)
    assert conv_row.key_value_pairs["initial_id"] == initial_id
    assert conv_row.key_value_pairs["legacy_gpaw"] == legacy_gpaw

    data = conv_row.data
    assert data["scf_iter"] == [1, 2, 3]
    assert len(data["scf_energies"]) == 3
    assert len(data["scf_log10_eigst"]) == 3
    assert len(data["scf_log10_dens"]) == 3

    calc_hash = conv_row.key_value_pairs["atoms_hash"]
    assert (gpw_logs / f"{calc_hash}.txt").exists()


# ---------------------------------------------------------------------------
# Scenario tests
# ---------------------------------------------------------------------------

@LEGACY
def test_nonmagnetic(fe_atom, pw_params, db, work_dirs, legacy_gpaw):
    _, gpw_logs = work_dirs
    atoms, initial_id, converged_id = _run(
        fe_atom, pw_params, db, work_dirs,
        n_spins=1, magmom_str=None, legacy_gpaw=legacy_gpaw,
    )
    _assert_common(atoms, initial_id, converged_id, db, gpw_logs, legacy_gpaw)

    assert all(m is None for m in db.get(id=converged_id).data["scf_magmoms"])


@LEGACY
def test_collinear_magnetic(fe_atom, pw_params, db, work_dirs, legacy_gpaw):
    _, gpw_logs = work_dirs
    atoms, initial_id, converged_id = _run(
        fe_atom, pw_params, db, work_dirs,
        n_spins=2, magmom_str="+1.0000", legacy_gpaw=legacy_gpaw,
    )
    _assert_common(atoms, initial_id, converged_id, db, gpw_logs, legacy_gpaw)

    magmoms = db.get(id=converged_id).data["scf_magmoms"]
    assert all(isinstance(m, float) for m in magmoms)
    assert all(m == pytest.approx(1.0) for m in magmoms)


@LEGACY
def test_noncollinear_magnetic(fe_atom, pw_params, db, work_dirs, legacy_gpaw):
    _, gpw_logs = work_dirs
    atoms, initial_id, converged_id = _run(
        fe_atom, pw_params, db, work_dirs,
        n_spins=2, magmom_str="+1.000,+0.000,-1.000", legacy_gpaw=legacy_gpaw,
    )
    _assert_common(atoms, initial_id, converged_id, db, gpw_logs, legacy_gpaw)

    magmoms = db.get(id=converged_id).data["scf_magmoms"]
    assert all(isinstance(m, list) and len(m) == 3 for m in magmoms)
    assert magmoms[0] == pytest.approx([1.0, 0.0, -1.0])


@LEGACY
def test_save_gpw_writes_file_and_updates_db(fe_atom, pw_params, db, work_dirs, legacy_gpaw):
    gpw_dir, gpw_logs = work_dirs
    atoms, initial_id, converged_id = _run(
        fe_atom, pw_params, db, work_dirs,
        n_spins=1, magmom_str=None, legacy_gpaw=legacy_gpaw, save_gpw=True,
    )
    conv_row = db.get(id=converged_id)
    calc_hash = conv_row.key_value_pairs["atoms_hash"]
    expected_gpw = gpw_dir / f"{calc_hash}.gpw"
    assert expected_gpw.exists()
    assert conv_row.key_value_pairs["gpw_file"] == str(expected_gpw)


@LEGACY
def test_load_calculation_roundtrip(fe_atom, pw_params, db, work_dirs, legacy_gpaw):
    """Run with save_gpw=True then load; verify atoms and calc are returned."""
    gpw_dir, gpw_logs = work_dirs
    FakeGPAW = make_fake_gpaw_class(
        n_spins=1,
        log_content=make_log(n_iters=3, magmom=None),
    )
    with patch("gpaw_weaver.calculations.GPAW", FakeGPAW), \
         patch("gpaw_weaver.calculations._NewGPAW", FakeGPAW):
        _, _, converged_id = run_and_store_gpaw_calculation(
            fe_atom, pw_params, db=db,
            legacy_gpaw=legacy_gpaw, save_gpw=True,
            gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )
        atoms_loaded, calc_loaded = load_gpaw_calculation(
            fe_atom, db=db, calc_params=pw_params, gpw_logs=gpw_logs,
        )

    assert atoms_loaded is not None
    assert calc_loaded is not None
    assert atoms_loaded.info["key_value_pairs"]["db_id"] == converged_id
    assert atoms_loaded.get_chemical_symbols() == ["Fe"]
    assert atoms_loaded.cell.lengths() == pytest.approx([2.87, 2.87, 2.87])


def test_load_distinguishes_different_structures(pw_params, db, work_dirs):
    """Two different structures stored under the same system are loaded independently."""
    from ase import Atoms
    gpw_dir, gpw_logs = work_dirs
    bcc = Atoms("Fe", positions=[(0, 0, 0)], cell=[2.87, 2.87, 2.87], pbc=True)
    fcc = Atoms("Fe", positions=[(0, 0, 0)], cell=[3.52, 3.52, 3.52], pbc=True)
    FakeGPAW = make_fake_gpaw_class(n_spins=1, log_content=make_log(n_iters=3))
    with patch("gpaw_weaver.calculations.GPAW", FakeGPAW), \
         patch("gpaw_weaver.calculations._NewGPAW", FakeGPAW):
        _, _, bcc_id = run_and_store_gpaw_calculation(
            bcc, pw_params, db=db, save_gpw=True,
            gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )
        run_and_store_gpaw_calculation(
            fcc, pw_params, db=db, save_gpw=True,
            gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )
        atoms_loaded, _ = load_gpaw_calculation(bcc, db=db, calc_params=pw_params, gpw_logs=gpw_logs)
    assert atoms_loaded.info["key_value_pairs"]["db_id"] == bcc_id


def test_load_distinguishes_magnetic_configurations(pw_params, db, work_dirs):
    """FM and AFM via atoms magmoms are loaded independently."""
    from ase import Atoms
    gpw_dir, gpw_logs = work_dirs
    fm  = Atoms("Fe2", positions=[(0,0,0),(1.435,1.435,1.435)],
                cell=[2.87, 2.87, 2.87], pbc=True)
    afm = fm.copy()
    fm.set_initial_magnetic_moments([2.0, 2.0])
    afm.set_initial_magnetic_moments([2.0, -2.0])
    FakeGPAW = make_fake_gpaw_class(n_spins=2, log_content=make_log(n_iters=3, magmom="+2.0000"))
    with patch("gpaw_weaver.calculations.GPAW", FakeGPAW), \
         patch("gpaw_weaver.calculations._NewGPAW", FakeGPAW):
        _, _, fm_id = run_and_store_gpaw_calculation(
            fm, pw_params, db=db, save_gpw=True,
            gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )
        run_and_store_gpaw_calculation(
            afm, pw_params, db=db, save_gpw=True,
            gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )
        atoms_loaded, _ = load_gpaw_calculation(fm, db=db, calc_params=pw_params, gpw_logs=gpw_logs)
    assert atoms_loaded.info["key_value_pairs"]["db_id"] == fm_id


def test_load_distinguishes_calc_params_magmoms(pw_params, db, work_dirs):
    """FM and AFM specified via calc_params magmoms are loaded independently."""
    from ase import Atoms
    gpw_dir, gpw_logs = work_dirs
    base = Atoms("Fe2", positions=[(0,0,0),(1.435,1.435,1.435)],
                 cell=[2.87, 2.87, 2.87], pbc=True)
    fm_params  = {**pw_params, 'magmoms': [2.0,  2.0]}
    afm_params = {**pw_params, 'magmoms': [2.0, -2.0]}
    FakeGPAW = make_fake_gpaw_class(n_spins=2, log_content=make_log(n_iters=3, magmom="+2.0000"))
    with patch("gpaw_weaver.calculations.GPAW", FakeGPAW), \
         patch("gpaw_weaver.calculations._NewGPAW", FakeGPAW):
        _, _, fm_id = run_and_store_gpaw_calculation(
            base, fm_params, db=db, save_gpw=True,
            gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )
        run_and_store_gpaw_calculation(
            base, afm_params, db=db, save_gpw=True,
            gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )
        atoms_loaded, _ = load_gpaw_calculation(
            base, db=db, calc_params=fm_params, gpw_logs=gpw_logs,
        )
    assert atoms_loaded.info["key_value_pairs"]["db_id"] == fm_id


def test_load_distinguishes_calc_params(fe_atom, pw_params, db, work_dirs):
    """Same atoms with different XC functionals are loaded independently."""
    gpw_dir, gpw_logs = work_dirs
    pbe_params  = {**pw_params, 'xc': 'PBE'}
    lda_params  = {**pw_params, 'xc': 'LDA'}
    FakeGPAW = make_fake_gpaw_class(n_spins=1, log_content=make_log(n_iters=3))
    with patch("gpaw_weaver.calculations.GPAW", FakeGPAW), \
         patch("gpaw_weaver.calculations._NewGPAW", FakeGPAW):
        _, _, pbe_id = run_and_store_gpaw_calculation(
            fe_atom, pbe_params, db=db, save_gpw=True,
            gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )
        run_and_store_gpaw_calculation(
            fe_atom, lda_params, db=db, save_gpw=True,
            gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )
        atoms_loaded, _ = load_gpaw_calculation(
            fe_atom, db=db, calc_params=pbe_params, gpw_logs=gpw_logs,
        )
    assert atoms_loaded.info["key_value_pairs"]["db_id"] == pbe_id


def test_load_raises_on_ambiguous_match(fe_atom, pw_params, db, work_dirs):
    """load_gpaw_calculation raises ValueError when multiple rows match."""
    gpw_dir, gpw_logs = work_dirs
    FakeGPAW = make_fake_gpaw_class(n_spins=1, log_content=make_log(n_iters=3))
    with patch("gpaw_weaver.calculations.GPAW", FakeGPAW), \
         patch("gpaw_weaver.calculations._NewGPAW", FakeGPAW):
        run_and_store_gpaw_calculation(
            fe_atom, pw_params, db=db, save_gpw=True,
            gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )
        run_and_store_gpaw_calculation(
            fe_atom, pw_params, db=db, save_gpw=True,
            gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )
        with pytest.raises(ValueError, match="2 converged rows match"):
            load_gpaw_calculation(fe_atom, db=db, calc_params=pw_params, gpw_logs=gpw_logs)


# ---------------------------------------------------------------------------
# _resolve_db tests
# ---------------------------------------------------------------------------

def test_db_accepts_path(fe_atom, pw_params, work_dirs, tmp_path):
    """Passing a Path connects automatically to that file."""
    gpw_dir, gpw_logs = work_dirs
    db_path = tmp_path / "custom.db"
    FakeGPAW = make_fake_gpaw_class(n_spins=1, log_content=make_log(n_iters=3))
    with patch("gpaw_weaver.calculations.GPAW", FakeGPAW), \
         patch("gpaw_weaver.calculations._NewGPAW", FakeGPAW):
        _, _, converged_id = run_and_store_gpaw_calculation(
            fe_atom, pw_params, db=db_path,
            gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )
    assert db_path.exists()
    from ase.db import connect
    assert connect(str(db_path)).get(id=converged_id) is not None


def test_db_accepts_str(fe_atom, pw_params, work_dirs, tmp_path):
    """Passing a str path connects automatically to that file."""
    gpw_dir, gpw_logs = work_dirs
    db_path = str(tmp_path / "custom.db")
    FakeGPAW = make_fake_gpaw_class(n_spins=1, log_content=make_log(n_iters=3))
    with patch("gpaw_weaver.calculations.GPAW", FakeGPAW), \
         patch("gpaw_weaver.calculations._NewGPAW", FakeGPAW):
        _, _, converged_id = run_and_store_gpaw_calculation(
            fe_atom, pw_params, db=db_path,
            gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )
    from ase.db import connect
    assert connect(db_path).get(id=converged_id) is not None


def test_db_path_no_extension(fe_atom, pw_params, work_dirs, tmp_path):
    """A path without an extension gets .db appended automatically."""
    gpw_dir, gpw_logs = work_dirs
    db_path = tmp_path / "myproject"
    FakeGPAW = make_fake_gpaw_class(n_spins=1, log_content=make_log(n_iters=3))
    with patch("gpaw_weaver.calculations.GPAW", FakeGPAW), \
         patch("gpaw_weaver.calculations._NewGPAW", FakeGPAW):
        run_and_store_gpaw_calculation(
            fe_atom, pw_params, db=db_path,
            gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )
    assert (tmp_path / "myproject.db").exists()


# ---------------------------------------------------------------------------
# query_gpaw_calculations tests
# ---------------------------------------------------------------------------

def test_query_returns_all_matching_rows(fe_atom, pw_params, db, work_dirs):
    """Returns all rows (initial + converged) whose structure and params match."""
    gpw_dir, gpw_logs = work_dirs
    FakeGPAW = make_fake_gpaw_class(n_spins=1, log_content=make_log(n_iters=3))
    with patch("gpaw_weaver.calculations.GPAW", FakeGPAW), \
         patch("gpaw_weaver.calculations._NewGPAW", FakeGPAW):
        run_and_store_gpaw_calculation(
            fe_atom, pw_params, db=db, gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )
    rows = query_gpaw_calculations(fe_atom, db=db, calc_params=pw_params)
    assert len(rows) == 2  # initial + converged


def test_query_partial_params_match(fe_atom, pw_params, db, work_dirs):
    """Partial calc_params matches rows regardless of unspecified params."""
    from ase import Atoms
    gpw_dir, gpw_logs = work_dirs
    pbe_params = {**pw_params, 'xc': 'PBE'}
    lda_params = {**pw_params, 'xc': 'LDA'}
    FakeGPAW = make_fake_gpaw_class(n_spins=1, log_content=make_log(n_iters=3))
    with patch("gpaw_weaver.calculations.GPAW", FakeGPAW), \
         patch("gpaw_weaver.calculations._NewGPAW", FakeGPAW):
        run_and_store_gpaw_calculation(
            fe_atom, pbe_params, db=db, gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )
        run_and_store_gpaw_calculation(
            fe_atom, lda_params, db=db, gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )
    # Querying with only xc='PBE' should return just the PBE rows
    rows = query_gpaw_calculations(fe_atom, db=db, calc_params={'xc': 'PBE'})
    assert len(rows) == 2
    assert all(r.key_value_pairs.get('xc') == 'PBE' for r in rows)


def test_query_excludes_different_structure(pw_params, db, work_dirs):
    """Rows for a different atoms object are not returned."""
    from ase import Atoms
    gpw_dir, gpw_logs = work_dirs
    bcc = Atoms("Fe", positions=[(0, 0, 0)], cell=[2.87, 2.87, 2.87], pbc=True)
    fcc = Atoms("Fe", positions=[(0, 0, 0)], cell=[3.52, 3.52, 3.52], pbc=True)
    FakeGPAW = make_fake_gpaw_class(n_spins=1, log_content=make_log(n_iters=3))
    with patch("gpaw_weaver.calculations.GPAW", FakeGPAW), \
         patch("gpaw_weaver.calculations._NewGPAW", FakeGPAW):
        run_and_store_gpaw_calculation(
            bcc, pw_params, db=db, gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )
        run_and_store_gpaw_calculation(
            fcc, pw_params, db=db, gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )
    rows = query_gpaw_calculations(bcc, db=db)
    assert len(rows) == 2
    assert all(r.toatoms().cell.lengths()[0] == pytest.approx(2.87) for r in rows)


def test_query_no_calc_params_matches_by_structure_only(fe_atom, pw_params, db, work_dirs):
    """Passing calc_params=None matches all rows for the given structure."""
    gpw_dir, gpw_logs = work_dirs
    pbe_params = {**pw_params, 'xc': 'PBE'}
    lda_params = {**pw_params, 'xc': 'LDA'}
    FakeGPAW = make_fake_gpaw_class(n_spins=1, log_content=make_log(n_iters=3))
    with patch("gpaw_weaver.calculations.GPAW", FakeGPAW), \
         patch("gpaw_weaver.calculations._NewGPAW", FakeGPAW):
        run_and_store_gpaw_calculation(
            fe_atom, pbe_params, db=db, gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )
        run_and_store_gpaw_calculation(
            fe_atom, lda_params, db=db, gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )
    rows = query_gpaw_calculations(fe_atom, db=db, calc_params=None)
    assert len(rows) == 4  # 2 runs × 2 rows each


def test_db_default(fe_atom, pw_params, work_dirs, tmp_path, monkeypatch):
    """Omitting db writes to calculations.db in the working directory."""
    monkeypatch.chdir(tmp_path)
    gpw_dir, gpw_logs = work_dirs
    FakeGPAW = make_fake_gpaw_class(n_spins=1, log_content=make_log(n_iters=3))
    with patch("gpaw_weaver.calculations.GPAW", FakeGPAW), \
         patch("gpaw_weaver.calculations._NewGPAW", FakeGPAW):
        run_and_store_gpaw_calculation(
            fe_atom, pw_params,
            gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )
    assert (tmp_path / "calculations.db").exists()


@LEGACY
def test_delete_calculation(fe_atom, pw_params, db, work_dirs, legacy_gpaw):
    """Test deleting a calculation from the database."""
    gpw_dir, gpw_logs = work_dirs
    FakeGPAW = make_fake_gpaw_class(n_spins=1, log_content=make_log(n_iters=3))
    with patch("gpaw_weaver.calculations.GPAW", FakeGPAW), \
         patch("gpaw_weaver.calculations._NewGPAW", FakeGPAW):
        run_and_store_gpaw_calculation(
            fe_atom, pw_params, db=db,
            legacy_gpaw=legacy_gpaw, save_gpw=True,
            gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )
    # Verify it exists
    rows_before = query_gpaw_calculations(fe_atom, db=db, calc_params=pw_params)
    assert len(rows_before) == 2  # initial and converged rows

    # Delete it
    deleted_count = delete_gpaw_calculation(fe_atom, pw_params, db=db, legacy_gpaw=legacy_gpaw)
    assert deleted_count == 2  # initial and converged rows

    # Verify it's gone
    rows_after = query_gpaw_calculations(fe_atom, db=db, calc_params=pw_params)
    assert len(rows_after) == 0


@LEGACY
def test_delete_calculation_by_row_id(fe_atom, pw_params, db, work_dirs, legacy_gpaw):
    """Test deleting a calculation by database row ID."""
    gpw_dir, gpw_logs = work_dirs
    FakeGPAW = make_fake_gpaw_class(n_spins=1, log_content=make_log(n_iters=3))
    with patch("gpaw_weaver.calculations.GPAW", FakeGPAW), \
         patch("gpaw_weaver.calculations._NewGPAW", FakeGPAW):
        run_and_store_gpaw_calculation(
            fe_atom, pw_params, db=db,
            legacy_gpaw=legacy_gpaw, save_gpw=True,
            gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )

    rows_before = query_gpaw_calculations(fe_atom, db=db, calc_params=pw_params)
    assert len(rows_before) == 2

    row_id = rows_before[0].id
    deleted_count = delete_gpaw_calculation(db=db, row_id=row_id)
    assert deleted_count == 2

    rows_after = query_gpaw_calculations(fe_atom, db=db, calc_params=pw_params)
    assert len(rows_after) == 0


@LEGACY
def test_list_calculations(fe_atom, pw_params, db, work_dirs, legacy_gpaw):
    """Test listing calculations with specified columns."""
    gpw_dir, gpw_logs = work_dirs
    FakeGPAW = make_fake_gpaw_class(n_spins=1, log_content=make_log(n_iters=3))
    with patch("gpaw_weaver.calculations.GPAW", FakeGPAW), \
         patch("gpaw_weaver.calculations._NewGPAW", FakeGPAW):
        run_and_store_gpaw_calculation(
            fe_atom, pw_params, db=db,
            legacy_gpaw=legacy_gpaw, save_gpw=True,
            gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )

    # List all calculations
    all_list = list_gpaw_calculations(db=db)
    assert len(all_list) == 2  # initial and converged

    # Check structure
    for entry in all_list:
        assert 'id' in entry
        assert isinstance(entry['id'], int)

    # List with specific columns
    columns = ['atoms_hash', 'legacy_gpaw']
    filtered_list = list_gpaw_calculations(atoms_initial=fe_atom, calc_params=pw_params, db=db, columns=columns)
    assert len(filtered_list) == 2
    for entry in filtered_list:
        assert 'id' in entry
        for col in columns:
            assert col in entry

    # List with no matching
    empty_list = list_gpaw_calculations(atoms_initial=fe_atom, calc_params={'xc': 'LDA'}, db=db, columns=columns)
    assert len(empty_list) == 0


@LEGACY
def test_parallel_forwarded_but_excluded_from_identity(
        fe_atom, pw_params, work_dirs, tmp_path, legacy_gpaw):
    """`parallel` reaches the GPAW constructor but is not part of identity.

    It must be forwarded to the constructor's ``parallel`` keyword, yet
    excluded from the ``atoms_hash`` and the stored DB key-value pairs, so
    the same physical calculation run with different parallelization maps to
    the same hash.
    """
    from ase.db import connect

    gpw_dir, gpw_logs = work_dirs
    FakeGPAW = make_fake_gpaw_class(n_spins=1, log_content=make_log())
    parallel = {"augment_grids": True}

    db_par = connect(str(tmp_path / "par.db"))
    with patch("gpaw_weaver.calculations.GPAW", FakeGPAW), \
         patch("gpaw_weaver.calculations._NewGPAW", FakeGPAW):
        atoms_par, _, conv_par = run_and_store_gpaw_calculation(
            fe_atom, pw_params, db=db_par, legacy_gpaw=legacy_gpaw,
            parallel=parallel, gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )

    # (1) forwarded to the constructor
    assert atoms_par.calc.init_kwargs.get("parallel") == parallel

    db_plain = connect(str(tmp_path / "plain.db"))
    with patch("gpaw_weaver.calculations.GPAW", FakeGPAW), \
         patch("gpaw_weaver.calculations._NewGPAW", FakeGPAW):
        atoms_plain, _, conv_plain = run_and_store_gpaw_calculation(
            fe_atom, pw_params, db=db_plain, legacy_gpaw=legacy_gpaw,
            gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )

    # (2) default None means the constructor is not given a parallel kwarg
    assert "parallel" not in atoms_plain.calc.init_kwargs

    # (3) identity hash is identical with vs without parallel
    hash_par = db_par.get(id=conv_par).key_value_pairs["atoms_hash"]
    hash_plain = db_plain.get(id=conv_plain).key_value_pairs["atoms_hash"]
    assert hash_par == hash_plain

    # (4) parallel is not persisted as a DB key-value pair
    assert "parallel" not in db_par.get(id=conv_par).key_value_pairs


# ---------------------------------------------------------------------------
# van der Waals wrapper-correction hook (vdw_factory)
# ---------------------------------------------------------------------------

class _FakeVdw(Calculator):
    """Minimal wrapper calculator: inner DFT energy + a fixed dispersion delta.

    Stands in for a real wrapper scheme (Tkatchenko-Scheffler, DFT-D3/D4) in
    tests without pulling in a dispersion backend.
    """

    implemented_properties = ["energy", "forces"]
    name = "fakevdw"

    def __init__(self, dft, delta):
        super().__init__()
        self.dft = dft
        self.delta = delta

    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        # Record the state (as real ASE calculators do) so ase.db can read the
        # stored energy back without recomputing.
        self.atoms = atoms.copy()
        self.dft.calculate(atoms, ["energy", "forces"], system_changes)
        r = self.dft.results
        self.results = {"energy": r["energy"] + self.delta, "forces": r["forces"]}


def _vdw_factory(dft_calc, atoms, descriptor):
    return _FakeVdw(dft_calc, descriptor["delta"])


@LEGACY
def test_vdw_wrapper_applied_and_stored(fe_atom, pw_params, db, work_dirs, legacy_gpaw):
    """A `vdw` descriptor wraps the DFT calc; corrected energy + identity stored."""
    gpw_dir, gpw_logs = work_dirs
    vdw_params = {**pw_params, "vdw": {"name": "fake", "delta": -5.0}}
    FakeGPAW = make_fake_gpaw_class(n_spins=1, log_content=make_log(n_iters=3))
    with patch("gpaw_weaver.calculations.GPAW", FakeGPAW), \
         patch("gpaw_weaver.calculations._NewGPAW", FakeGPAW):
        atoms, initial_id, converged_id = run_and_store_gpaw_calculation(
            fe_atom, vdw_params, db=db, legacy_gpaw=legacy_gpaw,
            vdw_factory=_vdw_factory, save_gpw=True,
            gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )
    # corrected energy (a bare FakeGPAW would give -100.0)
    assert atoms.get_potential_energy() == pytest.approx(-105.0)
    conv_row = db.get(id=converged_id)
    assert conv_row.energy == pytest.approx(-105.0)
    # the vdw descriptor is part of the stored identity
    assert conv_row.key_value_pairs["vdw"] == json.dumps(
        {"name": "fake", "delta": -5.0}, sort_keys=True)
    # save_gpw wrote via the inner GPAW (the wrapper has no .write)
    calc_hash = conv_row.key_value_pairs["atoms_hash"]
    assert (gpw_dir / f"{calc_hash}.gpw").exists()


@LEGACY
def test_vdw_changes_identity_hash(fe_atom, pw_params, db, work_dirs, legacy_gpaw):
    """Same structure/params hash differently with vs without a vdw descriptor,
    and a plain (no-vdw) call keeps the exact pre-existing hash (back-compat)."""
    gpw_dir, gpw_logs = work_dirs
    vdw_params = {**pw_params, "vdw": {"name": "fake", "delta": -5.0}}
    FakeGPAW = make_fake_gpaw_class(n_spins=1, log_content=make_log(n_iters=3))
    with patch("gpaw_weaver.calculations.GPAW", FakeGPAW), \
         patch("gpaw_weaver.calculations._NewGPAW", FakeGPAW):
        _, _, plain_id = run_and_store_gpaw_calculation(
            fe_atom, pw_params, db=db, legacy_gpaw=legacy_gpaw,
            gpw_dir=gpw_dir, gpw_logs=gpw_logs)
        _, _, vdw_id = run_and_store_gpaw_calculation(
            fe_atom, vdw_params, db=db, legacy_gpaw=legacy_gpaw,
            vdw_factory=_vdw_factory, gpw_dir=gpw_dir, gpw_logs=gpw_logs)
    h_plain = db.get(id=plain_id).key_value_pairs["atoms_hash"]
    h_vdw = db.get(id=vdw_id).key_value_pairs["atoms_hash"]
    assert h_plain != h_vdw
    # a no-vdw calc_params must hash exactly as a bare make_pw_params would,
    # i.e. the vdw feature does not perturb existing entries' identity.
    from gpaw_weaver.calculations import _calculation_hash
    assert h_plain == _calculation_hash(fe_atom, pw_params)


@LEGACY
def test_vdw_without_factory_raises(fe_atom, pw_params, db, work_dirs, legacy_gpaw):
    """A `vdw` descriptor without a vdw_factory is a hard error, not silent."""
    gpw_dir, gpw_logs = work_dirs
    vdw_params = {**pw_params, "vdw": {"name": "fake", "delta": -5.0}}
    FakeGPAW = make_fake_gpaw_class(n_spins=1, log_content=make_log(n_iters=3))
    with patch("gpaw_weaver.calculations.GPAW", FakeGPAW), \
         patch("gpaw_weaver.calculations._NewGPAW", FakeGPAW):
        with pytest.raises(ValueError, match="vdw_factory"):
            run_and_store_gpaw_calculation(
                fe_atom, vdw_params, db=db, legacy_gpaw=legacy_gpaw,
                gpw_dir=gpw_dir, gpw_logs=gpw_logs)
