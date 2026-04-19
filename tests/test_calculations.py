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
from unittest.mock import patch

import pytest

from gpaw_weaver.calculations import load_gpaw_calculation, run_and_store_gpaw_calculation
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
            fe_atom, pw_params, system="Fe", db=db,
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
    assert conv_row.key_value_pairs["system"] == "Fe"
    assert conv_row.key_value_pairs["legacy_gpaw"] == legacy_gpaw

    data = conv_row.data
    assert data["scf_iter"] == [1, 2, 3]
    assert len(data["scf_energies"]) == 3
    assert len(data["scf_log10_eigst"]) == 3
    assert len(data["scf_log10_dens"]) == 3

    assert (gpw_logs / f"{converged_id}.txt").exists()
    assert not (gpw_logs / f"{initial_id}.txt").exists()


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
    expected_gpw = gpw_dir / f"{converged_id}.gpw"
    assert expected_gpw.exists()
    assert db.get(id=converged_id).key_value_pairs["gpw_file"] == str(expected_gpw)


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
            fe_atom, pw_params, system="Fe", db=db,
            legacy_gpaw=legacy_gpaw, save_gpw=True,
            gpw_dir=gpw_dir, gpw_logs=gpw_logs,
        )
        atoms_loaded, calc_loaded = load_gpaw_calculation(
            pw_params, db, system="Fe", gpw_logs=gpw_logs,
        )

    assert atoms_loaded is not None
    assert calc_loaded is not None
    assert atoms_loaded.info["key_value_pairs"]["db_id"] == converged_id
    assert atoms_loaded.get_chemical_symbols() == ["Fe"]
    assert atoms_loaded.cell.lengths() == pytest.approx([2.87, 2.87, 2.87])
