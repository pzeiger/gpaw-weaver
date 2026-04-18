import pytest
from ase import Atoms
from ase.db import connect

from gpaw_weaver import make_pw_params


@pytest.fixture
def fe_atom():
    return Atoms("Fe", positions=[(0, 0, 0)], cell=[2.87, 2.87, 2.87], pbc=True)


@pytest.fixture
def pw_params():
    return make_pw_params(ecut=400, kpts={"size": (4, 4, 4), "gamma": True})


@pytest.fixture
def db(tmp_path):
    return connect(str(tmp_path / "test.db"))


@pytest.fixture
def work_dirs(tmp_path):
    return tmp_path / "gpw_files", tmp_path / "gpw_logs"
