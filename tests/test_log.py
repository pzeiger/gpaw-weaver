import pytest

from gpaw_weaver.log import extract_scf_convergence
from helpers import make_log


def test_nonmagnetic(tmp_path):
    log = tmp_path / "calc.txt"
    log.write_text(make_log(n_iters=3, magmom=None))

    iters = extract_scf_convergence(log)

    assert len(iters) == 3
    assert [d["iter"] for d in iters] == [1, 2, 3]
    assert all(d["magmom"] is None for d in iters)
    assert iters[0]["energy"] == pytest.approx(-100.01)
    assert iters[2]["energy"] == pytest.approx(-100.03)
    assert iters[0]["log10_eigst"] == pytest.approx(-2.0)
    assert iters[0]["log10_dens"] == pytest.approx(-3.0)


def test_collinear_magnetic(tmp_path):
    log = tmp_path / "calc.txt"
    log.write_text(make_log(n_iters=2, magmom="+1.0000"))

    iters = extract_scf_convergence(log)

    assert len(iters) == 2
    assert all(isinstance(d["magmom"], float) for d in iters)
    assert all(d["magmom"] == pytest.approx(1.0) for d in iters)


def test_noncollinear_magnetic(tmp_path):
    log = tmp_path / "calc.txt"
    log.write_text(make_log(n_iters=2, magmom="+1.000,+0.000,-1.000"))

    iters = extract_scf_convergence(log)

    assert len(iters) == 2
    assert all(isinstance(d["magmom"], list) for d in iters)
    assert iters[0]["magmom"] == pytest.approx([1.0, 0.0, -1.0])


def test_empty_scf_block(tmp_path):
    log = tmp_path / "calc.txt"
    log.write_text(" iter  time      total     log10-change:\nConverged\n")

    assert extract_scf_convergence(log) == []


def test_stops_at_blank_line(tmp_path):
    """A blank line terminates the SCF block even without 'Converged'."""
    log = tmp_path / "calc.txt"
    log.write_text(
        " iter  time      total     log10-change:\n"
        "iter:   1  10:30:01  -100.0100   -2.000   -3.000\n"
        "\n"
        "iter:   2  10:30:02  -100.0200   -3.000   -4.000\n"
    )
    iters = extract_scf_convergence(log)
    assert len(iters) == 1
