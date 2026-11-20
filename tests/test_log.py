from pathlib import Path

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


SCF_LOGS = Path(__file__).parent / "data" / "scf_logs"

# Logs written by GPAW @ 45d3025 (2-atom H2 and an O atom, PW 300 eV, Gamma only),
# trimmed to the SCF table: (file, iterations, first row, last row), each row as
# (energy, log10_eigst, log10_dens, magmom).
REAL_LOGS = [
    ("nonmag_legacy_True.txt", 8,
     (-7.522752, None, None, None), (-7.022996, -8.94, -4.28, None)),
    ("nonmag_legacy_False.txt", 10,
     (-7.522752, 1.03, None, None), (-7.022996, -10.08, -6.11, None)),
    ("spin_legacy_True.txt", 12,
     (5.744739, None, None, 1.4997), (5.938483, -8.37, -4.92, 2.0)),
    ("spin_legacy_False.txt", 11,
     (5.742937, -0.88, None, 1.4997), (5.938483, -9.15, -6.11, 2.0)),
    ("ncol_legacy_True.txt", 11,
     (5.865159, None, None, [0.5, 0.5, 1.0]), (5.938483, -7.92, -4.54, [0.8, 0.8, 1.6])),
    ("ncol_legacy_False.txt", 12,
     (5.863355, -0.88, None, [0.5, 0.5, 1.0]), (5.938483, -9.47, -6.14, [0.8, 0.8, 1.6])),
]


@pytest.mark.parametrize("name, n_iters, first, last", REAL_LOGS,
                         ids=[r[0][:-4] for r in REAL_LOGS])
def test_real_gpaw_logs(name, n_iters, first, last):
    iters = extract_scf_convergence(SCF_LOGS / name)

    assert [d["iter"] for d in iters] == list(range(1, n_iters + 1))
    for row, expected in ((iters[0], first), (iters[-1], last)):
        got = (row["energy"], row["log10_eigst"], row["log10_dens"], row["magmom"])
        for g, e in zip(got, expected):
            assert g == (None if e is None else pytest.approx(e))
