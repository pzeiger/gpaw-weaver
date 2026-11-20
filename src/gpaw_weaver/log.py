import re

# A row of GPAW's SCF table: ``iter:   2 ...`` (legacy implementation) or
# ``|iter:   2| ...`` (new implementation, cells separated by ``|``).
_ROW = re.compile(r'\|?iter:\s*(\d+)\|? \d{2}:\d{2}:\d{2}(.*)')
# GPAW writes the collinear magmom as ``%+.4f`` and the noncollinear one as a
# ``%+.1f`` triplet; eigst, dens and custom criteria have fewer decimals.
_MAGMOM = re.compile(
    r'[-+]\d+\.\d{4}|[-+]?\d+\.\d+,[-+]?\d+\.\d+,[-+]?\d+\.\d+')


def _number(cell):
    cell = cell.strip().rstrip('c')
    return float(cell) if cell else None


def _magmom(cell):
    if ',' in cell:
        return [float(v) for v in cell.split(',')]
    return float(cell)


def _parse_row(rest):
    """Return (energy, log10_eigst, log10_dens, magmom) from the part of a row
    after its time stamp."""
    if '|' in rest:
        # New layout: cells separated by ``|``.
        cells = rest.rstrip().rstrip('|').split('|')[1:]
        head, tail = (cells + ['', ''])[:3], cells[3:]
    else:
        # Legacy layout: after the time stamp and a space, energy
        # ``{:>12s}``, eigst ``{:>6s}`` and dens ``{:>5s}``, each followed by
        # ``c `` or two spaces, then optional custom criteria and the magmom.
        # A value wider than its column pushes the later columns right.
        head, pos = [], 1
        for width in (12, 6, 5):
            end = pos + width
            while end < len(rest) and rest[end] not in ' c':
                end += 1
            head.append(rest[pos:end])
            pos = end + 2
        tail = rest[pos:].split()
    energy, eigst, dens = (_number(c) for c in head)
    magmoms = [c.strip() for c in tail if _MAGMOM.fullmatch(c.strip())]
    magmom = _magmom(magmoms[-1]) if magmoms else None
    return energy, eigst, dens, magmom


def extract_scf_convergence(logfile):
    """Extract SCF iteration data from a GPAW txt log file.

    Reads the ``iter:`` rows of the SCF table written by GPAW's legacy
    implementation (fixed-width columns, no header) and by its new one
    (``|``-separated cells), from the first row to the next blank line or
    ``Converged`` line. A value GPAW leaves blank (eigst and dens in the first
    iteration) is ``None``. ``'magmom'`` is a float for collinear
    spin-polarized output, a list of three floats for noncollinear output,
    and ``None`` for non-magnetic calculations.

    Every line shaped like an ``iter:`` row counts, so other ``iter:`` tables
    GPAW writes (the PZ-SIC localization loop) are not told apart.
    """
    iterations = []
    with open(logfile, 'r') as f:
        in_scf = False
        for line in f:
            if in_scf and (line.strip() == '' or 'Converged' in line):
                break
            match = _ROW.match(line.rstrip())
            if not match:
                continue
            in_scf = True
            energy, eigst, dens, magmom = _parse_row(match.group(2))
            iterations.append({
                'iter': int(match.group(1)),
                'energy': energy,
                'log10_eigst': eigst,
                'log10_dens': dens,
                'magmom': magmom,
            })
    return iterations
