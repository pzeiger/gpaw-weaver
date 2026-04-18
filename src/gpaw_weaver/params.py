import warnings


def make_pw_params(ecut, kpts, xc='PBE', convergence=None,
                   mixer=None, setups=None, **kwargs):
    """Return a GPAW calc_params dict for plane-wave (PW) mode.

    Parameters
    ----------
    ecut : float
        Plane-wave energy cutoff in eV.
    kpts : dict
        k-point specification, e.g. ``{'size': (4, 4, 4), 'gamma': True}``.
    xc : str
        Exchange-correlation functional (default 'PBE').
    convergence : dict, optional
        Convergence criteria, e.g. ``{'density': 1e-6}``.
    mixer : dict, optional
        Density mixer settings.
    setups : dict, optional
        PAW setup overrides, e.g. ``{'Mn': ':d,4.0'}``.
    **kwargs
        Any additional GPAW keyword arguments.
    """
    params = {
        'mode': {'name': 'pw', 'ecut': ecut},
        'kpts': kpts,
        'xc': xc,
    }
    if convergence is not None:
        params['convergence'] = convergence
    if mixer is not None:
        params['mixer'] = mixer
    if setups is not None:
        params['setups'] = setups
    params.update(kwargs)
    return params


def make_fd_params(h=None, kpts=None, xc='PBE', convergence=None,
                   mixer=None, setups=None, gpts=None, **kwargs):
    """Return a GPAW calc_params dict for finite-difference (FD) mode.

    Parameters
    ----------
    h : float, optional
        Target real-space grid spacing in Ångström. Used when *gpts* is not given.
    gpts : tuple of int, optional
        Explicit grid-point counts ``(na, nb, nc)``. Takes precedence over *h*.
        Compute a commensurate grid via ``gpaw.utilities.h2gpts(h, cell, idiv=8)``.
        A warning is issued when both *h* and *gpts* are supplied.
    kpts : dict
        k-point specification, e.g. ``{'size': (4, 4, 4), 'gamma': True}``.
    xc : str
        Exchange-correlation functional (default 'PBE').
    convergence : dict, optional
        Convergence criteria, e.g. ``{'density': 1e-6}``.
    mixer : dict, optional
        Density mixer settings.
    setups : dict, optional
        PAW setup overrides, e.g. ``{'Mn': ':d,4.0'}``.
    **kwargs
        Any additional GPAW keyword arguments.
    """
    if gpts is None and h is None:
        raise ValueError('make_fd_params requires either h or gpts.')

    params = {'mode': 'fd', 'kpts': kpts, 'xc': xc}

    if gpts is not None:
        if h is not None:
            warnings.warn(
                'Both h and gpts were supplied to make_fd_params; '
                'gpts takes precedence — h will be ignored.',
                UserWarning, stacklevel=2,
            )
        params['gpts'] = tuple(int(g) for g in gpts)
    else:
        params['h'] = h

    if convergence is not None:
        params['convergence'] = convergence
    if mixer is not None:
        params['mixer'] = mixer
    if setups is not None:
        params['setups'] = setups
    params.update(kwargs)
    return params


def get_mode_filestr(calc_params):
    """Return a human-readable string encoding the calc mode and k-grid.

    Examples
    --------
    PW 500 eV, 3×3×1 Γ-centred        →  ``'pw_ecut500_3x3x1_gamma'``
    FD h=0.15 Å, 3×3×1 Γ-centred      →  ``'fd_h0.15_3x3x1_gamma'``
    FD gpts=(36,36,1), 3×3×1 Γ-centred →  ``'fd_gpts36x36x1_3x3x1_gamma'``
    """
    mode = calc_params.get('mode', 'fd')
    kpts = calc_params.get('kpts', {})

    if isinstance(kpts, dict):
        size = kpts.get('size', [1, 1, 1])
        gamma = kpts.get('gamma', False)
        kstr = 'x'.join(str(k) for k in size)
        if gamma:
            kstr += '_gamma'
    else:
        kstr = str(kpts)

    if isinstance(mode, dict) and mode.get('name') == 'pw':
        ecut = mode.get('ecut', 400)
        return f'pw_ecut{ecut:.0f}_{kstr}'
    else:
        if 'gpts' in calc_params:
            gpts = calc_params['gpts']
            gstr = 'x'.join(str(g) for g in gpts)
            return f'fd_gpts{gstr}_{kstr}'
        else:
            h = calc_params.get('h', 0.2)
            return f'fd_h{h:.2f}_{kstr}'
