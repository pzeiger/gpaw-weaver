import pytest

from gpaw_weaver import make_fd_params, make_pw_params, get_mode_filestr


class TestMakePwParams:
    def test_basic_fields(self):
        p = make_pw_params(400, {"size": (4, 4, 4), "gamma": True})
        assert p["mode"] == {"name": "pw", "ecut": 400}
        assert p["kpts"] == {"size": (4, 4, 4), "gamma": True}
        assert p["xc"] == "PBE"

    def test_xc_override(self):
        assert make_pw_params(400, {}, xc="LDA")["xc"] == "LDA"

    def test_optional_keys_absent_when_not_given(self):
        p = make_pw_params(400, {})
        assert "convergence" not in p
        assert "mixer" not in p
        assert "setups" not in p

    def test_optional_keys_present_when_given(self):
        conv = {"density": 1e-6}
        mixer = {"beta": 0.1}
        setups = {"Fe": ":d,4.0"}
        p = make_pw_params(400, {}, convergence=conv, mixer=mixer, setups=setups)
        assert p["convergence"] == conv
        assert p["mixer"] == mixer
        assert p["setups"] == setups

    def test_extra_kwargs_forwarded(self):
        p = make_pw_params(400, {}, nbands=10, charge=-1)
        assert p["nbands"] == 10
        assert p["charge"] == -1


class TestMakeFdParams:
    def test_with_h(self):
        p = make_fd_params(h=0.2, kpts={"size": (4, 4, 4)})
        assert p["mode"] == "fd"
        assert p["h"] == pytest.approx(0.2)
        assert "gpts" not in p

    def test_with_gpts(self):
        p = make_fd_params(gpts=(36, 36, 1), kpts={})
        assert p["gpts"] == (36, 36, 1)
        assert "h" not in p

    def test_gpts_coerced_to_int(self):
        p = make_fd_params(gpts=(36.0, 36.0, 1.0), kpts={})
        assert all(isinstance(g, int) for g in p["gpts"])

    def test_gpts_takes_precedence_over_h_with_warning(self):
        with pytest.warns(UserWarning, match="gpts takes precedence"):
            p = make_fd_params(h=0.2, gpts=(36, 36, 1), kpts={})
        assert "gpts" in p
        assert "h" not in p

    def test_neither_h_nor_gpts_raises(self):
        with pytest.raises(ValueError, match="requires either h or gpts"):
            make_fd_params(kpts={})

    def test_xc_default(self):
        assert make_fd_params(h=0.2, kpts={})["xc"] == "PBE"

    def test_optional_keys_absent_when_not_given(self):
        p = make_fd_params(h=0.2, kpts={})
        assert "convergence" not in p
        assert "mixer" not in p
        assert "setups" not in p


class TestGetModeFilestr:
    def test_pw_with_gamma(self):
        p = make_pw_params(500, {"size": (3, 3, 1), "gamma": True})
        assert get_mode_filestr(p) == "pw_ecut500_3x3x1_gamma"

    def test_pw_without_gamma(self):
        p = make_pw_params(400, {"size": (4, 4, 4)})
        assert get_mode_filestr(p) == "pw_ecut400_4x4x4"

    def test_fd_with_h(self):
        p = make_fd_params(h=0.15, kpts={"size": (3, 3, 1), "gamma": True})
        assert get_mode_filestr(p) == "fd_h0.15_3x3x1_gamma"

    def test_fd_with_gpts(self):
        p = make_fd_params(gpts=(36, 36, 1), kpts={"size": (3, 3, 1), "gamma": True})
        assert get_mode_filestr(p) == "fd_gpts36x36x1_3x3x1_gamma"
