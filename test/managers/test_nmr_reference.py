"""NMR reference runner: optimisation level, caching and reuse (calculations mocked)."""

import json
from unittest.mock import MagicMock

import numpy as np
import pytest

from ensemble_analyzer._managers.nmr_reference import NMRReferenceRunner
from ensemble_analyzer.conformer.nmr_data import NMRRecord
from ensemble_analyzer.protocol.nmr import check_protocols
from ensemble_analyzer.protocol.protocol import Protocol

TMS = {"name": "TMS", "nuclei": ["1H", "13C"]}


class FakeExecutor:
    """Records each calculation; optimisations shift the geometry, NMR steps return shieldings."""

    def __init__(self):
        self.calls = []

    def execute(self, idx, conf, protocol, cpu=None):
        self.calls.append((protocol.number, protocol.opt, protocol.functional, conf.folder,
                           protocol.charge, protocol.mult, conf.last_geometry.copy()))
        if protocol.opt:
            conf.last_geometry = conf.last_geometry + 0.01
        else:
            sh = {i: (180.0 + i if s == "C" else 31.0 + i * 0.01) for i, s in enumerate(conf.atoms)
                  if s in ("C", "H")}
            conf.nmr_data.add(protocol.number, NMRRecord(shielding=sh))
        return True


def _p(number, **kw):
    base = dict(functional="PBE0", basis="def2-SVP", calculator="orca")
    base.update(kw)
    return Protocol(number=number, **base)


def _nmr(**kw):
    base = {"shielding": True, "nuclei": ["1H", "13C"], "references": [TMS]}
    base.update(kw)
    return base


@pytest.fixture
def workdir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    return tmp_path


def _runner(fake):
    r = NMRReferenceRunner(MagicMock(cpu=4), MagicMock())
    r.executor = fake
    return r


def test_optimised_at_last_opt_level_then_shielding(workdir):
    steps = [_p(0, functional="r2SCAN-3c", basis="", opt=True), _p(1, functional="PBE0", nmr=_nmr())]
    check_protocols(steps, MagicMock())
    fake = FakeExecutor()
    sigma = _runner(fake).ensure(steps[1], steps)

    (n_opt, opt_flag, func, folder, *_), (n_nmr, nmr_opt, func2, folder2, *_) = fake.calls
    assert (n_opt, opt_flag, func) == (0, True, "r2SCAN-3c") and folder.endswith("TMS/opt")
    assert (n_nmr, nmr_opt, func2) == (1, False, "PBE0") and folder2.endswith("TMS/nmr")
    # TMS: atoms 0,2,3,4 are C, 5..16 are H
    assert sigma["13C"]["sigma"] == pytest.approx(np.mean([180.0, 182.0, 183.0, 184.0]))
    assert sigma["1H"]["sigma"] == pytest.approx(np.mean([31.0 + i * 0.01 for i in range(5, 17)]))
    assert sigma["1H"]["reference"] == "TMS" and sigma["1H"]["opt_step"] == 0
    cache = json.loads((workdir / "nmr_references" / "TMS" / "reference.json").read_text())
    assert set(cache["geometries"]) == {"0"} and set(cache["shielding"]) == {"1"}


def test_second_nmr_step_reuses_geometry(workdir):
    steps = [_p(0, opt=True), _p(1, nmr=_nmr()), _p(2, functional="B3LYP", nmr=_nmr())]
    check_protocols(steps, MagicMock())
    fake = FakeExecutor()
    runner = _runner(fake)
    runner.ensure(steps[1], steps)
    runner.ensure(steps[2], steps)
    assert [(c[0], c[1]) for c in fake.calls] == [(0, True), (1, False), (2, False)]
    assert np.allclose(fake.calls[2][6], fake.calls[1][6])          # same optimised geometry


def test_new_opt_level_reoptimises_from_cached_geometry(workdir):
    steps = [_p(0, opt=True), _p(1, nmr=_nmr()), _p(2, functional="B3LYP", opt=True), _p(3, nmr=_nmr())]
    check_protocols(steps, MagicMock())
    fake = FakeExecutor()
    runner = _runner(fake)
    runner.ensure(steps[1], steps)
    runner.ensure(steps[3], steps)
    assert [(c[0], c[1]) for c in fake.calls] == [(0, True), (1, False), (2, True), (3, False)]
    assert np.allclose(fake.calls[2][6], fake.calls[0][6] + 0.01)  # starts from the step-0 geometry


def test_restart_uses_cache(workdir):
    steps = [_p(0, opt=True), _p(1, nmr=_nmr())]
    check_protocols(steps, MagicMock())
    _runner(FakeExecutor()).ensure(steps[1], steps)
    fake = FakeExecutor()
    sigma = _runner(fake).ensure(steps[1], steps)
    assert fake.calls == [] and "1H" in sigma


def test_changed_level_invalidates_cache(workdir):
    steps = [_p(0, opt=True), _p(1, nmr=_nmr())]
    check_protocols(steps, MagicMock())
    _runner(FakeExecutor()).ensure(steps[1], steps)
    steps[1].basis = "def2-TZVP"
    fake = FakeExecutor()
    _runner(fake).ensure(steps[1], steps)
    assert [(c[0], c[1]) for c in fake.calls] == [(1, False)]


def test_no_opt_step_uses_nmr_level(workdir):
    steps = [_p(0), _p(1, nmr=_nmr())]
    check_protocols(steps, MagicMock())
    fake = FakeExecutor()
    _runner(fake).ensure(steps[1], steps)
    assert [(c[0], c[1]) for c in fake.calls] == [(1, True), (1, False)]


def test_reference_atoms_charge_and_variant_cleanup(workdir):
    (workdir / "ref.xyz").write_text("3\n\nC 0 0 0\nH 0 0 1.1\nH 0 1.1 0\n")
    ref = {"name": "mine", "xyz": "ref.xyz", "nuclei": ["1H"], "atoms": [2], "charge": 1, "mult": 2}
    steps = [_p(0, opt=True, constraints=[[0, 1]], add_input="TightSCF"),
             _p(1, nmr=_nmr(nuclei=["1H"], references=[ref]))]
    check_protocols(steps, MagicMock())
    fake = FakeExecutor()
    sigma = _runner(fake).ensure(steps[1], steps)
    assert sigma["1H"]["sigma"] == pytest.approx(31.02) and sigma["1H"]["n_atoms"] == 1
    assert all(c[4:6] == (1, 2) for c in fake.calls)


def test_failed_reference_raises(workdir):
    steps = [_p(0, opt=True), _p(1, nmr=_nmr())]
    check_protocols(steps, MagicMock())
    fake = MagicMock()
    fake.execute.return_value = False
    with pytest.raises(RuntimeError, match="TMS"):
        _runner(fake).ensure(steps[1], steps)


def test_coupling_only_step_needs_nothing(workdir):
    steps = [_p(0, nmr={"coupling": True, "nuclei": ["1H"]})]
    fake = FakeExecutor()
    assert _runner(fake).ensure(steps[0], steps) == {} and fake.calls == []


def test_cached_sigmas_without_running(workdir):
    from ensemble_analyzer._managers.nmr_reference import cached_reference_sigmas
    steps = [_p(0, opt=True), _p(1, nmr=_nmr()), _p(2, nmr=_nmr())]
    check_protocols(steps, MagicMock())
    sigma = _runner(FakeExecutor()).ensure(steps[1], steps)
    assert cached_reference_sigmas(steps[1], steps, MagicMock(cpu=1), MagicMock()) == sigma
    with pytest.raises(RuntimeError, match="never computed"):
        cached_reference_sigmas(steps[2], steps, MagicMock(cpu=1), MagicMock())
