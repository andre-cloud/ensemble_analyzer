"""ORCA %pal block: automatic nprocs_group."""

from unittest.mock import patch

import pytest

from enan_calculators._orca import (
    create_orca_calc,
    orca_major_version,
    orca_nprocs_group,
    orca_pal_block,
)


@pytest.mark.parametrize(
    "cpu, group",
    [(0, 1), (1, 1), (2, 1), (3, 1), (4, 2), (6, 3), (7, 1), (8, 4), (9, 3), (12, 6), (16, 8), (49, 7)],
)
def test_nprocs_group_is_largest_proper_divisor(cpu, group):
    assert orca_nprocs_group(cpu) == group


@pytest.mark.parametrize(
    "cpu, expected",
    [
        (1, "%pal nprocs 1 end\n"),
        (2, "%pal nprocs 2 end\n"),
        (7, "%pal nprocs 7 end\n"),
        (4, "%pal nprocs 4 nprocs_group 2 end\n"),
        (8, "%pal nprocs 8 nprocs_group 4 end\n"),
        (12, "%pal nprocs 12 nprocs_group 6 end\n"),
    ],
)
def test_pal_block(cpu, expected):
    assert orca_pal_block(cpu) == expected


@pytest.mark.parametrize("add_input", ["%pal nprocs_group 4 end", "%PAL NPROCS_GROUP 2 END"])
def test_user_nprocs_group_wins(add_input):
    assert orca_pal_block(8, add_input) == "%pal nprocs 8 end\n"


@pytest.mark.parametrize("version, has_group", [(None, True), (4, False), (5, True), (6, True)])
def test_version_gate(version, has_group):
    assert ("nprocs_group" in orca_pal_block(8, "", version)) is has_group


@pytest.mark.parametrize("env, major", [(None, None), ("", None), ("6.1.0", 6), ("4.2.1", 4), ("10.0", 10)])
def test_orca_major_version(monkeypatch, env, major):
    if env is None:
        monkeypatch.delenv("ORCAVERSION", raising=False)
    else:
        monkeypatch.setenv("ORCAVERSION", env)
    assert orca_major_version() == major


def test_create_orca_calc_writes_group(monkeypatch, tmp_path):
    monkeypatch.setenv("ORCAVERSION", "6.1.0")
    with patch("enan_calculators._orca.shutil.which", return_value="orca"):
        calc = create_orca_calc(0, 1, "B3LYP", "def2-SVP", cpu=8, directory=str(tmp_path))
    assert calc.parameters["orcablocks"].startswith("%pal nprocs 8 nprocs_group 4 end\n")


def test_orca_calc_common_str_respects_version(setup_calc_orca):
    from ensemble_analyzer.calculators.orca import OrcaCalc

    conf, proto = setup_calc_orca
    with patch("ensemble_analyzer.calculators.orca.orca_profile"), \
         patch("ensemble_analyzer.calculators.orca.OrcaCalc.VERSION", 6):
        _, ob, _ = OrcaCalc(proto, 8, conf).common_str()
    assert "%pal nprocs 8 nprocs_group 4 end" in ob

    with patch("ensemble_analyzer.calculators.orca.orca_profile"), \
         patch("ensemble_analyzer.calculators.orca.OrcaCalc.VERSION", 4):
        _, ob, _ = OrcaCalc(proto, 8, conf).common_str()
    assert "%pal nprocs 8 end" in ob


@pytest.fixture
def setup_calc_orca(mock_conformer, mock_protocol):
    mock_protocol.add_input = ""
    mock_protocol.solvent = None
    mock_protocol.dispersion = None
    return mock_conformer, mock_protocol
