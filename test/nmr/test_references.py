import numpy as np
import pytest

from ensemble_analyzer._nmr.references import (
    bundled_reference,
    needs_rdkit,
    resolve_reference_geometry,
)


@pytest.mark.parametrize("name", ["TMS", "tms", "Tetramethylsilane", "SiMe4"])
def test_bundled_lookup_is_case_insensitive_with_aliases(name):
    assert bundled_reference(name)["canonical"] == "TMS"


def test_bundled_unknown():
    assert bundled_reference("my_ref") is None


def test_bundled_geometry_used_without_smiles():
    g = resolve_reference_geometry({"name": "tms"})
    assert g.source == "bundled:TMS"
    assert g.symbols.count("Si") == 1 and g.symbols.count("H") == 12
    assert g.coords.shape == (17, 3) and g.charge == 0 and g.mult == 1
    assert not needs_rdkit({"name": "tms"})


def test_every_bundled_geometry_loads():
    from ensemble_analyzer._nmr.references import _INDEX
    for name in _INDEX:
        g = resolve_reference_geometry({"name": name})
        assert len(g.symbols) == len(g.coords) > 0


def test_user_xyz_wins(tmp_path):
    f = tmp_path / "ref.xyz"
    f.write_text("2\nc\nH 0 0 0\nF 0 0 0.92\n")
    g = resolve_reference_geometry({"name": "TMS", "xyz": str(f), "charge": 0, "mult": 1})
    assert g.symbols == ["H", "F"] and g.source == f"xyz:{f}"
    assert np.allclose(g.coords[1], [0, 0, 0.92])


def test_relative_xyz_resolved_from_base_dir(tmp_path):
    (tmp_path / "ref.xyz").write_text("1\n\nNe 0 0 0\n")
    g = resolve_reference_geometry({"name": "x", "xyz": "ref.xyz"}, base_dir=tmp_path)
    assert g.symbols == ["Ne"]


def test_smiles_needs_rdkit_only_when_not_bundled():
    assert needs_rdkit({"name": "PhF", "smiles": "Fc1ccccc1"})
    assert not needs_rdkit({"name": "TMS", "smiles": "C[Si](C)(C)C"})


def test_smiles_geometry():
    pytest.importorskip("rdkit")
    g = resolve_reference_geometry({"name": "PhF", "smiles": "Fc1ccccc1"})
    assert g.symbols.count("F") == 1 and len(g.symbols) == 12 and g.source == "smiles:Fc1ccccc1"
    g = resolve_reference_geometry({"name": "ion", "smiles": "C[N+](C)(C)C"})
    assert g.charge == 1


def test_nothing_to_build_from():
    with pytest.raises(ValueError, match="xyz.*smiles"):
        resolve_reference_geometry({"name": "unknown_ref"})
