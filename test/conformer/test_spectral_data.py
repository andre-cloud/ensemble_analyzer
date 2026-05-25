import pytest
import numpy as np
from ensemble_analyzer.conformer.spectral_data import SpectralRecord, SpectralStore


class TestSpectralRecord:

    def test_creation(self):
        rec = SpectralRecord(X=np.array([1000.0, 2000.0]), Y=np.array([0.5, 1.0]))
        assert len(rec) == 2
        assert not rec.is_empty

    def test_converts_to_array(self):
        rec = SpectralRecord(X=[1000.0, 2000.0], Y=[0.5, 1.0])
        assert isinstance(rec.X, np.ndarray)
        assert isinstance(rec.Y, np.ndarray)

    def test_shape_mismatch(self):
        with pytest.raises(ValueError, match="same shape"):
            SpectralRecord(X=np.array([1.0, 2.0]), Y=np.array([0.5]))

    def test_not_1d(self):
        with pytest.raises(ValueError):
            SpectralRecord(X=np.array([[1.0]]), Y=np.array([[0.5]]))

    def test_is_empty(self):
        rec = SpectralRecord(X=np.array([]), Y=np.array([]))
        assert rec.is_empty

    def test_as_dict_roundtrip(self):
        rec = SpectralRecord(X=np.array([1000.0]), Y=np.array([0.5]))
        d = rec.as_dict()
        assert d["X"] == [1000.0]
        restored = SpectralRecord.from_dict(d)
        assert np.allclose(restored.X, [1000.0])
        assert np.allclose(restored.Y, [0.5])


class TestSpectralStore:

    @pytest.fixture
    def store(self):
        s = SpectralStore()
        rec = SpectralRecord(X=np.array([1000.0]), Y=np.array([0.5]))
        s.add(1, "IR", rec)
        return s

    def test_add_and_getitem(self, store):
        rec = store[1, "IR"]
        assert np.allclose(rec.X, [1000.0])

    def test_contains(self, store):
        assert 1 in store
        assert 2 not in store

    def test_has_graph_type(self, store):
        assert store.has_graph_type(1, "IR")
        assert not store.has_graph_type(1, "VCD")

    def test_as_dict_roundtrip(self, store):
        d = store.as_dict()
        assert 1 in d
        assert "IR" in d[1]
        restored = SpectralStore()
        restored.load({"data": d})
        assert "IR" in restored.data[1]
        rec = restored[1, "IR"]
        assert np.allclose(rec.X, [1000.0])
