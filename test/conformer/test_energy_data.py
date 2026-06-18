import pytest
import numpy as np
from ensemble_analyzer.conformer.energy_data import EnergyRecord, EnergyStore


class TestEnergyRecord:

    def test_default_values(self):
        rec = EnergyRecord()
        assert rec.E == 0.0
        assert np.isnan(rec.G)
        assert np.isnan(rec.Pop)

    def test_as_dict_roundtrip(self):
        rec = EnergyRecord(
            E=-100.0, G=-99.5, Pop=25.0,
            B_vec=np.array([1.0, 2.0, 3.0]),
            Freq=np.array([100.0, 200.0]),
            NormalModes=np.ones((2, 5, 3)),
        )
        d = rec.as_dict()
        assert d["E"] == -100.0
        assert d["Pop"] == 25.0
        assert d["B_vec"] == [1.0, 2.0, 3.0]
        assert d["Freq"] == [100.0, 200.0]
        assert d["NormalModes"] == []  # no negative freqs → no modes saved
        restored = EnergyRecord.from_dict(d)
        assert restored.E == -100.0
        assert np.allclose(restored.B_vec, [1.0, 2.0, 3.0])
        assert np.allclose(restored.Freq, [100.0, 200.0])
        assert restored.NormalModes.shape == (0,)

    def test_from_dict_with_none_arrays(self):
        d = {"E": -50.0, "G": np.nan, "H": np.nan, "S": np.nan, "G_E": np.nan,
             "zpve": np.nan, "B": None, "B_vec": None, "m": None, "m_vec": None,
             "Pop": np.nan, "time": None, "Erel": np.nan, "Freq": None}
        rec = EnergyRecord.from_dict(d)
        assert rec.E == -50.0
        assert rec.B_vec is None


class TestEnergyStore:

    @pytest.fixture
    def store(self):
        s = EnergyStore()
        s.add(1, EnergyRecord(E=-100.0, G=-99.0, Pop=50.0, B=1.5, time=5.0,
                              Freq=np.array([100.0, 200.0]),
                              NormalModes=np.ones((2, 5, 3))))
        s.add(2, EnergyRecord(E=-200.0, G=-198.5, Pop=50.0, B=1.5, time=5.0))
        return s

    def test_add_and_retrieve(self, store):
        assert 1 in store
        rec = store[1]
        assert rec.E == -100.0

    def test_getitem_missing(self, store):
        rec = store[99]
        assert rec.E == 0.0

    def test_contains(self, store):
        assert 1 in store
        assert 2 in store
        assert 99 not in store

    def test_last(self, store):
        rec = store.last()
        assert rec.E == -200.0

    def test_last_empty(self):
        s = EnergyStore()
        rec = s.last()
        assert rec.E == 0.0

    def test_get_energy_uses_g(self, store):
        assert store.get_energy() == -198.5

    def test_get_energy_fallback(self):
        s = EnergyStore()
        s.add(1, EnergyRecord(E=-100.0))
        assert s.get_energy() == -100.0

    def test_set_valid_attribute(self, store):
        store.set(1, "Pop", 75.0)
        assert store[1].Pop == 75.0

    def test_set_missing_protocol(self, store):
        with pytest.raises(KeyError):
            store.set(99, "Pop", 50.0)

    def test_set_invalid_attribute(self, store):
        with pytest.raises(AttributeError):
            store.set(1, "nonexistent", 0.0)

    def test_log_info(self, store):
        data = store.log_info(1)
        assert len(data) == 7

    def test_as_dict_roundtrip(self, store):
        d = store.as_dict()
        assert "1" in d
        assert "2" in d
        restored = EnergyStore()
        restored.load({"data": d})
        assert 1 in restored
        assert restored[1].E == -100.0

    def test_get_last_freq_found(self, store):
        freq = store.get_last_freq(1)
        assert np.allclose(freq, [100.0, 200.0])

    def test_get_last_freq_fallback(self, store):
        freq = store.get_last_freq(2)
        assert np.allclose(freq, [100.0, 200.0])

    def test_get_last_freq_empty(self):
        s = EnergyStore()
        assert len(s.get_last_freq(1)) == 0
