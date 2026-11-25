import json

import numpy as np

from ensemble_analyzer.conformer.conformer import Conformer
from ensemble_analyzer.conformer.nmr_data import NMRRecord, NMRStore
from ensemble_analyzer.io_utils import write_json


def _record():
    return NMRRecord(shielding={0: 150.2, 3: 30.1}, aniso={0: 20.0, 3: 5.0},
                     coupling={(3, 0): 125.0, (0, 1): 35.5})


def test_record_orders_pairs_and_reads_both_ways():
    r = _record()
    assert set(r.coupling) == {(0, 3), (0, 1)}
    assert r.j(3, 0) == r.j(0, 3) == 125.0
    assert r.j(1, 3) is None
    assert not NMRRecord().has_data() and r.has_data()


def test_store_round_trip_through_json(tmp_path):
    store = NMRStore()
    store.add(2, _record())
    assert 2 in store and store.get(2).shielding[3] == 30.1 and store.get(5) is None
    with open(tmp_path / "x.json", "w") as fp:
        write_json({"nmr": store.as_dict()}, fp)
    loaded = NMRStore()
    loaded.load(json.loads((tmp_path / "x.json").read_text())["nmr"])
    r = loaded.get(2)
    assert r.shielding == {0: 150.2, 3: 30.1} and r.coupling == {(0, 1): 35.5, (0, 3): 125.0}


def test_conformer_checkpoint_round_trip_and_old_checkpoints():
    c = Conformer(number=1, geom=np.zeros((2, 3)), atoms=["H", "H"], raw=True)
    c.nmr_data.add(1, _record())
    data = {"number": 1, "last_geometry": c.last_geometry, "atoms": ["H", "H"],
            "energies": c.energies.as_dict(), "graphs_data": c.graphs_data.as_dict(),
            "active": True, "nmr_data": c.nmr_data.as_dict()}
    assert Conformer.load_raw(data).nmr_data.get(1).shielding[0] == 150.2
    data.pop("nmr_data")
    assert 1 not in Conformer.load_raw(data).nmr_data
