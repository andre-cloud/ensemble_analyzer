import json
import numpy as np

from ase.calculators import orca

from src.protocol import Protocol, Solvent
from src.IOsystem import SerialiseEncoder


def test_protocol():
    inp = {
        "number": "2",
        "functional": "R2SCAN-3C",
        "basis": "DEF2-MTZVPP",
        "solvent": None,
        "opt": True,
        "freq": False,
        "add_input": "",
        "thrG": 0.5,
        "thrB": 5e-05,
        "thrGMAX": 4,
        "calculator": "orca",
        "freq_fact": 1,
        "graph": False,
        "constrains": [],
        "maxstep": 0.1,
        "fmax": 0.01,
        "cluster": False,
    }

    p = Protocol.load_raw(inp)
    print(
        json.dumps(
            p,
            cls=SerialiseEncoder,
        )
    )
    print(json.dumps(inp))
    assert json.dumps(
        p,
        cls=SerialiseEncoder,
    ) == json.dumps(inp)
    assert p.calculation_level == "OPT"
    assert p.level == "R2SCAN-3C/DEF2-MTZVPP"

    assert type(p.get_calculator(1)) is tuple


def test_solvent():
    inp = {"solvent": "acetonitrile", "smd": True}
    assert json.dumps(
        Solvent(inp),
        cls=SerialiseEncoder,
    ) == json.dumps(inp)

    s = Solvent(inp)
    assert s.orca_input_smd() == '%cpcm smd true smdsolvent "acetonitrile" end'
    s.smd = False
    assert s.orca_input_smd() == ""
