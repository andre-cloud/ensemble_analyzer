#!/opt/miniconda3/bin/python

from src.launch import restart
from src.logger import create_log
from src.protocol import Protocol
from src.graph import main_graph, Compared
from src.pruning import calculate_rel_energies
from src.parser_parameter import get_data_for_graph

import json, logging
import argparse, os


parser = argparse.ArgumentParser()

parser.add_argument('idx', nargs='+', help="Protocol's number to (re-)generate the graphs")
parser.add_argument('-rb', '--read-boltz', help='Read Boltzmann population from a specific protocol', type=int)
args = parser.parse_args()


def load_protocol(p: dict, log: logging): 
    protocol = [Protocol.load_raw(p[d]) for d in p]

    return protocol


ensemble, protocol, _ = restart()

settings = json.load(open("settings.json"))
output = 'regraph.log'
temperature = settings.get("temperature", 298.15)

invert = settings.get("invert", False)

calculate_rel_energies(ensemble, 298.15)

# initiate the log
log = create_log(output)
# deactivate the log of matplotlib
logging.getLogger("matplotlib").disabled = False

protocol = load_protocol(json.load(open('protocol_dump.json')), log)


if args.read_boltz: 
    assert args.read_boltz in [p.number for p in protocol], f"{args.read_boltz} is not a specified step in the protocol file"
    for conf in ensemble:
        if not conf.active: continue
        for p in protocol: 
            conf.energies[str(p.number)]["Pop"] = conf.energies[str(args.read_boltz)]["Pop"]

for i in protocol:
    print(f"Retriving datas for Protocol {i.number}")
    get_data_for_graph(conformers=ensemble, protocol=i, log=log)
    print(f"Creating graphs for Protocol {i.number}")
    main_graph(ensemble, i, log=log, invert=invert)

for j in ["IR", "VCD", "UV", "ECD"]:
    g = Compared(protocol, graph_type=j, log=log)
    g.save_graph()
