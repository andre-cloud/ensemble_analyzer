#!/opt/miniconda3/bin/python

from src.launch import restart
from src.logger import create_log
from src.protocol import Protocol
from src.graph import main_spectra, plot_comparative_graphs
from src.pruning import calculate_rel_energies

import json, logging
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('idx', nargs='+', help="Protocol's number to (re-)generate the graphs", type=int)
parser.add_argument('-rb', '--read-boltz', help='Read Boltzmann population from a specific protocol', type=int)
parser.add_argument('-no-nm', '--no-nm', help='Do not save the nm graphs', action='store_false')
args = parser.parse_args()


def load_protocol(p: dict, log: logging): 
    protocol = [Protocol.load_raw(p[d]) for d in p]

    return protocol


ensemble, protocol, _ = restart()

settings = json.load(open("settings.json"))
output = 'regraph.log'
temperature = settings.get("temperature", 298.15)

invert = settings.get("invert", False)

FWHM={'vibro':settings.get("fwhm_vibro", None), "electro":settings.get("fwhm_electro", None)}
SHIFT={'vibro':settings.get("shift_vibro", None), "electro":settings.get("shift_electro", None)}

INTERESTED_AREA = {'vibro': settings.get('interested_vibro', None), "electro": settings.get('interested_electro', None)}

calculate_rel_energies(ensemble, 298.15)

# initiate the log
log = create_log(output)
# deactivate the log of matplotlib
logging.getLogger("matplotlib").disabled = False

protocol = load_protocol(json.load(open('protocol_dump.json')), log)


if args.read_boltz: 

    assert str(args.read_boltz) in [p.number for p in protocol], f"{args.read_boltz} is not a specified step in the protocol file"
    for conf in ensemble:
        if not conf.active: continue
        for p in protocol: 
            conf.energies[str(p.number)]["Pop"] = conf.energies[str(args.read_boltz)]["Pop"]


for i in args.idx:
    prot_obj = protocol[i]
    log.info(f'{"-"*30}\nProtocol {prot_obj.number}\n{"-"*30}')
    main_spectra(ensemble, prot_obj, log=log, invert=invert, shift=SHIFT, fwhm=FWHM, interested_area=INTERESTED_AREA)


plot_comparative_graphs(log, args.idx, show=False, nm=args.no_nm)

