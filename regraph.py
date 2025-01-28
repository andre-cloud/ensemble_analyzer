from src.launch import restart
from src.logger import create_log
from src.protocol import Protocol
from src.grapher import Graph, plot_conv_graph
import json, logging
import argparse, os


parser = argparse.ArgumentParser()

parser.add_argument('idx', nargs='+', help="Protocol's number to (re-)generate the graphs")

args = parser.parse_args()


def load_protocol(p: dict, log: logging): 
    protocol = [Protocol.load_raw(p[d]) for d in p]

    return protocol


conformers, protocol, start_from = restart()

settings = json.load(open("settings.json"))
output = 'regraph.log'
temperature = settings.get("temperature", 298.15)
final_lambda = settings.get("final_lambda", 800)
definition = settings.get("definition", 4)
fwhm = settings.get("fwhm", None)
shift = settings.get("shift", None)
invert = settings.get("invert", False)

# initiate the log
log = create_log(output)
# deactivate the log of matplotlib
logging.getLogger("matplotlib").disabled = False

protocol = load_protocol(json.load(open('protocol_dump.json')), log)



for p in args.idx: 
    log.info(f'\n\nProtocol number: {p}')
    Graph(
        confs=conformers,
        protocol=protocol[int(p)],
        log=log,
        T=temperature,
        final_lambda=final_lambda,
        definition=definition,
        FWHM=fwhm,
        shift=shift,
        invert=invert,
        regraph = True
    )


plot_conv_graph(args.idx, protocol)
