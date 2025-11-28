#!/opt/miniconda3/bin/python

from src._logger.create_log import create_logger
from src.graph import main_spectra, plot_comparative_graphs

from src._managers.checkpoint_manager import CheckpointManager
from src._managers.protocol_manager import ProtocolManager
from src._managers.calculation_config import CalculationConfig

import argparse

parser = argparse.ArgumentParser()

parser.add_argument('idx', nargs='+', help="Protocol's number to (re-)generate the graphs", type=int)
parser.add_argument('-rb', '--read-boltz', help='Read Boltzmann population from a specific protocol', type=int)
parser.add_argument('-no-nm', '--no-nm', help='Do not save the nm graphs', action='store_false')
parser.add_argument('-w', '--weight', help='Show Weighting function', action='store_true')
args = parser.parse_args()


fname_out = 'regraph.log'
log = create_logger(fname_out, logger_name="enan_regraphy") # logger

checkpoint_mgr = CheckpointManager()
protocol_mgr = ProtocolManager()
config_mgr = CalculationConfig().load() # settings

ensemble = checkpoint_mgr.load() # ensemble
protocol = protocol_mgr.load() # protocol

if args.read_boltz: 
    assert str(args.read_boltz) in [p.number for p in protocol], f"{args.read_boltz} is not a specified step in the protocol file"
    for conf in ensemble:
        if not conf.active: continue
        for p in args.idx: 
            conf.energies.__getitem__(p).Pop = conf.energies.__getitem__(args.read_boltz).Pop


for i in args.idx:
    prot_obj = protocol[i]
    log.spectra_start(protocol_number=i)
    main_spectra(ensemble, prot_obj, log=log, invert=config_mgr.invert, shift=config_mgr.shift, fwhm=config_mgr.fwhm, interested_area=config_mgr.interested)


plot_comparative_graphs(log, args.idx, show=False, nm=args.no_nm, show_ref_weight=args.weight)

