#!/opt/miniconda3/bin/python

from src._logger.create_log import create_logger
from src.graph import main_spectra, plot_comparative_graphs

from src._managers.checkpoint_manager import CheckpointManager
from src._managers.protocol_manager import ProtocolManager
from src._managers.calculation_config import CalculationConfig

from src._conformer.conformer import Conformer

from src.constants import *

from typing import List

import argparse

parser = argparse.ArgumentParser()

parser.add_argument('idx', nargs='+', help="Protocol's number to (re-)generate the graphs", type=int)
parser.add_argument('-rb', '--read-boltz', help='Read Boltzmann population from a specific protocol', type=int)
parser.add_argument('-no-nm', '--no-nm', help='Do not save the nm graphs', action='store_false')
parser.add_argument('-w', '--weight', help='Show Weighting function', action='store_true')
args = parser.parse_args()


fname_out = 'regraph.log'
log = create_logger(fname_out, logger_name="enan_regraphy", debug=True, ) # logger

checkpoint_mgr = CheckpointManager()
protocol_mgr = ProtocolManager()
config_mgr = CalculationConfig().load() # settings

ensemble = checkpoint_mgr.load() # ensemble
protocol = protocol_mgr.load() # protocol

protocol_number = args.idx[-1]
log.debug(f'E({protocol_number}): {[conf.energies.__getitem__(protocol_number=protocol_number).E for conf in ensemble if conf.active]}')
log.debug(f'G({protocol_number}): {[conf.energies.__getitem__(protocol_number=protocol_number).G for conf in ensemble if conf.active]}')

def calc_boltzmann(confs: List[Conformer], temperature: float, protocol_number:int) -> None: 

    active: List[Conformer] = [conf for conf in confs if conf.active]
    energy = np.array([conf.get_energy(protocol_number=protocol_number) for conf in active])
    log.debug(energy)
    rel_en = energy - np.min(energy)

    log.debug(rel_en)

    exponent = np.exp(-rel_en * CAL_TO_J * 1000 * EH_TO_KCAL /(R*temperature))
    populations = exponent/exponent.sum()
    for idx, conf in enumerate(active):
        conf.energies.__getitem__(protocol_number=protocol_number).Pop = populations[idx] * 100
        conf.energies.__getitem__(protocol_number=protocol_number).Erel = rel_en[idx] * EH_TO_KCAL

    return None

if args.read_boltz: 
    assert str(args.read_boltz) in [p.number for p in protocol], f"{args.read_boltz} is not a specified step in the protocol file"
    for conf in ensemble:
        if not conf.active: continue
        for p in args.idx: 
            conf.energies.__getitem__(p).Pop = conf.energies.__getitem__(args.read_boltz).Pop
else:
    for protocol_number in args.idx:
        calc_boltzmann(ensemble, protocol_number=protocol_number, temperature=config_mgr.temperature)
        log.debug(f'{[conf.energies.__getitem__(protocol_number=protocol_number).Erel for conf in ensemble if conf.active]}')

for i in args.idx:
    prot_obj = protocol[i]
    log.spectra_start(protocol_number=i)
    main_spectra(ensemble, prot_obj, log=log, invert=config_mgr.invert, shift=config_mgr.shift, fwhm=config_mgr.fwhm, interested_area=config_mgr.interested)


plot_comparative_graphs(log, args.idx, show=False, nm=args.no_nm, show_ref_weight=args.weight)

