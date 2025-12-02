#!/opt/miniconda3/bin/python

from ensemble_analyzer._logger.create_log import create_logger
from ensemble_analyzer.graph import main_spectra, plot_comparative_graphs

from ensemble_analyzer._managers.checkpoint_manager import CheckpointManager
from ensemble_analyzer._managers.protocol_manager import ProtocolManager
from ensemble_analyzer._managers.calculation_config import CalculationConfig

from ensemble_analyzer._conformer.conformer import Conformer

from ensemble_analyzer.constants import *
from ensemble_analyzer.title import title
from typing import List

import argparse
from datetime import datetime

def parse_argument():
    parser = argparse.ArgumentParser()

    parser.add_argument('idx', nargs='+', help="Protocol's number to (re-)generate the graphs", type=int)
    parser.add_argument('-rb', '--read-boltz', help='Read Boltzmann population from a specific protocol', type=int)
    parser.add_argument('-no-nm', '--no-nm', help='Do not save the nm graphs', action='store_false')
    parser.add_argument('-w', '--weight', help='Show Weighting function', action='store_true')
    parser.add_argument('--disable-color', action='store_false', help='Disable colored output' )
    args = parser.parse_args()
    return args

def main():
    args = parse_argument()
    fname_out = 'regraph.log'
    log = create_logger(fname_out, logger_name="enan_regraphy", debug=True, disable_color=False if not args.disable_color else True) # logger
    log.info(title)
    log._separator("Regraphing computed spectrum")

    checkpoint_mgr = CheckpointManager()
    protocol_mgr = ProtocolManager()
    config_mgr = CalculationConfig().load() # settings

    ensemble = checkpoint_mgr.load() # ensemble
    protocol = protocol_mgr.load() # protocol

    start = datetime.now()

    if args.read_boltz: 
        assert str(args.read_boltz) in [p.number for p in protocol], f"{args.read_boltz} is not a specified step in the protocol file"
        for conf in ensemble:
            if not conf.active: continue
            for p in args.idx: 
                conf.energies.set(p,'Pop',conf.energies.__getitem__(args.read_boltz).Pop)
    else:
        for protocol_number in args.idx:
            calc_boltzmann(confs=ensemble, protocol_number=protocol_number, temperature=config_mgr.temperature)

    for i in args.idx:
        prot_obj = protocol[i]
        main_spectra(ensemble, prot_obj, log=log, invert=config_mgr.invert, shift=config_mgr.shift, fwhm=config_mgr.fwhm, interested_area=config_mgr.interested)


    plot_comparative_graphs(log, args.idx, show=False, nm=args.no_nm, show_ref_weight=args.weight)


    log.application_correct_end(total_conformers=len([c for c in ensemble if c.active]), total_time=datetime.now()-start)

def calc_boltzmann(confs: List[Conformer], temperature: float, protocol_number:int) -> None: 

    active: List[Conformer] = [conf for conf in confs if conf.active]
    energy = np.array([conf.get_energy(protocol_number=protocol_number) for conf in active])
    rel_en = energy - np.min(energy)

    exponent = np.exp(-rel_en * CAL_TO_J * 1000 * EH_TO_KCAL /(R*temperature))
    populations = exponent/exponent.sum()
    for idx, conf in enumerate(active):
        conf.energies.set(protocol_number=protocol_number, property="Pop", value=populations[idx] * 100)
        conf.energies.set(protocol_number=protocol_number, property='Erel', value=rel_en[idx] * EH_TO_KCAL)

    return None


if __name__ == '__main__':

    main()