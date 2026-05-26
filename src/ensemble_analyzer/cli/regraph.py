#!/opt/miniconda3/bin/python

from ensemble_analyzer._logger.create_log import create_logger
from ensemble_analyzer.graph import main_spectra, plot_comparative_graphs

from ensemble_analyzer.ensemble_io import load_workflow_data
from ensemble_analyzer._managers.calculation_config import CalculationConfig

from ensemble_analyzer.conformer.conformer import Conformer

from ensemble_analyzer.constants import compute_boltzmann_populations
from ensemble_analyzer._title import title
from typing import List

import argparse
import numpy as np
from datetime import datetime

def parse_argument() -> argparse.Namespace:
    """
    Parse command line arguments for the regrapher tool.

    Returns:
        argparse.Namespace: Parsed arguments including protocol indices and flags.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument('idx', nargs='+', help="Protocol's number to (re-)generate the graphs", type=int)
    parser.add_argument('-rb', '--read-boltz', help='Read Boltzmann population from a specific protocol', type=int)
    parser.add_argument('-no-nm', '--no-nm', help='Do not save the nm graphs', action='store_false')
    parser.add_argument('-w', '--weight', help='Show Weighting function', action='store_true')
    parser.add_argument('--disable-color', action='store_false', help='Disable colored output' )
    args = parser.parse_args()
    return args

def main() -> None:
    """
    Main execution flow for enan_regraph.
    
    Loads checkpoint and protocol data, recalculates Boltzmann populations 
    (unless overridden), and triggers spectrum generation and plotting.
    """

    args = parse_argument()
    fname_out = 'regraph.log'
    log = create_logger(fname_out, logger_name="enan_regraphy", debug=True, disable_color=False if not args.disable_color else True) # logger
    log.info(title)
    log._separator("Regraphing computed spectrum")

    config_mgr = CalculationConfig().load() # settings

    ensemble, protocol = load_workflow_data()

    start = datetime.now()

    if args.read_boltz: 
        assert str(args.read_boltz) in [p.number for p in protocol], f"{args.read_boltz} is not a specified step in the protocol file"
        for conf in ensemble:
            if not conf.active: continue
            for p in args.idx: 
                conf.energies.set(p, 'Pop', conf.energies[args.read_boltz].Pop)
    else:
        for protocol_number in args.idx:
            compute_boltzmann_populations(confs=ensemble, protocol_number=protocol_number, temperature=config_mgr.temperature)

    for i in args.idx:
        prot_obj = protocol[i]
        main_spectra(ensemble, prot_obj, log=log, invert=config_mgr.invert, shift=config_mgr.shift, fwhm=config_mgr.fwhm, interested_area=config_mgr.interested)


    plot_comparative_graphs(log, args.idx, show=False, nm=args.no_nm, show_ref_weight=args.weight)


    log.application_correct_end(total_conformers=len([c for c in ensemble if c.active]), total_time=datetime.now()-start)

if __name__ == '__main__':

    main()