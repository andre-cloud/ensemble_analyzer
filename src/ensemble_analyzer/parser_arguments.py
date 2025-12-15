import argparse
import os
import sys
import json

from importlib.resources import files


def print_help_protocol():
    example = json.dumps(
        {
            "0": {
                "functional": "str: DEFINE THE DFT FUNCTIONAL",
                "basis": "str: DEFINE THE BASIS SET FOR THE CALCULATION. DEFAULT: def2-svp",
                "opt": "bool: TRUE IF WANT TO OPTIMIZE. DEFAULT: False",
                "freq": "bool: TRUE IF WANT ANALYTICAL FREQUENCY CALCULATION. DEFAULT: False",
                "freq_fact": "float: FREQUENCY SCALE FACTOR",
                "graph": "bool : TRUE IF WANT SIMULATION OF ELECTRONIC GRAPH",
                "mult": "int: DEFINE THE MULTIPLICITY OF THE SYSTEM. DEFAULT: 1",
                "charge": "int: DEFINE THE CHARGE OF THE SYSTEM. DEFAULT: 0",
                "solv": {
                    "solvent": "str|null: NAME OF THE SOLVENT. IF GAS PHASE DEFINE AS NULL",
                    "smd": "bool: TRUE IF SMD MODEL EMPLOYED FOR IMPLICIT CALCULATION, ELSE CPCM USE",
                    "_comment": "SOLV KEYWORD CAN BE OMITTED, AND A GAS PHASE CALCULATION WILL BE CARRIED OUT.",
                },
                "add_input": "str: ADDITIONAL INPUT TO BE PASSED TO THE CALCULATOR. NO SANITY CHECK ON THIS BLOCK! USE WITH CARE. DEFAULT null",
                "thrG": "float: ADD SPECIAL THRESHOLD FOR THIS PECULIAR STEP OF THE PROTOL FOR thrG, IF NOT PRESENT, THE DEFAULT ONE WILL BE USED",
                "thrB": "float: ADD SPECIAL THRESHOLD FOR THIS PECULIAR STEP OF THE PROTOL FOR thrB, IF NOT PRESENT, THE DEFAULT ONE WILL BE USED",
                "thrGMAX": "float: ADD SPECIAL THRESHOLD FOR THIS PECULIAR STEP OF THE PROTOL FOR thrGMAX, IF NOT PRESENT, THE DEFAULT ONE WILL BE USED",
            }
        },
        indent=4,
    )

    txt = f'The protocol file must be a JSON file and it can contain as many calculation parts as desired. These calculations can be single-point, optimizations, solo frequency, or optimization-frequency calculations. The JSON file is formatted as follows. Note that only "func" and "solvent" (if the "solv" keyword is used) are mandatory:\n{example}'

    print(txt)
    sys.exit()


def print_help_threshold():
    with open(
        str(
                files("ensemble_analyzer").joinpath("parameters_file/default_threshold.json")
        ),
        "r",
    ) as j:
        contents = json.loads(j.read())
    example = json.dumps(contents, indent=4)

    txt = f"With this JSON file, different pruning parameters can be changed. USE WITH CARE!\nHere, different types of calculations have different thresholds:\n\t- thrG: Refers to the energy (G or E) threshold to consider two conformers equivalent together with thrB.\n\t- thrB: Refers to the rotary constant threshold to consider two conformers equivalent together with thrG.\n\t- thrGMAX: Refers to the maximum energy window considered. Conformers lying above it will be sorted out immediately.\n\nTwo conformers are considered equivalent if G_(CONFi)-G(CONFi-1) < thrG AND B_(CONFi)-B_(CONFi-1) < thrB.\n\nThe default parametes are:\n{example}\n\nIf you want to change only some of the parameters, ALL parameters must be passed.\nDO NOT CHANGE THE KEYWORDS"

    print(txt)
    sys.exit()


def parser_arguments():
    """
    Argparse argument parse
    """

    parser = argparse.ArgumentParser(
        prog="Ensemble Analyser",
        description="Quickly and easily identify the most important conformers and reduce the computational time",
        add_help=False,
    )

    input_group = parser.add_argument_group("Input Files")
    input_group.add_argument(
        "-e",
        "--ensemble",
        help="The ensemble file. Could be an xyz file (preferably) or other type parsable by OpenBabel",
        required=("--restart" not in sys.argv),
    )
    input_group.add_argument(
        "--restart", help="Restart the calculation", action="store_true"
    )
    input_group.add_argument(
        "-p",
        "--protocol",
        help="JSON file contains the computational protocol. Default: %(default)s",
        default=os.path.join(
            os.path.dirname(__file__), "parameters_file", "default_protocol.json"
        ),
    )
    input_group.add_argument(
        "-t",
        "--threshold",
        help="JSON file contains the threshold divided by calculation type. Default: %(default)s",
        default=os.path.join(
            os.path.dirname(__file__), "parameters_file", "default_threshold.json"
        ),
    )

    molecule_group = parser.add_argument_group("Molecule information")
    molecule_group.add_argument(
        "-T",
        "--temperature",
        help="Define the temperature in Kelvin. Default %(default)s",
        default=298.15,
        type=float,
    )
    molecule_group.add_argument(
        "-P",
        "--pressure",
        help="Define the pressure in kPa. Default %(default)s",
        default=101.325,
        type=float,
    )
    molecule_group.add_argument(
        
        "--linear",
        help="Define if molecule is linear.",
        action="store_true",
    )
    molecule_group.add_argument(
        "-alpha",
        "--alpha",
        help="Define alpha value for qRRHO dumping. Default %(default)s",
        default=4,
        type=float,
    )
    molecule_group.add_argument(
        "-cut-off",
        "--cut-off",
        help="Define cut-off frequency value for qRRHO dumping. Default %(default)s",
        default=100,
        type=float,
    )

    molecule_group.add_argument(
        "-no-H",
        "--exclude-H",
        help="Exclude hydrogen atoms for the PCA analysis and for the RSMD calculation. Default %(default)s: All atoms are considered for PCA and RMSD",
        action="store_false",
    )
    # molecule_group.add_argument('-no-enantio', '--exclude-enantiomers', help='Exclude the enantiomeric conformation. Default %(default)s: no check with others conformation\'s mirror image is performed', action='store_true')

    system_group = parser.add_argument_group("System Parameters")
    system_group.add_argument(
        "-cpu",
        type=int,
        help="Define the number of CPU used by the calculations",
        default=1,
    )
    system_group.add_argument(
        "-calc",
        "--calculator",
        help="Define the calculator to use. Default %(default)s",
        choices=["orca", "gaussian"],
        default="orca",
    )

    graph_group = parser.add_argument_group("Graph Parameters")

    graph_group.add_argument(
        "--definition",
        help="Define the Number of samples to generate the linear space of the X axis: 10^definition. USE WITH CAUTION. Default %(default)s",
        default=4,
    )

    graph_group.add_argument(
        "--fwhm-vibro",
        help="Define the Full Width @ Half Maximum (FWHM) for the lorentzian convolution. This can be a single value, a list, or None.",
        default=None,
        type=mix_type, 
        nargs='+'
    )
    graph_group.add_argument(
        "--fwhm-electro",
        help="Define the Full Width @ Half Maximum (FWHM) for the gaussian convolution. This can be a single value, a list, or None.",
        default=None,
        type=mix_type, 
        nargs='+'
    )

    graph_group.add_argument(
        "--shift-vibro",
        help="Define the multiplier of the vibronic graph after the convolution. This can be a single value, a list, or None.",
        default=None,
        type=mix_type, 
        nargs='+'
    )
    graph_group.add_argument(
        "--shift-electro",
        help="Define the shift of the electronic graph after the convolution. This can be a single value, a list, or None.",
        default=None,
        type=mix_type, 
        nargs='+'
    )

    graph_group.add_argument(
        "--interest-vibro",
        help="Define the interested area of the vibronic graph, so that these area can can be centered. This can be a single value, a list, or None.",
        default=None,
        type=mix_type, 
        nargs='+'
    )
    graph_group.add_argument(
        "--interest-electro",
        help="Define the interested area of the electronic graph, so that these area can can be centered. This can be a single value, a list, or None.",
        default=None,
        type=mix_type, 
        nargs='+'
    )

    graph_group.add_argument(
        "--invert",
        help="Invert the CHIRAL computed graph",
        action="store_true",
        default=False,
    )

    log_group = parser.add_argument_group("Parameters for the log")
    log_group.add_argument(
        "-o",
        "--output",
        help="Define the output filename. Default: %(default)s",
        default="output.out",
    )

    log_group.add_argument(
        '--disable-color',
        action='store_false',
        help='Disable colored output'
    )

    help_group = parser.add_argument_group("Get help")
    help_group.add_argument(
        "-h", "--help", action="help", help="Show this help message"
    )
    help_group.add_argument(
        "-h-p",
        "--help-protocol",
        help="Get help to format correctly the protocol JSON-file",
        action="store_true",
    )
    help_group.add_argument(
        "-h-t",
        "--help-threshold",
        help="Get help to format correctly the threshold JSON-file",
        action="store_true",
    )

    if __name__ == "__main__":  # pragma: no cover:
        parser.print_help()

    a = sys.argv
    if "-h-p" in a or "--help-protocol" in a:
        return print_help_protocol()

    if "-h-t" in a or "--help-threshold" in a:
        return print_help_threshold()

    return parser.parse_args()


def mix_type(values): 

    if values is None: 
        return None 
    if len(values) == 0: 
        return True

    if len(values) == 1:
        v = values[0].lower()
        if v in ["true", "t", "yes", "y", "1"]:
            return True
        if v in ["false", "f", "no", "n", "0"]:
            return False

        # Float singolo
        try:
            return float(values[0])
        except ValueError:
            raise ValueError(f"Cannot convert '{values[0]}' to float or boolean.")
    if len(values) == 2: 
        try: 
            return [float(v) for v in values]
        except ValueError: 
            raise ValueError(f'Cannot convert values {values=} to a float list.')
    if len(values) >2 : 
        raise ValueError(f'Too many arguments given')


if __name__ == "__main__":  # pragma: no cover:
    parser_arguments()
