import argparse
import sys
import logging
from pathlib import Path
from typing import Dict

from .core import MatplotlibPickleEditor, PickleSecurityError
from .tui import InteractiveTUI


logging.basicConfig(
    level=logging.WARNING,
    format='%(levelname)s: %(message)s'
)
logger = logging.getLogger(__name__)


def parse_mapping_file(filepath: Path) -> Dict[str, str]:
    mapping = {}
    with open(filepath, 'r', encoding='utf-8') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '=' not in line:
                logger.warning(f"Line {line_num} ignored: {line}")
                continue
            old, new = line.split('=', 1)
            mapping[old.strip()] = new.strip()
    return mapping


def _floats(pairs, what: str) -> Dict[str, float]:
    out = {}
    for label, value in pairs or []:
        try:
            out[label] = float(value)
        except ValueError:
            logger.error(f"Invalid {what} for '{label}': {value}")
    return out


def batch_mode(args: argparse.Namespace) -> int:
    try:
        editor = MatplotlibPickleEditor(args.pickle_file)
        editor.load()

        if args.list:
            labels = editor.get_legend_labels()
            if not labels:
                print("No labels found")
                return 0

            print(f"\nLabels in {args.pickle_file}:")
            print("-" * 50)
            colors = editor.get_line_colors()
            for idx, label in labels.items():
                print(f"  [{idx}] {label} (color: {colors.get(label, 'N/A')})")
            print()
            return 0

        rename_map = dict(args.rename or [])
        if args.rename_file:
            rename_map.update(parse_mapping_file(args.rename_file))

        visibility_map = {}
        for label, val in args.visibility or []:
            v = val.lower()
            if v in ('true', '1', 't', 'yes', 'on', 'false', '0', 'f', 'no', 'off'):
                visibility_map[label] = v in ('true', '1', 't', 'yes', 'on')
            else:
                logger.error(f"Invalid boolean value for '{label}': {val}")

        edits = [
            (rename_map, editor.rename_legend_labels, "labels renamed"),
            (dict(args.color or []), editor.change_line_colors, "colors changed"),
            (dict(args.linestyle or []), editor.change_line_linestyle, "linestyles changed"),
            (_floats(args.linewidth, "linewidth"), editor.change_line_linewidth, "line widths changed"),
            (_floats(args.alpha, "alpha"), editor.change_line_alpha, "alpha values changed"),
            (visibility_map, editor.change_line_visibility, "visibilities changed"),
        ]
        for mapping, apply, what in edits:
            if mapping:
                logger.info(f"{apply(mapping)} {what}")

        if args.xlim:
            editor.set_xlim(args.xlim[0], args.xlim[1])
            logger.info(f"Set x-limits: {args.xlim[0]} – {args.xlim[1]}")
        if args.ylim:
            editor.set_ylim(args.ylim[0], args.ylim[1])
            logger.info(f"Set y-limits: {args.ylim[0]} – {args.ylim[1]}")

        if any(m for m, _, _ in edits) or args.xlim or args.ylim:
            output_path = editor.save(args.output, args.format)
            print(f"✓ Saved: {output_path}")
        else:
            logger.warning("No changes specified")

        return 0

    except PickleSecurityError as e:
        logger.error(f"✗ Security error: {e}")
        return 1
    except FileNotFoundError as e:
        logger.error(f"✗ File not found: {e}")
        return 1
    except Exception as e:
        logger.error(f"✗ {e}")
        if args.verbose:
            logger.exception("Full traceback:")
        return 1


def main() -> int:
    parser = argparse.ArgumentParser(
        description='Interactive/batch editor for matplotlib pickle files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
DEFAULT MODE: Interactive TUI
  python -m pickle_editor.cli plot.pkl

BATCH MODE (examples):
  python -m pickle_editor.cli plot.pkl --batch --list
  python -m pickle_editor.cli plot.pkl --batch --rename "Protocol 1" "Proto A"
  python -m pickle_editor.cli plot.pkl --batch --color "Experimental" red --output new.pkl
        """
    )

    parser.add_argument('pickle_file', type=Path,
                       help='Matplotlib pickle file')

    parser.add_argument('--batch', '-b', action='store_true',
                       help='Batch mode (non-interactive)')

    batch_group = parser.add_argument_group('batch mode options')
    batch_group.add_argument('--list', '-l', action='store_true',
                            help='List labels and exit')
    batch_group.add_argument('--rename', '-r', nargs=2,
                            metavar=('OLD', 'NEW'), action='append',
                            help='Rename a label')
    batch_group.add_argument('--rename-file', '-rf', type=Path,
                            help='Mapping file OLD=NEW')
    batch_group.add_argument('--color', '-c', nargs=2,
                            metavar=('LABEL', 'COLOR'), action='append',
                            help='Change color')

    batch_group.add_argument('--linestyle', '-ls', nargs=2,
                        metavar=('LABEL', 'STYLE'), action='append',
                        help='Change line style (e.g. -, --, :, -.)')
    batch_group.add_argument('--linewidth', '-lw', nargs=2,
                            metavar=('LABEL', 'WIDTH'), action='append',
                            help='Change line width (float)')
    batch_group.add_argument('--alpha', '-a', nargs=2,
                            metavar=('LABEL', 'ALPHA'), action='append',
                            help='Change line transparency (0-1)')
    batch_group.add_argument('--visibility', '-vis', nargs=2,
                            metavar=('LABEL', 'bool'), action='append',
                            help='Change line visibility (true/false)')

    batch_group.add_argument('--xlim', nargs=2, type=float,
                             metavar=('XMIN', 'XMAX'),
                             help='Set x-axis limits (e.g. --xlim 150 500)')
    batch_group.add_argument('--ylim', nargs=2, type=float,
                             metavar=('YMIN', 'YMAX'),
                             help='Set y-axis limits (e.g. --ylim -1 1)')

    batch_group.add_argument('--output', '-o', type=Path,
                            help='Output file')
    batch_group.add_argument('--format', '-f', default='pickle',
                            choices=['pickle', 'png', 'pdf', 'svg'],
                            help='Output format')

    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Verbose output')

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    if args.batch:
        return batch_mode(args)

    try:
        editor = MatplotlibPickleEditor(args.pickle_file)
        editor.load()

        tui = InteractiveTUI(editor)
        tui.run()

        return 0

    except KeyboardInterrupt:
        print("\n\nInterrupted by user")
        return 130
    except PickleSecurityError as e:
        logger.error(f"✗ Security error: {e}")
        return 1
    except FileNotFoundError as e:
        logger.error(f"✗ File not found: {e}")
        return 1
    except Exception as e:
        logger.error(f"✗ Error: {e}")
        if args.verbose:
            logger.exception("Full traceback:")
        return 1


if __name__ == '__main__':
    sys.exit(main())
