#!/usr/bin/env python3
"""
CLI entry point for Matplotlib Pickle Editor.

Supports interactive (TUI) and batch modes.
"""

import argparse
import sys
import logging
from pathlib import Path
from typing import Dict

from .core import MatplotlibPickleEditor, PickleSecurityError
from .tui import InteractiveTUI, INQUIRER_AVAILABLE


# Setup logging
logging.basicConfig(
    level=logging.WARNING,
    format='%(levelname)s: %(message)s'
)
logger = logging.getLogger(__name__)


def parse_mapping_file(filepath: Path) -> Dict[str, str]:
    """
    Parse mapping file OLD=NEW.
    
    Args:
        filepath: Path to the mapping file
    
    Returns:
        Dictionary {old: new}
    """
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


def batch_mode(args):
    """
    Run in batch (non-interactive) mode.
    
    Args:
        args: Arguments parsed by argparse
    
    Returns:
        Exit code (0 = success, 1 = error)
    """
    try:
        editor = MatplotlibPickleEditor(args.pickle_file,
                                       strict_validation=not args.no_strict)
        editor.load()

        # Query mode: list and exit
        if args.list:
            labels = editor.get_legend_labels()
            if not labels:
                print("No labels found")
                return 0

            print(f"\nLabels in {args.pickle_file}:")
            print("-" * 50)
            for idx, label in labels.items():
                colors = editor.get_line_colors()
                color = colors.get(label, "N/A")
                print(f"  [{idx}] {label} (color: {color})")
            print()
            return 0

        # Rename
        rename_map = {}
        if args.rename:
            for old, new in args.rename:
                rename_map[old] = new

        if args.rename_file:
            file_map = parse_mapping_file(args.rename_file)
            rename_map.update(file_map)

        if rename_map:
            count = editor.rename_legend_labels(rename_map)
            logger.info(f"Renamed {count} labels")

        # Colors
        color_map = {}
        if args.color:
            for label, color in args.color:
                color_map[label] = color

        if color_map:
            count = editor.change_line_colors(color_map)
            logger.info(f"Changed {count} colors")

        # Linestyle
        linestyle_map = {}
        if args.linestyle:
            for label, style in args.linestyle:
                linestyle_map[label] = style

        if linestyle_map:
            count = editor.change_line_linestyle(linestyle_map)
            logger.info(f"Changed {count} linestyles")

        # Linewidth
        linewidth_map = {}
        if args.linewidth:
            for label, width in args.linewidth:
                try:
                    linewidth_map[label] = float(width)
                except ValueError:
                    logger.error(f"Invalid linewidth for '{label}': {width}")

        if linewidth_map:
            count = editor.change_line_linewidth(linewidth_map)
            logger.info(f"Changed {count} line widths")

        # Alpha channel
        alpha_map = {}
        if args.alpha:
            for label, alpha in args.alpha:
                try:
                    alpha_map[label] = float(alpha)
                except ValueError:
                    logger.error(f"Invalid alpha for '{label}': {alpha}")

        if alpha_map:
            count = editor.change_line_alpha(alpha_map)
            logger.info(f"Changed {count} alpha values")

        # # Preview
        # if args.preview:
        #     editor.preview()

        # Save
        if rename_map or color_map or linestyle_map or linewidth_map or alpha_map:
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


def main():
    """CLI entry point."""
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

    # Batch mode options
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
    
    batch_group.add_argument('--output', '-o', type=Path,
                            help='Output file')
    batch_group.add_argument('--format', '-f', default='pickle',
                            choices=['pickle', 'png', 'pdf', 'svg'],
                            help='Output format')
    # batch_group.add_argument('--preview', '-p', action='store_true',
    #                         help='Preview before saving')

    parser.add_argument('--no-strict', action='store_true',
                       help='Disable strict validation')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Verbose output')

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Batch mode
    if args.batch:
        return batch_mode(args)

    # Interactive TUI mode (default)
    if not INQUIRER_AVAILABLE:
        print("ERROR: Interactive mode requires InquirerPy", 
              file=sys.stderr)
        print("Install: pip install InquirerPy rich", file=sys.stderr)
        print("\nUse --batch for non-interactive mode", file=sys.stderr)
        return 1

    try:
        editor = MatplotlibPickleEditor(args.pickle_file,
                                       strict_validation=not args.no_strict)
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