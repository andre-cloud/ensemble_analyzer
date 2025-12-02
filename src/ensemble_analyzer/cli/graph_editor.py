#!/usr/bin/env python3
"""
Interactive TUI for editing legends and colours in matplotlib figure pickles.

Default mode: Interactive TUI
Batch mode available via --batch
"""

import argparse
import pickle
import sys
import logging
from pathlib import Path
from typing import Dict, Optional, List
import warnings

try:
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    from matplotlib.axes import Axes
except ImportError:
    print("ERROR: matplotlib not installed. Run: pip install matplotlib",
          file=sys.stderr)
    sys.exit(1)

# InquirerPy (required for TUI)
try:
    from InquirerPy import inquirer
    from InquirerPy.base.control import Choice
    from InquirerPy.separator import Separator
    INQUIRER_AVAILABLE = True
except ImportError:
    INQUIRER_AVAILABLE = False

#
# Rich (required for coloured output)
try:
    from rich.console import Console
    from rich.table import Table
    from rich.panel import Panel
    from rich.prompt import Confirm
    RICH_AVAILABLE = True
    console = Console()
except ImportError:
    RICH_AVAILABLE = False
    console = None


# Setup logging
logging.basicConfig(
    level=logging.WARNING,  # Meno verboso in TUI mode
    format='%(levelname)s: %(message)s'
)
logger = logging.getLogger(__name__)


class PickleSecurityError(Exception):
    """Exception for security issues when loading pickle."""
    pass


class MatplotlibPickleEditor:
    """Modify colors and labels figure matplotlib pickled."""
    
    COMMON_COLORS = [
        'red', 'blue', 'green', 'black', 'orange', 'purple', 'brown',
        'pink', 'gray', 'cyan', 'magenta', 'yellow',
        '#2E86AB', '#A23B72', '#F18F01', '#C73E1D', '#6A994E',
        '#BC4B51', '#5B8E7D', '#8B5A3C', '#264653', '#E76F51'
    ]
    
    def __init__(self, pickle_path: Path, strict_validation: bool = True):
        self.pickle_path = pickle_path
        self.strict_validation = strict_validation
        self.figure: Optional[Figure] = None
        self.axes: Optional[Axes] = None
        self._modifications_made = False
        
        if not self.pickle_path.exists():
            raise FileNotFoundError(f"File not found: {self.pickle_path}")
    
    def load(self) -> None:
        """Load and validate the pickle."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                with open(self.pickle_path, 'rb') as f:
                    obj = pickle.load(f)
            except pickle.UnpicklingError as e:
                raise PickleSecurityError(
                    f"Corrupted or invalid pickle file: {e}"
                ) from e

        if not isinstance(obj, Figure):
            if self.strict_validation:
                raise PickleSecurityError(
                    f"Object is not matplotlib.figure.Figure, but {type(obj)}"
                )
            else:
                logger.warning(f"WARNING: Unexpected type {type(obj)}")

        self.figure = obj

        if self.figure.axes:
            self.axes = self.figure.axes[0]
        else:
            raise PickleSecurityError("No axes found in figure")
    
    def get_legend_labels(self) -> Dict[int, str]:
        """Get the current legend labels."""
        if not self.axes:
            raise RuntimeError("You must call load() first")

        legend = self.axes.get_legend()
        if not legend:
            return {}

        labels = {}
        for idx, text in enumerate(legend.get_texts()):
            labels[idx] = text.get_text()

        return labels
    
    def get_line_colors(self) -> Dict[str, str]:
        """Get the current line colours."""
        if not self.axes:
            raise RuntimeError("You must call load() first")

        legend = self.axes.get_legend()
        if not legend:
            return {}

        lines = self.axes.get_lines()
        colors = {}

        for line, text in zip(lines, legend.get_texts()):
            label = text.get_text()
            color = matplotlib.colors.to_hex(line.get_color())
            colors[label] = color

        return colors
    
    def rename_legend_labels(self, mapping: Dict[str, str]) -> int:
        """Rename legend labels."""
        if not self.axes:
            raise RuntimeError("You must call load() first")

        legend = self.axes.get_legend()
        if not legend:
            return 0

        renamed = 0
        for text in legend.get_texts():
            current_label = text.get_text()
            if current_label in mapping:
                text.set_text(mapping[current_label])
                renamed += 1
                self._modifications_made = True

        return renamed
    
    def change_line_colors(self, label_color_map: Dict[str, str]) -> int:
        """Change line colours."""
        if not self.axes:
            raise RuntimeError("You must call load() first")

        lines = self.axes.get_lines()
        legend = self.axes.get_legend()

        if not legend:
            return 0

        label_to_lines = {}
        for line, text in zip(lines, legend.get_texts()):
            label = text.get_text()
            label_to_lines[label] = line

        changed = 0
        for label, color in label_color_map.items():
            if label in label_to_lines:
                try:
                    label_to_lines[label].set_color(color)
                    changed += 1
                    self._modifications_made = True
                    # Also update legend handles so the legend shows the new colour
                    
                    legend = self.axes.get_legend()
                    if legend:
                        legend_lines = legend.get_lines()
                        legend_texts = legend.get_texts()
                        for leg_line, leg_text in zip(legend_lines, legend_texts):
                            if leg_text.get_text() == label:
                                try:
                                    leg_line.set_color(color)
                                except Exception:
                                    pass
                except ValueError:
                    pass

        return changed
    
    def save(self, output_path: Optional[Path] = None,
             format: str = 'pickle') -> Path:
        """Save the modified figure."""
        if not self.figure:
            raise RuntimeError("You must call load() first")

        if output_path is None:
            if format == 'pickle':
                output_path = self.pickle_path
            else:
                output_path = self.pickle_path.with_suffix(f'.{format}')

        if format == 'pickle':
            with open(output_path, 'wb') as f:
                pickle.dump(self.figure, f, protocol=pickle.HIGHEST_PROTOCOL)
        else:
            self.figure.savefig(output_path, format=format, dpi=300,
                               bbox_inches='tight')

        self._modifications_made = False
        return output_path
    
    def preview(self) -> None:
        """Show the figure for preview."""
        if not self.figure:
            raise RuntimeError("You must call load() first")
        plt.show()
    
    def has_modifications(self) -> bool:
        """Check if there are unsaved modifications."""
        return self._modifications_made


class InteractiveTUI:
    """Interactive Terminal UI using InquirerPy."""

    def __init__(self, editor: MatplotlibPickleEditor):
        if not INQUIRER_AVAILABLE:
            raise RuntimeError(
                "InquirerPy not installed. Run: pip install InquirerPy"
            )

        self.editor = editor
        self.console = console if RICH_AVAILABLE else None

    def print_panel(self, message: str, title: str = "Info", style: str = "cyan") -> None:
        """Print a formatted panel."""
        if self.console:
            self.console.print(Panel(message, title=title, border_style=style))
        else:
            print(f"\n{title}: {message}\n")

    def show_current_state(self) -> None:
        """Show current figure state."""
        labels = self.editor.get_legend_labels()
        colors = self.editor.get_line_colors()

        if not labels:
            self.print_panel("No legend found", "Warning", "yellow")
            return

        if self.console:
            table = Table(title="Current Legend State", show_header=True,
                         header_style="bold cyan")
            table.add_column("Index", style="dim", width=8)
            table.add_column("Label", style="bold")
            table.add_column("Colour", style="magenta")

            for idx, label in labels.items():
                color = colors.get(label, "N/A")
                table.add_row(str(idx), label, color)

            self.console.print(table)
        else:
            print("\n=== Current State ===")
            for idx, label in labels.items():
                color = colors.get(label, "N/A")
                print(f"  [{idx}] {label} (colour: {color})")
            print()

    def rename_labels_flow(self) -> None:
        """Interactive flow to rename labels."""
        labels = self.editor.get_legend_labels()

        if not labels:
            self.print_panel("No labels to rename", "Error", "red")
            return

        # Select label to rename
        choices = [
            Choice(value=label, name=f"{label} [{idx}]")
            for idx, label in labels.items()
        ]
        choices.append(Separator())
        choices.append(Choice(value=None, name="← Back"))

        selected = inquirer.select(
            message="Select label to rename:",
            choices=choices,
            default=None
        ).execute()

        if selected is None:
            return

        # Input new name
        new_name = inquirer.text(
            message=f"New name for '{selected}':",
            default=selected
        ).execute()

        if new_name and new_name != selected:
            self.editor.rename_legend_labels({selected: new_name})
            self.print_panel(f"'{selected}' → '{new_name}'", "Success", "green")

    def change_colors_flow(self) -> None:
        """Interactive flow to change colours."""
        labels = self.editor.get_legend_labels()
        colors = self.editor.get_line_colors()

        if not labels:
            self.print_panel("No labels found", "Error", "red")
            return

        # Select label
        choices = [
            Choice(value=label,
                   name=f"{label} (current: {colors.get(label, 'N/A')})")
            for label in labels.values()
        ]
        choices.append(Separator())
        choices.append(Choice(value=None, name="← Back"))

        selected = inquirer.select(
            message="Select label to change colour:",
            choices=choices,
            default=None
        ).execute()

        if selected is None:
            return

        # Choose colour input method
        color_method = inquirer.select(
            message="How would you like to specify the colour?",
            choices=[
                Choice(value="preset", name="Choose from predefined palette"),
                Choice(value="custom", name="Enter manually (name or hex)"),
                Choice(value=None, name="← Cancel")
            ],
            default="preset"
        ).execute()

        if color_method is None:
            return

        if color_method == "preset":
            # Predefined palette with preview
            color_choices = [
                Choice(value=c, name=f"{c}")
                for c in self.editor.COMMON_COLORS
            ]
            color_choices.append(Separator())
            color_choices.append(Choice(value=None, name="← Cancel"))

            new_color = inquirer.select(
                message="Select colour:",
                choices=color_choices,
                default=None
            ).execute()
        else:
            # Manual input
            new_color = inquirer.text(
                message="Colour (name or hex #RRGGBB):",
                validate=lambda x: len(x) > 0
            ).execute()

        if new_color:
            self.editor.change_line_colors({selected: new_color})
            self.print_panel(
                f"Colour of '{selected}' changed to {new_color}",
                "Success", "green"
            )

    def save_flow(self) -> None:
        """Interactive flow to save."""
        # Output format
        format_choice = inquirer.select(
            message="Save format:",
            choices=[
                Choice(value="pickle", name="Pickle (editable later)"),
                Choice(value="png", name="PNG (image)"),
                Choice(value="pdf", name="PDF (vector)"),
                Choice(value="svg", name="SVG (vector)"),
                Choice(value=None, name="← Cancel")
            ],
            default="pickle"
        ).execute()

        if format_choice is None:
            return

        # Output path
        default_name = self.editor.pickle_path.stem
        if format_choice != "pickle":
            default_name = f"{default_name}_modified"

        custom_path = inquirer.confirm(
            message="Would you like to specify a custom path?",
            default=False
        ).execute()

        output_path = None
        if custom_path:
            path_str = inquirer.text(
                message="Output path:",
                default=f"{default_name}.{format_choice}"
            ).execute()
            output_path = Path(path_str)

        try:
            saved_path = self.editor.save(output_path, format_choice)
            self.print_panel(f"Saved: {saved_path}", "Success", "green")
        except Exception as e:
            self.print_panel(f"Save error: {e}", "Error", "red")

    def run(self) -> None:
        """Main loop TUI."""
        if self.console:
            self.console.clear()
            self.console.print(
                Panel.fit(
                    "[bold cyan]Interactive Matplotlib Pickle Editor[/]\n"
                    f"File: {self.editor.pickle_path}",
                    border_style="cyan"
                )
            )

        while True:
            # Show state
            self.show_current_state()

            # Main menu
            action = inquirer.select(
                message="What would you like to do?",
                choices=[
                    Choice(value="rename", name="📝 Rename legend label"),
                    Choice(value="color", name="🎨 Change line colours"),
                    Separator(),
                    Choice(value="preview", name="👁️  Preview figure"),
                    Choice(value="save", name="💾 Save modifications"),
                    Separator(),
                    Choice(value="reload", name="🔄 Reload original file"),
                    Choice(value="exit", name="🚪 Exit"),
                ],
                default="rename"
            ).execute()

            if action == "rename":
                self.rename_labels_flow()
            elif action == "color":
                self.change_colors_flow()
            elif action == "preview":
                self.print_panel("Close the matplotlib window to continue", "Info")
                self.editor.preview()
            elif action == "save":
                self.save_flow()
            elif action == "reload":
                if self.editor.has_modifications():
                    confirm = inquirer.confirm(
                        message="You have unsaved modifications. Reload anyway?",
                        default=False
                    ).execute()
                    if not confirm:
                        continue

                self.editor.load()
                self.print_panel("File reloaded", "Success", "green")
            elif action == "exit":
                if self.editor.has_modifications():
                    confirm = inquirer.confirm(
                        message="You have unsaved modifications. Exit anyway?",
                        default=False
                    ).execute()
                    if not confirm:
                        continue

                if self.console:
                    self.console.print("\n[cyan]Goodbye! 👋[/]\n")
                else:
                    print("\nGoodbye!\n")
                break


def parse_mapping_file(filepath: Path) -> Dict[str, str]:
    """Parsa file di mapping OLD=NEW."""
    mapping = {}
    with open(filepath, 'r', encoding='utf-8') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '=' not in line:
                logger.warning(f"Linea {line_num} ignorata: {line}")
                continue
            old, new = line.split('=', 1)
            mapping[old.strip()] = new.strip()
    return mapping


def batch_mode(args):
    """Run in batch mode (non-interactive)."""
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
                print(f"  [{idx}] {label} (colour: {color})")
            print()
            return 0

        # Renames
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

        # Colours
        color_map = {}
        if args.color:
            for label, color in args.color:
                color_map[label] = color

        if color_map:
            count = editor.change_line_colors(color_map)
            logger.info(f"Changed {count} colours")

        # Preview
        if args.preview:
            editor.preview()

        # Save
        if rename_map or color_map:
            output_path = editor.save(args.output, args.format)
            print(f"✓ Saved: {output_path}")
        else:
            logger.warning("No modifications specified")

        return 0

    except Exception as e:
        logger.error(f"✗ {e}")
        return 1


def main():
    """Entry point CLI."""
    parser = argparse.ArgumentParser(
        description='Interactive/batch editor for matplotlib pickles',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
DEFAULT MODE: Interactive TUI
  python %(prog)s plot.pkl

BATCH MODE (examples):
  python %(prog)s plot.pkl --batch --list
  python %(prog)s plot.pkl --batch --rename "Protocol 1" "Proto A"
  python %(prog)s plot.pkl --batch --color "Experimental" red --output new.pkl
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
                            help='Rename label')
    batch_group.add_argument('--rename-file', '-rf', type=Path,
                            help='Mapping file OLD=NEW')
    batch_group.add_argument('--color', '-c', nargs=2,
                            metavar=('LABEL', 'COLOR'), action='append',
                            help='Change colour')
    batch_group.add_argument('--output', '-o', type=Path,
                            help='Output file')
    batch_group.add_argument('--format', '-f', default='pickle',
                            choices=['pickle', 'png', 'pdf', 'svg'],
                            help='Output format')
    batch_group.add_argument('--preview', '-p', action='store_true',
                            help='Preview before saving')

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
        print("ERROR: Interactive mode requires InquirerPy", file=sys.stderr)
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
    except Exception as e:
        logger.error(f"✗ Error: {e}", exc_info=args.verbose)
        return 1


if __name__ == '__main__':
    sys.exit(main())