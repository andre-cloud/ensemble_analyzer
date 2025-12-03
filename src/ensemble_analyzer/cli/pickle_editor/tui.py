"""
Interactive Terminal User Interface for MatplotlibPickleEditor.

Interactive interface based on InquirerPy and Rich.
"""

import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .core import MatplotlibPickleEditor

# InquirerPy (required for TUI)
try:
    from InquirerPy import inquirer
    from InquirerPy.base.control import Choice
    from InquirerPy.separator import Separator
    INQUIRER_AVAILABLE = True
except ImportError:
    INQUIRER_AVAILABLE = False

# Rich (required for colored output)
try:
    from rich.console import Console
    from rich.table import Table
    from rich.panel import Panel
    RICH_AVAILABLE = True
    console = Console()
except ImportError:
    RICH_AVAILABLE = False
    console = None


class InteractiveTUI:
    """
    Interactive terminal user interface using InquirerPy.
    
    Provides an interactive menu for modifying matplotlib figures in a guided way.
    """

    def __init__(self, editor: 'MatplotlibPickleEditor'):
        """
        Initialize the TUI.
        
        Args:
            editor: Instance of MatplotlibPickleEditor
        
        Raises:
            RuntimeError: If InquirerPy is not installed
        """
        if not INQUIRER_AVAILABLE:
            raise RuntimeError(
                "InquirerPy not installed. Run: pip install InquirerPy"
            )

        self.editor = editor
        self.console = console if RICH_AVAILABLE else None

    def print_panel(self, message: str, title: str = "Info", 
                    style: str = "cyan") -> None:
        """
        Print a formatted panel.
        
        Args:
            message: Message to display
            title: Panel title
            style: Border color style
        """
        if self.console:
            self.console.print(Panel(message, title=title, border_style=style))
        else:
            print(f"\n{title}: {message}\n")

    def show_current_state(self) -> None:
        """Show the current state of the figure."""
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
            table.add_column("Color", style="magenta")

            for idx, label in labels.items():
                color = colors.get(label, "N/A")
                table.add_row(str(idx), label, color)

            self.console.print(table)
        else:
            print("\n=== Current State ===")
            for idx, label in labels.items():
                color = colors.get(label, "N/A")
                print(f"  [{idx}] {label} (color: {color})")
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

        # Enter new name
        new_name = inquirer.text(
            message=f"New name for '{selected}':",
            default=selected
        ).execute()

        if new_name and new_name != selected:
            self.editor.rename_legend_labels({selected: new_name})
            self.print_panel(f"'{selected}' → '{new_name}'", "Success", "green")

    def change_colors_flow(self) -> None:
        """Interactive flow to change colors."""
        labels = self.editor.get_legend_labels()
        colors = self.editor.get_line_colors()

        if not labels:
            self.print_panel("No label found", "Error", "red")
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
            message="Select label to change color:",
            choices=choices,
            default=None
        ).execute()

        if selected is None:
            return

        # Choose color input method
        color_method = inquirer.select(
            message="How do you want to specify the color?",
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
                message="Select color:",
                choices=color_choices,
                default=None
            ).execute()
        else:
            # Manual input
            new_color = inquirer.text(
                message="Color (name or hex #RRGGBB):",
                validate=lambda x: len(x) > 0
            ).execute()

        if new_color:
            changed = self.editor.change_line_colors({selected: new_color})
            if changed > 0:
                self.print_panel(
                    f"Color of '{selected}' changed to {new_color}",
                    "Success", "green"
                )
            else:
                self.print_panel(
                    f"Unable to change color (invalid color?)",
                    "Error", "red"
                )

    def change_linestyle_flow(self) -> None:
        """Interactive flow to change line style."""
        labels = self.editor.get_legend_labels()
        if not labels:
            self.print_panel("No label found", "Error", "red")
            return

        choices = [Choice(value=label, name=label) for label in labels.values()]
        choices.append(Separator())
        choices.append(Choice(value=None, name="← Back"))

        selected = inquirer.select(
            message="Select label to change line style:",
            choices=choices,
            default=None
        ).execute()

        if selected is None:
            return

        # Common styles
        style_choices = [
            Choice(value='-', name="Solid (-)"),
            Choice(value='--', name="Dashed (--)"),
            Choice(value=':', name="Dotted (:)"),
            Choice(value='-.', name="Dash-dot (-.)"),
            Choice(value='custom', name="Enter manually"),
            Choice(value=None, name="← Cancel")
        ]

        style_choice = inquirer.select(
            message="Select style:",
            choices=style_choices,
            default='-'
        ).execute()

        if style_choice is None:
            return

        if style_choice == 'custom':
            new_ls = inquirer.text(
                message="New line style (e.g. '-', '--', ':', '-.'):",
                default='-'
            ).execute()
        else:
            new_ls = style_choice

        changed = self.editor.change_line_linestyle({selected: new_ls})
        if changed > 0:
            self.print_panel(
                f"Line style of '{selected}' changed to {new_ls}", 
                "Success", "green"
            )
        else:
            self.print_panel("Unable to change style", "Error", "red")

    def change_linewidth_flow(self) -> None:
        """Interactive flow to change line width."""
        labels = self.editor.get_legend_labels()
        if not labels:
            self.print_panel("No label found", "Error", "red")
            return

        choices = [Choice(value=label, name=label) for label in labels.values()]
        choices.append(Separator())
        choices.append(Choice(value=None, name="← Back"))

        selected = inquirer.select(
            message="Select label to change line width:",
            choices=choices,
            default=None
        ).execute()

        if selected is None:
            return

        new_width = inquirer.text(
            message="New width (float):",
            default="1.5"
        ).execute()

        try:
            width_val = float(new_width)
            changed = self.editor.change_line_linewidth({selected: width_val})
            if changed > 0:
                self.print_panel(
                    f"Width of '{selected}' changed to {width_val}", 
                    "Success", "green"
                )
            else:
                self.print_panel("Unable to change width", "Error", "red")
        except ValueError:
            self.print_panel("Invalid width value", "Error", "red")

    def change_alpha_flow(self) -> None:
        """Interactive flow to change transparency (alpha)."""
        labels = self.editor.get_legend_labels()
        if not labels:
            self.print_panel("No label found", "Error", "red")
            return

        choices = [Choice(value=label, name=label) for label in labels.values()]
        choices.append(Separator())
        choices.append(Choice(value=None, name="← Back"))

        selected = inquirer.select(
            message="Select label to change transparency (alpha):",
            choices=choices,
            default=None
        ).execute()

        if selected is None:
            return

        new_alpha = inquirer.text(
            message="New alpha (0–1):",
            default="1.0"
        ).execute()

        try:
            alpha_val = float(new_alpha)
            changed = self.editor.change_line_alpha({selected: alpha_val})
            if changed > 0:
                self.print_panel(
                    f"Alpha of '{selected}' changed to {alpha_val}", 
                    "Success", "green"
                )
            else:
                self.print_panel("Unable to change alpha", "Error", "red")
        except ValueError:
            self.print_panel("Invalid alpha value", "Error", "red")

    def save_flow(self) -> None:
        """Interactive flow to save."""
        from pathlib import Path
        
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
            message="Do you want to specify a custom path?",
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
            self.print_panel(f"Error saving: {e}", "Error", "red")

    def run(self) -> None:
        """Main TUI loop."""
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
                message="What do you want to do?",
                choices=[
                    Choice(value="rename", name="📝 Rename legend labels"),
                    Choice(value="color", name="🎨 Change line colors"),
                    Choice(value="linestyle", name="⎯⎯ Change line style"),
                    Choice(value="linewidth", name="➖ Change line width"),
                    Choice(value="alpha", name="☰ Change transparency"),
                    Separator(),
                    Choice(value="preview", name="👁️  Preview figure"),
                    Choice(value="save", name="💾 Save changes"),
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
            elif action == "linestyle":
                self.change_linestyle_flow()
            elif action == "linewidth":
                self.change_linewidth_flow()
            elif action == "alpha":
                self.change_alpha_flow()
            elif action == "preview":
                self.print_panel(
                    "Not implemented yet, without having a freeze of the TUI", 
                    "Error", "red"
                )
                # self.editor.preview()
            elif action == "save":
                self.save_flow()
            elif action == "reload":
                if self.editor.has_modifications():
                    confirm = inquirer.confirm(
                        message="You have unsaved changes. Reload anyway?",
                        default=False
                    ).execute()
                    if not confirm:
                        continue

                self.editor.load()
                self.print_panel("File reloaded", "Success", "green")
            elif action == "exit":
                if self.editor.has_modifications():
                    confirm = inquirer.confirm(
                        message="You have unsaved changes. Exit anyway?",
                        default=False
                    ).execute()
                    if not confirm:
                        continue

                if self.console:
                    self.console.print("\n[cyan]Goodbye! 👋[/]\n")
                else:
                    print("\nGoodbye!\n")
                break