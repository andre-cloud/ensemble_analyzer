import sys
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .core import MatplotlibPickleEditor

try:
    from InquirerPy import inquirer
    from InquirerPy.base.control import Choice
    from InquirerPy.separator import Separator
    INQUIRER_AVAILABLE = True
except ImportError:
    INQUIRER_AVAILABLE = False

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

    def __init__(self, editor: 'MatplotlibPickleEditor'):
        if not INQUIRER_AVAILABLE:
            raise RuntimeError(
                "InquirerPy not installed. Run: pip install InquirerPy"
            )

        self.editor = editor
        self.console = console if RICH_AVAILABLE else None

    def print_panel(self, message: str, title: str = "Info",
                    style: str = "cyan") -> None:
        if self.console:
            self.console.print(Panel(message, title=title, border_style=style))
        else:
            print(f"\n{title}: {message}\n")

    def show_current_state(self) -> None:

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
        labels = self.editor.get_legend_labels()

        if not labels:
            self.print_panel("No labels to rename", "Error", "red")
            return

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

        new_name = inquirer.text(
            message=f"New name for '{selected}':",
            default=selected
        ).execute()

        if new_name and new_name != selected:
            self.editor.rename_legend_labels({selected: new_name})
            self.print_panel(f"'{selected}' → '{new_name}'", "Success", "green")

    def change_colors_flow(self) -> None:
        labels = self.editor.get_legend_labels()
        colors = self.editor.get_line_colors()

        if not labels:
            self.print_panel("No label found", "Error", "red")
            return

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

    def set_limits_flow(self) -> None:
        if not self.editor.axes:
            self.print_panel("No axes found", "Error", "red")
            return
        cur_x = self.editor.axes.get_xlim()
        cur_y = self.editor.axes.get_ylim()
        self.print_panel(
            f"Current X: [{cur_x[0]:.1f}, {cur_x[1]:.1f}]\n"
            f"Current Y: [{cur_y[0]:.3f}, {cur_y[1]:.3f}]",
            "Current Limits", "cyan"
        )
        axis_choice = inquirer.select(
            message="Which axis?",
            choices=[
                Choice(value="x", name="X axis"),
                Choice(value="y", name="Y axis"),
                Choice(value="xy", name="Both"),
                Choice(value=None, name="← Cancel"),
            ],
            default="x"
        ).execute()
        if axis_choice is None:
            return
        axes_to_set = ["x", "y"] if axis_choice == "xy" else [axis_choice]
        for a in axes_to_set:
            cur = cur_x if a == "x" else cur_y
            lo = inquirer.text(
                message=f"{a.upper()} min (current: {cur[0]:.3f}):",
                default=str(cur[0])
            ).execute()
            hi = inquirer.text(
                message=f"{a.upper()} max (current: {cur[1]:.3f}):",
                default=str(cur[1])
            ).execute()
            try:
                lo_f, hi_f = float(lo), float(hi)
                if a == "x":
                    self.editor.set_xlim(lo_f, hi_f)
                else:
                    self.editor.set_ylim(lo_f, hi_f)
                self.print_panel(
                    f"{a.upper()} set to [{lo_f:.3f}, {hi_f:.3f}]",
                    "Success", "green"
                )
            except ValueError:
                self.print_panel(f"Invalid number for {a.upper()} limits", "Error", "red")

    def change_visibility_flow(self) -> None:
        labels = self.editor.get_legend_labels()
        if not labels:
            self.print_panel("No label found", "Error", "red")
            return

        choices = [Choice(value=label, name=label) for label in labels.values()]
        choices.append(Separator())
        choices.append(Choice(value=None, name="← Back"))

        selected = inquirer.select(
            message="Select label to toggle visibility:",
            choices=choices,
            default=None
        ).execute()

        if selected is None:
            return

        visible = inquirer.select(
            message=f"Set visibility for '{selected}':",
            choices=[
                Choice(value=True, name="Visible (Show)"),
                Choice(value=False, name="Hidden (Hide)"),
                Choice(value=None, name="← Cancel")
            ],
            default=True
        ).execute()

        if visible is not None:
            changed = self.editor.change_line_visibility({selected: visible})
            status = "Visible" if visible else "Hidden"
            if changed > 0:
                self.print_panel(
                    f"'{selected}' is now {status}",
                    "Success", "green"
                )
            else:
                self.print_panel("Unable to change visibility", "Error", "red")

    def save_flow(self) -> None:

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
            self.show_current_state()

            action = inquirer.select(
                message="What do you want to do?",
                choices=[
                    Choice(value="rename", name="📝 Rename legend labels"),
                    Choice(value="color", name="🎨 Change line colors"),
                    Choice(value="linestyle", name="⎯⎯ Change line style"),
                    Choice(value="linewidth", name="➖ Change line width"),
                    Choice(value="alpha", name="☰ Change transparency"),
                    Choice(value="visibility", name="👁️  Change line visibility"),
                    Choice(value="limits", name="⟷  Set axis limits"),
                    Separator(),
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
            elif action == 'visibility':
                self.change_visibility_flow()
            elif action == "limits":
                self.set_limits_flow()
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
