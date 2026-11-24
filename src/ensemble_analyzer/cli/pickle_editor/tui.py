from pathlib import Path
from typing import TYPE_CHECKING, Optional

from InquirerPy import inquirer
from InquirerPy.base.control import Choice
from InquirerPy.separator import Separator
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

if TYPE_CHECKING:
    from .core import MatplotlibPickleEditor


class InteractiveTUI:

    def __init__(self, editor: 'MatplotlibPickleEditor'):
        self.editor = editor
        self.console = Console()

    def print_panel(self, message: str, title: str = "Info",
                    style: str = "cyan") -> None:
        self.console.print(Panel(message, title=title, border_style=style))

    def _report(self, changed: int, ok: str, fail: str) -> None:
        if changed > 0:
            self.print_panel(ok, "Success", "green")
        else:
            self.print_panel(fail, "Error", "red")

    def _pick_label(self, message: str, describe=None) -> Optional[str]:
        """Ask the user to pick a legend label (None = back / no labels)."""
        labels = self.editor.get_legend_labels()
        if not labels:
            self.print_panel("No label found", "Error", "red")
            return None
        choices = [Choice(value=label, name=describe(idx, label) if describe else label)
                   for idx, label in labels.items()]
        choices += [Separator(), Choice(value=None, name="← Back")]
        return inquirer.select(message=message, choices=choices, default=None).execute()

    def show_current_state(self) -> None:
        labels = self.editor.get_legend_labels()
        colors = self.editor.get_line_colors()

        if not labels:
            self.print_panel("No legend found", "Warning", "yellow")
            return

        table = Table(title="Current Legend State", show_header=True,
                      header_style="bold cyan")
        table.add_column("Index", style="dim", width=8)
        table.add_column("Label", style="bold")
        table.add_column("Color", style="magenta")
        for idx, label in labels.items():
            table.add_row(str(idx), label, colors.get(label, "N/A"))
        self.console.print(table)

    def rename_labels_flow(self) -> None:
        selected = self._pick_label("Select label to rename:",
                                    lambda idx, label: f"{label} [{idx}]")
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
        colors = self.editor.get_line_colors()
        selected = self._pick_label(
            "Select label to change color:",
            lambda idx, label: f"{label} (current: {colors.get(label, 'N/A')})")
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
            color_choices = [Choice(value=c, name=c) for c in self.editor.COMMON_COLORS]
            color_choices += [Separator(), Choice(value=None, name="← Cancel")]
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
            self._report(self.editor.change_line_colors({selected: new_color}),
                         f"Color of '{selected}' changed to {new_color}",
                         "Unable to change color (invalid color?)")

    def change_linestyle_flow(self) -> None:
        selected = self._pick_label("Select label to change line style:")
        if selected is None:
            return

        style_choice = inquirer.select(
            message="Select style:",
            choices=[
                Choice(value='-', name="Solid (-)"),
                Choice(value='--', name="Dashed (--)"),
                Choice(value=':', name="Dotted (:)"),
                Choice(value='-.', name="Dash-dot (-.)"),
                Choice(value='custom', name="Enter manually"),
                Choice(value=None, name="← Cancel")
            ],
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

        self._report(self.editor.change_line_linestyle({selected: new_ls}),
                     f"Line style of '{selected}' changed to {new_ls}",
                     "Unable to change style")

    def _float_flow(self, message: str, prompt: str, default: str,
                    what: str, setter) -> None:
        selected = self._pick_label(message)
        if selected is None:
            return
        raw = inquirer.text(message=prompt, default=default).execute()
        try:
            value = float(raw)
        except ValueError:
            self.print_panel(f"Invalid {what.lower()} value", "Error", "red")
            return
        self._report(setter({selected: value}),
                     f"{what} of '{selected}' changed to {value}",
                     f"Unable to change {what.lower()}")

    def change_linewidth_flow(self) -> None:
        self._float_flow("Select label to change line width:", "New width (float):",
                         "1.5", "Width", self.editor.change_line_linewidth)

    def change_alpha_flow(self) -> None:
        self._float_flow("Select label to change transparency (alpha):",
                         "New alpha (0–1):", "1.0", "Alpha",
                         self.editor.change_line_alpha)

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
        for a in ["x", "y"] if axis_choice == "xy" else [axis_choice]:
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
            except ValueError:
                self.print_panel(f"Invalid number for {a.upper()} limits", "Error", "red")
                continue
            (self.editor.set_xlim if a == "x" else self.editor.set_ylim)(lo_f, hi_f)
            self.print_panel(f"{a.upper()} set to [{lo_f:.3f}, {hi_f:.3f}]",
                             "Success", "green")

    def change_visibility_flow(self) -> None:
        selected = self._pick_label("Select label to toggle visibility:")
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
            status = "Visible" if visible else "Hidden"
            self._report(self.editor.change_line_visibility({selected: visible}),
                         f"'{selected}' is now {status}",
                         "Unable to change visibility")

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

        output_path = None
        if inquirer.confirm(message="Do you want to specify a custom path?",
                            default=False).execute():
            output_path = Path(inquirer.text(
                message="Output path:",
                default=f"{default_name}.{format_choice}"
            ).execute())

        try:
            saved_path = self.editor.save(output_path, format_choice)
            self.print_panel(f"Saved: {saved_path}", "Success", "green")
        except Exception as e:
            self.print_panel(f"Error saving: {e}", "Error", "red")

    def _confirm_discard(self, message: str) -> bool:
        return not self.editor.has_modifications() or inquirer.confirm(
            message=message, default=False).execute()

    def run(self) -> None:
        self.console.clear()
        self.console.print(
            Panel.fit(
                "[bold cyan]Interactive Matplotlib Pickle Editor[/]\n"
                f"File: {self.editor.pickle_path}",
                border_style="cyan"
            )
        )

        flows = {
            "rename": self.rename_labels_flow,
            "color": self.change_colors_flow,
            "linestyle": self.change_linestyle_flow,
            "linewidth": self.change_linewidth_flow,
            "alpha": self.change_alpha_flow,
            "visibility": self.change_visibility_flow,
            "limits": self.set_limits_flow,
            "save": self.save_flow,
        }

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

            if action in flows:
                flows[action]()
            elif action == "reload":
                if self._confirm_discard("You have unsaved changes. Reload anyway?"):
                    self.editor.load()
                    self.print_panel("File reloaded", "Success", "green")
            elif action == "exit":
                if self._confirm_discard("You have unsaved changes. Exit anyway?"):
                    self.console.print("\n[cyan]Goodbye! 👋[/]\n")
                    break
