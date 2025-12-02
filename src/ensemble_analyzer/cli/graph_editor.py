#!/usr/bin/env python3
"""
TUI interattiva per modificare legende e colori in pickle di figure matplotlib.

Mode predefinito: Interactive TUI
Mode batch disponibile con flag --batch

Sicurezza: Carica SOLO pickle da fonti fidate. Il modulo pickle può eseguire
codice arbitrario durante la deserializzazione.

Dipendenze:
    pip install matplotlib InquirerPy rich

Autore: Senior Python Developer
Versione: 2.0.0
"""

import argparse
import pickle
import sys
import logging
from pathlib import Path
from typing import Dict, Optional, List
import warnings

# Matplotlib (required)
try:
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    from matplotlib.axes import Axes
except ImportError:
    print("ERRORE: matplotlib non installato. Esegui: pip install matplotlib", 
          file=sys.stderr)
    sys.exit(1)

# InquirerPy (required per TUI)
try:
    from InquirerPy import inquirer
    from InquirerPy.base.control import Choice
    from InquirerPy.separator import Separator
    INQUIRER_AVAILABLE = True
except ImportError:
    INQUIRER_AVAILABLE = False

# Rich (required per output colorato)
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
    """Eccezione per problemi di sicurezza nel caricamento pickle."""
    pass


class MatplotlibPickleEditor:
    """Editor per modificare legende e colori in figure matplotlib pickled."""
    
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
            raise FileNotFoundError(f"File non trovato: {self.pickle_path}")
    
    def load(self) -> None:
        """Carica e valida il pickle."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            try:
                with open(self.pickle_path, 'rb') as f:
                    obj = pickle.load(f)
            except pickle.UnpicklingError as e:
                raise PickleSecurityError(
                    f"File pickle corrotto o non valido: {e}"
                ) from e
        
        if not isinstance(obj, Figure):
            if self.strict_validation:
                raise PickleSecurityError(
                    f"Oggetto non è matplotlib.figure.Figure, ma {type(obj)}"
                )
            else:
                logger.warning(f"ATTENZIONE: Tipo inatteso {type(obj)}")
        
        self.figure = obj
        
        if self.figure.axes:
            self.axes = self.figure.axes[0]
        else:
            raise PickleSecurityError("Nessun axes trovato nella figure")
    
    def get_legend_labels(self) -> Dict[int, str]:
        """Ottiene le label correnti della legenda."""
        if not self.axes:
            raise RuntimeError("Devi chiamare load() prima")
        
        legend = self.axes.get_legend()
        if not legend:
            return {}
        
        labels = {}
        for idx, text in enumerate(legend.get_texts()):
            labels[idx] = text.get_text()
        
        return labels
    
    def get_line_colors(self) -> Dict[str, str]:
        """Ottiene i colori correnti delle linee."""
        if not self.axes:
            raise RuntimeError("Devi chiamare load() prima")
        
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
        """Rinomina le label della legenda."""
        if not self.axes:
            raise RuntimeError("Devi chiamare load() prima")
        
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
        """Cambia colori delle linee."""
        if not self.axes:
            raise RuntimeError("Devi chiamare load() prima")
        
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
                except ValueError:
                    pass
        
        return changed
    
    def save(self, output_path: Optional[Path] = None, 
             format: str = 'pickle') -> Path:
        """Salva la figure modificata."""
        if not self.figure:
            raise RuntimeError("Devi chiamare load() prima")
        
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
        """Mostra la figure per preview."""
        if not self.figure:
            raise RuntimeError("Devi chiamare load() prima")
        plt.show()
    
    def has_modifications(self) -> bool:
        """Controlla se ci sono modifiche non salvate."""
        return self._modifications_made


class InteractiveTUI:
    """Terminal User Interface interattiva usando InquirerPy."""
    
    def __init__(self, editor: MatplotlibPickleEditor):
        if not INQUIRER_AVAILABLE:
            raise RuntimeError(
                "InquirerPy non installato. Esegui: pip install InquirerPy"
            )
        
        self.editor = editor
        self.console = console if RICH_AVAILABLE else None
    
    def print_panel(self, message: str, title: str = "Info", style: str = "cyan") -> None:
        """Stampa un pannello formattato."""
        if self.console:
            self.console.print(Panel(message, title=title, border_style=style))
        else:
            print(f"\n{title}: {message}\n")
    
    def show_current_state(self) -> None:
        """Mostra stato corrente della figura."""
        labels = self.editor.get_legend_labels()
        colors = self.editor.get_line_colors()
        
        if not labels:
            self.print_panel("Nessuna legenda trovata", "Attenzione", "yellow")
            return
        
        if self.console:
            table = Table(title="Stato Corrente Legenda", show_header=True, 
                         header_style="bold cyan")
            table.add_column("Indice", style="dim", width=8)
            table.add_column("Label", style="bold")
            table.add_column("Colore", style="magenta")
            
            for idx, label in labels.items():
                color = colors.get(label, "N/A")
                table.add_row(str(idx), label, color)
            
            self.console.print(table)
        else:
            print("\n=== Stato Corrente ===")
            for idx, label in labels.items():
                color = colors.get(label, "N/A")
                print(f"  [{idx}] {label} (colore: {color})")
            print()
    
    def rename_labels_flow(self) -> None:
        """Flow interattivo per rinominare label."""
        labels = self.editor.get_legend_labels()
        
        if not labels:
            self.print_panel("Nessuna label da rinominare", "Errore", "red")
            return
        
        # Selezione label da rinominare
        choices = [
            Choice(value=label, name=f"{label} [{idx}]") 
            for idx, label in labels.items()
        ]
        choices.append(Separator())
        choices.append(Choice(value=None, name="← Torna indietro"))
        
        selected = inquirer.select(
            message="Seleziona label da rinominare:",
            choices=choices,
            default=None
        ).execute()
        
        if selected is None:
            return
        
        # Input nuovo nome
        new_name = inquirer.text(
            message=f"Nuovo nome per '{selected}':",
            default=selected
        ).execute()
        
        if new_name and new_name != selected:
            self.editor.rename_legend_labels({selected: new_name})
            self.print_panel(f"'{selected}' → '{new_name}'", "Successo", "green")
    
    def change_colors_flow(self) -> None:
        """Flow interattivo per cambiare colori."""
        labels = self.editor.get_legend_labels()
        colors = self.editor.get_line_colors()
        
        if not labels:
            self.print_panel("Nessuna label trovata", "Errore", "red")
            return
        
        # Selezione label
        choices = [
            Choice(value=label, 
                   name=f"{label} (attuale: {colors.get(label, 'N/A')})") 
            for label in labels.values()
        ]
        choices.append(Separator())
        choices.append(Choice(value=None, name="← Torna indietro"))
        
        selected = inquirer.select(
            message="Seleziona label per cambiare colore:",
            choices=choices,
            default=None
        ).execute()
        
        if selected is None:
            return
        
        # Scelta metodo input colore
        color_method = inquirer.select(
            message="Come vuoi specificare il colore?",
            choices=[
                Choice(value="preset", name="Scegli da palette predefinita"),
                Choice(value="custom", name="Inserisci manualmente (nome o hex)"),
                Choice(value=None, name="← Annulla")
            ],
            default="preset"
        ).execute()
        
        if color_method is None:
            return
        
        if color_method == "preset":
            # Palette predefinita con preview
            color_choices = [
                Choice(value=c, name=f"{c}") 
                for c in self.editor.COMMON_COLORS
            ]
            color_choices.append(Separator())
            color_choices.append(Choice(value=None, name="← Annulla"))
            
            new_color = inquirer.select(
                message="Seleziona colore:",
                choices=color_choices,
                default=None
            ).execute()
        else:
            # Input manuale
            new_color = inquirer.text(
                message="Colore (nome o hex #RRGGBB):",
                validate=lambda x: len(x) > 0
            ).execute()
        
        if new_color:
            self.editor.change_line_colors({selected: new_color})
            self.print_panel(
                f"Colore di '{selected}' cambiato in {new_color}", 
                "Successo", "green"
            )
    
    def save_flow(self) -> None:
        """Flow interattivo per salvare."""
        if not self.editor.has_modifications():
            self.print_panel("Nessuna modifica da salvare", "Info", "yellow")
            return
        
        # Formato output
        format_choice = inquirer.select(
            message="Formato di salvataggio:",
            choices=[
                Choice(value="pickle", name="Pickle (modificabile successivamente)"),
                Choice(value="png", name="PNG (immagine)"),
                Choice(value="pdf", name="PDF (vettoriale)"),
                Choice(value="svg", name="SVG (vettoriale)"),
                Choice(value=None, name="← Annulla")
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
            message="Vuoi specificare un path custom?",
            default=False
        ).execute()
        
        output_path = None
        if custom_path:
            path_str = inquirer.text(
                message="Path di output:",
                default=f"{default_name}.{format_choice}"
            ).execute()
            output_path = Path(path_str)
        
        try:
            saved_path = self.editor.save(output_path, format_choice)
            self.print_panel(f"Salvato: {saved_path}", "Successo", "green")
        except Exception as e:
            self.print_panel(f"Errore salvataggio: {e}", "Errore", "red")
    
    def run(self) -> None:
        """Main loop TUI."""
        if self.console:
            self.console.clear()
            self.console.print(
                Panel.fit(
                    "[bold cyan]Editor Interattivo Matplotlib Pickle[/]\n"
                    f"File: {self.editor.pickle_path}",
                    border_style="cyan"
                )
            )
        
        while True:
            # Mostra stato
            self.show_current_state()
            
            # Menu principale
            action = inquirer.select(
                message="Cosa vuoi fare?",
                choices=[
                    Choice(value="rename", name="📝 Rinomina label legenda"),
                    Choice(value="color", name="🎨 Cambia colori linee"),
                    Separator(),
                    Choice(value="preview", name="👁️  Preview figura"),
                    Choice(value="save", name="💾 Salva modifiche"),
                    Separator(),
                    Choice(value="reload", name="🔄 Ricarica file originale"),
                    Choice(value="exit", name="🚪 Esci"),
                ],
                default="rename"
            ).execute()
            
            if action == "rename":
                self.rename_labels_flow()
            elif action == "color":
                self.change_colors_flow()
            elif action == "preview":
                self.print_panel("Chiudi la finestra matplotlib per continuare", "Info")
                self.editor.preview()
            elif action == "save":
                self.save_flow()
            elif action == "reload":
                if self.editor.has_modifications():
                    confirm = inquirer.confirm(
                        message="Hai modifiche non salvate. Ricaricare comunque?",
                        default=False
                    ).execute()
                    if not confirm:
                        continue
                
                self.editor.load()
                self.print_panel("File ricaricato", "Successo", "green")
            elif action == "exit":
                if self.editor.has_modifications():
                    confirm = inquirer.confirm(
                        message="Hai modifiche non salvate. Uscire comunque?",
                        default=False
                    ).execute()
                    if not confirm:
                        continue
                
                if self.console:
                    self.console.print("\n[cyan]Arrivederci! 👋[/]\n")
                else:
                    print("\nArrivederci!\n")
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
    """Esegue in batch mode (non interattivo)."""
    try:
        editor = MatplotlibPickleEditor(args.pickle_file, 
                                       strict_validation=not args.no_strict)
        editor.load()
        
        # Query mode: lista e termina
        if args.list:
            labels = editor.get_legend_labels()
            if not labels:
                print("Nessuna label trovata")
                return 0
            
            print(f"\nLabel in {args.pickle_file}:")
            print("-" * 50)
            for idx, label in labels.items():
                colors = editor.get_line_colors()
                color = colors.get(label, "N/A")
                print(f"  [{idx}] {label} (colore: {color})")
            print()
            return 0
        
        # Rinominazioni
        rename_map = {}
        if args.rename:
            for old, new in args.rename:
                rename_map[old] = new
        
        if args.rename_file:
            file_map = parse_mapping_file(args.rename_file)
            rename_map.update(file_map)
        
        if rename_map:
            count = editor.rename_legend_labels(rename_map)
            logger.info(f"Rinominate {count} label")
        
        # Colori
        color_map = {}
        if args.color:
            for label, color in args.color:
                color_map[label] = color
        
        if color_map:
            count = editor.change_line_colors(color_map)
            logger.info(f"Cambiati {count} colori")
        
        # Preview
        if args.preview:
            editor.preview()
        
        # Salva
        if rename_map or color_map:
            output_path = editor.save(args.output, args.format)
            print(f"✓ Salvato: {output_path}")
        else:
            logger.warning("Nessuna modifica specificata")
        
        return 0
        
    except Exception as e:
        logger.error(f"✗ {e}")
        return 1


def main():
    """Entry point CLI."""
    parser = argparse.ArgumentParser(
        description='Editor interattivo/batch per pickle matplotlib',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
MODE PREDEFINITO: Interactive TUI
  python %(prog)s plot.pkl

MODE BATCH (esempi):
  python %(prog)s plot.pkl --batch --list
  python %(prog)s plot.pkl --batch --rename "Protocol 1" "Proto A"
  python %(prog)s plot.pkl --batch --color "Experimental" red --output new.pkl
        """
    )
    
    parser.add_argument('pickle_file', type=Path, 
                       help='File pickle matplotlib')
    
    parser.add_argument('--batch', '-b', action='store_true',
                       help='Mode batch (non interattivo)')
    
    # Opzioni batch mode
    batch_group = parser.add_argument_group('batch mode options')
    batch_group.add_argument('--list', '-l', action='store_true',
                            help='Lista label e termina')
    batch_group.add_argument('--rename', '-r', nargs=2, 
                            metavar=('OLD', 'NEW'), action='append',
                            help='Rinomina label')
    batch_group.add_argument('--rename-file', '-rf', type=Path,
                            help='File mappature OLD=NEW')
    batch_group.add_argument('--color', '-c', nargs=2,
                            metavar=('LABEL', 'COLOR'), action='append',
                            help='Cambia colore')
    batch_group.add_argument('--output', '-o', type=Path,
                            help='File output')
    batch_group.add_argument('--format', '-f', default='pickle',
                            choices=['pickle', 'png', 'pdf', 'svg'],
                            help='Formato output')
    batch_group.add_argument('--preview', '-p', action='store_true',
                            help='Preview prima di salvare')
    
    parser.add_argument('--no-strict', action='store_true',
                       help='Disabilita validazione rigorosa')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Output verboso')
    
    args = parser.parse_args()
    
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Batch mode
    if args.batch:
        return batch_mode(args)
    
    # Interactive TUI mode (default)
    if not INQUIRER_AVAILABLE:
        print("ERRORE: Mode interattivo richiede InquirerPy", file=sys.stderr)
        print("Installa: pip install InquirerPy rich", file=sys.stderr)
        print("\nUsa --batch per mode non interattivo", file=sys.stderr)
        return 1
    
    try:
        editor = MatplotlibPickleEditor(args.pickle_file,
                                       strict_validation=not args.no_strict)
        editor.load()
        
        tui = InteractiveTUI(editor)
        tui.run()
        
        return 0
        
    except KeyboardInterrupt:
        print("\n\nInterrotto dall'utente")
        return 130
    except Exception as e:
        logger.error(f"✗ Errore: {e}", exc_info=args.verbose)
        return 1


if __name__ == '__main__':
    sys.exit(main())