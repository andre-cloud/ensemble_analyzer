#!/usr/bin/env python3

import argparse
import pickle
import sys
import logging
from pathlib import Path
from typing import Dict, Optional, Tuple
import warnings

# Matplotlib potrebbe non essere presente
try:
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    from matplotlib.axes import Axes
except ImportError:
    print("ERRORE: matplotlib non installato. Esegui: pip install matplotlib", 
          file=sys.stderr)
    sys.exit(1)


# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)s: %(message)s'
)
logger = logging.getLogger(__name__)


class PickleSecurityError(Exception):
    """Eccezione per problemi di sicurezza nel caricamento pickle."""
    pass


class MatplotlibPickleEditor:
    """Editor per modificare legende e colori in figure matplotlib pickled."""
    
    def __init__(self, pickle_path: Path, strict_validation: bool = True):
        """
        Args:
            pickle_path: Path al file pickle
            strict_validation: Se True, valida rigorosamente il tipo dell'oggetto
        
        Raises:
            FileNotFoundError: Se il file non esiste
            PickleSecurityError: Se la validazione fallisce
        """
        self.pickle_path = pickle_path
        self.strict_validation = strict_validation
        self.figure: Optional[Figure] = None
        self.axes: Optional[Axes] = None
        
        if not self.pickle_path.exists():
            raise FileNotFoundError(f"File non trovato: {self.pickle_path}")
    
    def load(self) -> None:
        """Carica e valida il pickle.
        
        Raises:
            PickleSecurityError: Se l'oggetto non è una Figure matplotlib valida
            pickle.UnpicklingError: Se il file è corrotto
        """
        logger.info(f"Caricamento pickle: {self.pickle_path}")
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            try:
                with open(self.pickle_path, 'rb') as f:
                    obj = pickle.load(f)
            except pickle.UnpicklingError as e:
                raise PickleSecurityError(
                    f"File pickle corrotto o non valido: {e}"
                ) from e
        
        # Validazione tipo
        if not isinstance(obj, Figure):
            if self.strict_validation:
                raise PickleSecurityError(
                    f"Oggetto non è matplotlib.figure.Figure, ma {type(obj)}"
                )
            else:
                logger.warning(
                    f"ATTENZIONE: Tipo inatteso {type(obj)}, procedo comunque"
                )
        
        self.figure = obj
        
        # Ottieni primo axes (assunzione: single subplot)
        if self.figure.axes:
            self.axes = self.figure.axes[0]
            logger.info(f"Trovati {len(self.figure.axes)} axes, uso il primo")
        else:
            raise PickleSecurityError("Nessun axes trovato nella figure")
    
    def get_legend_labels(self) -> Dict[int, str]:
        """Ottiene le label correnti della legenda.
        
        Returns:
            Dict mappando indice -> label testuale
        """
        if not self.axes:
            raise RuntimeError("Devi chiamare load() prima")
        
        legend = self.axes.get_legend()
        if not legend:
            logger.warning("Nessuna legenda presente")
            return {}
        
        labels = {}
        for idx, text in enumerate(legend.get_texts()):
            labels[idx] = text.get_text()
        
        return labels
    
    def rename_legend_labels(self, mapping: Dict[str, str]) -> int:
        """Rinomina le label della legenda.
        
        Args:
            mapping: Dict {vecchio_nome: nuovo_nome}
        
        Returns:
            Numero di label rinominate
        """
        if not self.axes:
            raise RuntimeError("Devi chiamare load() prima")
        
        legend = self.axes.get_legend()
        if not legend:
            logger.warning("Nessuna legenda da rinominare")
            return 0
        
        renamed = 0
        for text in legend.get_texts():
            current_label = text.get_text()
            if current_label in mapping:
                text.set_text(mapping[current_label])
                logger.info(f"'{current_label}' -> '{mapping[current_label]}'")
                renamed += 1
        
        return renamed
    
    def change_line_colors(self, label_color_map: Dict[str, str]) -> int:
        """Cambia colori delle linee basandosi sulla label.
        
        Args:
            label_color_map: Dict {label: colore}
                colore può essere: 'red', '#FF0000', (1,0,0), etc.
        
        Returns:
            Numero di linee modificate
        """
        if not self.axes:
            raise RuntimeError("Devi chiamare load() prima")
        
        lines = self.axes.get_lines()
        legend = self.axes.get_legend()
        
        if not legend:
            logger.warning("Nessuna legenda per mappare label->linee")
            return 0
        
        # Mappa label -> line objects
        label_to_lines = {}
        for line, text in zip(lines, legend.get_texts()):
            label = text.get_text()
            label_to_lines[label] = line
        
        changed = 0
        for label, color in label_color_map.items():
            if label in label_to_lines:
                try:
                    label_to_lines[label].set_color(color)
                    logger.info(f"Colore '{label}' -> {color}")
                    changed += 1
                except ValueError as e:
                    logger.error(f"Colore invalido '{color}': {e}")
            else:
                logger.warning(f"Label '{label}' non trovata")
        
        return changed
    
    def save(self, output_path: Optional[Path] = None, 
             format: str = 'pickle') -> Path:
        """Salva la figure modificata.
        
        Args:
            output_path: Path di output (None = overwrite input)
            format: 'pickle' o formato immagine ('png', 'pdf', etc.)
        
        Returns:
            Path del file salvato
        """
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
            logger.info(f"Pickle salvato: {output_path}")
        else:
            self.figure.savefig(output_path, format=format, dpi=300, 
                               bbox_inches='tight')
            logger.info(f"Immagine salvata: {output_path}")
        
        return output_path
    
    def preview(self) -> None:
        """Mostra la figure per preview (blocca l'esecuzione)."""
        if not self.figure:
            raise RuntimeError("Devi chiamare load() prima")
        
        plt.show()


def parse_mapping_file(filepath: Path) -> Dict[str, str]:
    """Parsa file di mapping nel formato: vecchio_nome=nuovo_nome
    
    Args:
        filepath: Path al file di configurazione
    
    Returns:
        Dict di mappatura
    
    Example file content:
        Protocol 1=Protocollo A
        Experimental=Sperimentale
        # commento ignorato
    """
    mapping = {}
    
    with open(filepath, 'r', encoding='utf-8') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            
            # Salta commenti e linee vuote
            if not line or line.startswith('#'):
                continue
            
            if '=' not in line:
                logger.warning(
                    f"Linea {line_num} ignorata (formato invalido): {line}"
                )
                continue
            
            old, new = line.split('=', 1)
            mapping[old.strip()] = new.strip()
    
    return mapping


def main():
    """Entry point CLI."""
    parser = argparse.ArgumentParser(
        description='Modifica legende e colori in pickle matplotlib',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Esempi d'uso:
  
  # Lista label correnti
  %(prog)s plot.pkl --list
  
  # Rinomina singola label
  %(prog)s plot.pkl --rename "Protocol 1" "Protocollo A"
  
  # Rinomina da file di configurazione
  %(prog)s plot.pkl --rename-file mappings.txt
  
  # Cambia colore
  %(prog)s plot.pkl --color "Experimental" red
  
  # Operazioni multiple + salva come nuova figura
  %(prog)s plot.pkl --rename "Protocol 1" "Proto A" \\
                    --color "Proto A" "#FF5733" \\
                    --output modified.pkl
  
  # Preview prima di salvare
  %(prog)s plot.pkl --rename "Protocol 1" "Proto A" --preview
  
  # Esporta come PNG invece di pickle
  %(prog)s plot.pkl --rename "Protocol 1" "Proto A" \\
                    --output result.png --format png

Formato file mapping (--rename-file):
  vecchio_nome=nuovo_nome
  Protocol 1=Protocollo A
  Experimental=Sperimentale
  # commenti iniziano con #
        """
    )
    
    parser.add_argument(
        'pickle_file',
        type=Path,
        help='File pickle contenente matplotlib Figure'
    )
    
    # Azioni di query
    query_group = parser.add_argument_group('query')
    query_group.add_argument(
        '--list', '-l',
        action='store_true',
        help='Lista le label correnti e termina'
    )
    
    # Azioni di modifica
    mod_group = parser.add_argument_group('modifiche')
    mod_group.add_argument(
        '--rename', '-r',
        nargs=2,
        metavar=('OLD', 'NEW'),
        action='append',
        help='Rinomina label (può essere specificato più volte)'
    )
    mod_group.add_argument(
        '--rename-file', '-rf',
        type=Path,
        metavar='FILE',
        help='File contenente mappature OLD=NEW (una per linea)'
    )
    mod_group.add_argument(
        '--color', '-c',
        nargs=2,
        metavar=('LABEL', 'COLOR'),
        action='append',
        help='Cambia colore linea per label (può essere specificato più volte)'
    )
    
    # Output
    output_group = parser.add_argument_group('output')
    output_group.add_argument(
        '--output', '-o',
        type=Path,
        help='File di output (default: sovrascrive input se pickle)'
    )
    output_group.add_argument(
        '--format', '-f',
        default='pickle',
        choices=['pickle', 'png', 'pdf', 'svg', 'jpg'],
        help='Formato output (default: pickle)'
    )
    output_group.add_argument(
        '--preview', '-p',
        action='store_true',
        help='Mostra preview prima di salvare'
    )
    
    # Opzioni avanzate
    adv_group = parser.add_argument_group('avanzate')
    adv_group.add_argument(
        '--no-strict',
        action='store_true',
        help='Disabilita validazione rigorosa del tipo'
    )
    adv_group.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Output verboso'
    )
    
    args = parser.parse_args()
    
    # Setup logging level
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    try:
        # Inizializza editor
        editor = MatplotlibPickleEditor(
            args.pickle_file,
            strict_validation=not args.no_strict
        )
        editor.load()
        
        # Query mode: lista e termina
        if args.list:
            labels = editor.get_legend_labels()
            if not labels:
                print("Nessuna label trovata nella legenda")
                return 0
            
            print(f"\nLabel correnti in {args.pickle_file}:")
            print("-" * 50)
            for idx, label in labels.items():
                print(f"  [{idx}] {label}")
            print()
            return 0
        
        # Costruisci mappature di rinominazione
        rename_map = {}
        
        if args.rename:
            for old, new in args.rename:
                rename_map[old] = new
        
        if args.rename_file:
            file_map = parse_mapping_file(args.rename_file)
            rename_map.update(file_map)
            logger.info(f"Caricate {len(file_map)} mappature da file")
        
        # Applica rinominazioni
        if rename_map:
            count = editor.rename_legend_labels(rename_map)
            logger.info(f"Rinominate {count}/{len(rename_map)} label")
        
        # Costruisci mappature colori
        color_map = {}
        if args.color:
            for label, color in args.color:
                color_map[label] = color
        
        # Applica colori
        if color_map:
            count = editor.change_line_colors(color_map)
            logger.info(f"Cambiati {count}/{len(color_map)} colori")
        
        # Preview
        if args.preview:
            logger.info("Apertura preview (chiudi finestra per continuare)")
            editor.preview()
        
        # Salva se ci sono state modifiche
        if rename_map or color_map:
            output_path = editor.save(args.output, args.format)
            logger.info(f"✓ Modifiche salvate: {output_path}")
        else:
            logger.warning("Nessuna modifica specificata. Usa --help per info.")
        
        return 0
        
    except (FileNotFoundError, PickleSecurityError, RuntimeError) as e:
        logger.error(f"✗ {e}")
        return 1
    except KeyboardInterrupt:
        logger.info("\nInterrotto dall'utente")
        return 130
    except Exception as e:
        logger.error(f"✗ Errore inatteso: {e}", exc_info=args.verbose)
        return 1


if __name__ == '__main__':
    sys.exit(main())