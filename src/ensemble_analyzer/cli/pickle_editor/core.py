import pickle
import logging
import warnings
from pathlib import Path
from typing import Dict, Optional
import copy

try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    from matplotlib.axes import Axes
except ImportError as e:
    raise ImportError(
        "matplotlib not installed. Run: pip install matplotlib"
    ) from e


logger = logging.getLogger(__name__)


class PickleSecurityError(Exception):
    pass


class MatplotlibPickleEditor:

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
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                with open(self.pickle_path, 'rb') as f:
                    obj = pickle.load(f)
            except pickle.UnpicklingError as e:
                raise PickleSecurityError(
                    f"Pickle file corrupted or invalid: {e}"
                ) from e

        if not isinstance(obj, Figure):
            if self.strict_validation:
                raise PickleSecurityError(
                    f"Object is not matplotlib.figure.Figure, but {type(obj)}"
                )
            else:
                logger.warning(f"WARNING: unexpected type {type(obj)}")

        self.figure = obj

        self._patch_stale_pickle_state()

        if self.figure.axes:
            self.axes = self.figure.axes[0]
        else:
            raise PickleSecurityError("No axes found in figure")

    def _patch_stale_pickle_state(self):
        from matplotlib.artist import Artist
        from matplotlib.font_manager import FontProperties
        from matplotlib.patches import FancyBboxPatch
        from matplotlib.text import Text
        from collections.abc import Iterable

        def _walk_artists(obj, seen):
            if id(obj) in seen:
                return
            seen.add(id(obj))
            yield obj
            for attr in ('axes', 'children', '_children', 'child_axes', 'lines',
                         'patches', 'texts', 'images', 'legends',
                         'tables'):
                children = getattr(obj, attr, [])
                if not isinstance(children, Iterable):
                    continue
                for child in children:
                    if isinstance(child, Artist):
                        yield from _walk_artists(child, seen)
            for attr in ('xaxis', 'yaxis'):
                child = getattr(obj, attr, None)
                if child is not None and isinstance(child, Artist):
                    yield from _walk_artists(child, seen)

        seen = set()
        artists = list(_walk_artists(self.figure, seen))

        for a in artists:
            if hasattr(a, 'xaxis') and hasattr(a, 'yaxis') and not hasattr(a, '_axis_map'):
                a._axis_map = {'x': a.xaxis, 'y': a.yaxis}

        self._fix_pickle_state()

    def get_legend_labels(self) -> Dict[int, str]:
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
        if not self.axes:
            raise RuntimeError("You must call load() first")

        legend = self.axes.get_legend()
        if not legend:
            return {}

        lines = self.axes.get_lines()
        colors = {}

        for line, text in zip(lines, legend.get_texts()):
            label = text.get_text()
            color = mpl.colors.to_hex(line.get_color())
            colors[label] = color

        return colors

    def rename_legend_labels(self, mapping: Dict[str, str]) -> int:
        if not self.axes:
            raise RuntimeError("You must call load() first")

        legend = self.axes.get_legend()
        if not legend:
            return 0

        changed = 0
        for text in legend.get_texts():
            current_label = text.get_text()
            if current_label in mapping:
                text.set_text(mapping[current_label])
                changed += 1
                self._modifications_made = True

        return changed

    def change_line_colors(self, label_color_map: Dict[str, str]) -> int:
        if not self.axes:
            raise RuntimeError("You must call load() first")

        legend = self.axes.get_legend()
        if not legend:
            return 0

        lines = self.axes.get_lines()
        legend_texts = legend.get_texts()
        legend_lines = legend.get_lines()

        changed = 0
        for line, leg_line, text in zip(lines, legend_lines, legend_texts):
            label = text.get_text()
            if label in label_color_map:
                color = label_color_map[label]
                try:
                    line.set_color(color)
                    leg_line.set_color(color)
                    changed += 1
                    self._modifications_made = True
                except ValueError as e:
                    logger.warning(f"Invalid color '{color}' for '{label}': {e}")

        return changed

    def change_line_linestyle(self, style_map: Dict[str, str]) -> int:
        if not self.axes:
            raise RuntimeError("You must call load() first")

        legend = self.axes.get_legend()
        if not legend:
            return 0

        lines = self.axes.get_lines()
        legend_texts = legend.get_texts()
        legend_lines = legend.get_lines()

        changed = 0
        for line, leg_line, text in zip(lines, legend_lines, legend_texts):
            label = text.get_text()
            if label in style_map:
                style = style_map[label]
                try:
                    line.set_linestyle(style)
                    leg_line.set_linestyle(style)
                    changed += 1
                    self._modifications_made = True
                except Exception as e:
                    logger.warning(f"Invalid style '{style}' for '{label}': {e}")

        return changed

    def change_line_linewidth(self, width_map: Dict[str, float]) -> int:
        if not self.axes:
            raise RuntimeError("You must call load() first")

        legend = self.axes.get_legend()
        if not legend:
            return 0

        lines = self.axes.get_lines()
        legend_texts = legend.get_texts()
        legend_lines = legend.get_lines()

        changed = 0
        for line, leg_line, text in zip(lines, legend_lines, legend_texts):
            label = text.get_text()
            if label in width_map:
                width = width_map[label]
                try:
                    line.set_linewidth(width)
                    leg_line.set_linewidth(width)
                    changed += 1
                    self._modifications_made = True
                except Exception as e:
                    logger.warning(f"Invalid width '{width}' for '{label}': {e}")

        return changed

    def change_line_alpha(self, alpha_map: Dict[str, float]) -> int:
        if not self.axes:
            raise RuntimeError("You must call load() first")

        legend = self.axes.get_legend()
        if not legend:
            return 0

        lines = self.axes.get_lines()
        legend_texts = legend.get_texts()
        legend_lines = legend.get_lines()

        changed = 0
        for line, leg_line, text in zip(lines, legend_lines, legend_texts):
            label = text.get_text()
            if label in alpha_map:
                alpha = alpha_map[label]
                try:
                    if not 0 <= alpha <= 1:
                        logger.warning(f"Alpha must be between 0 and 1, received {alpha}")
                        continue
                    line.set_alpha(alpha)
                    leg_line.set_alpha(alpha)
                    changed += 1
                    self._modifications_made = True
                except Exception as e:
                    logger.warning(f"Invalid alpha '{alpha}' for '{label}': {e}")

        return changed

    def change_line_visibility(self, visibility_map: Dict[str, bool]) -> int:
        if not self.axes:
            raise RuntimeError("You must call load() first")

        legend = self.axes.get_legend()
        if not legend:
            return 0

        lines = self.axes.get_lines()
        legend_texts = legend.get_texts()
        legend_lines = legend.get_lines()

        changed = 0
        for line, leg_line, text in zip(lines, legend_lines, legend_texts):
            label = text.get_text()
            if label in visibility_map:
                visible = visibility_map[label]
                try:
                    line.set_visible(visible)
                    leg_line.set_visible(visible)
                    text.set_alpha(1.0 if visible else 0.5)
                    changed += 1
                    self._modifications_made = True
                except Exception as e:
                    logger.warning(f"Invalid visibility '{visible}' for '{label}': {e}")

        return changed

    @staticmethod
    def _fix_pickle_state():
        import gc
        from matplotlib.font_manager import FontProperties
        from matplotlib.patches import FancyBboxPatch
        for obj in gc.get_objects():
            if isinstance(obj, FontProperties):
                for key, val in list(obj.__dict__.items()):
                    if isinstance(val, list):
                        obj.__dict__[key] = tuple(val)
            if isinstance(obj, FancyBboxPatch) and not hasattr(obj, '_original_hatchcolor'):
                obj._original_hatchcolor = getattr(obj, '_hatch_color', None)

    def save(self, output_path: Optional[Path] = None,
             format: str = 'pickle') -> Path:
        if not self.figure:
            raise RuntimeError("You must call load() first")

        self._fix_pickle_state()

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
        if not self.figure:
            raise RuntimeError("You must call load() first")

        plt.show()

    def has_modifications(self) -> bool:
        return self._modifications_made
