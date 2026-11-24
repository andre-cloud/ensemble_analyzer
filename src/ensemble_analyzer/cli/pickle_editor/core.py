import pickle
import logging
import warnings
from pathlib import Path
from typing import Dict, Optional

import matplotlib as mpl
from matplotlib.figure import Figure
from matplotlib.axes import Axes


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

    def __init__(self, pickle_path: Path):
        self.pickle_path = pickle_path
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
            raise PickleSecurityError(
                f"Object is not matplotlib.figure.Figure, but {type(obj)}"
            )

        old_fig = obj
        if not old_fig.axes:
            raise PickleSecurityError("No axes found in figure")
        old_ax = old_fig.axes[0]

        lines_data = self._extract_lines(old_ax)
        legend_texts = self._extract_legend_texts(old_ax)

        new_fig = Figure(figsize=old_fig.get_size_inches(), dpi=old_fig.get_dpi())
        new_ax = new_fig.add_subplot(111)

        self._copy_axes_props(old_ax, new_ax)

        new_lines = self._replot_lines(new_ax, lines_data)

        # Restore legend labels — prefer saved per-line labels, fall back to extracted
        saved = getattr(old_ax, '_ea_labels', None)
        if saved is None:
            saved = list(legend_texts)
        if saved:
            if len(saved) == len(new_lines):
                texts = [
                    saved[i] if saved[i] is not None else l.get_label()
                    for i, l in enumerate(new_lines)
                ]
            else:
                texts = list(saved) + [
                    l.get_label() for l in new_lines[len(saved):]
                ]
            new_ax.legend(new_lines, texts)

        self.figure = new_fig
        self.axes = new_ax

    def _extract_lines(self, ax: Axes) -> list:
        data = []
        for line in ax.get_lines():
            ld = {
                'xdata': line.get_xdata(),
                'ydata': line.get_ydata(),
                'color': line.get_color(),
                'linestyle': line.get_linestyle(),
                'linewidth': line.get_linewidth(),
                'alpha': line.get_alpha(),
                'visible': line.get_visible(),
                'label': line.get_label(),
                'marker': line.get_marker(),
                'markersize': line.get_markersize(),
                'markerfacecolor': line.get_markerfacecolor(),
                'markeredgecolor': line.get_markeredgecolor(),
                'markevery': line.get_markevery(),
                'zorder': line.get_zorder(),
                'drawstyle': line.get_drawstyle(),
                'dash_capstyle': line.get_dash_capstyle(),
                'dash_joinstyle': line.get_dash_joinstyle(),
                'solid_capstyle': line.get_solid_capstyle(),
                'solid_joinstyle': line.get_solid_joinstyle(),
            }
            data.append(ld)
        return data

    def _extract_legend_texts(self, ax: Axes) -> list:
        legend = ax.get_legend()
        if not legend:
            return []
        return [t.get_text() for t in legend.get_texts()]

    def _copy_axes_props(self, old: Axes, new: Axes):
        xl = old.get_xlabel()
        if xl:
            new.set_xlabel(xl)
        yl = old.get_ylabel()
        if yl:
            new.set_ylabel(yl)
        t = old.get_title()
        if t:
            new.set_title(t)
        new.set_xscale(old.get_xscale())
        new.set_yscale(old.get_yscale())
        new.xaxis.set_ticks(old.get_xticks())
        new.yaxis.set_ticks(old.get_yticks())
        new.set_xlim(old.get_xlim())
        new.set_ylim(old.get_ylim())
        new.xaxis.set_ticklabels([t.get_text() for t in old.get_xticklabels()])
        new.yaxis.set_ticklabels([t.get_text() for t in old.get_yticklabels()])
        self._copy_grid(old, new)
        self._add_secondary_xaxis(old, new)

    def _copy_grid(self, old: Axes, new: Axes):
        x_lines = old.get_xgridlines()
        y_lines = old.get_ygridlines()
        x_on = any(l.get_visible() for l in x_lines) if x_lines else False
        y_on = any(l.get_visible() for l in y_lines) if y_lines else False
        if not x_on and not y_on:
            return
        kw = {}
        for line in x_lines + y_lines:
            if line.get_visible():
                kw['linestyle'] = line.get_linestyle()
                kw['linewidth'] = line.get_linewidth()
                kw['alpha'] = line.get_alpha()
                c = line.get_color()
                if c:
                    kw['color'] = c
                break
        new.grid(x_on or y_on, **kw)
        if not y_on:
            new.yaxis.grid(False)
        if not x_on:
            new.xaxis.grid(False)

    def _add_secondary_xaxis(self, old: Axes, new: Axes):
        from ensemble_analyzer.constants import eV_to_nm
        xl = old.get_xlabel().lower()
        has_nm = 'nm' in xl or 'wavelength' in xl
        has_ev = 'ev' in xl or 'energy' in xl
        if has_nm:
            secax = new.secondary_xaxis("top", functions=(eV_to_nm, eV_to_nm))
            secax.set_xlabel("Energy [eV]")
        elif has_ev:
            secax = new.secondary_xaxis("top", functions=(eV_to_nm, eV_to_nm))
            secax.set_xlabel(r"Wavelength $\lambda$ [nm]")

    def _replot_lines(self, ax: Axes, lines_data: list) -> list:
        new_lines = []
        for ld in lines_data:
            kwargs = {}
            for k in ('color', 'linestyle', 'linewidth', 'visible',
                      'marker', 'markersize', 'markerfacecolor', 'markeredgecolor',
                      'zorder', 'drawstyle', 'dash_capstyle', 'dash_joinstyle',
                      'solid_capstyle', 'solid_joinstyle', 'label'):
                v = ld.get(k)
                if v is not None and v != 'None':
                    kwargs[k] = v
            alpha = ld.get('alpha')
            if alpha is not None:
                kwargs['alpha'] = alpha
            markevery = ld.get('markevery')
            if markevery is not None:
                kwargs['markevery'] = markevery
            line, = ax.plot(ld['xdata'], ld['ydata'], **kwargs)
            new_lines.append(line)
        return new_lines

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
            colors[text.get_text()] = mpl.colors.to_hex(line.get_color())
        return colors

    def rename_legend_labels(self, mapping: Dict[str, str]) -> int:
        if not self.axes:
            raise RuntimeError("You must call load() first")
        legend = self.axes.get_legend()
        if not legend:
            return 0
        changed = 0
        for text in legend.get_texts():
            current = text.get_text()
            if current in mapping:
                text.set_text(mapping[current])
                changed += 1
                self._modifications_made = True
        if changed:
            texts = [t.get_text() for t in legend.get_texts()]
            txt_iter = iter(texts)
            self.axes._ea_labels = [
                next(txt_iter) if l.get_visible() else None
                for l in self.axes.get_lines()
            ]
        return changed

    def _apply(self, mapping: Dict[str, object], setter: str, what: str,
               check=None) -> int:
        """Call ``setter`` on each plotted line (and its legend handle) whose
        legend label is in ``mapping``; return how many were changed."""
        if not self.axes:
            raise RuntimeError("You must call load() first")
        legend = self.axes.get_legend()
        if not legend:
            return 0
        changed = 0
        for line, leg_line, text in zip(self.axes.get_lines(),
                                        legend.get_lines(), legend.get_texts()):
            label = text.get_text()
            if label not in mapping:
                continue
            value = mapping[label]
            if check and not check(value):
                logger.warning(f"Invalid {what} '{value}' for '{label}'")
                continue
            try:
                getattr(line, setter)(value)
                getattr(leg_line, setter)(value)
                changed += 1
                self._modifications_made = True
            except Exception as e:
                logger.warning(f"Invalid {what} '{value}' for '{label}': {e}")
        return changed

    def change_line_colors(self, label_color_map: Dict[str, str]) -> int:
        return self._apply(label_color_map, "set_color", "color")

    def change_line_linestyle(self, style_map: Dict[str, str]) -> int:
        return self._apply(style_map, "set_linestyle", "style")

    def change_line_linewidth(self, width_map: Dict[str, float]) -> int:
        return self._apply(width_map, "set_linewidth", "width")

    def change_line_alpha(self, alpha_map: Dict[str, float]) -> int:
        return self._apply(alpha_map, "set_alpha", "alpha (must be 0-1)",
                           check=lambda a: 0 <= a <= 1)

    def change_line_visibility(self, visibility_map: Dict[str, bool]) -> int:
        if not self.axes:
            raise RuntimeError("You must call load() first")
        legend = self.axes.get_legend()
        if not legend:
            return 0
        lines = self.axes.get_lines()
        legend_texts = legend.get_texts()
        changed = 0
        for line, text in zip(lines, legend_texts):
            label = text.get_text()
            if label in visibility_map:
                visible = visibility_map[label]
                try:
                    line.set_visible(visible)
                    changed += 1
                    self._modifications_made = True
                except Exception as e:
                    logger.warning(f"Invalid visibility '{visible}' for '{label}': {e}")
        if changed:
            self._rebuild_legend()
        return changed

    def _rebuild_legend(self):
        legend = self.axes.get_legend()
        if not legend:
            return
        lines = self.axes.get_lines()
        labels = [t.get_text() for t in legend.get_texts()]
        visible = [(l, lb) for l, lb in zip(lines, labels) if l.get_visible()]
        if visible:
            self.axes.legend([v[0] for v in visible], [v[1] for v in visible])
            new_leg = self.axes.get_legend()
            new_texts = iter(t.get_text() for t in new_leg.get_texts())
            self.axes._ea_labels = [
                next(new_texts) if l.get_visible() else None
                for l in lines
            ]
        else:
            legend.remove()
            self.axes._ea_labels = [None] * len(lines)

    def set_xlim(self, xmin: float, xmax: float) -> None:
        if not self.axes:
            raise RuntimeError("You must call load() first")
        self.axes.set_xlim(xmin, xmax)
        self._modifications_made = True

    def set_ylim(self, ymin: float, ymax: float) -> None:
        if not self.axes:
            raise RuntimeError("You must call load() first")
        self.axes.set_ylim(ymin, ymax)
        self._modifications_made = True

    def save(self, output_path: Optional[Path] = None,
             format: str = 'pickle') -> Path:
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

    def has_modifications(self) -> bool:
        return self._modifications_made
