import pickle
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Optional, List, Tuple, Dict
from dataclasses import dataclass
from ensemble_analyzer.constants import eV_to_nm, CHIRALS, FACTOR_EV_NM
from ensemble_analyzer._spectral.graph_default import GraphDefault


@dataclass
class ComparedGraph:

    graph_type: str
    experimental_file: Optional[str] = None
    log: Optional[any] = None
    protocol_index: Optional[List[int]] = None
    nm: bool = True

    def __post_init__(self):
        self._validate_graph_type()
        self.Xr, self.Yr, self.bounders, self.weighted = self._load_experimental()
        self.data = self._load_computed()        
        self.defaults = GraphDefault(self.graph_type)
        
        
    def _validate_graph_type(self) -> None:

        valid_types = ["UV", "IR", "ECD", "VCD"]
        if self.graph_type.upper() not in valid_types:
            raise ValueError(f"Graph_Type ({self.graph_type.upper()}) not included: {valid_types}")
    
    def _load_computed(self) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
        data = {}
        pattern = f"{self.graph_type.upper()}_p"
        files = sorted(
            Path('.').glob('*.xy'),
            key=lambda p: int(self._extract_protocol_number(p.stem)) if pattern in p.name else float('inf')
        )

        for filepath in files:
            if pattern not in filepath.name:
                continue
                
            proto = self._extract_protocol_number(filepath.name)
            
            if not self._is_protocol_included(proto):
                if self.log:
                    self.log.debug(f'Protocol {proto} skipped')
                continue
            
            X, Y = np.loadtxt(filepath, unpack=True, dtype=np.float64)
            Y = self._normalize_spectrum(Y)
            data[proto] = (X, Y)

        if self.log:
            self.log.debug(f"Loaded {len(data)} computed {self.graph_type} spectra.")
        
        return data
    
    def _extract_protocol_number(self, filename: str) -> str:
        return filename.split("_p")[1].split("_")[0]
    
    def _is_protocol_included(self, proto: str) -> bool:
        if self.protocol_index is None:
            return True
        return int(proto) in self.protocol_index
    
    def _normalize_spectrum(self, Y: np.ndarray) -> np.ndarray:
        if not isinstance(self.bounders, np.ndarray):
            return Y
            
        idx_start, idx_end = int(self.bounders[0]), int(self.bounders[1])
        max_y = np.max(np.abs(Y[idx_start:idx_end]))
        
        if max_y > 1:
            return Y / max_y
        return Y

    def _load_experimental(self) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray]]:
        if not self.experimental_file:
            return None, None, None, None
        
        try:
            X, Y = np.loadtxt(self.experimental_file, unpack=True)
            bounders = np.loadtxt(f'{self.graph_type.upper()}_index_lim.xy')
            _, weighted = np.loadtxt(f'{self.graph_type}_weighted.xy', unpack=True)
            return X, Y, bounders, weighted
        except FileNotFoundError as e:
            if self.log:
                self.log.error(f"{self.log.FAIL} File not found: {e}")
            return None, None, None, None

    def plot(self, save: bool = True, show: bool = False, show_ref_weight: bool = False) -> None:
        self._plot_spectrum(save, show, in_nm=False, show_ref_weight=show_ref_weight)
        
        if self.graph_type.upper() in ["UV", "ECD"] and self.nm:
            self._plot_spectrum(save, show, in_nm=True, show_ref_weight=show_ref_weight)
    
    def _plot_spectrum(self, save: bool, show: bool, in_nm: bool = False, show_ref_weight:bool = False) -> None:

        plt.style.use("seaborn-v0_8-paper")
        fig, ax = plt.subplots()
        
        self._plot_computed_data(ax, in_nm)
        
        self._plot_experimental_data(ax, in_nm, show_ref_weight)
        
        self._configure_axes(ax, in_nm)
        self._configure_limits(ax, in_nm)
        
        ax.legend(fancybox=True, shadow=True)
        ax.set_title(f'{self.graph_type.upper()} spectra comparison')
        ax.grid(linestyle="-", linewidth=0.2, alpha=0.5)
        ax.yaxis.grid(False)
        
        plt.tight_layout()
        
        self._save_or_show(fig, save, show, in_nm)
    
    def _plot_computed_data(self, ax: plt.Axes, in_nm: bool) -> None:

        for proto, (X, Y) in self.data.items():
            if Y[~np.isnan(Y)].size > 0:
                x_values = FACTOR_EV_NM / X if in_nm else X
                ax.plot(x_values, Y, lw=1, label=f"Protocol {proto}", alpha=.75)
    
    def _plot_experimental_data(self, ax: plt.Axes, in_nm: bool, show_ref_weight:bool = False) -> None:

        if self.Xr is None or self.bounders is None:
            return
            
        idx_start, idx_end = int(self.bounders[0]), int(self.bounders[1])
        x_exp = self.Xr[idx_start:idx_end]
        y_exp = self.Yr[idx_start:idx_end]
        x_exp_weigh = self.Xr
        
        if in_nm:
            x_exp = FACTOR_EV_NM / x_exp
            x_exp_weigh = FACTOR_EV_NM / x_exp_weigh
            
        ax.plot(x_exp, y_exp, color='black', lw=1.5, label='Experimental')
        if show_ref_weight: 
            ax.plot(x_exp_weigh, self.weighted, color='black', lw=.4, label='Weighting function', alpha=0.5)

    
    def _configure_axes(self, ax: plt.Axes, in_nm: bool) -> None:

        if in_nm:
            ax.set_xlabel(r"Wavelength $\lambda$ [nm]")
            secax = ax.secondary_xaxis("top", functions=(eV_to_nm, eV_to_nm))
            secax.set_xlabel(self.defaults.axis_label['x'])
        else:
            ax.set_xlabel(self.defaults.axis_label['x'])
            if self.graph_type.upper() in ["UV", "ECD"]:
                secax = ax.secondary_xaxis("top", functions=(eV_to_nm, eV_to_nm))
                secax.set_xlabel(r"Wavelength $\lambda$ [nm]")
        
        ax.set_ylabel(self.defaults.axis_label['y'])
    
    def _configure_limits(self, ax: plt.Axes, in_nm: bool) -> None:

        y_lim = [-1.05, 1.05] if self.graph_type.upper() in CHIRALS else [-0.05, 1.05]
        ax.set_ylim(y_lim)
        
        if self.Xr is not None and self.bounders is not None:
            idx_start, idx_end = int(self.bounders[0]), int(self.bounders[1])
            buffer = self.defaults.X_buffer
            
            if in_nm:
                x_min = FACTOR_EV_NM / (self.Xr[idx_end] + buffer)
                x_max = FACTOR_EV_NM / (self.Xr[idx_start] - buffer)
            else:
                x_min = self.Xr[idx_start] - buffer
                x_max = self.Xr[idx_end] + buffer
            
            ax.set_xlim(x_min, x_max)
        
        if self.graph_type.upper() in ["IR", "VCD"] and not in_nm:
            ax.invert_xaxis()
    
    def _save_or_show(self, fig: plt.Figure, save: bool, show: bool, in_nm: bool) -> None:

        if save:
            suffix = "_nm" if in_nm else ""
            fname = f"{self.graph_type.upper()}_comparison{suffix}.png"
            pickle.dump(fig, open(f"{self.graph_type.upper()}_comparison{suffix}.pickle", 'wb'))
            plt.savefig(fname, dpi=300)
            if self.log:
                self.log.info(f"Saved {fname}")
        
        if show:
            plt.show()
        else:
            plt.close(fig)