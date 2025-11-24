import os
import numpy as np
import matplotlib.pyplot as plt
import pickle as pl
from typing import Optional, Union, List
from dataclasses import dataclass, field
from src.constants import eV_to_nm, CHIRALS, FACTOR_EV_NM

from src._spectral.graph_default import GraphDefault

@dataclass
class ComparedGraph:
    graph_type: str
    experimental_file: str = None
    log: any = None
    protocol_index: Optional[Union[List[str], None]] = None
    nm : bool = True

    def __post_init__(self):
        self.Xr, self.Yr, self.bounders = self._load_experimental()

        self.data = self._load_computed()

        self.defaults = GraphDefault(self.graph_type)
        
    def _load_computed(self):
        data = {}
        for f in os.listdir():
            if f.endswith('.xy') and f"{self.graph_type.upper()}_p" in f:
                proto = f.split("_p")[1].split("_")[0]

                if self.protocol_index is not None: 
                    if int(proto) not in self.protocol_index:
                        self.log.info(f'Protocol {proto} skipped')
                        continue
                    
                X, Y = np.loadtxt(f, unpack=True, dtype=np.float64)

                if isinstance(self.bounders, np.ndarray):
                    max_y = np.max(np.abs(Y[int(self.bounders[0]):int(self.bounders[1])]))
                    Y /= max_y

                data[proto] = (X, Y)

        if self.log:
            self.log.info(f"Loaded {len(data)} computed {self.graph_type} spectra.")
        return data

    def _load_experimental(self):
        if not self.experimental_file:
            return None, None, None
        X, Y = np.loadtxt(self.experimental_file, unpack=True)
        bounders = np.loadtxt(f'{self.graph_type.upper()}_index_lim.xy')

        return X, Y, bounders

    def plot(self, save=True, show=False):
        plt.style.use("seaborn-v0_8-paper")

        fig, ax = plt.subplots()

        if self.Xr is not None:
            ax.plot(self.Xr[int(self.bounders[0]):int(self.bounders[1])], self.Yr[int(self.bounders[0]):int(self.bounders[1])], color='black', lw=1.5, label='Experimental')
            ax.set_xlim(self.Xr[int(self.bounders[0])]-self.defaults.X_buffer,self.Xr[int(self.bounders[1])]+self.defaults.X_buffer)
        
        for proto, (X, Y) in self.data.items():
            if Y[~np.isnan(Y)].size > 0:
                ax.plot(X, Y, lw=1, label=f"Protocol {proto}", alpha=.75)

        ax.legend(fancybox=True, shadow=True)

        ax.set_xlabel(self.defaults.axis_label['x'])
        ax.set_ylabel(self.defaults.axis_label['y'])

        ax.grid(linestyle="-", linewidth=0.2, alpha=0.5)
        plt.gca().yaxis.grid(False)


        if self.graph_type in ["IR", "VCD"]:
            ax.xaxis.set_inverted(True)

        elif self.graph_type in ["UV", "ECD"]:
            secax = ax.secondary_xaxis("top", functions=(eV_to_nm, eV_to_nm))
            secax.set_xlabel(r"Wavelength $\lambda$ [$nm$]")

            if hasattr(self, "ref"):
                secax.set_xlim(eV_to_nm(self.ref.x_max-self.defaults.X_buffer), eV_to_nm(self.ref.x_min+self.defaults.X_buffer))

        if self.graph_type in CHIRALS:
            ax.set_ylim([-1.05,1.05])
        else:
            ax.set_ylim([-.05,1.05])

        ax.set_title(f'{self.graph_type.upper()} spectra comparison')
        
        plt.tight_layout()

        if save:
            fname = f"{self.graph_type.upper()}_comparison.png"
            plt.savefig(fname, dpi=300)
            if self.log:
                self.log.info(f"Saved {fname}")
        if show:
            plt.show()
        else:
            plt.close(fig)

        if self.graph_type in ["UV", "ECD"] and self.nm:  

            fig, ax = plt.subplots()


            ax.set_xlabel(r"Wavelength $\lambda$ [$nm$]")
            ax.set_ylabel(self.defaults.axis_label['y'])
            
            secax = ax.secondary_xaxis("top", functions=(eV_to_nm, eV_to_nm))
            secax.set_xlabel(self.defaults.axis_label['x'])
            
            for proto, (X, Y) in self.data.items():
                if Y[~np.isnan(Y)].size > 0:
                    ax.plot(FACTOR_EV_NM/X, Y, lw=1, label=f"Protocol {proto}", alpha=.75)

            # plot REF
            if self.Xr is not None:
                ax.plot(FACTOR_EV_NM/self.Xr[int(self.bounders[0]):int(self.bounders[1])], self.Yr[int(self.bounders[0]):int(self.bounders[1])], color='black', lw=1.5, label='Experimental')
                ax.set_xlim(FACTOR_EV_NM/(self.Xr[int(self.bounders[1])]+self.defaults.X_buffer), FACTOR_EV_NM/(self.Xr[int(self.bounders[0])]-self.defaults.X_buffer))

                secax.set_xlim(eV_to_nm(self.ref.x_max-self.defaults.X_buffer), eV_to_nm(self.ref.x_min+self.defaults.X_buffer))
    
            if self.graph_type in CHIRALS:
                ax.set_ylim([-1.05,1.05])
            else:
                ax.set_ylim([-.05,1.05])

            if save:
                fname = f"{self.graph_type.upper()}_comparison_nm.png"
                plt.savefig(fname, dpi=300)
                if self.log:
                    self.log.info(f"Saved {fname}")
            if show:
                plt.show()
            else:
                plt.close(fig)
