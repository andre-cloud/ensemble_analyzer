import numpy as np
import os
import pickle as pl

import matplotlib.pyplot as plt
from scipy.optimize import minimize

try:
    from src.constants import *
except ModuleNotFoundError: 
    from constants import *


def eV_to_nm(eV):
    return FACTOR_EV_NM / eV


class Graph:
    START = {"IR": 0.1, "VCD": 0.1, "UV": 0.1, "ECD": 0.1}  # cm-1  # cm-1  # eV  # eV
    END = {"IR": 6000, "VCD": 6000, "UV": 9, "ECD": 9}  # cm-1  # cm-1  # eV  # eV

    def __init__(
        self,
        graph_type,
        log=None,
        protocol=None,
        definition=4,
        norm=1,
        read_pop=None,
        **kwargs,
    ):
        assert graph_type in list(self.START.keys())

        self.protocol = protocol
        self.read_pop = read_pop

        self.graph_type = graph_type
        self.x = np.linspace(
            self.START[graph_type], self.END[graph_type], num=10**definition
        )
        self.y = np.zeros(self.x.shape)
        self.norm = norm
        self.log = log

        self.CONV = {
            "IR": self.lorentzian,
            "VCD": self.lorentzian,
            "UV": self.gaussian,
            "ECD": self.gaussian,
            "lorentzian": self.lorentzian,
            "gaussian": self.gaussian,
        }

    def normalize(self):
        return self.y / np.max(np.abs(self.y)) * self.norm

    def gaussian(self, x0, I, fwhm):
        sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
        y = np.sum(
            I[:, np.newaxis]
            / (sigma * np.sqrt(2 * np.pi))
            * np.exp(-0.5 * ((self.x - x0[:, np.newaxis]) / sigma) ** 2),
            axis=0,
        )
        return y

    def lorentzian(self, x0, I, fwhm):
        y = np.sum(
            (I[:, np.newaxis] * fwhm**2)
            / (fwhm**2 + 4 * (self.x - x0[:, np.newaxis]) ** 2),
            axis=0,
        )
        return y

    def dump_graph(self, fname):
        with open(fname, "w") as f:
            for x, y in zip(self.x, self.y):
                f.write(f"{x:.6f} {y:.6f}\n")

class Computed(Graph):
    """
    Computed spectrum from conformers with optional convolution and (auto-)optimization.
    Handles shift and fwhm for UV/ECD as additive (shift in eV), and for IR/VCD as multiplicative (shift as scaling factor).
    Provides grid search + local optimization (L-BFGS-B) for best shift/fwhm against experimental data.
    """
    DEFs = {
        "IR": 50,    # σ in cm-1
        "VCD": 50,   # σ in cm-1
        "UV": 0.25,  # FWHM in eV
        "ECD": 0.25, # FWHM in eV
    }

    def __init__(
        self,
        conf,
        invert,
        convolution=None,
        shift=None,
        fwhm=None,
        read_pop=None,
        **kwargs,
    ):
        """
        Parameters
        ----------
        conf : list
            List of conformer objects.
        invert : bool
            Whether to invert the spectrum (for chiral).
        convolution : str or None
            Convolution type (defaults to graph_type).
        shift : None, float, or list
            For UV/ECD: additive shift in eV (None, float, or [min,max]).
            For IR/VCD: ignored.
        fwhm : None, float, or list
            FWHM for convolution (None, float, or [min,max]).
        read_pop : str or None
            Population protocol to read.
        kwargs : dict
            Passed to Graph and Experimental.
        """
        super().__init__(**kwargs)
        self.invert = invert
        self.g = convolution if convolution else self.graph_type
        self.conformers = conf
        self.read_pop = read_pop
        print = self.log

        # Validate at least one conformer is active and has required data
        idx_first_valid_conf = None
        for idx, i in enumerate(conf):
            if i.active:
                idx_first_valid_conf = idx
                break
        if idx_first_valid_conf is None:
            print("No active conformers found.")
            return
        conf0 = conf[idx_first_valid_conf]
        if conf0.energies[self.protocol].get("graph", None) is None:
            print("No 'graph' key in conformer energies.")
            return
        if conf0.energies[self.protocol]["graph"].get(self.graph_type, None) is None:
            print(f"No '{self.graph_type}' data in conformer energies.")
            return

        # Only allow shift/fwhm for UV/ECD as input
        if self.graph_type in ["UV", "ECD"]:
            self.shift = shift if (shift is None or isinstance(shift, (float, list))) else None
            self.fwhm = fwhm if (fwhm is None or isinstance(fwhm, (float, list))) else None
        else:
            self.shift = None
            self.fwhm = fwhm if (fwhm is None or isinstance(fwhm, (float, list))) else None

        # Retrieve conformer data (no convolution, no file output)
        self.retrieve_data()

        # Check if experimental data exists for autoconvolution
        self.auto = os.path.exists(
            os.path.join(os.getcwd(), f"{self.graph_type.lower()}_ref.dat")
        )
        if not np.all(self.y_comp == 0):
            if (isinstance(self.fwhm, float) or isinstance(self.fwhm, int)) and (isinstance(self.shift, float) or isinstance(self.shift, int)):
                self.y = self.compute_convolution(
                    shift=self.shift,
                    fwhm=self.fwhm
                )
                
                if self.invert:
                    self.y *= -1
                self.y = self.normalize()
            elif self.auto:
                kwargs = dict(kwargs)
                kwargs.update({"log": self.log})
                self.ref = Experimental(**kwargs)
                self.autoconvolute()
            else:
                # Use default convolution
                self.y = self.compute_convolution(
                    shift=None,
                    fwhm=self.fwhm if self.fwhm is not None else self.DEFs[self.graph_type]
                )
                if self.invert:
                    self.y *= -1
                self.y = self.normalize()

    def retrieve_data(self):
        """
        Collects and sums the raw data from all active conformers, weighted by population.
        No convolution or file output is performed here.
        """
        self.y_comp = np.zeros(self.x.shape)
        for i in self.conformers:
            if not i.active:
                continue
            try:
                xvals = np.array(i.energies[self.protocol]["graph"][self.graph_type]["x"])
                yvals = np.array(i.energies[self.protocol]["graph"][self.graph_type]["y"])
                pop = (
                    i.energies[self.protocol]["Pop"]
                    if not self.read_pop
                    else i.energies[self.read_pop]["Pop"]
                )
            except Exception as e:
                print(f"Skipping conformer due to error: {e}")
                continue
            # Interpolate conformer contribution to self.x grid
            y_interp = np.interp(self.x, xvals, yvals, left=0, right=0)
            self.y_comp += y_interp * float(pop)
        print(f"Retrieved and summed conformer data for {self.graph_type}.")

    def compute_convolution(self, shift=None, fwhm=None):
        """
        Apply convolution to the composite spectrum, with shift/fwhm.
        For UV/ECD: x + shift (additive).
        For IR/VCD: x * shift (multiplicative), but shift is always 1 unless optimized.
        """
        conv_func = self.CONV[self.g]

        if fwhm:
            if isinstance(fwhm, list):
                fwhm_val = fwhm[0]
            else:
                fwhm_val = fwhm
        else: 
            fwhm_val = self.DEFs[self.graph_type]

        # For shift: only relevant for UV/ECD, else use 1.0
        if self.graph_type in ["UV", "ECD"]:
            shift_val = shift if shift is not None else 0.0
            x_conv = self.x + shift_val
        else:
            shift_val = shift if shift is not None else 1.0
            x_conv = self.x * shift_val
        y = conv_func(x_conv, self.y_comp, fwhm_val)
        return y

    def autoconvolute(self):
        """
        Optimize shift/fwhm to best match experimental spectrum.
        Uses grid search + local optimization (L-BFGS-B).
        """
        print(f"Starting autoconvolution for {self.graph_type}.")
        conv_func = self.CONV[self.graph_type]
        # Define bounds and initial guess for shift/fwhm
        # UV/ECD: shift is additive, IR/VCD: shift is multiplicative
        if self.graph_type in ["UV", "ECD"]:
            # Shift: None/float/list
            if isinstance(self.shift, list) and len(self.shift) == 2:
                shift_bounds = (self.shift[0], self.shift[1])
                shift_init = np.mean(self.shift)
            elif isinstance(self.shift, float):
                shift_bounds = (self.shift, self.shift)
                shift_init = self.shift
            else:
                shift_bounds = (-0.2, 0.2)
                shift_init = 0.0
        else:
            # For IR/VCD, shift is scaling factor (multiplicative)
            shift_bounds = (0.95, 1.05)
            shift_init = 1.0

        if isinstance(self.fwhm, list) and len(self.fwhm) == 2:
            fwhm_bounds = (self.fwhm[0], self.fwhm[1])
            fwhm_init = np.mean(self.fwhm)
        elif isinstance(self.fwhm, float):
            fwhm_bounds = (self.fwhm, self.fwhm)
            fwhm_init = self.fwhm
        else:
            fwhm_bounds = (0.2, self.DEFs[self.graph_type] * 2)
            fwhm_init = self.DEFs[self.graph_type]

        # Grid search (coarse) + local optimization
        best_rmsd = np.inf
        best_params = (shift_init, fwhm_init)
        grid_n = 5
        shift_grid = np.linspace(shift_bounds[0], shift_bounds[1], grid_n)
        fwhm_grid = np.linspace(fwhm_bounds[0], fwhm_bounds[1], grid_n)
        for s in shift_grid:
            for f in fwhm_grid:
                y = self.compute_convolution(shift=s, fwhm=f)
                if self.invert:
                    y = -y
                # Normalize in the experimental window
                y_norm = y / np.max(np.abs(y[self.ref.x_min_idx:self.ref.x_max_idx]))
                ref_norm = self.ref.y / np.max(np.abs(self.ref.y[self.ref.x_min_idx:self.ref.x_max_idx]))
                rmsd = self.diversity_function(
                    y_norm[self.ref.x_min_idx:self.ref.x_max_idx],
                    ref_norm[self.ref.x_min_idx:self.ref.x_max_idx],
                )
                if rmsd < best_rmsd:
                    best_rmsd = rmsd
                    best_params = (s, f)
        # Now local optimization from best grid point
        def opt_fun(params):
            s, f = params
            y = self.compute_convolution(shift=s, fwhm=f)
            if self.invert:
                y = -y
            y_norm = y / np.max(np.abs(y[self.ref.x_min_idx:self.ref.x_max_idx]))
            ref_norm = self.ref.y / np.max(np.abs(self.ref.y[self.ref.x_min_idx:self.ref.x_max_idx]))
            rmsd = self.diversity_function(
                y_norm[self.ref.x_min_idx:self.ref.x_max_idx],
                ref_norm[self.ref.x_min_idx:self.ref.x_max_idx],
            )
            print(f"Shift: {s:.4f}, FWHM: {f:.4f}, RMSD: {rmsd:.6f}")
            return rmsd
        res = minimize(
            opt_fun,
            best_params,
            bounds=[shift_bounds, fwhm_bounds],
            options={"maxiter": 500},
            method="L-BFGS-B",
        )
        if res.success:
            opt_shift, opt_fwhm = res.x
            self.shift = opt_shift
            self.fwhm = opt_fwhm
            self.y = self.compute_convolution(shift=opt_shift, fwhm=opt_fwhm)
            if self.invert:
                self.y *= -1
            self.y = self.normalize()
            similarity = ((1 if self.graph_type not in CHIRALS else 2) - res.fun) / (1 if self.graph_type not in CHIRALS else 2) * 100
            print(
                f"{'='*20}\nResults for {self.graph_type} graph calculated in Protocol {self.protocol}\n{'='*20}\n"
                + f"Optimal shift: {self.shift:.4f}, Optimal FWHM: {self.fwhm:.4f}, Similarity: {similarity:.2f}%\n"
                + "="*20
            )
        else:
            print("Optimization failed. Using default convolution.")
            self.shift = 0.0 if self.graph_type in ["UV", "ECD"] else 1.0
            self.fwhm = self.DEFs[self.graph_type]
            self.y = self.compute_convolution(shift=self.shift, fwhm=self.fwhm)
            if self.invert:
                self.y *= -1
            self.y = self.normalize()

    def diversity_function(self, a, b):
        """
        Compute normalized RMSD between two spectra arrays.
        """
        return (np.sqrt(np.mean((a - b) ** 2))) / (
            1 if self.graph_type not in CHIRALS else 2
        )


class Experimental(Graph):
    DEFs = {
        "IR": 3,
        "VCD": 3,
        "UV": 100,
        "ECD": 100,
    }

    FACTOR = {
        "IR": FACTOR_EV_CM_1,
        "VCD": FACTOR_EV_CM_1,
        "UV": FACTOR_EV_NM,
        "ECD": FACTOR_EV_NM,
    }

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # Read the reference data

        self.filename = os.path.join(os.getcwd(), f"{self.graph_type.lower()}_ref.dat")

        self.log.debug("Reading the reference data from " + self.filename)

        self.data = np.loadtxt(self.filename)

        self.retrive_data()

        self.y = self.normalize()
        self.dump_graph(f"{self.filename[:-4]}_norm.xy")
        self.log.debug(
            "Reference data normalized and saved in " + f"{self.filename[:-4]}_norm.xy"
        )

    def retrive_data(self):
        # Sort by the Xs
        self.data = self.data[np.argsort(self.data[:, 0])]

        self.x_min = float(min(self.data[:, 0]))
        self.x_max = float(max(self.data[:, 0]))
        self.x_max_idx = np.argmax(self.data[:, 0])

        self.convert = (self.graph_type in ["UV", "ECD"] and np.mean(self.data[:,0]) > 20)

        if self.convert:
            if self.graph_type not in ["IR", "VCD"]:
                self.x_exp = self.FACTOR[self.graph_type] / self.data[:, 0]
                self.x_max, self.x_min = (
                    self.FACTOR[self.graph_type] / self.x_min,
                    self.FACTOR[self.graph_type] / self.x_max,
                )

        else:
            self.x_exp = self.data[:, 0]

        self.reverse = self.x_exp[0] > self.x_exp[-1]

        self.y_exp = self.data[:, 1]
        if self.graph_type in ["UV", "IR"]:
            self.y_exp -= np.min(self.y_exp)

        if self.reverse:
            self.x_exp = self.x_exp[::-1]
            self.y_exp = self.y_exp[::-1]

        self.y = self.interpolate()
        tmp = (self.x > self.x_min) & (self.x < self.x_max)

        self.x_min_idx = np.where(self.x == np.min(self.x[tmp]))[0][0]
        self.x_max_idx = np.where(self.x == np.max(self.x[tmp]))[0][0]

    def interpolate(self):
        return np.interp(self.x, self.x_exp, self.y_exp, left=0, right=0)


class Compared(Graph):
    BUFFER = {"IR": 0.1, "VCD": 0.1, "UV": 0.5, "ECD": 0.5}

    axis_label = {
        "IR": {
            "x": r"Wavenumber $\widetilde \nu$ [$cm^{-1}$]",
            "y": r"Intensity [$a.u.$]",
        },
        "VCD": {
            "x": r"Wavenumber $\widetilde \nu$ [$cm^{-1}$]",
            "y": r"Intensity [$a.u.$]",
        },
        "UV": {"x": r"Energy [$eV$]", "y": r"Intensity [$a.u.$]"},
        "ECD": {"x": r"Energy [$eV$]", "y": r"Intensity [$a.u.$]"},
    }

    def __init__(self, protocol, **kwargs):
        super().__init__(**kwargs)

        self.labels = []
        if os.path.exists(
            os.path.join(os.getcwd(), f"{self.graph_type.lower()}_ref.dat")
        ):
            self.ref = Experimental(**kwargs)

        comps = {
            int(os.path.basename(i).split("_")[1]): np.loadtxt(i)
            for i in os.listdir(os.getcwd())
            if self.graph_type.lower() in str(i).lower()
            and "comp" in str(i).lower()
            and "xy" in str(i).lower()
        }
        self.comps = {}
        
        for i in comps:
            if np.any(np.abs(comps[i][:, 1]) > 0):
                self.comps[i] = comps[i].copy()


        self.comps = dict(sorted(self.comps.items()))
        self.labels += [f"Protocol {i} @ {protocol[i].functional}" for i in self.comps.keys()]

        self.title = f"{self.graph_type.upper()} graph comparison"

    def save_graph(self):
        plt.style.use("seaborn-v0_8-paper")
        self.fig, self.axis = plt.subplots()

        if hasattr(self, "ref"):
            self.axis.plot(
                self.ref.x[self.ref.x_min_idx : self.ref.x_max_idx],
                self.ref.y[self.ref.x_min_idx : self.ref.x_max_idx],
                label="Experimental",
                lw=2,
            )
            self.axis.set_xlim(
                self.ref.x_min - self.BUFFER[self.graph_type],
                self.ref.x_max + self.BUFFER[self.graph_type],
            )

            m, M = (
                -np.max(np.abs(self.ref.y))
                if self.graph_type
                in [
                    "ECD",
                    "VCD",
                ]
                else 0
            ), np.max(np.abs(self.ref.y))
            self.axis.vlines(
                self.ref.x[self.ref.x_min_idx], m, M, color="grey", ls="dashed", lw=1
            )
            self.axis.vlines(
                self.ref.x[self.ref.x_max_idx], m, M, color="grey", ls="dashed", lw=1
            )

        for idx, (i, j) in enumerate(list(self.comps.items())):
            self.axis.plot(j[:, 0], j[:, 1], label=self.labels[idx], lw=1)

        self.axis.set_xlabel(self.axis_label[self.graph_type]["x"])
        self.axis.set_ylabel(self.axis_label[self.graph_type]["y"])

        # creating secondary x-xis for nm only for UV and ECD
        if self.graph_type in ["UV", "ECD"]:
            secax = self.axis.secondary_xaxis("top", functions=(eV_to_nm, eV_to_nm))
            secax.set_xlabel(r"Wavelength $\lambda$ [$nm$]")

            if hasattr(self, "ref"):
                secax.set_xlim(eV_to_nm(self.ref.x_max), eV_to_nm(self.ref.x_min))

        if self.graph_type in ["IR", "VCD"]:
            self.axis.xaxis.set_inverted(True)

        self.axis.grid(linestyle="-", linewidth=0.2, alpha=0.5)
        plt.gca().yaxis.grid(False)

        self.axis.legend(ncols=2, shadow=True, fancybox=True)
        plt.title(self.title)

        plt.tight_layout()

        with open(f"{self.graph_type.lower()}_comparison.pickle", "wb") as f:
            pl.dump(self.fig, f)

        plt.savefig(f"{self.graph_type.lower()}_comparison.png", dpi=300)

    def eV_to_nm(eV):
        return FACTOR_EV_NM / eV

    def nm_to_eV(nm):
        return FACTOR_EV_NM / nm


def main_graph(ensemble, p, log, invert, shift=None, fwhm=None, read_pop=None):
    for j in ["IR", "VCD", "UV", "ECD"]:
        graph = Computed(
            ensemble,
            invert=(j in CHIRALS) and invert,
            graph_type=j,
            protocol=str(p.number),
            shift=shift,
            fwhm=fwhm,
            log=log,
            read_pop=read_pop,
        )

        if not hasattr(graph, "y_comp"):
            continue

        graph.dump_graph(f"protocol_{p.number}_{j.lower()}_comp.xy")
        shift, fwhm = graph.shift, graph.fwhm
        plt.plot(graph.x, graph.y)
        if hasattr(graph, "ref"):
            plt.plot(
                graph.ref.x[graph.ref.x_min_idx : graph.ref.x_max_idx],
                graph.ref.y[graph.ref.x_min_idx : graph.ref.x_max_idx],
            )
            plt.xlim(graph.ref.x_min - 0.5, graph.ref.x_max + 0.5)

        plt.title(f"{j} of protocol {p.number}")
        plt.show()


if __name__ == "__main__":
    # from launch import restart
    # from pruning import calculate_rel_energies
    # from logger import create_log

    # ensemble, protocol, _ = restart()
    # calculate_rel_energies(ensemble, 298.15)

    # log = create_log("test.out")
    # for i in protocol:
    #     main_graph(ensemble, i, log=log)
    # for j in ["IR", "VCD", "UV", "ECD"]:
    #     g = Compared(protocol, graph_type=j, log=log)
    #     g.save_graph()

    import mock 
    Experimental(graph_type='UV', log=mock.MagicMock(), convert=True )