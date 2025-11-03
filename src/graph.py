import numpy as np
from types import UnionType
import os
import pickle as pl

import matplotlib.pyplot as plt
from scipy.optimize import minimize, differential_evolution

from src.constants import * 



def eV_to_nm(eV):
    return FACTOR_EV_NM / eV


class Graph:
    START = {"IR": 0.1, "VCD": 0.1, "UV": 0.1, "ECD": 0.1}  # cm-1  # cm-1  # eV  # eV
    END = {"IR": 6000, "VCD": 6000, "UV": 9, "ECD": 9}  # cm-1  # cm-1  # eV  # eV

    def __init__(self, graph_type, log=None, protocol=None, definition=4, norm=1, read_pop = None, **kwargs):
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
    DEFs = {
        "IR": 50,  # σ in cm-1
        "VCD": 50,  # σ in cm-1
        "UV": 0.25,  # FWHM in eV
        "ECD": 0.25,  # FWHM in eV
    }

    def __init__(self, conf, invert, convolution=None, shift=None, fwhm=None, read_pop=None, **kwargs):
        
        # kwargs.update({'shift': shift, 'fwmh':fwhm})
        super().__init__(**kwargs)

        self.invert = invert
        self.g = convolution if convolution else self.graph_type
        self.conformers = conf
        idx_first_valid_conf = [idx for idx, i in enumerate(conf) if i.active][0]
        if conf[idx_first_valid_conf].energies[self.protocol].get("graph", None) is None:
            return
        elif (
            conf[idx_first_valid_conf].energies[self.protocol].get("graph").get(self.graph_type, None)
            is None
        ):
            return

        self.shift = shift or None
        self.fwhm = fwhm or None

        self.retrive_data()

        self.auto = os.path.exists(
            os.path.join(os.getcwd(), f"{self.graph_type.lower()}_ref.dat")
        )
        
        if self.auto:
            kwargs.update({'log': self.log})
            self.ref = Experimental(**kwargs)
            self.autoconvolute()
        else:
            self.y = self.CONV[self.g](self.x, self.y_comp, self.DEFs[self.graph_type])
            if self.invert: 
                self.y *= -1
            self.y = self.normalize()
            # self.shift = 0
            # self.fwhm = self.DEFs[self.graph_type]

    def retrive_data(self):
        self.y_comp = np.zeros(self.x.shape)
        conv_func = self.CONV[self.graph_type]
        for i in self.conformers:
            if not i.active:
                continue
            x = np.array(i.energies[self.protocol]["graph"][self.graph_type]["x"])
            y = np.array(i.energies[self.protocol]["graph"][self.graph_type]["y"])
            
            pop = i.energies[self.protocol]["Pop"] if not self.read_pop else i.energies[self.read_pop]["Pop"]

            conv_graph = conv_func(x, y, self.DEFs[self.graph_type])

            with open(os.path.join(i.folder,f'p{self.protocol}_{self.graph_type}.xy'), "w") as f:
                for x, y in zip(self.x, self.y_comp):
                    f.write(f"{x:.6f} {y:.6f}\n")

            self.y_comp += conv_graph * float(pop)

    def autoconvolute(self):
        f = self.CONV[self.graph_type]

        def wrapper(variables):
            shift, fwhm = variables
            x = self.x + shift if self.graph_type in ["UV", "ECD"] else self.x * shift

            y = f(x, self.y_comp, fwhm)
            if self.invert: 
                y *= -1

            y_norm = y / np.max(np.abs(y[self.ref.x_min_idx : self.ref.x_max_idx]))
            ref_norm = self.ref.y / np.max(
                np.abs(self.ref.y[self.ref.x_min_idx : self.ref.x_max_idx])
            )
            d = self.diversity_function(
                y_norm[self.ref.x_min_idx : self.ref.x_max_idx],
                ref_norm[self.ref.x_min_idx : self.ref.x_max_idx],
            )
            print(f"Shift: {shift:.4f}, FWHM: {fwhm:.4f}, RMSD: {d:.6f}")
            return d  # Minimizzare questa funzione

        sb, fb = 0.2, 0.3
        ss, sf = 0, self.DEFs[self.graph_type]

        if hasattr(self, "shift"):
            ss = self.shift
        if hasattr(self, "fwhm"):
            sf = self.fwhm


        # Shift
        if type(self.shift) == list:
            shift_bounds = (self.shift[0], self.shift[1])
            ss = ss[0]
        elif type(self.shift) == float:
            shift_bounds = (self.shift, self.shift)
        elif self.graph_type in ["UV", "ECD"]:
            shift_bounds = (-sb + ss, sb + ss)
        elif self.graph_type in ["IR", "VCD"]:
            shift_bounds = (0.5, 1)

        # FWHM
        if type(self.fwhm) == list:
            fwhm_bounds = (self.fwhm[0], self.fwhm[1])
            sf = sf[0]
        elif type(self.fwhm) == float:
            fwhm_bounds = (self.fwhm, self.fwhm)
        else:
            fwhm_bounds = (0.2, sf + fb)

        initial_guess = [ss, sf]  # RED SHIFT NEGATIVE
        bounds = [shift_bounds, fwhm_bounds]

        result = minimize(
            wrapper,
            initial_guess,
            bounds=bounds,
            options={"maxiter": 1000},
            method="L-BFGS-B",
        )
        # result = minimize(wrapper, initial_guess, bounds=bounds, options={'maxiter': 1000}, method='Nelder-Mead')

        if result.success:
            self.shift, self.fwhm = result.x
            x = (
                self.x + self.shift
                if self.graph_type in ["UV", "ECD"]
                else self.x * self.shift
            )
            self.y = self.CONV[self.graph_type](x, self.y_comp, self.fwhm)
            self.y = self.normalize()
            self.log.info(
                f"{'='*20}\nResults for {self.graph_type} graph calcualted in Protocol {self.protocol}\n{'='*20}\n" +
                f"Optimal shift: {self.shift:.3f}, Optimal FWHM: {self.fwhm:.3f}, Similarity: {((1 if self.graph_type not in CHIRALS else 2)-result.fun)/(1 if self.graph_type not in CHIRALS else 2)*100:.2f}%\n"+
                "="*20+'\n'
            )
        else:
            self.log.error("Optimization failed. Convolution with defaul values")
            self.shift, self.fwhm = 0, self.DEFs[self.graph_type]
            x = self.x + self.shift if self.graph_type in ["UV", "ECD"] else self.x
            self.y = self.CONV[self.graph_type](
                self.x + self.shift, self.y_comp, self.fwhm
            )
            self.y = self.normalize()

    def diversity_function(self, a, b):
        return (np.sqrt(np.mean((a - b) ** 2))) / (1 if self.graph_type not in CHIRALS else 2) # RMSD


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

        self.convert = (
            len(self.data[self.data[:, 0] > self.DEFs[self.graph_type], :]) != 0
        )

        if self.convert:
            assert self.graph_type not in ["IR", "VCD"]
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

        self.comps = {
            int(os.path.basename(i).split("_")[1]): np.loadtxt(i)
            for i in os.listdir(os.getcwd())
            if self.graph_type.lower() in str(i).lower()
            and "comp" in str(i).lower()
            and "xy" in str(i).lower()
        }
        self.comps = dict(sorted(self.comps.items()))
        self.labels += [f"{protocol[i].functional}" for i in self.comps.keys()]

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

            m, M = -np.max(np.abs(self.ref.y)) if self.graph_type in [
                "ECD",
                "VCD",
            ] else 0, np.max(np.abs(self.ref.y))
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
        
        with open(f"{self.graph_type.lower()}_comparison.pickle", 'wb') as f:
            pl.dump(self.fig, f)
        
        plt.savefig(f"{self.graph_type.lower()}_comparison.png", dpi=300)

        

    def eV_to_nm(eV):
        return FACTOR_EV_NM / eV

    def nm_to_eV(nm):
        return FACTOR_EV_NM / nm


def main_graph(ensemble, p, log, invert, shift = None, fwhm = None, read_pop = None):
    for j in ["IR", "VCD", "UV", "ECD"]:
        graph = Computed(
            ensemble,
            invert=(j in CHIRALS) and invert,
            graph_type=j,
            protocol=str(p.number),
            shift=shift,
            fwhm=fwhm,
            log=log,
            read_pop = read_pop
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
    from launch import restart
    from pruning import calculate_rel_energies
    from logger import create_log

    ensemble, protocol, _ = restart()
    calculate_rel_energies(ensemble, 298.15)


    log = create_log("test.out")
    for i in protocol:
        main_graph(ensemble, i, log=log)
    for j in ["IR", "VCD", "UV", "ECD"]:
        g = Compared(protocol, graph_type=j, log=log)
        g.save_graph()
