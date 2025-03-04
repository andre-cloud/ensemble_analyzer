import numpy as np

import os
import scipy.optimize as opt
from scipy.signal import argrelextrema
from scipy.constants import c, h, electron_volt, R
import matplotlib.pyplot as plt
from typing import Union

plt.set_loglevel("error")

try:
    from src.regex_parsing import regex_parsing
except ImportError as e:  # pragma: no cover
    print(e)
    from regex_parsing import regex_parsing

FACTOR_EV_NM = h * c / (10**-9 * electron_volt)


class Graph:
    """
    Generate electronic graphs (UV and ECD) and, if present "uv_ref.dat" or "ecd_ref.dat" file as mere xy file, it will automatically align the calculated graph with the reference.
    """

    def __init__(
        self,
        confs,
        protocol,
        log,
        T,
        final_lambda=800.0,
        definition=4,
        FWHM: Union[None, float] = None,
        shift: Union[None, float] = None,
        invert: bool = False,
        regraph: bool = False,
    ):
        self.confs = [i for i in confs if i.active]
        self.protocol = protocol
        self.log = log
        self.pop = self.calc_pop(T, regraph=regraph)

        # self.log.debug(self.pop)

        self.x = np.linspace(
            FACTOR_EV_NM / (final_lambda), FACTOR_EV_NM / 150, 10**definition
        )  # eV x axis

        self.spectra = []
        self.filter_outputs()

        self.uv_impulses = np.array(
            [self.get_uv(i, pop) for i, pop in zip(self.spectra, self.pop)]
        )
        self.ecd_impulses = np.array(
            [self.get_ecd(i, pop) for i, pop in zip(self.spectra, self.pop)]
        )

        shift_from_uv = None

        # creating the ECD and UV graph. If uv_ref.dat and/or ecd_ref.dat auto-convolution of the reference is performed
        if os.path.exists(os.path.join(os.getcwd(), "uv_ref.dat")):
            uv, _ = self.auto_convolution(
                os.path.join(os.getcwd(), "uv_ref.dat"),
                fname_ref_damp_norm=os.path.join(os.getcwd(), "uv_ref_norm_eV.dat"),
                impulses=self.uv_impulses,
                fname=f"uv_protocol_{self.protocol.number}_auto_conv.dat",
                user_sigma=FWHM / (2 * np.sqrt(2 * np.log(2))) if FWHM else None,
                user_shift=shift,
                title="UV autoconvolution",
            )
        else:
            uv = self.calc_graph(
                impulses=self.uv_impulses,
                sigma=FWHM / (2 * np.sqrt(2 * np.log(2))) if FWHM else 1 / 3,
                fname=f"uv_protocol_{self.protocol.number}.dat",
                save=True,
            )

        if os.path.exists(os.path.join(os.getcwd(), "ecd_ref.dat")):
            ecd, _ = self.auto_convolution(
                os.path.join(os.getcwd(), "ecd_ref.dat"),
                fname_ref_damp_norm=os.path.join(os.getcwd(), "ecd_ref_norm_eV.dat"),
                impulses=self.ecd_impulses,
                fname=f"ecd_protocol_{self.protocol.number}_auto_conv.dat",
                user_sigma=FWHM / (2 * np.sqrt(2 * np.log(2))) if FWHM else None,
                user_shift=shift,
                invert=invert,
                is_ecd=True,
                title="ECD autoconvolution",
            )
        else:
            ecd = self.calc_graph(
                impulses=self.ecd_impulses,
                sigma=FWHM / (2 * np.sqrt(2 * np.log(2))) if FWHM else 1 / 3,
                fname=f"ecd_protocol_{self.protocol.number}.dat",
                save=True,
            )

        Graph.damp_graph(f"ecd_protocol_{self.protocol.number}.dat", self.x, ecd)
        Graph.damp_graph(f"uv_protocol_{self.protocol.number}.dat", self.x, uv)

    def filter_outputs(self) -> None:
        """
        Read once all conformers output, keeping only graph information

        :rtype: None
        """

        st, en = (
            regex_parsing[self.protocol.calculator]["start_spec"],
            regex_parsing[self.protocol.calculator]["end_spec"],
        )

        for i in self.confs:
            with open(
                os.path.join(
                    os.getcwd(), i.folder, f"protocol_{self.protocol.number}.out"
                )
            ) as f:
                fl = f.read()

            sp = fl.split(st)[-1].split(en)[0]
            self.spectra.append(sp)

        return None

    def get_uv(self, spectra, pop):
        """
        Get the impulses for the UV spectra calculation

        :param spectra: parsed output file
        :type spectra: str
        :param pop: bolzmann population
        :type pop: float

        :return: energy and impulse tuple
        :rtype: tuple(float, float)
        """
        graph = spectra.split(regex_parsing[self.protocol.calculator]["s_UV"])[
            -1
        ].split(regex_parsing[self.protocol.calculator]["break"])[0]

        return [
            (
                FACTOR_EV_NM
                / float(
                    i.strip().split()[
                        regex_parsing[self.protocol.calculator]["idx_en_tddft"]
                    ]
                ),
                float(
                    i.strip().split()[
                        regex_parsing[self.protocol.calculator]["idx_imp_tddft"]
                    ]
                )
                * pop,
            )
            for i in graph.splitlines()
            if i
        ]

    def get_ecd(self, spectra, pop):
        """
        Get the impulses for the ECD spectra calculation

        :param spectra: parse output file
        :type spectra: str
        :param pop: bolzmann population
        :type pop: float


        :return: energy and impulse tuple
        :rtype: tuple(float, float)
        """
        graph = spectra.split(regex_parsing[self.protocol.calculator]["s_ECD"])[
            -1
        ].split(regex_parsing[self.protocol.calculator]["break"])[0]

        return [
            (
                FACTOR_EV_NM
                / float(
                    i.strip().split()[
                        regex_parsing[self.protocol.calculator]["idx_en_tddft"]
                    ]
                ),
                float(
                    i.strip().split()[
                        regex_parsing[self.protocol.calculator]["idx_imp_tddft"]
                    ]
                )
                * pop,
            )
            for i in graph.splitlines()
            if i
        ]

    @staticmethod
    def gaussian(x, ev, I, sigma) -> np.ndarray:
        r"""
        Create a gaussian convolution for each impulse

        :param ev: energy of the impulse
        :type ev: float
        :param I: intensity of the impulse (Fosc for UV, R(vel) for ECD)
        :type I: float
        :param sigma: sigma of the gaussian distribution
        :type sigma: float

        .. math::
            f(x) = \frac {I}{σ*\sqrt{2π}} e^{-\frac 12(\frac {x-eV}{σ})^2}

        :return: Gaussian of the impulse
        :rtype: np.array
        """
        return I / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x - ev) / sigma) ** 2)

    def calc_pop(self, T, regraph) -> np.ndarray:
        """
        Boltzmann population, if necessary with correction

        :param T: temperature [K]
        :type T: float
        :param regraph: the routine ragraph is running
        :type regraph: bool

        :return: population distribution
        :rtype: np.array (1D array)
        """

        if not regraph:
            n, n_1 = self.protocol.number, str(int(self.protocol.number) - 1)

            # CONF do have frequency calculation before
            if self.protocol.freq:
                pass
            elif self.confs[0].energies[n_1]["G"]:
                # So energy is corrected: if functionals are the same, nothing change; else energy of the new function is corrected with lower frequency correction
                for i in self.confs:
                    i.energies[n]["G"] = i.energies[n]["E"] + (
                        i.energies[n_1]["G"] - i.energies[n_1]["E"]
                    )
            else:
                raise IOError("No frequency calculation")

            ens = np.array([i.get_energy for i in self.confs])
            ens_rel = ens - min(ens)
            bolz = np.exp((-ens_rel * 4186) / (R * T))
            pop = bolz / np.sum(bolz)
            for idx, i in enumerate(list(ens)):
                self.confs[idx]._last_energy["G"] = ens[idx]
                self.confs[idx]._last_energy["Erel"] = i
                self.confs[idx]._last_energy["Pop"] = pop[idx] * 100

        else:
            if self.confs[0].energies.get("Pop"):
                pop = np.array(
                    [i.energies[self.protocol.number]["Pop"] for i in self.confs]
                )
            else:
                pop = self.calc_pop(T, False)

        return pop

    def auto_convolution(
        self,
        fname_ref,
        fname_ref_damp_norm,
        impulses,
        fname,
        title,
        norm=1,
        user_sigma=None,
        user_shift=None,
        invert=False,
        is_ecd=False,
    ) -> np.ndarray:
        """
        Optimization to find the best fitting values for the Gaussian convolution.
        Optimization "Fitness Function" is the sum of the absolute value of the differences between the computed and the experimental graph that lay above the threshold.

        :param fname_ref: filename of the reference file
        :type fname_ref: str
        :param fname_ref_damp_norm: filename to save the autoconvoluted graph
        :type fname_ref_damp_norm: str
        :param impulses: list of single excitation [eV, I] where I can be a UV or ECD excitation
        :type impulses: np.array((eV, I))
        :param fname: filename to save the final convoluted graph
        :type fname: str
        :param user_sigma: user defined sigma for gaussian convolution
        :type user_sigma: float
        :param user_shift: user defined shift after gaussian convolution
        :type user_shift: float
        :param invert: invert the ECD graph when True
        :type invert: bool
        :param is_ecd: is convoluting an ecd graph
        :type is_ecd: bool
        :param title: title of the section for logging purposes
        :type title: str

        :return: normalized graph
        :rtype: np.array (1D array)
        """

        X = self.x.copy()
        ref = Ref_graph(
            fname_ref,
            None,
            invert=invert,
            X=X,
            fname_ref_damp_norm=fname_ref_damp_norm,
            norm=norm,
        )

        x_min_index, x_max_index = ref.x_min_index, ref.x_max_index
        if x_min_index > x_max_index:
            x_max_index, x_min_index = x_min_index, x_max_index

        m = np.max(np.max(impulses, axis=2))
        cut = np.where(impulses[:, :, 0] == m)
        max_impulse_x = impulses[cut[0], cut[1], 0]
        user_shift = list(user_shift) if user_shift is not None else user_shift
        user_sigma = list(user_sigma) if user_sigma is not None else user_sigma

        # dx = X[0]-X[1]

        def optimiser(variables):
            """
            Callback for the scipy.optimize.minimize
            """
            (sigma, shift) = variables

            Y_comp = Graph.normalise(
                self.calc_graph(
                    impulses=impulses, shift=shift, sigma=sigma, save=False
                ),
                norm=norm,
                x_min=x_min_index,
                x_max=x_max_index,
            )

            diversity = self.diversity_function(
                Y_comp[x_min_index:x_max_index], ref.y[x_min_index:x_max_index]
            )

            return diversity

        in_delta = 1

        initial_guess = [0.1415, in_delta]
        if user_shift is not None:
            if len(user_shift) > 1:
                shift = user_shift
            else:
                shift = user_shift * 2
        else:
            shift = (-1.5, 1.5)

        if user_sigma is not None:
            if len(user_sigma) > 1:
                sigma = user_sigma
            else:
                sigma = user_sigma * 2
        else:
            sigma = (0.08, 0.27)  # corresponding to FWMH = (0.19, 0.64) eV

        default_guess = [0.1415, 0]  # the σ correspond to a FWHM of 0.33 eV
        result = opt.minimize(
            optimiser,
            initial_guess,
            bounds=[sigma, shift],  # FWHM is between 0.2 and 0.5 eV
            options={"maxiter": 10000},
            method="Powell",
        )

        if is_ecd:
            result.fun = result.fun / 2

        if result.success:
            sigma, shift = result.x
            self.log.info(
                title + "\n"
                f"Convergence of parameters succeeded in {result.nfev} steps.\n"
                f"Confidence level: {(1-result.fun)*100:.2f}%. Parameters obtained\n"
                f"\t- σ = {sigma:.4f} eV (that correspond to a FWHM = {(sigma*np.sqrt(2*np.log(2))*2):.4f} eV)\n"
                f"\t- Δ = {shift:.4f} eV (in this case, a negative shift corresponds to a RED-shift)"
            )
            Y_COMP = Graph.normalise(
                self.calc_graph(
                    impulses=impulses, shift=shift, sigma=sigma, save=True, fname=fname
                ),
                norm=norm,
                x_min=x_min_index,
                x_max=x_max_index,
            )

        else:
            self.log.info(
                title + "\n"
                f"Convergence of parameters NOT succeeded.\n"
                "Parameters used to convolute the saved graph\n"
                f"\t- σ = {default_guess[0]:.4f} eV (that correspond to a FWHM = {(default_guess[0]*np.sqrt(2*np.log(2))*2):.4f} eV)\n"
                f"\t- Δ = 0.0000 eV"
            )
            Y_COMP = Graph.normalise(
                self.calc_graph(
                    impulses=impulses, shift=0, sigma=1 / 3, save=True, fname=fname
                ),
                norm=norm,
                x_min=x_min_index,
                x_max=x_max_index,
            )

        return Y_COMP, shift

    def diversity_function(self, comp, ref, **kwargs) -> float:
        dx = kwargs.get("dx", 1e-5)
        if set(ref == 0) != set([False, True]):
            # mape
            # diff = ref - comp
            # diff_p = np.abs(diff)/ref
            # diversity = np.mean(diff_p)/2

            # # RMSD
            diversity = np.sqrt(np.mean((comp - ref) ** 2))

            # integrals
            # i_ref = np.trapz(np.flip(ref)**2, dx=dx)
            # i_comp = np.trapz(np.flip(comp)**2, dx=dx)
            # diversity = np.sqrt(np.abs(i_comp-i_ref))

        else:
            diversity = 1

        return diversity

    def calc_graph(self, impulses, sigma, shift=0, fname="", save=False) -> np.ndarray:
        """
        Build the Spectra

        :param impulses: list of the (eV, I)s for each root of each conformer
        :type impulses: list
        :param sigma: dispersion for the gaussian convolution [eV]
        :type sigma: float
        :param fname: the name of the file to store the graph
        :type fname: str
        :param save: save the .dat file for each conformer
        :type save: bool

        :return: Convoluted and weighted spectra
        :rtype: np.array (1D array)
        """

        x = self.x.copy()
        y = np.zeros(x.shape)

        for idx in range(len(self.confs)):
            y_ = np.zeros(x.shape)
            for ev, I in impulses[idx]:
                y_ += Graph.gaussian(x + shift, ev, I, sigma)

            if save:
                Graph.damp_graph(
                    fname=os.path.join(os.getcwd(), self.confs[idx].folder, fname),
                    x=x,
                    y=y_,
                )

            y += y_
        return y

    @staticmethod
    def damp_graph(fname: str, x: np.ndarray, y: np.ndarray) -> None:
        """
        Damp an xy graph into file

        :param fname: filename to store the graph
        :type fname: str
        :param x: x axis
        :type x: np.array (1D array)
        :param y: y axis
        :type y: np.array (1D array)

        :return: None
        """
        data = np.array([(xi, yi) for xi, yi in zip(x, y)])
        np.savetxt(fname, data)
        return

    @staticmethod
    def load_graph(fname, is_ev=False):
        """
        Load an already saved graph from a file XY

        :param fname: filename
        :type fname: str
        :param is_ev: if the Y is already converted in eV, defaults to False
        :type is_ev: bool, optional
        :return: the graph it self divided in X and Y
        :rtype: tuple (X, Y)
        """
        arr = np.loadtxt(fname)
        if not is_ev:
            return FACTOR_EV_NM / arr[:, 0], arr[:, 1]
        return arr[:, 0], arr[:, 1]

    @staticmethod
    def normalise(y: np.ndarray, norm=1, x_min=0, x_max=-1) -> np.array:
        """
        Normalize an ensemble between 1 and -1, if not set otherwise.

        :param y: 1D array
        :type y: np.array
        :param norm: max value to normalize at
        :type norm: float
        :param x_min: index of the minimum value of the experimental value
        :type x_min: int
        :param x_max: index of the maximum value of the experimental value
        :type x_max: int


        :return: 1D normalized array
        :rtype: np.array
        """
        if x_max != -1:
            if x_min > x_max:
                x_min, x_max = x_max, x_min

        min_neg = np.min(y[x_min:x_max]) < 0
        y_max = np.max(
            [np.max(y[x_min:x_max]), np.min(y[x_min:x_max]) * (-1 if min_neg else 1)]
        )

        if (
            abs(y_max / np.max([np.max(y), np.min(y) * (-1 if np.min(y) < 0 else 1)]))
            < 1e-20
        ):
            y_max = np.max([np.max(y), np.min(y) * (-1 if np.min(y) < 0 else 1)])
        y_out = (y / y_max) * norm
        return y_out


class Ref_graph:
    """
    Load the reference graph in order to shift and convolute properly the calculated one.
    """

    def __init__(
        self,
        fname: str,
        log,
        X,
        fname_ref_damp_norm: str,
        norm: int,
        is_ev: bool = False,
        invert: bool = False,
    ):
        self.norm = norm

        data = np.loadtxt(fname, dtype=float)
        self.data = data[np.argsort(data[:, 0])]
        self.x_min = float(min(data[:, 0]))
        self.x_max = float(max(data[:, 0]))

        self.x = data[:, 0] if is_ev else FACTOR_EV_NM / data[:, 0]
        self.y = data[:, 1] * (1 if not invert else -1)

        Y_exp_interp, self.win_high, self.win_low = self.interpolate(
            X, fname_ref_damp_norm
        )

        self.log = log

    @staticmethod
    def get_maximum(x, y):
        """Get the maximum between the maxima.

        :return: X,Y of the maximum of the maxima. This point must lay between two minima to be considered
        :rtype: tuple
        """

        y = np.abs(y)

        max_indices = argrelextrema(y, np.greater, order=10)[0]
        min_indices = argrelextrema(y, np.less, order=10)[0]

        first_nonzero_index = np.argmax(y != 0)
        last_nonzero_index = len(y) - np.argmax(y[::-1] != 0) - 1

        if y[first_nonzero_index] >= y[first_nonzero_index + 1]:
            max_indices = (
                np.concatenate((max_indices, [first_nonzero_index]))
                if len(max_indices) != 0
                else np.array([y[0]])
            )
        elif y[first_nonzero_index] < y[first_nonzero_index + 1]:
            min_indices = (
                np.concatenate((min_indices, [first_nonzero_index]))
                if len(min_indices) != 0
                else np.array([y[0]])
            )

        if y[last_nonzero_index] >= y[last_nonzero_index - 1]:
            max_indices = (
                np.concatenate((max_indices, [last_nonzero_index]))
                if len(max_indices) != 0
                else np.array([last_nonzero_index])
            )
        elif y[last_nonzero_index] < y[last_nonzero_index - 1]:
            min_indices = (
                np.concatenate((min_indices, [last_nonzero_index]))
                if len(min_indices) != 0
                else np.array([last_nonzero_index])
            )

        massimi = np.array([(x[i], y[i]) for i in max_indices])

        minimi = np.array([(x[i], y[i]) for i in min_indices])

        max_with_minima = []
        for i in max_indices:
            left_min_index = np.argmax(min_indices < i)
            right_min_index = np.argmax(min_indices > i)
            if left_min_index >= 0 and right_min_index <= len(min_indices) - 1:
                max_with_minima.append((x[i], y[i]))

        max_with_minima = np.array(max_with_minima)
        tmp = np.argmax(max_with_minima[:, 1])
        max_with_minima = list(max_with_minima[tmp])

        min_min = minimi[:, 0] < max_with_minima[0]
        min_max = minimi[:, 0] > max_with_minima[0]

        low = (
            minimi[min_min, 0][-1]
            if len(minimi[min_min, 0]) != 0
            else max_with_minima[0] - 0.3
        )
        up = (
            minimi[min_max, 0][0]
            if len(minimi[min_max, 0]) != 0
            else max_with_minima[0] + 0.3
        )
        return max_with_minima, [low, up]

    def get_ref_limits_index(self):
        lim = FACTOR_EV_NM / self.x_min
        lower_limit = self.x[self.x <= lim]
        lower_limit = list(self.x).index(lower_limit[-1])

        lim = FACTOR_EV_NM / self.x_max
        higher_limit = self.x[self.x >= lim]
        higher_limit = list(self.x).index(higher_limit[0])

        self.x_max_index = higher_limit
        self.x_min_index = lower_limit

    def interpolate(self, X, fname_ref_damp):
        Y_exp_interp = np.interp(X, self.x, self.y, left=0, right=0)
        if set(list(Y_exp_interp)) == set([0.0]):
            Y_exp_interp = np.interp(X, self.x[::-1], self.y[::-1], left=0, right=0)

        self.x = X
        self.max_exp, [low, up] = Ref_graph.get_maximum(self.x, self.y)
        window_low = np.where(self.x >= low)[0][-1]
        window_high = np.where(self.x <= up)[0][0]

        self.get_ref_limits_index()

        self.y = Graph.normalise(Y_exp_interp, norm=self.norm)

        Graph.damp_graph(fname_ref_damp, self.x, self.y)

        return Y_exp_interp, window_low, window_high


def eV_to_nm(eV):
    return FACTOR_EV_NM / eV


def main_graph(graphs, protocol, fname, output, title):
    fig, ax = plt.subplots(1, 1, dpi=300)
    ref = None
    if os.path.exists(os.path.join(os.getcwd(), fname)):
        ref = np.loadtxt(fname)
        ax.plot(ref[:, 0], ref[:, 1], label="Sperimentale", lw=1.5)
    for p in graphs:
        data_uv = np.loadtxt(f"{output}_protocol_{p}.dat")
        ax.plot(data_uv[:, 0], data_uv[:, 1], label=protocol[int(p)].functional, lw=1)

    if not (ref is None):
        if output == "uv":
            ax.set_ylim(0, 1.2)
        elif output == "ecd":
            ax.set_ylim(-1.2, 1.2)
        # pass

    plt.legend(
        loc="upper left",
        # bbox_to_anchor=(0.5, -0.25),
        fancybox=True,
        shadow=True,
        ncol=2,
    )

    ax.set_xlabel("Energia [eV]")

    ax2 = ax.secondary_xaxis("bottom", functions=(eV_to_nm, eV_to_nm))
    p = ax.get_position()
    # p = [p.x0, p.y0-0.5]
    ax2.spines["bottom"].set_position(("outward", p.y0 + 0.073 * fig.dpi))

    ax2.set_xlabel("Lunghezza d'onda (nm)")

    ax.set_xlim(1.8)
    plt.title(title)
    plt.ylabel("Intensità [u.a.]")
    plt.tight_layout()
    plt.savefig(f"regraphed_{output}.png", dpi=300)


def plot_conv_graph(graphs, protocol):
    """Create the electronic graphs calculated

    :param graphs: Indexes of the "graph" protocols
    :type graphs: list[int]
    :param protocol: The protocols list
    :type protocol: list[Protocol]
    """

    main_graph(graphs, protocol, "uv_ref_norm_eV.dat", "uv", "Comparazione grafico UV")
    main_graph(graphs, protocol, "ecd_ref_norm_eV.dat", "ecd", "ECD comparison graph")


if __name__ == "__main__":
    X = np.linspace(FACTOR_EV_NM / 150, FACTOR_EV_NM / 800, 10**4)
    Ref_graph("_ref.dat", None, X=X, fname_ref_damp_norm="None", norm=1)
