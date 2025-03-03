
import numpy as np
from types import UnionType
import os
from scipy.constants import c, h, electron_volt, R
import matplotlib.pyplot as plt 
from scipy.optimize import minimize, differential_evolution



FACTOR_EV_NM = h * c / (10**-9 * electron_volt)
FACTOR_EV_CM_1 = 1/8065.544 # to yield eV

def eV_to_nm(eV):
    return FACTOR_EV_NM / eV



class Graph: 

    START = {
        'IR':  0.1, # cm-1
        'VCD': 0.1, # cm-1
        'UV':  0.1, # eV
        'ECD': 0.1  # eV
    }
    END = {
        'IR': 6000,  # cm-1
        'VCD': 6000, # cm-1
        'UV': 9,     # eV
        'ECD': 9     # eV
    }


    def __init__(self, graph_type, protocol, definition=4, norm=1):
        
        assert graph_type in list(self.START.keys()) 

        self.protocol = protocol

        self.graph_type = graph_type
        self.x = np.linspace(self.START[graph_type], self.END[graph_type], num=10**definition)
        self.y = np.zeros(self.x.shape)
        self.norm = norm

        self.CONV = {
            'IR': self.lorentzian,
            'VCD': self.lorentzian,
            'UV': self.gaussian,
            'ECD': self.gaussian,
            'lorentzian':self.lorentzian,
            'gaussian':self.gaussian,
        }


    def normalize(self): 
        return self.y / np.max(np.abs(self.y)) * self.norm

    def gaussian(self, x0, I, fwhm):
        sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
        y = np.sum(I[:, np.newaxis] / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((self.x - x0[:, np.newaxis]) / sigma) ** 2), axis=0)
        return y

    def lorentzian(self, x0, I, fwhm):
        y = np.sum((I[:, np.newaxis] * fwhm**2) / (fwhm**2 + 4 * (self.x - x0[:, np.newaxis]) ** 2), axis=0)
        return y



class Computed(Graph): 

    DEFs = {
        'IR': 50,      # σ in cm-1
        'VCD': 50,     # σ in cm-1
        'UV': 0.4,     # FWHM in eV
        'ECD': 0.4     # FWHM in eV
    }

    def __init__(self, conf, invert, convolution = None, shift = None, fwhm = None, **kwargs):
        super().__init__(**kwargs)

        self.invert = invert
        self.g = convolution if convolution else self.graph_type
        self.conformers = conf
        if conf[0].energies[self.protocol].get("graph", None) is None:
            return
        elif conf[0].energies[self.protocol].get("graph").get(self.graph_type, None) is None:
            return
        
        if shift:
            self.shift = shift
        if fwhm:        
            self.fwhm = fwhm

        self.retrive_data()

        self.auto = os.path.exists(os.path.join(os.getcwd(), f'{self.graph_type.lower()}_ref.dat'))

        if self.auto: 
            self.ref = Experimental(**kwargs)
            self.autoconvolute()
        else:
            self.y = self.CONV[self.g](self.x, self.y_comp, self.DEFs[self.graph_type])
            self.y = self.normalize()
            self.shift = 0
            self.fwhm = self.DEFs[self.graph_type]


    def retrive_data(self):
        self.y_comp = np.zeros(self.x.shape)
        f = self.CONV[self.graph_type]
        for i in self.conformers:
            if not i.active:
                continue
            x = np.array(i.energies[self.protocol]["graph"][self.graph_type]['x'])
            y = np.array(i.energies[self.protocol]["graph"][self.graph_type]['y'])
            self.y_comp += f(x, y, self.DEFs[self.graph_type]) * float(i.energies[self.protocol]['Pop'])


    def autoconvolute(self):
        f = self.CONV[self.graph_type]
        
        def wrapper(variables):
            shift, fwhm = variables
            y = f(self.x + shift, self.y_comp, fwhm)
            y_norm = y / np.max(np.abs(y[self.ref.x_min_idx:self.ref.x_max_idx]))
            ref_norm = self.ref.y / np.max(np.abs(self.ref.y[self.ref.x_min_idx:self.ref.x_max_idx]))
            d = self.diversity_function(y_norm[self.ref.x_min_idx:self.ref.x_max_idx],
                                        ref_norm[self.ref.x_min_idx:self.ref.x_max_idx])
            print(f"Shift: {shift:.4f}, FWHM: {fwhm:.4f}, RMSD: {d:.6f}")
            return d  # Minimizzare questa funzione

        sb, fb = 0.6, 0.3
        ss, sf = 0, self.DEFs[self.graph_type]

        if hasattr(self, 'shift'):
            ss = self.shift
        if hasattr(self, 'fwhm'):
            sf = self.fwhm

        initial_guess = [ss, sf]  # BLUE SHIFT NEGATIVE
        bounds = [(-sb, sb), (0.2, sf + fb)]

        result = minimize(wrapper, initial_guess, bounds=bounds, options={'maxiter': 1000}, method='L-BFGS-B')
        # result = minimize(wrapper, initial_guess, bounds=bounds, options={'maxiter': 1000}, method='Nelder-Mead')

        if result.success:
            self.shift, self.fwhm = result.x
            self.y = self.CONV[self.graph_type](self.x + self.shift, self.y_comp, self.fwhm)
            self.y = self.normalize()
            print(f"Optimal shift: {self.shift:.3f}, Optimal FWHM: {self.fwhm:.3f}, Similarity: {(1-result.fun)*100:.2f}%")
        else:
            print("Optimization failed")

    def diversity_function(self, a, b):
        return np.sqrt(np.mean((a - b) ** 2))  # Calcolo corretto di RMSD

    

class Experimental(Graph): 

    DEFs = {
        'IR': 3,
        'VCD': 3,
        'UV': 100,
        'ECD': 100,
    }

    FACTOR = {
        'IR': FACTOR_EV_CM_1, 
        'VCD': FACTOR_EV_CM_1, 
        'UV': FACTOR_EV_NM,
        'ECD': FACTOR_EV_NM,
    }

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # Read the reference data
        self.data = np.loadtxt(os.path.join(os.getcwd(), f'{self.graph_type.lower()}_ref.dat'))

        # Sort by the Xs
        self.data = self.data[np.argsort(self.data[:, 0])]

        self.x_min = float(min(self.data[:, 0]))
        self.x_max = float(max(self.data[:, 0]))
        self.x_max_idx = np.argmax(self.data[:, 0])

        self.convert = len(self.data[self.data[:, 0]>self.DEFs[self.graph_type], :]) != 0 

        if self.convert:
            assert self.graph_type not in ['IR', 'VCD']
            self.x_exp = self.FACTOR[self.graph_type] / self.data[:, 0]  
            self.x_max, self.x_min = self.FACTOR[self.graph_type]/self.x_min, self.FACTOR[self.graph_type]/self.x_max

        else:
            self.x_exp = self.data[:, 0]



        self.reverse = self.x_exp[0] > self.x_exp[-1]

        self.y_exp = self.data[:, 1]
        if self.graph_type in ['UV', 'IR']:
            self.y_exp -= np.min(self.y_exp)

        if self.reverse: 
            self.x_exp = self.x_exp[::-1]
            self.y_exp = self.y_exp[::-1]

        self.y = self.interpolate()
        tmp = (self.x > self.x_min) & (self.x < self.x_max)
        print(tmp)
        print(self.x[tmp])

        self.x_min_idx = np.where(self.x == np.min(self.x[tmp]))[0][0]
        self.x_max_idx = np.where(self.x == np.max(self.x[tmp]))[0][0]


        self.y = self.normalize()

    def interpolate(self):
        return np.interp(self.x, self.x_exp, self.y_exp, left=0, right=0)






if __name__ == '__main__':

    from launch import restart
    from pruning import calculate_rel_energies

    ensemble, protocol, _= restart()
    calculate_rel_energies(ensemble, 298.15)

    for i in  protocol: 
        shift, fwhm = None, None
        for j in ['IR', 'UV', 'ECD']:
            graph = Computed(ensemble, False, graph_type=j, protocol=str(i.number), shift=shift, fwhm=fwhm)

            if not hasattr(graph, "y_comp"):
                continue

            shift, fwhm = graph.shift, graph.fwhm
            plt.plot(graph.x, graph.y)
            if hasattr(graph, 'ref'):
                plt.plot(graph.ref.x[graph.ref.x_min_idx:graph.ref.x_max_idx], graph.ref.y[graph.ref.x_min_idx:graph.ref.x_max_idx])
                plt.xlim(graph.ref.x_min-.5, graph.ref.x_max+.5)

            plt.title(f'{j} of protocol {i.number}')
            plt.show()
