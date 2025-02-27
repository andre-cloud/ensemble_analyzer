
import numpy as np
from types import UnionType
import os
from scipy.constants import c, h, electron_volt, R
import matplotlib.pyplot as plt 
from scipy.optimize import minimize



FACTOR_EV_NM = h * c / (10**-9 * electron_volt)
FACTOR_EV_CM_1 = 1/8065.544 # to yield eV

def eV_to_nm(eV):
    return FACTOR_EV_NM / eV



class Graph: 

    START = {
        'IR': 300,  # cm-1
        'VCD': 300, # cm-1
        'UV': 0.1,  # eV
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

    def gaussian(self, x0, I, fwmh):
        sigma = fwmh / (2 * np.sqrt(2 * np.log(2)))
        y = np.zeros(self.x.shape)
        for x, i in zip(x0, I):
            y += i / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((self.x - x) / sigma) ** 2)
        return y
    
    def lorentzian(self, x0, I, fwhm):
        return (I*fwhm**2)/(fwhm**2+4*(self.x-x0)**2)



class Computed(Graph): 

    DEFs = {
        'IR': 50,      # σ in cm-1
        'VCD': 50,     # σ in cm-1
        'UV': 1/3,     # FWHM in eV
        'ECD': 1/3     # FWHM in eV
    }

    def __init__(self, conf, invert, convolution = None, **kwargs):
        super().__init__(**kwargs)

        self.invert = invert
        self.g = convolution if convolution else self.graph_type
        self.conformers = conf


        self.convolution_model = {
            'IR': self.lorentzian,
            'VCD': self.lorentzian,
            'UV': self.gaussian,
            'ECD': self.gaussian,
        }


        self.retrive_data()

        self.auto = os.path.exists(os.path.join(os.getcwd(), f'{self.graph_type.lower()}_ref.dat'))

        if self.auto: 
            self.ref = Experimental(**kwargs)
            self.autoconvolute()


    def retrive_data(self):
        self.y_comp = np.zeros(self.x.shape)
        f = self.convolution_model[self.graph_type]
        for i in self.conformers:
            if not i.active:
                continue
            x = np.array(i.energies[self.protocol]["graph"][self.graph_type]['x'])
            y = np.array(i.energies[self.protocol]["graph"][self.graph_type]['y'])

            print(f(
                x, y, 
                self.DEFs[self.graph_type]))
            self.y_comp += f(
                x, y, 
                self.DEFs[self.graph_type]) \
                    * float(i.energies[self.protocol]['Pop'])


    def autoconvolute(self):
        
        f = self.convolution_model[self.graph_type]

        def wrapper(variables):
            shift, fwhm = variables
            y = f(self.x-shift, self.y_comp,fwhm)
            n = np.max(np.abs(y[self.ref.x_min_idx:self.ref.x_max_idx]))
            y /= n
            d = self.diversity_function(y[self.ref.x_min_idx:self.ref.x_max_idx]/n
                                        , self.ref.y[self.ref.x_min_idx:self.ref.x_max_idx])
            print(d)
            return d
        
        # Initial guess for the parameters
        initial_guess = [0, self.DEFs[self.graph_type]]

        # Bounds for the parameters
        bounds = [(-.75, .75), (0.1, 0.5)]  # Example bounds

        # Fit the convolution model to the experimental data

        result = minimize(wrapper, initial_guess, bounds=bounds, method='L-BFGS-B', options={'maxiter': 1000})

        # Extract the optimal parameters
        if result.success:
            shift_opt, fwhm_opt = result.x
            self.y = self.convolution_model[self.graph_type](self.x - shift_opt, self.y_comp, fwhm_opt)
            self.y = self.normalize()
            print(f"Optimal shift: {shift_opt}, Optimal FWHM: {fwhm_opt}")
        else:
            print("Optimization failed")

    def diversity_function(self, y_comp, y_ref):
        return np.sqrt(np.mean((y_comp - y_ref) ** 2))

    

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
        self.x_min_idx = np.argmin(self.data[:, 0])
        self.x_max = float(max(self.data[:, 0]))
        self.x_max_idx = np.argmax(self.data[:, 0])

        self.convert = len(self.data[self.data[:, 0]>self.DEFs[self.graph_type], :]) != 0 

        if self.convert:
            assert self.graph_type not in ['IR', 'VCD']
            self.x_exp = self.FACTOR[self.graph_type] / self.data[:, 0]  
            self.x_max, self.x_min = self.FACTOR[self.graph_type]/self.x_min, self.FACTOR[self.graph_type]/self.x_max

        else:
            self.x_exp = self.data[:, 0]


        print(self.x_exp)
        self.reverse = self.x_exp[0] > self.x_exp[-1]

        self.y_exp = self.data[:, 1]
        if self.graph_type in ['UV', 'IR']:
            self.y_exp -= np.min(self.y_exp)

        if self.reverse: 
            self.x_exp = self.x_exp[::-1]
            self.y_exp = self.y_exp[::-1]

        self.y = self.interpolate()
        self.y = self.normalize()

    def interpolate(self):
        return np.interp(self.x, self.x_exp, self.y_exp, left=0, right=0)






if __name__ == '__main__':

    from conformer import Conformer
    import json
    
    j = json.loads(open('checkpoint.json').read())
    c = Conformer.load_raw(j['1'])

    graph = Computed([c], False, graph_type='ECD', protocol="6")

    plt.plot(graph.x, graph.y)
    plt.plot(graph.ref.x, graph.ref.y)
    # plt.xlim(graph.ref.x_min, graph.ref.x_max)
    plt.show()
