from abc import ABC, abstractmethod


CALCULATOR_REGISTRY = {}

def register_calculator(name):
    """Decorator to register each calculator """
    def decorator(cls):
        CALCULATOR_REGISTRY[name.lower()] = cls
        return cls
    return decorator

class BaseCalc(ABC):


    def __init__(self, protocol, cpu: int, conf=None):
        self.protocol = protocol
        self.cpu = cpu
        self.conf = conf
        self.constrains = protocol.constrains

    @abstractmethod
    def common_str(self):
        """String for the solvent definition"""
        pass

    @abstractmethod
    def single_point(self):
        """Create the common ASE calculator"""
        pass

    @abstractmethod
    def optimisation(self):
        """Add to the calculator the opt flag"""
        pass

    @abstractmethod
    def frequency(self):
        """Add to the calculator the opt flag"""
        pass