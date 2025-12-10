from abc import ABC, abstractmethod
from typing import Dict

def register_parser(name):
    """Decorator to register each parser"""

    def decorator(cls):
        PARSER_REGISTRY[name.lower()] = cls
        return cls

    return decorator

class BaseParser(ABC):

    def __init__(self, output_name, log): 

        with open(output_name) as f:
            self.fl = f.read()

        self.log = log

        self.skip_message = "ATTENTION: Calculation CRASHED, impossible parsing. Conformer will be deactivated and no longer considered"
    
    @abstractmethod
    def parse_geom(self):
        pass

    @abstractmethod
    def parse_B_m(self):
        pass

    @abstractmethod
    def parse_energy(self):
        pass
    
    @abstractmethod
    def parse_freq(self):
        "Obtain: frequencies, IR XY graph and VCD XY graph"
        pass

    @abstractmethod
    def parse_tddft(self):
        "Obtain: absortion energies, f intensities and Rot Strength"
        pass

    @abstractmethod
    def opt_done(self) -> bool:
        pass

    @abstractmethod
    def normal_termination(self) -> bool:
        pass
    
    def get_filtered_text(self, start:str, end:str) -> str:
        return self.fl.split(start)[-1].split(end)[0]
    

    def parse_table(self, table:list, list_index:list):
        data = []
        for line in table: 
            if not line: 
                continue
            if '---' in line: 
                continue
            line_splitted = line.split()
            data.append([line_splitted[i] for i in list_index])
        
        return data



PARSER_REGISTRY : Dict[str, BaseParser]= {}