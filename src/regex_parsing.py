regex_parsing = {
    "orca": {
        "B": "Rotational constants in cm-1",
        'units_B': 'cm-1',
        "m": r"Total Dipole Moment\s*:\s*([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)",
        "E": r"FINAL SINGLE POINT ENERGY\s*(-?\d*.\d*)",
        "start_spec": "SPECTRA",
        "end_spec": "***",
        "s_UV": """ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS    
----------------------------------------------------------------------------------------------------
     Transition      Energy     Energy  Wavelength fosc(D2)      D2        DX        DY        DZ   
                      (eV)      (cm-1)    (nm)                 (au**2)    (au)      (au)      (au)  
----------------------------------------------------------------------------------------------------""",
        "s_ECD": """CD SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS    
------------------------------------------------------------------------------------------
     Transition      Energy     Energy  Wavelength    R        MX        MY        MZ   
                      (eV)      (cm-1)    (nm)   (1e40*cgs)   (au)      (au)      (au)  
------------------------------------------------------------------------------------------""",
        "s_IR": """Mode   freq       eps      Int      T**2         TX        TY        TZ
       cm**-1   L/(mol*cm) km/mol    a.u.
----------------------------------------------------------------------------
""",
        "s_VCD": """Mode   Freq    VCD-Intensity    
       (1/cm) (1E-44*esu^2*cm^2) 
---------------------------------""",
        "break": "\n\n",
        "idx_en_tddft": 3,  # index for energy in the UV & ECD table in eV
        "idx_imp_tddft": 6,  # index for oscillator strength in the UV table
        "idx_en_ir": 1,  # index for energy in the IR table in cm**-1
        "idx_imp_ir": 3,  # index for oscillator strength in the IR table
        "idx_imp_vcd": 2,  # index for oscillator strength in the VCD table
        "s_freq": "VIBRATIONAL FREQUENCIES",
        "e_freq": "------------",
        "idx_freq": 1,  # index for frequency in frequency table
        "opt_done": "THE OPTIMIZATION HAS CONVERGED",
        "geom_start": """CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------""",
        "finish": "ORCA TERMINATED NORMALLY",
        "ext": "out",
    },
    "gaussian": {
        "B": "Rotational constants",
        'units_B': 'GHz',
        "m": r"Dipole moment \(field-independent basis, Debye\):[\s\S]*?X=\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s+Y=\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s+Z=\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)",
        "E": r" SCF Done:.+ (-?\d+.\d+)",
        "start_spec": "",
        "end_spec": "",
        "s_UV": "",
        "s_ECD": "",
        "s_IR": "",
        "s_VCD": "",
        "break": "",
        "idx_en_tddft": 3,  # index for energy in the UV & ECD table in eV
        "idx_imp_tddft": 6,  # index for oscillator strength in the UV table
        "idx_en_ir": 1,  # index for energy in the IR table in cm**-1
        "idx_imp_ir": 3,  # index for oscillator strength in the IR table
        "idx_imp_vcd": 2,  # index for oscillator strength in the VCD table
        "s_freq": "",
        "e_freq": "",
        "idx_freq": 1,  # index for frequency in frequency table
        "opt_done": "",
        "geom_start": "",
        "finish": "",
        "ext": "log",
    },
}
