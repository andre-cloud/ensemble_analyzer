regex_parsing = {
    "orca": {
        "B": "Rotational constants in cm-1",
        "m": "Total Dipole Moment",
        "E": "FINAL SINGLE POINT ENERGY",
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
        "break": "\n\n",
        "idx_en_UV": 5,  # index for energy in the UV table in nm
        "idx_imp_UV": 6,  # index for oscillator strength in the UV table
        "idx_en_ECD": 5,  # index for energy in the ECD table in nm
        "idx_imp_ECD": 6,  # index for rotatory strength in the ECD table
        "s_freq": "VIBRATIONAL FREQUENCIES",
        "e_freq": "------------",
        "idx_freq": 1,  # index for frequency in frequency table
        "opt_done": "                    ***        THE OPTIMIZATION HAS CONVERGED     ***",
        "geom_start": """CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------""",
        "finish": "ORCA TERMINATED NORMALLY",
        "ext": "out",
    }
}
