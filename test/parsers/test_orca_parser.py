import pytest
import numpy as np
from unittest.mock import MagicMock, mock_open, patch
from ensemble_analyzer.parsers.orca import OrcaParser


ORCA5_OPT = """
Program Version 5

ORCA TERMINATED NORMALLY
"""

ORCA6_OPT = """
Program Version 6

ORCA TERMINATED NORMALLY
"""

ORCA6_OPT_CONVERGED = """
Program Version 6

CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  C    0.000000    0.000000    0.000000
  O    0.000000    0.000000    1.200000

THE OPTIMIZATION HAS CONVERGED

ORCA TERMINATED NORMALLY
"""

ORCA6_ENERGY = """
Program Version 6

FINAL SINGLE POINT ENERGY    -113.45678900

ORCA TERMINATED NORMALLY
"""

ORCA6_BCONST = """
Program Version 6

Rotational constants in cm-1:       1.50000       0.20000       0.10000
Total Dipole Moment   :     0.50000     0.30000     0.10000

ORCA TERMINATED NORMALLY
"""

ORCA6_FREQ = """
Program Version 6

VIBRATIONAL FREQUENCIES
------------
   0:      500.00 cm**-1
   1:     1000.00 cm**-1
   2:       -50.00 cm**-1



Mode   freq       eps      Int      T**2         TX        TY        TZ
       cm**-1   L/(mol*cm) km/mol    a.u.
----------------------------------------------------------------------------
   0    500.00    1.00    10.00    0.0100    0.1000    0.0000    0.0000
   1   1000.00    2.00    20.00    0.0200    0.2000    0.0000    0.0000



Mode   Freq    VCD-Intensity    
       (1/cm) (1E-44*esu^2*cm^2) 
---------------------------------
   0    500.00    0.50
   1   1000.00   -0.30



ORCA TERMINATED NORMALLY
"""

ORCA5_FREQ = """
Program Version 5

VIBRATIONAL FREQUENCIES
------------
   0:      500.00 cm**-1
   1:     1000.00 cm**-1

Mode   freq       eps      Int      T**2         TX        TY        TZ
       cm**-1   L/(mol*cm) km/mol    a.u.
----------------------------------------------------------------------------
   0    500.00    1.00    10.00    0.0100    0.1000    0.0000    0.0000
   1   1000.00    2.00    20.00    0.0200    0.2000    0.0000    0.0000

ORCA TERMINATED NORMALLY
"""

# Layout of ORCA 6.1 TD-DFT outputs (tables, energy summary, dispersion, final energy).
_ORCA6_SPECTRA = """
----------------------------------------------------------------------------------------------------
                     ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
----------------------------------------------------------------------------------------------------
     Transition      Energy     Energy  Wavelength fosc(D2)      D2        DX        DY        DZ
                      (eV)      (cm-1)    (nm)                 (au**2)    (au)      (au)      (au)
----------------------------------------------------------------------------------------------------
  0-1A  ->  1-1A    3.500000   28229.4   354.2   0.123400000   0.50000   0.10000  -0.20000   0.30000
  0-1A  ->  2-1A    4.000000   32262.1   310.0   0.567800000   0.80000   0.40000   0.50000   0.60000

----------------------------------------------------------------------------------------------------
                     ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS
----------------------------------------------------------------------------------------------------
     Transition      Energy     Energy  Wavelength fosc(P2)      P2        PX        PY        PZ
                      (eV)      (cm-1)    (nm)                 (au**2)    (au)      (au)      (au)
----------------------------------------------------------------------------------------------------
  0-1A  ->  1-1A    3.500000   28229.4   354.2   0.999900000   0.50000   0.10000  -0.20000   0.30000
  0-1A  ->  2-1A    4.000000   32262.1   310.0   0.999900000   0.80000   0.40000   0.50000   0.60000

------------------------------------------------------------------------------------------
                     CD SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
------------------------------------------------------------------------------------------
     Transition      Energy     Energy  Wavelength    R        MX        MY        MZ
                      (eV)      (cm-1)    (nm)   (1e40*cgs)   (au)      (au)      (au)
------------------------------------------------------------------------------------------
  0-1A  ->  1-1A    3.500000   28229.4   354.2   12.34560   0.10000   0.20000   0.30000
  0-1A  ->  2-1A    4.000000   32262.1   310.0  -78.90120   0.40000   0.50000   0.60000

------------------------------------------------------------------------------------------
                     CD SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS
------------------------------------------------------------------------------------------
     Transition      Energy     Energy  Wavelength    R        MX        MY        MZ
                      (eV)      (cm-1)    (nm)   (1e40*cgs)   (au)      (au)      (au)
------------------------------------------------------------------------------------------
  0-1A  ->  1-1A    3.500000   28229.4   354.2   99.99990   0.10000   0.20000   0.30000
  0-1A  ->  2-1A    4.000000   32262.1   310.0   99.99990   0.40000   0.50000   0.60000

"""

_ORCA6_SOC_SPECTRA = """
--------------------------------------------------------------------------------------------------------
      SOC CORRECTED ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
--------------------------------------------------------------------------------------------------------
      Transition         Energy     Energy  Wavelength fosc(D2)      D2       |DX|      |DY|      |DZ|
                          (eV)      (cm-1)    (nm)                 (au**2)    (au)      (au)      (au)
--------------------------------------------------------------------------------------------------------
  0-1.0A  ->  1-3.0A    3.100000   25003.2   400.0   0.000000369   0.00000   0.00000   0.00150   0.00000
  0-1.0A  ->  2-3.0A    3.100000   25003.2   400.0   0.000000000   0.00000   0.00000   0.00000   0.00000
  0-1.0A  ->  3-3.0A    3.100010   25003.3   400.0   0.000001102   0.00001   0.00000   0.00000   0.00260
  0-1.0A  ->  4-1.0A    3.499970   28229.2   354.2   0.123390000   0.50000   0.10000   0.20000   0.30000

----------------------------------------------------------------------------------------------
      SOC CORRECTED CD SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
----------------------------------------------------------------------------------------------
      Transition         Energy     Energy  Wavelength    R       |MX|      |MY|      |MZ|
                          (eV)      (cm-1)    (nm)    (1e40*cgs)  (au)      (au)      (au)
----------------------------------------------------------------------------------------------
  0-1.0A  ->  1-3.0A    3.100000   25003.2   400.0   -0.00000   0.00008   0.00000   0.00000
  0-1.0A  ->  2-3.0A    3.100000   25003.2   400.0   -0.00000   0.00000   0.00000   0.00005
  0-1.0A  ->  3-3.0A    3.100010   25003.3   400.0    0.00000   0.00000   0.00000   0.00000
  0-1.0A  ->  4-1.0A    3.499970   28229.2   354.2   12.34500   0.10000   0.20000   0.30000

"""

_ORCA6_TDDFT_ENERGY = """
-----------------------
CIS/TD-DFT TOTAL ENERGY
-----------------------

    E(SCF)  =   -914.136983097 Eh
    DE(CIS) =      0.065909138 Eh (Root  1)
    ----------------------------- ---------
    E(tot)  =   -914.071073959 Eh

-------------------------   ----------------
Dispersion correction           -0.061064496
-------------------------   ----------------

-------------------------   --------------------
FINAL SINGLE POINT ENERGY      -914.132138455
-------------------------   --------------------
"""

_ORCA6_SOC_ENERGY = """
---------------------------
SOC CIS/TD-DFT TOTAL ENERGY
---------------------------

    E(SCF)        =  -3734.600483805 Eh
    E(SOC CIS)    =      0.092111444 Eh (SOC Root  1)
    ----------------------------- ---------
    E(tot)        =  -3734.508372362 Eh

-------------------------   ----------------
Dispersion correction           -0.086497399
-------------------------   ----------------

-------------------------   --------------------
FINAL SINGLE POINT ENERGY     -3734.594869760291
-------------------------   --------------------
"""

ORCA6_TDDFT = (
    "\n                         Program Version 6.1.1  -  RELEASE   -\n"
    + _ORCA6_SPECTRA + _ORCA6_TDDFT_ENERGY + "\nORCA TERMINATED NORMALLY\n"
)

# `iroot 0`: no CIS/TD-DFT TOTAL ENERGY block, FINAL is E(SCF) + dispersion
ORCA6_TDDFT_IROOT0 = (
    "\n                         Program Version 6.1.1  -  RELEASE   -\n"
    + _ORCA6_SPECTRA
    + "\nFINAL SINGLE POINT ENERGY      -914.198047593\n\nORCA TERMINATED NORMALLY\n"
)

ORCA6_TDDFT_SOC = (
    "\n                         Program Version 6.1.0  -  RELEASE   -\n"
    + _ORCA6_SPECTRA + _ORCA6_SOC_SPECTRA + _ORCA6_SOC_ENERGY + "\nORCA TERMINATED NORMALLY\n"
)

ORCA6_TDDFT_SOC_ROOT0 = (
    "\n                         Program Version 6.1.0  -  RELEASE   -\n"
    + _ORCA6_SPECTRA + _ORCA6_SOC_SPECTRA
    + _ORCA6_SOC_ENERGY.replace("0.092111444 Eh (SOC Root  1)", "-0.000000248 Eh (SOC Root  0)")
    + "\nORCA TERMINATED NORMALLY\n"
)

# Layout of ORCA 5.0 TD-DFT outputs: energies in cm-1 only, "CD SPECTRUM" title
ORCA5_TDDFT = """
                         Program Version 5.0.4 -  RELEASE  -

-----------------------------------------------------------------------------
         ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
-----------------------------------------------------------------------------
State   Energy    Wavelength  fosc         T2        TX        TY        TZ
        (cm-1)      (nm)                 (au**2)    (au)      (au)      (au)
-----------------------------------------------------------------------------
   1   28233.0    354.2   0.123400000   0.50000   0.10000  -0.20000   0.30000
   2   32266.0    309.9   0.567800000   0.80000   0.40000   0.50000   0.60000

-----------------------------------------------------------------------------
         ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS
-----------------------------------------------------------------------------
State   Energy    Wavelength  fosc         P2        PX        PY        PZ
        (cm-1)      (nm)                 (au**2)    (au)      (au)      (au)
-----------------------------------------------------------------------------
   1   28233.0    354.2   0.999900000   0.50000   0.10000  -0.20000   0.30000
   2   32266.0    309.9   0.999900000   0.80000   0.40000   0.50000   0.60000

-------------------------------------------------------------------
                             CD SPECTRUM
-------------------------------------------------------------------
State   Energy    Wavelength       R         MX        MY        MZ
        (cm-1)      (nm)       (1e40*cgs)   (au)      (au)      (au)
-------------------------------------------------------------------
   1   28233.0    354.2      12.34560   0.10000   0.20000   0.30000
   2   32266.0    309.9     -78.90120   0.40000   0.50000   0.60000

-------------------------------------------------------------------
         CD SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS
-------------------------------------------------------------------
State   Energy    Wavelength       R         MX        MY        MZ
        (cm-1)      (nm)       (1e40*cgs)   (au)      (au)      (au)
-------------------------------------------------------------------
   1   28233.0    354.2      99.99990   0.10000   0.20000   0.30000
   2   32266.0    309.9      99.99990   0.40000   0.50000   0.60000

-----------------------
CIS/TD-DFT TOTAL ENERGY
-----------------------

    E(SCF)  =  -3479.038876466 Eh
    DE(CIS) =      0.197863977 Eh (Root  1)
    ----------------------------- ---------
    E(tot)  =  -3478.841012489 Eh

-------------------------   ----------------
Dispersion correction           -0.061064496
-------------------------   ----------------

-------------------------   --------------------
FINAL SINGLE POINT ENERGY     -3478.902076985325
-------------------------   --------------------

ORCA TERMINATED NORMALLY
"""

ORCA_CRASHED = """
Program Version 6

Some partial output but no termination
"""

ORCA6_NORMAL_MODES = """
Program Version 6

NORMAL MODES
------------

These modes are the Cartesian displacements weighted by the diagonal matrix
M(i,i)=1/sqrt(m[i]) where m[i] is the mass of the displaced atom
Thus, these vectors are normalized but *not* orthogonal

                  0          1          2          3          4          5    
      0       1.000000   0.000000   0.000000   0.000000   0.000000   0.000000
      1       0.000000   1.000000   0.000000   0.000000   0.000000   0.000000
      2       0.000000   0.000000   1.000000   0.000000   0.000000   0.000000
      3       0.000000   0.000000   0.000000   1.000000   0.000000   0.000000
      4       0.000000   0.000000   0.000000   0.000000   1.000000   0.000000
      5       0.000000   0.000000   0.000000   0.000000   0.000000   1.000000
                  6    
      0       0.010000
      1       0.020000
      2       0.030000
      3       0.040000
      4       0.050000
      5       0.060000


-----------
IR SPECTRUM
-----------

ORCA TERMINATED NORMALLY
"""


class TestOrcaParser:

    @pytest.fixture
    def parser_v6(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=ORCA6_OPT)):
            p = OrcaParser("dummy.out", mock_logger)
            p.fl = ""
            return p

    def test_version_detection_v6(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=ORCA6_OPT)):
            p = OrcaParser("dummy.out", mock_logger)
            assert p.version == "6"

    def test_version_detection_v5(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=ORCA5_OPT)):
            p = OrcaParser("dummy.out", mock_logger)
            assert p.version == "5"

    def test_normal_termination(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=ORCA6_OPT)):
            p = OrcaParser("dummy.out", mock_logger)
            assert p.normal_termination() is True

    def test_crashed_calculation(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=ORCA_CRASHED)):
            p = OrcaParser("dummy.out", mock_logger)
            assert p.correct_exiting is False
            assert p.normal_termination() is False
            p.log.warning.assert_called()

    def test_parse_geom(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=ORCA6_OPT_CONVERGED)):
            p = OrcaParser("dummy.out", mock_logger)
            geom = p.parse_geom()
            assert geom.shape == (2, 3)
            assert np.allclose(geom[0], [0.0, 0.0, 0.0])
            assert np.allclose(geom[1], [0.0, 0.0, 1.2])

    def test_opt_done(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=ORCA6_OPT_CONVERGED)):
            p = OrcaParser("dummy.out", mock_logger)
            assert p.opt_done() is True

    def test_opt_not_done(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=ORCA5_OPT)):
            p = OrcaParser("dummy.out", mock_logger)
            assert p.opt_done() is False

    def test_parse_energy_v6(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=ORCA6_ENERGY)):
            p = OrcaParser("dummy.out", mock_logger)
            assert p.parse_energy() == -113.45678900

    def test_parse_energy_not_found(self, parser_v6):
        parser_v6.fl = "no energy here"
        assert parser_v6.parse_energy() == 0.0

    def test_parse_B_m(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=ORCA6_BCONST)):
            p = OrcaParser("dummy.out", mock_logger)
            B, M = p.parse_B_m()
            assert np.allclose(B, [1.5, 0.2, 0.1])
            assert np.allclose(M, [0.5, 0.3, 0.1])

    def test_parse_B_m_not_found(self, parser_v6):
        parser_v6.fl = "no B or M here"
        B, M = parser_v6.parse_B_m()
        assert np.allclose(B, [1, 0, 0])
        assert np.allclose(M, [1, 0, 0])

    def test_parse_freq_v6(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=ORCA6_FREQ)):
            p = OrcaParser("dummy.out", mock_logger)
            freq, ir, vcd = p.parse_freq()
            assert len(freq) == 3
            assert np.allclose(freq, [500.0, 1000.0, -50.0])
            assert ir.shape[0] == 2
            assert vcd.shape[0] == 2

    def test_parse_freq_no_vcd_in_v5(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=ORCA5_FREQ)):
            p = OrcaParser("dummy.out", mock_logger)
            freq, ir, vcd = p.parse_freq()
            assert vcd.shape == (1, 2)

    def test_parse_freq_not_found(self, parser_v6):
        parser_v6.fl = "no frequencies"
        freq, ir, vcd = parser_v6.parse_freq()
        assert len(freq) == 0
        assert ir.shape == (1, 2)
        assert vcd.shape == (1, 2)

    @staticmethod
    def _parser(mock_logger, text):
        with patch("builtins.open", mock_open(read_data=text)):
            return OrcaParser("dummy.out", mock_logger)

    # --- TD-DFT spectra -------------------------------------------------

    def test_parse_tddft_v6_energy_ev_and_fosc(self, mock_logger):
        uv, ecd = self._parser(mock_logger, ORCA6_TDDFT).parse_tddft()
        assert np.allclose(uv, [[3.5, 0.1234], [4.0, 0.5678]])
        assert np.allclose(ecd, [[3.5, 12.3456], [4.0, -78.9012]])

    def test_parse_tddft_v6_ignores_velocity_tables(self, mock_logger):
        uv, ecd = self._parser(mock_logger, ORCA6_TDDFT).parse_tddft()
        assert not np.any(np.isclose(uv[:, 1], 0.9999))
        assert not np.any(np.isclose(ecd[:, 1], 99.9999))

    def test_parse_tddft_v5_converts_cm1_to_ev(self, mock_logger):
        uv, ecd = self._parser(mock_logger, ORCA5_TDDFT).parse_tddft()
        assert uv.shape == (2, 2)
        assert np.allclose(uv[:, 0], [28233.0 / 8065.544, 32266.0 / 8065.544])
        assert np.isclose(uv[0, 0], 3.5004, atol=1e-4)
        assert np.allclose(uv[:, 1], [0.1234, 0.5678])
        # "CD SPECTRUM" (length gauge), not the velocity table
        assert np.allclose(ecd[:, 1], [12.3456, -78.9012])

    def test_parse_tddft_soc_corrected_preferred(self, mock_logger):
        uv, ecd = self._parser(mock_logger, ORCA6_TDDFT_SOC).parse_tddft()
        assert uv.shape == (4, 2)
        assert np.isclose(uv[0, 0], 3.1)
        assert np.isclose(uv[3, 1], 0.12339)
        assert ecd.shape == (4, 2)
        assert np.isclose(ecd[3, 1], 12.345)

    def test_parse_tddft_not_found(self, parser_v6):
        parser_v6.fl = "no spectra"
        uv, ecd = parser_v6.parse_tddft()
        assert uv.shape == (1, 2)
        assert ecd.shape == (1, 2)

    def test_parse_tddft_missing_table_does_not_parse_whole_file(self, mock_logger):
        text = "Program Version 6\nSPECTRA\n" + ORCA6_FREQ
        uv, ecd = self._parser(mock_logger, text).parse_tddft()
        assert np.allclose(uv, 0.0) and uv.shape == (1, 2)
        assert np.allclose(ecd, 0.0) and ecd.shape == (1, 2)

    # --- TD-DFT energies ------------------------------------------------

    def test_tddft_energy_is_ground_state(self, mock_logger):
        p = self._parser(mock_logger, ORCA6_TDDFT)
        # FINAL (-914.132138455) - DE(CIS, root 1) = E(SCF) + dispersion
        assert np.isclose(p.parse_energy(), -914.136983097 - 0.061064496)

    def test_tddft_energy_excited_state_kept_on_request(self, mock_logger):
        p = self._parser(mock_logger, ORCA6_TDDFT)
        assert np.isclose(p.parse_energy(excited_state=True), -914.132138455)

    def test_tddft_energy_v5(self, mock_logger):
        p = self._parser(mock_logger, ORCA5_TDDFT)
        assert np.isclose(p.parse_energy(), -3479.038876466 - 0.061064496)

    def test_tddft_energy_iroot0_unchanged(self, mock_logger):
        p = self._parser(mock_logger, ORCA6_TDDFT_IROOT0)
        assert np.isclose(p.parse_energy(), -914.198047593)

    def test_tddft_energy_soc_root1(self, mock_logger):
        p = self._parser(mock_logger, ORCA6_TDDFT_SOC)
        assert np.isclose(p.parse_energy(), -3734.600483805 - 0.086497399)

    def test_tddft_energy_soc_root0_unchanged(self, mock_logger):
        p = self._parser(mock_logger, ORCA6_TDDFT_SOC_ROOT0)
        assert np.isclose(p.parse_energy(), -3734.594869760291)

    def test_tddft_energy_ignores_excitation_of_previous_job(self, mock_logger):
        # TD-DFT job followed by a plain SCF: the excitation belongs to the first FINAL
        text = ORCA6_TDDFT + "\nFINAL SINGLE POINT ENERGY      -914.300000000\n"
        assert np.isclose(self._parser(mock_logger, text).parse_energy(), -914.3)

    # --- version handling -----------------------------------------------

    def test_version_detection_full_string(self, mock_logger):
        assert self._parser(mock_logger, ORCA6_TDDFT).version == "6"
        assert self._parser(mock_logger, ORCA5_TDDFT).version == "5"

    def test_version_4_uses_orca5_layout(self, mock_logger):
        p = self._parser(mock_logger, "Program Version 4.2.1\nORCA TERMINATED NORMALLY\n")
        assert p.version == "4"
        assert p.regex is OrcaParser.REGEX["5"]

    def test_version_missing_does_not_raise(self, mock_logger):
        p = self._parser(mock_logger, "ORCA TERMINATED NORMALLY\n")
        assert p.version == "0"
        assert p.regex is OrcaParser.REGEX["6"]
        mock_logger.warning.assert_called()

    def test_parse_normal_modes(self, mock_logger):
        with patch("builtins.open", mock_open(read_data=ORCA6_NORMAL_MODES)):
            p = OrcaParser("dummy.out", mock_logger)
            modes = p.parse_normal_modes(n_atoms=2)
        assert modes.shape == (7, 2, 3)
        assert np.allclose(modes[0, 0], [1.0, 0.0, 0.0])
        assert np.allclose(modes[0, 1], [0.0, 0.0, 0.0])
        assert np.allclose(modes[1, 0], [0.0, 1.0, 0.0])
        assert np.allclose(modes[1, 1], [0.0, 0.0, 0.0])
        assert np.allclose(modes[2, 0], [0.0, 0.0, 1.0])
        assert np.allclose(modes[2, 1], [0.0, 0.0, 0.0])
        assert np.allclose(modes[3, 0], [0.0, 0.0, 0.0])
        assert np.allclose(modes[3, 1], [1.0, 0.0, 0.0])
        assert np.allclose(modes[4, 0], [0.0, 0.0, 0.0])
        assert np.allclose(modes[4, 1], [0.0, 1.0, 0.0])
        assert np.allclose(modes[5, 0], [0.0, 0.0, 0.0])
        assert np.allclose(modes[5, 1], [0.0, 0.0, 1.0])
        assert np.allclose(modes[6, 0], [0.01, 0.02, 0.03])
        assert np.allclose(modes[6, 1], [0.04, 0.05, 0.06])

    def test_parse_normal_modes_not_found(self, parser_v6):
        parser_v6.fl = "no normal modes"
        modes = parser_v6.parse_normal_modes(n_atoms=2)
        assert modes.shape == (0, 2, 3)

from ensemble_analyzer.calculators.orca import _split_post_blocks

class TestSplitPostBlocks:
    """Tests for _split_post_blocks edge cases (indented vs compact)."""

    def test_indented_frag(self):
        t = "%frag\n Definition\n  1 {18:27} end\n  2 {0:17} end\n end\nend"
        pre, post = _split_post_blocks(t)
        assert pre == ""
        assert post == "%frag\n Definition\n  1 {18:27} end\n  2 {0:17} end\n end\nend"

    def test_compact_frag(self):
        t = "%frag\nDefinition\n1 {18:27} end\n2 {0:17} end\nend\nend"
        pre, post = _split_post_blocks(t)
        assert pre == ""
        assert post == "%frag\nDefinition\n1 {18:27} end\n2 {0:17} end\nend\nend"

    def test_newline_prefix_frag(self):
        t = "\n%frag\nDefinition\n1 {18:27} end\n2 {0:17} end\nend\nend"
        pre, post = _split_post_blocks(t)
        assert post == "%frag\nDefinition\n1 {18:27} end\n2 {0:17} end\nend\nend"

    def test_mixed_pre_and_post(self):
        t = "%scf maxiter 500 end\n%frag\nDefinition\n1 {18:27} end\nend\nend"
        pre, post = _split_post_blocks(t)
        assert pre == "%scf maxiter 500 end"
        assert post == "%frag\nDefinition\n1 {18:27} end\nend\nend"

    def test_multiple_post_blocks(self):
        t = "%frag\n1 {0:5} end\nend\n%eprnmr\ngtensor 1\nend"
        pre, post = _split_post_blocks(t)
        assert pre == ""
        assert post == "%frag\n1 {0:5} end\nend\n%eprnmr\ngtensor 1\nend"

    def test_pre_non_post_block_and_post(self):
        t = "%output\n print[ P_Mulliken 1 ] end\n%frag\n1 {0:5} end\nend"
        pre, post = _split_post_blocks(t)
        assert pre == "%output\n print[ P_Mulliken 1 ] end"
        assert post == "%frag\n1 {0:5} end\nend"

    def test_empty(self):
        pre, post = _split_post_blocks("")
        assert pre == ""
        assert post == ""

    def test_no_post_blocks(self):
        t = "%scf maxiter 500 end\n%maxcore 8000"
        pre, post = _split_post_blocks(t)
        assert "%scf maxiter 500 end" in pre
        assert "%maxcore 8000" in pre
        assert post == ""

    def test_missing_newlines(self):
        t = "! optts %maxcore 4000 %eprnmr nuclei H = all {shift} end end \n! Freq \n%pal nprocs_group 4 end %scf MaxIter 200 end"
        pre, post = _split_post_blocks(t)
        assert "! optts" in pre
        assert "%maxcore 4000" in pre
        assert "! Freq" in pre
        assert "%pal nprocs_group 4 end" in pre
        assert "%scf MaxIter 200 end" in pre
        assert "%eprnmr nuclei H = all {shift} end end" in post
        # Ensure that no newlines are wrongly stripped internally
        assert pre.count("\n") >= 4
