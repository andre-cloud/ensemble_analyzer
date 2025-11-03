try:
    from src.IOsystem import mkdir
except ImportError as e:  # pragma: no cover
    print(e)
    from IOsystem import mkdir

import numpy as np
import random
from ase.atoms import Atoms

EH_TO_KCAL = 627.5096080305927


class Conformer:
    """
    Storing all the information on each conformer for all the parts of the protocol
    """

    def __init__(
        self,
        number: int,
        geom: np.ndarray,
        atoms: np.ndarray,
        raw=False,
    ) -> None:
        self.number = number
        self._initial_geometry = geom.copy()

        self.last_geometry = geom.copy()
        self.atoms = atoms.copy()
        self.energies = {}
        self.active = True
        self.color = "#%06x" % random.randint(0, 0xFFFFFF)

        self.cluster = None
        # IO
        self.folder = f"conf_{self.number}"
        if not raw:
            mkdir(self.folder)

    def get_ase_atoms(self, calc=None):
        """Return the atoms needed for the calculation

        :param calc: the type of calculator to use
        :type calc: ase.calculator
        :return: the ase instance ready to start the calculation
        :rtype: ase.Atoms
        """
        return Atoms(
            symbols="".join(list(self.atoms)),
            positions=self.last_geometry,
            calculator=calc,
        )

    @property
    def weight_mass(self):
        return np.sum(
            Atoms(
                symbols="".join(list(self.atoms)),
                positions=self.last_geometry,
            ).get_masses()
        )

    @property
    def rotatory(self):
        return self.energies[list(self.energies.keys())[-1]]["B"]

    @property
    def moment(self):
        return self.energies[list(self.energies.keys())[-1]]["m"]

    @property
    def get_energy(self):
        en = self._last_energy
        if en["G"] not in [np.nan, 0]:
            return en["G"]
        return en["E"]

    # def distance_matrix(self, exclude_H=False):
    #     if exclude_H:
    #         n = len(self.atoms[self.atoms != 'H'])
    #         dm = np.zeros((n, n))
    #         geo = self.last_geometry[self.atoms != 'H']
    #         for i, pos1 in enumerate(geo):
    #             for j, pos2 in enumerate(geo):
    #                 if i==j: continue
    #                 if dm[i,j]: continue
    #                 dist = np.linalg.norm(pos1 - pos2)
    #                 dm[i, j] = dist
    #                 dm[j, i] = dist
    #     else:
    #         dm = np.linalg.norm(self.last_geometry[:, None, :] - self.last_geometry[None, :, :], axis=-1)

    #     return dm

    def distance_matrix(self, include_H, geom=None):
        if geom:
            geo = geom
        else:
            geo = self.last_geometry

        if include_H:
            geo = np.array(geo)
        else:
            mask = self.atoms != "H"
            geo = np.array(geo)[mask]

        dm = np.linalg.norm(geo[:, None, :] - geo[None, :, :], axis=-1)
        return dm

    @property
    def _last_energy(self):
        if self.energies:
            return self.energies[list(self.energies.keys())[-1]]
        return {"E": 0, "G": None}

    def write_xyz(self):
        """Write the XYZ string to be stored in a file

        :return: the string in the XYZ formatting
        :rtype: str
        """
        if not self.active:
            return ""

        # Header
        header = f'{len(self.atoms)}\nCONFORMER {self.number} {"  {:.10f}".format(self._last_energy["G"]/EH_TO_KCAL) if self._last_energy["G"] else "  {:.10f}".format(self._last_energy["E"]/EH_TO_KCAL)}'

        # Atoms and positions
        atom_lines = [
            f"{a}  {x:14f}  {y:14f}  {z:14f}"
            for a, (x, y, z) in zip(self.atoms, self.last_geometry)
        ]

        # Combine the header and atom lines
        txt = "\n".join([header] + atom_lines)

        return txt

    def create_log(self, monitor_internals):
        """Generate all the information needed for the tabulation

        :return: a long tuple with all the information. (Number, E, G, B, Erel, Pop, Elapsed Time)
        :rtype: tuple
        """
        en = self._last_energy
        number, e, g_e, g, b, erel, time, pop, cluster = (
            self.number,
            en.get("E", float(0)),
            en.get("G-E", float(0)),
            en.get("G", float(0)),
            en.get("B", float(0)),
            en.get("Erel", float(0)),
            en.get("time"),
            en.get("Pop", float(0)),
            self.cluster,
        )

        monitor = None
        if len(monitor_internals) > 0:
            monitor = []
            atoms = Atoms(
                symbols="".join(list(self.atoms)),
                positions=self.last_geometry,
            )
            for internal in monitor_internals:
                if len(internal) == 2:
                    monitor.append(float(atoms.get_distance(*internal)))
                if len(internal) == 3:
                    monitor.append(float(atoms.get_angle(*internal)))
                if len(internal) == 4:
                    monitor.append(float(atoms.get_dihedral(*internal)))

        if g:
            g /= 627.51
        return number, e / 627.51, g_e, g, b, erel, pop, time, cluster, *monitor

    @staticmethod
    def load_raw(json):
        """Load raw configuration. Restart purposes.

        :param json: All the information of the conformer
        :type json: dict
        :return: the conformer instance
        :rtype: Conformer
        """
        a = Conformer(
            number=json["number"],
            geom=json["last_geometry"],
            atoms=json["atoms"],
            raw=True,
        )
        a.energies = json["energies"]
        a.active = json["active"]
        # a.cluster = json["cluster"]
        return a

    def __str__(self) -> str:
        return str(self.number)

    def __repr__(self) -> str:
        return str(self.number)

    # Functions needed for sorting the conformers' ensemble

    def __lt__(self, other):
        if not self.active:
            return 0 < other.get_energy
        return self.get_energy < other.get_energy

    def __gt__(self, other):
        if not self.active:
            return 0 > other.get_energy
        return self.get_energy > other.get_energy

    def __eq__(self, other):
        if not self.active:
            return 0 == other.get_energy
        return self.get_energy == other.get_energy
