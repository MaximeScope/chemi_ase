import numpy as np

import ase
from ase.calculators.lammpslib import LAMMPSlib
from ase.calculators import lj, morse

def LJ(frame, sigma, epsilon):
    print(frame)

    frame.calc = lj.LennardJones(sigma=sigma, epsilon=epsilon)

    data = frame.get_potential_energies()
    globalData = frame.get_potential_energy()

    return {
        "LJTempData": {
            "target": "atom",
            "values": data,
            "units": "eV",
        },
        "LJTempData" + "_global": {
            "target": "structure",
            "values": np.array([globalData]),
            "units": "eV",
        },
    }

def LAMMPS(frame, potential, cutoff_radius, parameters, atoms):
    print(frame)
        
    atom_types = {}
    for i, atom in enumerate(atoms):
        atom_types[atom] = i + 1

    cmds = [ f"pair_style {potential} {cutoff_radius}" ]
    index = 0
    for x in range(int((((1 + 8 * len(parameters)) ** 0.5) - 1)/2)):
        for y in range(x + 1):
            params = " ".join(map(str, parameters[index])) if len(parameters[index]) > 1 else parameters[index][0]
            cmds.append(f"pair_coeff {x + 1} {y + 1} {params}")
            print(f"pair_coeff {x + 1} {y + 1} {params}")
            index += 1
    
    print(cmds)
    print(atom_types)
    
    frame.calc = LAMMPSlib(lmpcmds=cmds, atom_types=atom_types, log_file='log.lammps')

    data = frame.get_potential_energies()
    globalData = frame.get_potential_energy()

    return {
        potential: {
            "target": "atom",
            "values": data,
            "units": "eV",
        },
        potential + "_global": {
            "target": "structure",
            "values": np.array([globalData]),
            "units": "eV",
        },
    }