import numpy as np

import ase
from ase.calculators.lammpslib import LAMMPSlib
from ase.calculators import lj, morse
from ase.optimize import LBFGS

def LJ(frame, sigma, epsilon):
    print(frame)

    frame.calc = lj.LennardJones(sigma=sigma, epsilon=epsilon)

    return store_potential(frame, "LJTempData")

def LAMMPS(potential, settings):
    print(settings['frame'])
    frame, parameters = settings['frame'], settings['parameters']
        
    atom_types = {}
    for i, atom in enumerate(settings['atoms']):
        atom_types[atom] = i + 1

    cmds = [ f"pair_style {potential} {settings['cutoff_radius']}" ]
    index = 0
    for x in range(int((((1 + 8 * len(parameters)) ** 0.5) - 1)/2)):
        for y in range(x + 1):
            params = " ".join(map(str, parameters[index])) if len(parameters[index]) > 1 else parameters[index][0]
            cmds.append(f"pair_coeff {x + 1} {y + 1} {params}")
            index += 1
    
    frame.calc = LAMMPSlib(lmpcmds=cmds, atom_types=atom_types, log_file='log.lammps')

    if settings['relax']:
        optimizer = LBFGS(frame)
        optimizer.run(fmax=settings['fmax'], steps=100)
        frame.calc = LAMMPSlib(lmpcmds=cmds, atom_types=atom_types, log_file='log.lammps')

    return store_potential(frame, potential)

def store_potential(frame, potential):

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