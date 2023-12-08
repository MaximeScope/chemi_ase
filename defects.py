import os
import numpy as np
import ase
from ase.build import bulk
from ase import Atoms

def vacancies(atom, structure, aseLatticeData, n_defects=1):
    # Duplicate the XYZ file
    os.system(f'cp {atom}/XYZ/{atom}_{structure}.xyz {atom}/XYZ/{atom}_{structure}_vacancy.xyz')

    # Generate random defects in the XYZ file:
    n_atoms = len(aseLatticeData)
    for i in range(n_defects):
        defect = np.random.randint(1, n_atoms) + 2
        # Remove the ith line of the XYZ file corresponding to the defect:
        with open(f'{atom}/XYZ/{atom}_{structure}_vacancy.xyz', 'r') as f:
            lines = f.readlines()
        with open(f'{atom}/XYZ/{atom}_{structure}_vacancy.xyz', 'w') as f:
            for j, line in enumerate(lines):
                if j == 0:
                    f.write(f'{n_atoms - n_defects}\n')
                elif j != defect:
                    f.write(line)

        print(f'Removed atom {defect} from {atom}_{structure}_vacancy.xyz')

def dislocations(atom, structure, aseLatticeData, threshold=[-10., 5., -10.], factor=[0., 0., 0.5]):
    """
    Generate a crystal structure with a slip dislocation.

    Parameters:
        atom (str): The atom name.
        structure (str): The structure type (fcc, bcc, and so on).
        aseLatticeData (Atoms): The lattice structure.
        line_direction (str): The line direction along which the dislocation will occur ('x', 'y', or 'z').
        line_position (tuple): The position of the dislocation line in fractional coordinates.

    Returns:
        Atoms: Crystal structure with a dislocation.
    """
    # Create a copy of the bulk structure
    dislocated_structure = aseLatticeData.copy()

    # Create a dislocation line
    for i in dislocated_structure:
        pos = i.position
        if (pos > np.array(threshold)).all():
            print(f'Atom {i.index} is above the threshold')
            new_pos = np.array(factor) + pos
            i.position = new_pos
        else:
            print(f'Atom {i.index} is below the threshold')

    # Write the new structure to the XYZ file
    dislocated_structure.write(f'{atom}/XYZ/{atom}_{structure}_dislocation.xyz')

def interstitials(atom, structure, aseLatticeData, interstitial):
    # Create a copy of the bulk structure
    data_copy = aseLatticeData.copy()

    # randomly select the defect position
    random_index = np.random.randint(len(data_copy))
    base_position = data_copy[random_index].position
    # faire une génération uniforme et pas gaussienne
    # random offset in the range [-0.1, 0.1] autour de
    # la position entre deux atomes les plus proches voisins
    offset = np.random.rand(3) * 0.2  # Random offset in the range [-0.1, 0.1]
    position = base_position + offset

    # Create a new atom with the same symbol and the new position
    new_atom = Atoms(interstitial, positions=[position])

    # Add the new atom to the existing structure
    data_copy += new_atom

    # Write the new structure to the XYZ file
    data_copy.write(f'{atom}/XYZ/{atom}_{structure}_interstitial.xyz')