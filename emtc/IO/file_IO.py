
'''
    QSKIT File IO
'''
from pathlib import Path

script_path = Path(__file__).resolve()
script_path = script_path.parent
current_path = Path.cwd()


from pathlib import Path
from typing import Tuple, List, Union, Dict
import numpy as np

def read_poscar(poscar_path: Union[str, Path]) -> Dict[str, Union[str, List[str], np.ndarray, List[int], float]]:
    """
    read VASP POSCAR/CONTCAR file
    
    parameters :
        poscar_path: POSCAR file path (str or Path objective)
    
    return :
        dict:
        - 'comment': comment line
        - 'scaling_factor': scaling factor (float)
        - 'lattice': 3x3 cell (numpy array)
        - 'element_symbols': element symbols list (List[str])
        - 'element_counts': element number list (List[int])
        - 'selective_dynamics': selective dynamics (bool)
        - 'coordinate_type': coordinate type ('Cartesian' or 'Direct')
        - 'positions': atom positions (numpy array, Nx3)
        - 'selective_flags': selective flags (None or Nx3 bool array)
    """
    poscar_path = Path(poscar_path)
    with poscar_path.open('r') as f:
        lines = [line.strip() for line in f.readlines() if line.strip()]
    
    data = {}
    data['comment'] = "QSKIT Structure"  # comment line
    
    data['scaling_factor'] = float(lines[1])
    
    lattice = []
    for line in lines[2:5]:
        lattice.append([float(x) for x in line.split()[:3]])
    data['lattice'] = np.array(lattice)
    
    elements_line = lines[5].split()
    counts_line = lines[6].split()
    
    # check element line or not 
    try:
        float(elements_line[0])
        data['element_symbols'] = [f'Element{i+1}' for i in range(len(counts_line))]
        data['element_counts'] = [int(x) for x in elements_line]
        remaining_lines = lines[6:]
    except ValueError:
        data['element_symbols'] = elements_line
        data['element_counts'] = [int(x) for x in counts_line]
        remaining_lines = lines[7:]
    
    # check selective dynamics or not 
    if 'Selective dynamics' in remaining_lines[0]:
        data['selective_dynamics'] = True
        coord_type_line = remaining_lines[1]
        positions_start = 2
    else:
        data['selective_dynamics'] = False
        coord_type_line = remaining_lines[0]
        positions_start = 1
    
    data['coordinate_type'] = 'Direct' if coord_type_line[0].lower() in 'd' else 'Cartesian'
    
    positions = []
    selective_flags = [] if data['selective_dynamics'] else None
    
    total_atoms = sum(data['element_counts'])
    pos_lines = remaining_lines[positions_start:positions_start + total_atoms]
    
    for line in pos_lines:
        parts = line.split()
        positions.append([float(x) for x in parts[:3]])
        if data['selective_dynamics']:
            flags = [s.upper() == 'T' for s in parts[3:6]]
            selective_flags.append(flags)
    
    data['positions'] = np.array(positions)
    if data['selective_dynamics']:
        data['selective_flags'] = np.array(selective_flags)
    else:
        data['selective_flags'] = None
    return data

def write_poscar(data, output_path):
    """writ to poscar"""
    with open(output_path, 'w') as f:
        f.write(f"{data['comment']}\n")
        f.write(f"  {data['scaling_factor']}\n")
        for vec in data['lattice']:
            f.write(f"  {vec[0]:.9f}  {vec[1]:.9f}  {vec[2]:.9f}\n")
        f.write("  " + "  ".join(data['element_symbols']) + "\n")
        f.write("  " + "  ".join(map(str, data['element_counts'])) + "\n")
        f.write(data['coordinate_type'] + "\n")
        for pos in data['positions']:
            f.write(f"  {pos[0]:.9f}  {pos[1]:.9f}  {pos[2]:.9f}\n")

if __name__ == "__main__":
    poscar_data = read_poscar('POSCAR')
    
    print("Comment:", poscar_data['comment'])
    print("Scaling factor:", poscar_data['scaling_factor'])
    print("Lattice vectors:\n", poscar_data['lattice'])
    print("Element symbols:", poscar_data['element_symbols'])
    print("Element counts:", poscar_data['element_counts'])
    print("Coordinate type:", poscar_data['coordinate_type'])
    print("Positions:\n", poscar_data['positions'])
    
    if poscar_data['selective_dynamics']:
        print("Selective dynamics flags:\n", poscar_data['selective_flags'])