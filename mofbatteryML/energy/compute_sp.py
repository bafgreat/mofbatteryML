import os
import shutil
import random
import subprocess
import tempfile
from mofbatteryml.io import coords_library
from mofbatteryml.io import filetyper


def compute_xtb_energy(ase_atoms):
    """
    Function to compute the xtb energy
    parameter
    ----------
    ase_atoms: ase.Atoms object
    """
    base_dir = os.getcwd()
    tmp_dir = 'tmp_folder'
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    os.chdir(tmp_dir)

    tmp_in = 'tmp_energy.gen'
    tmp_out = 'tmp_energy.out'

    # with open(tmp_in, 'w') as f:
    #     coords_library.write_ase_atoms(ase_atoms, f)
    ase_atoms.write(tmp_in)
    command = f'xtb --sp --gfn 2 --tblite --spinpol {tmp_in} > {tmp_out}'
    os.system(command)
    energy = read_xtb_energy(tmp_out)
    # try:
    #     completed_process = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)
    #     energy = read_xtb_energy(tmp_out)
    # except subprocess.CalledProcessError as e:
    #     print(f"Error executing xtb: {e.stderr}")
    #     return None

    os.chdir(base_dir)
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

    return energy



def read_xtb_energy(filename):
    """
    Function to read the xtb energy
    parameter
    ----------
    filename: string
    """
    contents = filetyper.get_contents(filename)
    energy = {}
    for line in contents:
        if 'TOTAL ENERGY' in line:
            data = line.split()
            energy['energy_kcal_mol'] = float(data[3])*627.5095
        if 'HOMO-LUMO GAP' in line:
            data = line.split()
            energy['homo_lumo_ev'] = float(data[3])
    return energy


def compute_energy_of_atom(folder_of_atoms, output_folder):
    """
    Function to compute the xtb energy
    parameter
    ----------
    ase_atoms: ase.Atoms object
    """
    all_energy = {}
    json_filename = f'{output_folder}/energy_of_atoms.json'
    for filename in folder_of_atoms:
        basename = filename.split('/')[-1].split('_')[0]
        ase_atoms = coords_library.read_and_return_ase_atoms(filename)
        energy = compute_xtb_energy(ase_atoms)
        all_energy[basename] = energy
        filetyper.append_json(all_energy, json_filename)


# def compute_xtb

# folder_of_atoms = glob.glob("../data/metal_atoms/*.xyz")
# output_folder = '../data/json_files'
# compute_energy_of_atom(folder_of_atoms, output_folder)
# compute_xtb_energy("../data/metal_atoms/Al_atom.xyz")
