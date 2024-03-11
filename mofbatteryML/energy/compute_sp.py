import os
from mofbatteryML.io import coords_library
from mofbatteryML.io import filetyper


def compute_xtb_energy(ase_atoms):
    """
    Function to compute the xtb energy
    parameter
    ----------
    ase_atoms: ase.Atoms object
    """
    base_dir = os.getcwd()
    result_folder = 'tmp_dir'
    if not os.path.exists(result_folder):
        os.makedirs(result_folder)
    os.chdir(result_folder)
    tmp_in = 'tmp_energy.gen'
    tmp_out = 'tmp_energy.out'
    coords_library.write_ase_atoms(ase_atoms, tmp_in)
    os.system(f'xtb --sp --gfn 2 --tblite --spinpol {tmp_in} > {tmp_out}')
    energy = read_xtb_energy(tmp_out)
    os.chdir(base_dir)
    # if os.path.exists(result_folder):
    #     shutil.rmtree(result_folder)
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

