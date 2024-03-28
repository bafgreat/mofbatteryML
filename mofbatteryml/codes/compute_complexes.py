import os
import glob
import argparse
import numpy as np
from ase.io import read
from mofbatteryml.io import filetyper
from mofbatteryml.energy import docker

# all_mofs = sorted(glob.glob('../data/selected_cifs/*cif'))
# path_to_monomer = '../data/metal_atoms'

all_mofs = sorted(glob.glob('../../battary_data/selected_cifs/*cif'))
path_to_monomer = '../../battary_data/metal_atoms'

def compute_xtb_for_complexes(list_of_hosts, list_of_monomers, results_folder, number_of_host=1, number_of_monomers=1, number_of_complexes=10):
    """
    A function to extract the energy of a given host-guest system.
    parameter
    ----------
    list_of_hosts: list
        List of host system files.
    list_of_monomers: list
        List of monomer files.
    results_folder: str
        Folder to store results.
    number_of_host: int, optional
        Number of host molecules.
    number_of_monomers: int, optional
        Number of monomer molecules.
    number_of_complexes: int, optional
        Number of complexes to generate.
    """
    seen = []
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)
    all_complexes = os.path.join(results_folder, 'complexes')
    if not os.path.exists(all_complexes):
        os.makedirs(all_complexes)

    json_energy_filename = os.path.join(results_folder, 'energy.json')

    if os.path.exists(json_energy_filename):
        json_energy_filename = filetyper.load_data(json_energy_filename)
        seen = list(json_energy_filename.keys())

    new_mol = {}
    new_energy = {}

    list_of_monomers = sorted([os.path.join(list_of_monomers, i) for i in os.listdir(list_of_monomers)])

    for host_system_file in list_of_hosts:
        print(host_system_file)
        for monomer_file in list_of_monomers:
            monomer = read(monomer_file)
            host_system = read(host_system_file)
            host_base_name = os.path.basename(host_system_file).split('.')[0]
            monomer_base_name = os.path.basename(monomer_file).split('.')[0]
            base_name = host_base_name + '_' + monomer_base_name
            if base_name not in seen:
                energy_dict, complex_molecules = docker.Dock(
                    host_system, monomer, number_of_host, number_of_monomers, number_of_complexes)
                new_mol[base_name] = complex_molecules
                # new_energy[base_name] = energy_dict
                mol_files = os.path.join(all_complexes, host_base_name+'.json')
                # filetyper.append_json(new_energy, json_energy_filename)
                filetyper.append_json_atom(new_mol, mol_files)
            else:
                print(f'{base_name} has already been computed')
    return

def run_xtb_for_complexes(all_mofs, path_to_monomer, result_folder, batch):
    if len(batch) == 1:
        list_of_mofs = [all_mofs[batch[0]]]
    elif len(batch) == 0:
        list_of_mofs = all_mofs
    else:
        list_of_mofs = all_mofs[batch[0]:batch[-1]+1]

    compute_xtb_for_complexes(list_of_mofs, path_to_monomer, results_folder=result_folder)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run work_flow function with optional verbose output')
    parser.add_argument('-rf', '--result_folder', type=str, default='Results_folder', help='Folder to store results')
    parser.add_argument('-bs', '--batch', type=int, nargs='+', help='batch size to run')
    args = parser.parse_args()
    run_xtb_for_complexes(all_mofs, path_to_monomer, result_folder=args.result_folder, batch=args.batch)
