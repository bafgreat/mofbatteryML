import os
import glob
import argparse
from ase.io import read
from mofbatteryML.io import filetyper
from mofbatteryML.energy import docker

def extract_energy_molecules_from_file(host_system_file, monomer_file, number_of_host, number_of_monomers, number_of_complexes, results_folder):
    """
    A fucntion to extract the energy of a given host-guest system.
    parameter
    ----------
    host_folder: str
    monomer_folder: str
    """
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)
    json_mol_filename = os.path.join(results_folder, 'complexes.json')
    json_energy_filename = os.path.join(results_folder, 'energy.json')
    new_mol = {}
    new_energy = {}
    monomer = read(monomer_file)
    host_system = read(host_system_file)
    host_base_name = os.path.basename(host_system_file).split('.')[0]
    monomer_base_name = os.path.basename(monomer_file).split('.')[0]
    base_name = host_base_name + '_' + monomer_base_name
    energy_dict, complex_molecules = docker. Dock(host_system, monomer, number_of_host, number_of_monomers, number_of_complexes)
    new_mol[ base_name]=complex_molecules
    new_energy[base_name] = energy_dict
    filetyper.append_json(new_energy, json_energy_filename)
    filetyper.append_json_atom(new_mol, json_mol_filename)
    return


def main():
    '''
    Command line interface for computing docker
    '''
    parser = argparse.ArgumentParser(
        description='Run work_flow function with optional verbose output')
    parser.add_argument('host_file', type=str,
                        help='A file containing host systems')

    parser.add_argument('monomer_file', type=str,
                        help='A file containing guest systems')

    parser.add_argument('-nh', '--number_of_host', type=int,
                        default=1, help='The number of host systems')

    parser.add_argument('-nm', '--number_of_monomers', type=int,
                        default=1, help='The number of monomer systems')

    parser.add_argument('-nc', '--number_of_complexes', type=int,
                        default=1, help='The number of monomer systems')

    parser.add_argument('-r', '--results_folder', type=str,
                        default='results_folders', help='directory to save output files')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='print verbose output')
    args = parser.parse_args()


    extract_energy_molecules_from_folder(args.host_file, args.monomer_folder, args.number_of_host, args.number_of_monomers, args.number_of_complexes, args.results_folder)
