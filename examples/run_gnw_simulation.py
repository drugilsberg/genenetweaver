#! /usr/bin/env python3
# -*- coding: utf-8 -*-

# Import
import os
from genenetweaver.gene_net_weaver import GeneNetWeaver
import numpy as np
import argparse


def argument_parser():
    parser = argparse.ArgumentParser(
        description='Run GeneNetWeaver (GNW) to simulate gene expression data '
                    'from which gene regulatory networks can be inferred. '
                    'Please provide the following arguments:')

    parser.add_argument('--path_to_jar',
                        required=True,
                        help='Enter the path to the GNW jar.')
    parser.add_argument('--network',
                        required=True,
                        type=str,
                        help=(
                            'Choose between "e.coli" and "yeast".' +
                            'If you prefer to use a custom network,' +
                            'please provide the path.'
                        ))
    parser.add_argument('--path_to_network',
                        help='Enter the path to the custom network file.',
                        default=None)
    parser.add_argument('--n_samples',
                        required=True,
                        type=int,
                        help='Enter how many samples GNW should simulate.')
    parser.add_argument('--n_entities',
                        required=True,
                        type=int,
                        help='Enter the size of the simulated network.')
    parser.add_argument('--output_directory',
                        required=True,
                        help='Please enter the path to the output directory.')
    parser.add_argument('--noise_type',
                        help=(
                            'Choose between: "normal",' +
                            '"lognormal" or "microarray".'
                        ),
                        default='normal')
    parser.add_argument('--simulation_type',
                        help='Choose between: "ode" and "sde"',
                        default='ode')

    return parser.parse_args()


def return_network_string(name, path):
    if name == 'yeast':
        # Returns Yeast network provided by GeneNetWeaver
        return 'yeast_transcriptional_network_Balaji2006.tsv'
    elif name == 'e.coli':
        # Returns E.coli network provided by GeneNetWeaver
        return 'ecoli_transcriptional_network_regulonDB_6_7.tsv'
    else:
        # Upload example network
        # TODO: Write function to include custom network into file.
        return path


def simulate(n_samples, n_entities, gnw_jar, noise, simulation_type,
             input_network, outdir):
    # Simulation
    print('Simulation of network with {} samples and {} entities.'
          .format(n_samples, n_entities))

    gnw = GeneNetWeaver(gnw_jar=gnw_jar, output_path=outdir)

    network, rna, protein = gnw.generate_data_set(
        n_samples,
        n_entities,
        noise_type=noise,
        simulation_type=simulation_type,
        input_network_string=input_network,
        remove_temporary_directory=True
    )

    # Write files
    directory = outdir + 'simulations/sim_net_{}samples_{}entities'.format(
        n_samples, n_entities
    )

    if not os.path.exists(directory):
        os.makedirs(directory)

    print('Write simulation results to {}.'.format(directory))
    with open(directory + '/true_network.tsv', 'w') as network_file:
        network.to_csv(network_file, sep='\t', index=True, header=True)

    with open(directory + '/rna_data.tsv', 'w') as rna_file:
        rna.to_csv(rna_file, sep='\t', index=True, header=True)

    with open(directory + '/protein_data.tsv', 'w') as protein_file:
        protein.to_csv(protein_file, sep='\t', index=True, header=True)


if __name__ == '__main__':
    # Get arguments
    args = argument_parser()

    # Simulation settings
    n_samples = args.n_samples
    n_entities = args.n_entities
    simulation_type = args.simulation_type
    noise_type = args.noise_type
    gnw_jar = os.path.abspath(args.path_to_jar)
    if not os.path.exists(gnw_jar):
        error = '{} does not exist.'.format(gnw_jar)
        raise OSError(error)
    input_network_string = return_network_string(
        args.network, args.path_to_network
    )

    # Run simulation
    simulate(n_samples, n_entities, gnw_jar, noise_type, simulation_type,
             input_network_string, args.output_directory)
