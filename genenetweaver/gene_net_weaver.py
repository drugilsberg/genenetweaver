import logging
import os.path
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
import xml.dom.minidom
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__.split('.')[-1])
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

"""
This a wrapper for geneNetWeaver (gnw.jar), see http://gnw.sourceforge.net/
It implements
- extraction of a subnetwork (there are sample networks in the jar,
  see ch/epfl/lis/inetworks). Per default it extracts the network
  from e.coli. the subnetwork is saved as an xml (kinetic model, required
  in order to simulate it)
- simulation of the subnetwork (there are sample networks in the jar,
  see ch/epfl/lis/inetworks)
- if no output_path is defined an temporary directory will be created
  (removed at the end)
"""


class GeneNetWeaver:
    subnetwork_file = None
    network_prefix = None
    rna_suffix = '_multifactorial.tsv'
    protein_suffix = '_proteins_multifactorial.tsv'
    gold_standard_network_suffix = '_goldstandard.tsv'
    output_path = None
    gnw_jar = None
    default_settings = {}
    noises = {
        'normal': {'addNormalNoise': 1, 'addLognormalNoise': 0, 'addMicroarrayNoise': 0},
        'lognormal': {'addNormalNoise': 0, 'addLognormalNoise': 1, 'addMicroarrayNoise': 0},
        'microarray': {'addNormalNoise': 0, 'addLognormalNoise': 0, 'addMicroarrayNoise': 1}
    }
    simulators = {
        'ode': {'simulateODE': 1, 'simulateSDE': 0},
        'sde': {'simulateODE': 0, 'simulateSDE': 1},
    }
    temporary_directory = False  # set to true if created temporary directory, remove it at the end

    def __init__(self, gnw_jar, output_path=None):

        # check jar location
        if not os.path.exists(gnw_jar):
            logger.error('{} does not exist.'.format(gnw_jar))
            raise RuntimeError
        else:
            self.gnw_jar = gnw_jar

        # get output_path
        if output_path:
            if not os.path.exists(output_path):
                logger.error('Directory {} does not exist.'.format(output_path))
                raise RuntimeError
        else:
            logger.info("No output directory specified. Create temporary directory.")
            output_path = create_temporary_dir()
            self.temporary_directory = True
        logger.info('Files will be stored in {}'.format(output_path))
        self.output_path = output_path

        # get default settings
        settings_file = convert_settings_to_dict(gnw_jar, 'settings.txt', output_path)
        with open(settings_file, 'r') as f:
            self.default_settings = dict([parse_line(line.strip()) for line in f if re.search('^[a-z]', line)])
        self.default_settings.update(
            {
                'ssKnockouts': 0,
                'ssKnockdowns': 0,
                'ssMultifactorial': 1,
                'ssDREAM4TimeSeries': 0,
                'ssDualKnockouts': 0,
                'tsDREAM4TimeSeries': 0,
                'numTimeSeries': 0,
                'maxtTimeSeries': 1000,
                'dt': 50,
                'simulateODE': 1,  # only choose SDE or ODE otherwise filenameing is confusing
                'simulateSDE': 0,
                'addNormalNoise': 0,  # choose one of the 3 noise
                'addLognormalNoise': 0,
                'addMicroarrayNoise': 0,
                'normalizeAfterAddingNoise': 0
            })


    def get_number_of_species(self, network_file):
        # parses the xml file, looks for listOfSpecies tag and counts
        # (omits the species '_void_')
        number_of_species = len(
            xml.dom.minidom
            .parse(network_file)
            .getElementsByTagName('listOfSpecies')[0]
            .getElementsByTagName('species')
        ) - 1
        logger.debug('number of species: {}'.format(number_of_species))
        return number_of_species


    def extract_subnetwork(self, input_network_string, subnet_size=None):

        commandline_args = ['--extract',
                            '--random-seed',
                            '--rat-selection=20',  # random vs greedy-selection
                            '--num-subnets=1',
                            '--input-net-format=0',  # (0=TSV, 1=GML, 2=DOT, 3=TSV DREAM, 4=SBML)
                            '--keep-self-interactions',
                            '--output-net-format=4'  # need to use 4=SBML in order to simulate it
                            ]

        if subnet_size is not None:
            commandline_args.append('--subnet-size={}'.format(subnet_size))

        # get input_network_file from input_network_string
        input_network_file = get_input_net_path(input_network_string, self.output_path, self.gnw_jar)

        # extract subnetwork from input_network_file
        if input_network_file:
            call_string = 'java -jar {} --input-net {} --output-path {} {}'.format(self.gnw_jar,
                                                                                   input_network_file,
                                                                                   self.output_path,
                                                                                   ' '.join(commandline_args))
            logger.info('extract subnetwork from {} with size {}'.format(input_network_file, subnet_size))
            logger.debug(call_string)
            try:
                retcode = subprocess.call(call_string, shell=True)
            except Exception as e:
                logger.error('Run of gnw failed: {}'.format(call_string))
        else:
            logger.error('Cannot find {}'.format(self.input_network_string))

        # check if subnetwork_file was created
        self.network_prefix = '{}-1'.format(filename_to_basename_without_extension(input_network_file))
        subnetwork_file = '{}/{}.xml'.format(self.output_path, self.network_prefix)

        if os.path.exists(subnetwork_file):
            logger.debug('{} exists'.format(subnetwork_file))
            self.subnetwork_file = subnetwork_file
            return self.get_number_of_species(subnetwork_file)


    def transform_network(self, input_network_string):

        commandline_args = ['--transform',
                            '--input-net-format=0',  # (0=TSV, 1=GML, 2=DOT, 3=TSV DREAM, 4=SBML)
                            '--keep-self-interactions',
                            '--output-net-format=4'  # need to use 4=SBML in order to simulate it
                            ]

        # get input_network_file from input_network_string
        input_network_file = get_input_net_path(input_network_string, self.output_path, self.gnw_jar)

        # extract subnetwork from input_network_file
        if input_network_file:
            call_string = 'java -jar {} --input-net {} --output-path {} {}'.format(self.gnw_jar,
                                                                                   input_network_file,
                                                                                   self.output_path,
                                                                                   ' '.join(commandline_args))
            logger.info('transform network to xml')
            logger.debug(call_string)
            try:
                retcode = subprocess.call(call_string, shell=True)
            except Exception as e:
                logger.error('Run of gnw failed: {}'.format(call_string))
        else:
            logger.error('Cannot find {}'.format(self.input_network_string))

        # check if subnetwork_file was created
        self.network_prefix = filename_to_basename_without_extension(input_network_file)
        network_file = '{}/{}.xml'.format(self.output_path, self.network_prefix)
        logger.debug('expected location of network file: {}'.format(network_file))

        if os.path.exists(network_file):
            logger.debug('{} exists'.format(network_file))
            self.subnetwork_file = network_file
            return self.get_number_of_species(network_file)


    def simulate(self, **settings):

        if self.subnetwork_file:
            # combine default and specific settings and write it to file
            combined_settings = self.default_settings.copy()
            combined_settings.update(settings)
            settings_filename = '{}/{}'.format(self.output_path, 'settings.txt')
            create_settings_file(combined_settings, settings_filename)
            logger.info('settings file {} was created'.format(settings_filename))

            # simulate
            call_string = 'cd {} ; java -jar {} -c {} --input-net {} --simulate'.format(self.output_path, self.gnw_jar,
                                                                                        settings_filename,
                                                                                        self.subnetwork_file)
            logger.info('simulate network {}'.format(self.subnetwork_file))
            logger.debug(call_string)
            try:
                retcode = subprocess.call(call_string, shell=True)
            except Exception as e:
                logger.error('Run of gnw failed: {}'.format(call_string))

        else:
            logger.error('No network file. Run extract first')

    def get_simulation_settings(self, noise_type, simulation_type):
        settings = {}
        if noise_type:
            settings.update(self.noises[noise_type])
        settings.update(self.simulators[simulation_type])
        return settings

    def generate_data_set(self, number_of_samples, number_of_features=None, noise_type=None, simulation_type='ode',
                          input_network_string='ecoli_transcriptional_network_regulonDB_6_7.tsv',
                          remove_temporary_directory=True):

        # extract subnetwork if number_of_features is specified, otherwise only transform
        if number_of_features is None:
            number_of_features = self.transform_network(input_network_string)
        else:
            number_of_features = self.extract_subnetwork(input_network_string, subnet_size=number_of_features)

        rnas = pd.DataFrame()
        proteins = pd.DataFrame()

        number_of_simulation_rounds = int(np.ceil(number_of_samples / float(number_of_features)))
        logger.info('number of simulation rounds: {}'.format(number_of_simulation_rounds))

        for _ in range(number_of_simulation_rounds):
            # simulation
            self.simulate(**self.get_simulation_settings(noise_type, simulation_type))

            # read output files
            rna, protein = self.read_current_states()
            rnas = pd.concat([rnas, rna])
            proteins = pd.concat([proteins, protein])

        # get true network
        network = self.read_true_network()

        # remove temporary directory
        if self.temporary_directory and remove_temporary_directory:
            remove_dir(self.output_path)
        else:
            logger.info('Files are stored in {}'.format(self.output_path))

        # return
        fix = lambda df: filter_dataframe_index(df, number_of_samples)
        return network, fix(rnas), fix(proteins),

    def read_current_states(self):
        return pd.read_csv('{}/{}{}'.format(self.output_path, self.network_prefix, self.rna_suffix), sep='\t'), \
               pd.read_csv('{}/{}{}'.format(self.output_path, self.network_prefix, self.protein_suffix), sep='\t')

    def read_true_network(self):
        df = pd.read_csv('{}/{}{}'.format(self.output_path, self.network_prefix,
                                          self.gold_standard_network_suffix), sep='\t', header=None)
        df.columns = ['g1', 'g2', 'interacts']
        return df


def create_temporary_dir():
    try:
        temporary_directory = tempfile.mkdtemp()
        return temporary_directory
    except Exception as e:
        logger.error("Cannot create temporary directory. Provide writeable directory.")


def remove_dir(directory_name):
    try:
        shutil.rmtree(directory_name)
        logger.info('Removed directory {}.'.format(directory_name))
    except Exception as e:
        logger.error('Cannot remove directory {}.'.format(directory_name))


def get_input_net_path(input_net, output_path, zipfile_name):
    # check if
    # - file exists (if it is a path to a file)
    # - it is in the output_path
    # - if it is in the jar
    input_net_path = None
    if os.path.exists(input_net):
        input_net_path = input_net
    elif os.path.exists('{}/{}'.format(output_path, input_net)):
        input_net_path = '{}/{}'.format(output_path, input_net)
    else:
        full_filename = file_in_zipfile(zipfile_name, 'networks/{}'.format(input_net))
        if full_filename:
            z = zipfile.ZipFile(zipfile_name)
            z.extract(full_filename, output_path)
            input_net_path = '{}/{}'.format(output_path, full_filename)
    return input_net_path


def file_in_zipfile(zipfile_name, search_path):
    filepath_in_zipfile = None
    z = zipfile.ZipFile(zipfile_name)
    for file in z.infolist():
        if search_path in file.filename:
            filepath_in_zipfile = file.filename
    return filepath_in_zipfile


def extract_file_from_zipfile(zipfile_name, filename, output_path):
    z = zipfile.ZipFile(zipfile_name)
    z.extract(filename, output_path)
    return '{}/{}'.format(output_path, filename)


def filename_to_basename_without_extension(filename):
    return os.path.splitext(os.path.basename(filename))[0]


def convert_settings_to_dict(zipfile_name, filename, output_path):
    filename_in_zip = file_in_zipfile(zipfile_name, filename)
    file_path = extract_file_from_zipfile(zipfile_name, filename_in_zip, output_path)
    return file_path


def parse_line(line):
    key, value = tuple(line.replace(' ', '').split('='))
    try:
        return key, int(value)
    except:
        try:
            return key, float(value)
        except:
            return key, value


def create_settings_file(settings, filename):
    with open(filename, 'w') as f:
        for key, value in settings.items():
            f.write('{} = {}\n'.format(key, value))


def filter_dataframe_index(df, n):
    df_filtered = df[:n]
    df_filtered.index = range(n)
    return df_filtered
