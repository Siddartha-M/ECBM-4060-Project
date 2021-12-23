#from functionLib import *
import functionLib
import os
import sys


class identifier(object):
    mode = None
    input_directory = None
    output_directory = None
    frequency_threshold = None
    sites_override = None
    reference = None
    quality = None  # a tuple of (missing_base_threshold, degenerate_base_threshold, minimal_length)

    # constructor for default configuration
    def __init__(self, function='all', input=functionLib.inputDirectory, output=functionLib.outDirectory, ref=functionLib.refsequence, sites=None,
                 freq=0.05):
        self.mode = function
        self.input_directory = input
        self.output_directory = output
        self.reference = ref
        self.sites_override = sites
        self.frequency_threshold = freq
        self.quality = (45, 120, 29000)
        print('Initialization finished.')
        print('Reconfiguring input and output directory is highly recommended.')

    def operate(self):
        print('To see all commands,type "help".')
        while True:
            arg = input(">>>: ")
            if arg == 'help':
                self.help()
            elif arg == 'exit':
                break
            elif arg == 'default':
                self.restore_default_param()
            elif arg == 'set all':
                self.set_all_param()
            elif arg == 'reference':
                self.set_reference()
            elif arg == 'input':
                self.set_input_directory()
            elif arg == 'output':
                self.set_output_directory()
            elif arg == 'frequency':
                self.set_freq()
            elif arg == 'sites':
                self.set_sites()
            elif arg == 'mode':
                self.set_mode()
            elif arg == 'start':
                functionLib.clear_text(functionLib.ploidy)
                print('It may take hours to execute. Do not terminate the process!!!')
                if self.mode == 'all':
                    functionLib.clear_directory(self.output_directory)
                    snpFile = functionLib.mutation_identification(self.input_directory, self.output_directory, quality=self.quality)
                    functionLib.haplotype_identification(snpFile, self.output_directory, sites=self.sites_override,
                                             frequency=self.frequency_threshold)
                    print('Done! Output can be found in', self.output_directory)
                elif self.mode == 'haplotype':
                    functionLib.clear_directory(self.output_directory)
                    count = 0
                    for f in os.listdir(self.input_directory):
                        if count == 0 and f.endswith('.tsv'):
                            self.input_directory = os.path.join(self.input_directory, f)
                            count = count + 1
                    if count == 0:
                        print('No input file is found. Exiting...')
                        sys.exit()
                    elif count > 1:
                        print('Warning: multiple input files are found. Only"' + f + '" is used.')
                    functionLib.haplotype_identification(self.input_directory, self.output_directory, sites=self.sites_override,
                                             frequency=self.frequency_threshold)
                    print('Done! Output can be found in', self.output_directory)
                elif self.mode == 'debug':
                    # TODO debug features to be implemented
                    flag = input("Enter debugging machine (Desktop/LVM/WLS): ")
                    if flag == 'Desktop':
                        self.input_directory = '/mnt/DATA/ECBM4060Project/testdata'
                        self.output_directory = '/mnt/DATA/ECBM4060Project/testout'
                        functionLib.clear_directory(self.output_directory)
                        snpFile = functionLib.mutation_identification(self.input_directory, self.output_directory,
                                                          quality=self.quality)
                        functionLib.haplotype_identification(snpFile, self.output_directory, sites=self.sites_override,
                                                 frequency=self.frequency_threshold)
                    elif flag == 'LVM':
                        self.input_directory = '/home/zhang/ECBM4060Project/testdata'
                        self.output_directory = '/home/zhang/ECBM4060Project/testout'
                        functionLib.clear_directory(self.output_directory)
                        snpFile = functionLib.mutation_identification(self.input_directory, self.output_directory,
                                                          quality=self.quality)
                        functionLib.haplotype_identification(snpFile, self.output_directory, sites=self.sites_override,
                                                 frequency=self.frequency_threshold)
                    elif flag == 'WLS':
                        # TODO implement WLS2 features
                        pass
                    else:
                        print('Not implemented yet.')
                    print('Done! Output can be found in', self.output_directory)
            elif arg == 'quality':
                self.set_quality_filter()
            elif arg == 'current':
                self.check_param()
            else:
                print('Unrecognized command! To see all commands, type "help"')
        print('Exiting...')
        sys.exit()

    def set_all_param(self):
        self.set_reference()
        self.set_input_directory()
        self.set_mode()
        self.set_quality_filter()
        self.set_output_directory()
        self.set_freq()
        self.set_sites()
        print('All parameters has been set!')

    def check_param(self):
        print('Reference sequence: ', self.reference)
        print('Input files location: ', self.input_directory)
        print('Output directory: ', self.output_directory)
        print('Operation mode: ', self.mode)
        print('Quality filter (missing bases, degenerate bases, minimal length): ' + str(self.quality))
        print('Mutation frequency threshold: ', self.frequency_threshold)
        print('Sites of special interest: ', self.sites_override)

    def restore_default_param(self):
        self.reference = functionLib.refsequence
        self.output_directory = functionLib.outDirectory
        self.input_directory = functionLib.inputDirectory
        self.mode = 'all'
        self.frequency_threshold = 0.05
        self.sites_override = None
        self.quality = (15, 50, 29000)
        print('All parameters restored to default!')

    def set_reference(self):
        ref = input('Enter the path to reference sequence file: ')
        if os.path.isfile(ref):
            self.reference = ref
            print('Success! Reference sequence set to:', ref)
        else:
            self.reference = functionLib.refsequence
            print('Invalid input! Reference restored to default.')

    def set_output_directory(self):
        out = input('Enter the path to output directory: ')
        if os.path.isdir(out):
            self.output_directory = out
            print('Success! Output directory set to:', out)
        else:
            self.output_directory = functionLib.outDirectory
            print('Invalid input! Output directory restored to default.')

    def set_input_directory(self):
        inputdir = input('Enter the path to folder where input sequences are located: ')
        if os.path.isdir(inputdir):
            self.input_directory = inputdir
            print('Success! Input directory set to:', inputdir)
        else:
            self.input_directory = functionLib.inputDirectory
            print('Invalid input! Input directory restored to default.')

    def set_mode(self):
        operation_mode = input('Enter mode selection (all/haplotype/debug): ')
        if operation_mode in {'all', 'haplotype', 'debug'}:
            self.mode = operation_mode
            print('Success! Operation mode set to:', operation_mode)
        else:
            self.mode = 'all'
            print('Unrecognized mode. Operation mode reset to all. ')

    def set_freq(self):
        freq = input('Enter frequency threshold: ')
        if functionLib.is_frequency(freq):
            self.frequency_threshold = float(freq)
            print('Success! Mutation frequency threshold set to:', freq)
        else:
            self.frequency_threshold = 0.05
            print('Invalid input! Frequency threshold restored to default.')

    def set_sites(self):
        sites = input('Enter the sites of interest (separated by space): ')
        self.sites_override = list()
        try:
            for site in sites.strip().split():
                self.sites_override.append(int(site))
            self.sites_override.sort()
            if len(self.sites_override) < 1:
                self.sites_override = None
            elif len(self.sites_override) == 1:
                raise self.SiteError('Sites number is 1')
            else:
                if (self.sites_override[0] < 0) or (self.sites_override[-1] > 29903):
                    raise self.SiteError('Sites out of range')
                    print('Sites out of range.')
                    self.sites_override = None
            print('Sites set to:', sites)
        except ValueError:
            self.sites_override = None
            print('Invalid input! No sites of interest is taken.')
        except self.SiteError:
            self.sites_override = None
            print('Invalid sites! No sites of interest is taken.')

    def set_quality_filter(self):
        mb = input('Enter the maximum missing bases(suggested value 15): ')
        db = input('Enter the maximum degenerate bases (suggested value 50): ')
        ml = input('Enter the minimal sequence length (suggested value 29000): ')
        try:
            self.quality = (int(mb), int(db), int(ml))
            print('Success! Quality filter set to: missing bases <', mb + ',', 'degenerate bases <', db + ',',
                  'length >', ml + '.')
        except ValueError:
            print('Invalid input! Suggested values used for quality filter.')
            self.quality = (15, 50, 29000)

    @staticmethod
    def help():
        print('***************************HELP***************************')
        print('*          To execute analysis, type "start"             *')
        print('*      To display current settings, type "current"       *')
        print('*        To set all parameters, type "set all".          *')
        print('*   To reset all parameters to default, type "default".  *')
        print('*      To set reference sequence, type "reference"       *')
        print('*        To set input directory, type "input".           *')
        print('*       To set output directory, type "output"           *')
        print('*         To set operation mode, type "mode".            *')
        print('*     To set sequence quality filter, type "quality"     *')
        print('* To set mutation frequency threshold, type "frequency". *')
        print('*   To indicate sites of special interest, type "sites"  *')
        print('*           To exit the program, type "exit"             *')
        print('**********************************************************')

    class SiteError(Exception):
        pass
