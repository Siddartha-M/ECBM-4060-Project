import identifier
import sys
import platform

identifier.functionLib.welcome()
print('')
if platform.system() != 'Linux':
    print('This project must be executed on Linux!')
    print('Exiting')
    sys.exit()
missing_dependencies = identifier.functionLib.check_dependencies(['bowtie2', 'vcftools', 'samtools', 'bcftools', 'java'])
if len(missing_dependencies) > 0:
    print('Required dependencies: Java 10, Bowtie2, VCFTools, SAMTools, BCFTools are needed for proper execution!')
    missing = 'Missing dependencies:'
    for dep in missing_dependencies:
        missing = missing + ' ' + dep
    print(missing)
    print('Exiting')
    sys.exit()
else:
    ID = identifier.identifier()
    ID.operate()
