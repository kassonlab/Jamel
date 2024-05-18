from ChimeraGenerator import chimeracls,general_attr_set,fasta_creation
from Chimeragenesis.AccessiontoAlignment import no_gap_sequence_from_alignment
from ColorCoding import create_dictionary_from_alignment
import argparse
from sys import exit
from json import load
from os import path
from setup import alphafold_submission_for_chimera_container
parser = argparse.ArgumentParser(
    description='')
parser.add_argument('-i', '--jsoninput', dest='arg_jsons', required=False, type=str,
                    help='Comma seperated json config inputs. ex: constant_json,variant_json')
args = parser.parse_args()
class PointNamingArguments:
    placeholder: str = "*"
    '''a unique character to be replaced by fasta identifiers or provided nicknames for naming convention'''
    monomer_naming_convention: str = "*"
    '''A combination of words and placeholders to create a naming that will persist through all monomer chimera 
    related data files e.g. (fasta,pdb,plddt)'''
    slurm_naming: str
    '''Slurm naming convention. Placeholders in the convention will be replaced by numbers in ascending order as needed'''
    multimer_naming_convention: str
    chimera_naming_convention: str
    gromacs_slurm_dir: str
    fasta_directory: str
    plddt_directory: str
    pdb_directory: str
    alphafold_outputs_dir: str
    fasta_extension: str = '.fasta'
    '''File extension for fasta files'''
    plddt_extension: str = '.plddt'
    pdb_extension: str = '.pdb'
    gmx_setup_extension: str = '.setup'
    gmx_production_extension: str = '.prod'


class PointFastaArguments:
    fasta_toggles: dict
    '''A dictionary that holds optional operations within the fasta operation, they are turned on by 
    placing true/false'''
    Make_a_lst_of_created_fasta_files: bool
    '''Create a file with line separted list of all fastas created'''
    Create_an_alignment: bool
    constant_or_variant: bool
    '''Distinguishes whether the config json is for the constant part of the chimera or the variant protein that will be swapped with all protein in the protein list'''

    protein_list: str
    '''Two column list of accession number and their associated protein file_stems'''
    fasta_identifier: str
    '''Fasta label for the constant protein or the reference protein for the homologs to be iterated through'''
    # TODO add defaults like this
    muscle_command_for_msa: str = 'module load gcc/9.2.0 && module load muscle/3.8.31 && muscle'
    '''The start of the muscle command to be run to create the alignment.'''
    number_of_subunits: int
    sequence_of_interest: str
    '''This is the path to the fasta containing the sequence of either the constant protein you want spliced into the 
    variant proteins, or the section of the variants you want replaced with the constant protein. Each sequence of 
    interest should be in their own fasta file. All non-reference variant proteins will have the homologous sequence 
    from an alignment replaced'''
    full_reference_fasta: str
    '''The path to full-length protein of the constant protein or the reference for the variants'''
    msa_file: str
    '''Path to the multiple sequence alignment to be created or a user-generated msa'''
    email_for_accession: str
    msa_fasta: str
    '''Name of the new fasta that will contain all variant protein sequences to be aligned'''
    constant_fasta_for_alphafold: str
    '''Path to the fasta to be predicted by Alphafold'''
    fasta_file_list_name: str
    '''filename to hold the complete list of fastas created in the fasta operation'''



class PointMutationContainer:
    argument_dict: dict
    '''Dictionary containing all arguments for running the chimera script'''
    operation_toggles: dict
    '''Dictionary containing boolean toggles for all major operation:fasta creation, alphafold submission, 
    analysis of alphafold results, and gromacs simulation submission.'''
    run_fasta_operation: bool
    fasta_args: dict
    naming_args: dict
    alphafold_submission: bool
    run_analysis_operation: bool
    run_gromacs_operation: bool

    def __init__(self, shifted_json_file):
        with open(shifted_json_file, 'rb') as jfile:
            self.argument_dict = load(jfile)
        self.operation_toggles = self.argument_dict['operation_toggles']
        self.chimeras = ()

    def get_dict_args(self, dict_class, attr_name, arg_key):
        dict_args = self.argument_dict[arg_key]
        setattr(self, attr_name, general_attr_set(dict_class,dict_args))

    def add_chimera(self, chimera):
        self.chimeras += (chimera,)


container = PointMutationContainer(args.arg_jsons)
container.get_dict_args(PointFastaArguments,'fasta_args','fasta_arguments')
container.get_dict_args(PointFastaArguments,'fasta_args','fasta_arguments')
container.get_dict_args(PointNamingArguments,'naming_args','naming_arguments')
naming_args=container.naming_args
placeholder=naming_args.placeholder
fasta_toggles=container.fasta_args.fasta_toggles
aln=container.fasta_args.msa_file
subunits=container.fasta_args.subunits
mutate_to=container.fasta_args.mutate_to
seq_dict=create_dictionary_from_alignment(aln)
mutation_index=container.fasta_args.mutation_index

AMINO_ACIDS = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
if mutate_to.upper() not in AMINO_ACIDS:
    print('You entered an invalid one-letter amino acid code')
    exit()
chimera_container=()
for id,seq in seq_dict.items():
    chimera=chimeracls()
    chimera_container+=(chimera,)
    chimera.nickname=id
    chimera.seq=seq
    chimera.WT_stem = naming_args.multimer_naming_convention.replace(placeholder, chimera.nickname)
    # TODO should i fix them for this, this is only mutating wildtype but not chimera
    if mutate_to==seq[mutation_index]:
        chimera.mutated=False
    else:
        chimera.mutated = True
    # is this the broken version
    chimera.mutated_seq=no_gap_sequence_from_alignment(seq[:mutation_index]+mutate_to.upper()+seq[mutation_index+1:])
    chimera.mutated_stem=naming_args.mutated_naming_convention.replace(placeholder, chimera.nickname)


if container.operation_toggles['run_fasta_operation']:
    fasta_creation(path.join(naming_args.fasta_directory,chimera.mutated_stem+naming_args.fasta_extension), [(chimera.mutated_seq, subunits, chimera.nickname)])
    if fasta_toggles['make_native_fastas']:
        for chimera in container.chimeras:
            fasta_creation(chimera.WT_stem, [(chimera.seq, subunits, chimera.nickname)])

if container.operation_toggles['alphafold_submission']:
    mutated_list=
    alphafold_submission_for_chimera_container(container,)
    
