import os
import pathlib
from json import load
from pathlib import Path

import AccessiontoAlignment, ESM, ChimeraGenerator, Analysis
from setup import create_alphafold_slurm


class HomomerNamingArguments:
    file_stem_placeholder: str
    WT_convention: str
    chimera_convention: str
    output_directory: str

    # All non-nested keys from the json are set as attributes to be called later
    def __init__(self, naming_dict):
        for key, value in naming_dict.items():
            setattr(self, key, value)


class HomomerFastaArguments:
    Make_a_lst_of_created_fasta_files: bool
    '''Create a file with line separted list of all fastas created'''
    Create_an_alignment: bool
    number_of_subunits: int
    sequence_of_interest: str
    '''This is the path to the fasta containing the sequence of either the constant protein you want spliced into the 
    variant proteins, or the section of the variants you want replaced with the constant protein. Each sequence of 
    interest should be in their own fasta file. All non-reference variant proteins will have the Homomer sequence 
    from an alignment replaced'''
    base_identifier: str
    msa_file_name: str
    '''Path to the multiple sequence alignment to be created or a user-generated msa'''
    email_for_accession: str
    fasta_file_list_name: str
    '''filename to hold the complete list of fastas created in the fasta operation'''

    def __init__(self, fasta_dict):
        for key, value in fasta_dict.items():
            setattr(self, key, value)


class HomomerSubmissionArguments:
    submission_toggles: dict
    '''A dictionary that holds optional operations within the submission operation, they are generally turned on by 
    placing true within the quotes of their value'''
    create_slurms: bool
    sbatch_slurms: bool
    stragglers_or_custom_or_all: bool = 'stragglers'
    '''If stragglers is placed in the quotes it will check to see if an alphafold prediction is done for all created 
    chimeras and making a new the list of incomplete predictions or 'stragglers' to either create a file with their 
    fastas or put them into a slurm to be submitted again. If custom is placed, all chimeras within the 
    custom_list_to_run file will be submitted as slurm jobs. If all is selected all chimeras will be run 
    indiscriminately'''
    create_file_of_stragglers: bool
    slurm_naming: str
    custom_list_to_run: str
    '''Path to file that for list that either be created by toggling on create_file_of_stragglers, or user-generated that should contain 
    line separated fastas of proteins you want predicted'''
    proteins_per_slurm: int
    template_slurm: str
    '''Path to the file that has the template settings for alphafold slurm jobs'''
    alphafold_shell_script: str
    '''Path to the shell script that runs alphafold, it must take a list of comma-separated fastas 
    and an output as its arguments. It's recommended you use the on provided in github'''

    def __init__(self, fasta_dict):
        for key, value in fasta_dict.items():
            setattr(self, key, value)


class HomomerAnalysisArguments:
    analysis_toggles: dict
    make_plddts: bool
    make_pdbs: bool
    make_emboss_files: bool

    emboss_command: str
    analysis_output_csv: str
    '''Path of the file you want created that will contain all analysis created.'''
    column_names: dict
    '''A dictionary containing a multiple keys and their list value that each correspond with a type of data that is 
    output by the script, a user-selected title for the data column and whether that data is included in the final 
    output'''
    file_stem: list
    '''An array [true,"Protein"] that contains a boolean that will control whether the data column containing 
    filestems/nicknames created from the naming convention previously established will be included in the data csv.  
    The second index controls the column name in the csv'''
    overall_native_stability: list
    '''See file_stem docstring'''
    relative_stability: list
    '''See file_stem docstring'''
    overall_chimera_stability: list
    '''See file_stem docstring'''

    def __init__(self, fasta_dict):
        for key, value in fasta_dict.items():
            setattr(self, key, value)


class HomomerChimeraArgs:
    argument_dict: dict
    '''Dictionary containing all arguments for running the chimera script'''
    operation_toggles: dict
    '''Dictionary containing boolean toggles for all major operation:fasta creation, alphafold submission, 
    analysis of alphafold results, and gromacs simulation submission.'''
    run_fasta_operation: bool
    alphafold_submission: bool
    run_analysis_operation: bool
    run_gromacs_operation: bool
    fasta_args: HomomerFastaArguments
    naming_args: HomomerNamingArguments
    submission_args: HomomerSubmissionArguments
    analysis_args: HomomerAnalysisArguments
    all_fastas: list[str]
    fasta_dir: pathlib.Path

    def __init__(self, arg_json):
        with open(arg_json, 'rb') as jfile:
            self.argument_dict = load(jfile)
        self.operation_toggles = self.argument_dict['operation_toggles']
        self.get_dict_args(HomomerFastaArguments, 'fasta_args', 'fasta_arguments')
        self.get_dict_args(HomomerNamingArguments, 'naming_args', 'naming_arguments')
        self.get_dict_args(HomomerSubmissionArguments, 'submission_args', 'submission_arguments')
        self.seq_df = ESM.SequenceDataframe(self.fasta_args.msa_file_name)
        self.base_splice = AccessiontoAlignment.extract_seq_from_fasta(self.fasta_args.sequence_of_interest)
        self.seq_df['chi_seq'] = self.seq_df.index.map(
            lambda label: ChimeraGenerator.get_chimera_sequence(self.fasta_args.msa_file_name,
                                                                self.fasta_args.base_identifier, label,
                                                                self.base_splice))
        self.seq_df['file_stem'] = self.seq_df.index.map(
            lambda label: self.naming_args.WT_convention.replace(self.naming_args.file_stem_placeholder, label))
        self.seq_df['chi_file_stem'] = self.seq_df.index.map(
            lambda label: self.naming_args.chimera_convention.replace(self.naming_args.file_stem_placeholder, label))
        self.fasta_dir = Path(self.naming_args.output_directory).joinpath('Fasta')
        self.pdb_dir = Path(self.naming_args.output_directory).joinpath('PDB')
        self.alphafold_dir = Path(self.naming_args.output_directory).joinpath('AlphaFold')
        fastas = self.seq_df['file_stem'].map(
            lambda stem: self.fasta_dir.joinpath(stem).with_suffix('.fa').__str__()).to_list()
        chi_fastas = self.seq_df['chi_file_stem'].map(
            lambda stem: self.fasta_dir.joinpath(stem).with_suffix('.fa').__str__()).to_list()
        self.all_fastas = fastas + chi_fastas

    def get_dict_args(self, dict_class, attr_name, arg_key):
        dict_args = self.argument_dict[arg_key]
        setattr(self, attr_name, dict_class(dict_args))

    def fasta_operations(self):
        num_of_subunits = self.fasta_args.number_of_subunits
        self.naming_args.output_directory = Path(self.naming_args.output_directory)
        fasta_dir = Path(self.naming_args.output_directory).joinpath('Fasta')
        os.makedirs(fasta_dir, exist_ok=True)
        for label, data in self.seq_df.iterrows():
            AccessiontoAlignment.fasta_creation(fasta_dir.joinpath(data['file_stem']).with_suffix('.fa'),
                                                AccessiontoAlignment.create_seq_records(label, data['sequence'],
                                                                                        subunit_count=num_of_subunits))
            AccessiontoAlignment.fasta_creation(fasta_dir.joinpath(data['chi_file_stem']).with_suffix('.fa'),
                                                AccessiontoAlignment.create_seq_records(label, data['chi_seq'],
                                                                                        subunit_count=num_of_subunits))

    def alphafold_submission(self):
        # fasta list is a separate input for control
        fasta_to_run = ()
        placeholder = self.naming_args.file_stem_placeholder
        submission_toggles = self.submission_args.submission_toggles
        proteins_per_slurm = self.submission_args.proteins_per_slurm
        naming_convention = self.submission_args.slurm_naming
        # Loops through all fastas created and checks if they are complete by looking for their ranking_debug file
        if submission_toggles['stragglers_or_custom_or_all'] == 'stragglers':
            for fasta in self.all_fastas:
                if not os.path.exists(Path(self.alphafold_dir).joinpath(Path(fasta).stem, '/ranking_debug.json')):
                    fasta_to_run += (fasta,)
            # Puts all fastas in a line separated file specified by custom_list_to_run
            if submission_toggles['create_file_of_stragglers']:
                with open(self.submission_args.custom_list_to_run, 'w') as run_list:
                    run_list.write('\n'.join(fasta for fasta in fasta_to_run))
            # if all of them are complete and fasta_to_run is empty then all slurm actions are toggled off
            if not fasta_to_run:
                self.submission_args.submission_toggles['create_slurms'] = self.submission_args.submission_toggles[
                    'sbatch slurms'] = False
                print('All Done')

        elif submission_toggles['stragglers_or_custom_or_all'] == 'custom':
            with open(self.submission_args.custom_list_to_run, 'r') as run_list:
                fasta_to_run = [x.split()[0] for x in run_list.readlines()]

        elif submission_toggles['stragglers_or_custom_or_all'] == 'all':
            fasta_to_run = self.all_fastas
        for slurm_index, file_index in enumerate(range(0, len(fasta_to_run), proteins_per_slurm)):
            current_slurm = naming_convention.replace(placeholder, str(slurm_index))
            if self.submission_args.submission_toggles['create_slurms']:
                create_alphafold_slurm(fasta_to_run[file_index:file_index + proteins_per_slurm], current_slurm,
                                       self.alphafold_dir, self.fasta_args.number_of_subunits)
            if self.submission_args.submission_toggles['sbatch_slurms']:
                os.system(f'sbatch {current_slurm}')

    def analysis_operations(self):
        Analysis.generate_alphafold_files(self.alphafold_dir, self.naming_args.output_directory)
        msa_dict = AccessiontoAlignment.create_dictionary_from_alignment(self.fasta_args.msa_file_name)
        base_plddt = Analysis.get_plddt_dict_from_pdb(
            self.pdb_dir.joinpath(self.fasta_args.base_identifier).with_suffix('.pdb'))
        print(base_plddt)
        # master_plddt_dict = {label: AccessiontoAlignment.map_plddt_to_aln(self.seq_df.get_aln(label), tuple(
        #
        #                      for label in msa_dict.keys()}
        # for label, data in self.seq_df.iterrows():

        # chimera.overall_native_stability = Analysis.overall_confidence(chimera.native_plddt[chimera.native_seq])
        # chimera.overall_chimera_stability = Analysis.overall_confidence(chimera.chi_plddt[chimera.chi_seq])
        # chimera.sequence_similarity = Analysis.get_sequence_similarity(
        #     self.naming_args.emboss_naming.replace(name_placeholder, chimera.file_stem))
        # Analysis.relative_stabilty()
        data_columns = Analysis.determine_columns_from_container(self)
        Analysis.convert_data_dict_to_csv(data_columns, self)


class ShiftedFastaArguments:
    fasta_toggles: dict
    '''A dictionary that holds optional operations within the fasta operation, they are generally turned on by 
    placing True within the quotes of their value'''
    Make_a_lst_of_created_fasta_files: bool
    Manually_control_number_of_scanner_movements: bool
    Create_an_alignment: bool
    Make_pair_or_combo_heteromers: bool
    '''This toggle either creates every combination between all protein list generated by each json file if "combo" 
    is entered in the json file, if "pair" is entered all chimera lists in the per json file will be paired one to 
    one to one etc. or single file list will be created for homomeric proteins. Switch toggle to pair if you one of 
    your unique proteins will not be a chimera'''
    fasta_file_list_name: str
    scanner_start: int
    scanner_length: int
    scanner_rate: int
    parent_aln_file: str
    number_of_subunits: int

    def __init__(self, fasta_dict):
        for key, value in fasta_dict.items():
            setattr(self, key, value)
        self.scanner_end = self.scanner_start + self.scanner_length
        if self.scanner_length < 10:
            raise RuntimeError("scanner_length must be at 10 units")


class ShiftedAnalysisArguments:
    analysis_toggles: dict
    make_plddts: str
    make_pdbs: str

    analysis_output_file: str
    '''Path of the file you want created that will contain all analysis created.'''
    column_names: dict
    '''A dictionary containing a multiple keys and their list value that each correspond with a type of data that is 
    output by the script, a user-selected title for the data column and whether that data is included in the final 
    output'''
    filename_stems: list
    '''A list ["True","Protein"] that contains quotes that when filled with True, include a data column in the final 
    output that contains the nicknames created from the naming convention previously established. The second index 
    controls the column name in the final output'''
    relative_stability: list
    overall_chimera_stability: list

    def __init__(self, fasta_dict):
        for key, value in fasta_dict.items():
            setattr(self, key, value)


class ShiftedChimeraArgs:
    argument_dict: dict
    operation_toggles: dict
    run_fasta_operation: str
    alphafold_submission: str
    run_analysis_operation: str
    fasta_args: ShiftedFastaArguments
    naming_args: HomomerNamingArguments
    submission_args: HomomerSubmissionArguments
    analysis_args: HomomerAnalysisArguments

    def __init__(self, shifted_json_file):
        with open(shifted_json_file, 'rb') as jfile:
            self.argument_dict = load(jfile)
        self.get_dict_args(HomomerFastaArguments, 'fasta_args', 'fasta_arguments')
        self.get_dict_args(HomomerNamingArguments, 'naming_args', 'naming_arguments')
        self.get_dict_args(HomomerSubmissionArguments, 'submission_args', 'submission_arguments')
        self.operation_toggles = self.argument_dict['operation_toggles']
        self.chimeras = ()

    def get_dict_args(self, dict_class, attr_name, arg_key):
        dict_args = self.argument_dict[arg_key]
        setattr(self, attr_name, dict_class(dict_args))

    def fasta_operations(self):
        chi_aln_df = ChimeraGenerator.create_chimera_combinations(self.fasta_args.parent_aln_file,
                                                                  self.fasta_args.scanner_length,
                                                                  self.fasta_args.scanner_start,
                                                                  self.fasta_args.scanner_rate)
        fasta_dir = Path(self.naming_args.output_directory).joinpath('Fasta')
        os.makedirs(fasta_dir, exist_ok=True)
        chi_aln_df.make_individual_fasta(fasta_dir, self.fasta_args.number_of_subunits)
