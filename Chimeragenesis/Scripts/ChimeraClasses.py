import os
import subprocess
import time
from json import load
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import AccessiontoAlignment, ESM, ChimeraGenerator, Analysis


class HomomerNamingArguments:
    file_stem_placeholder: str
    WT_convention: str
    chimera_convention: str

    # All non-nested keys from the json are set as attributes to be called later
    def __init__(self, naming_dict):
        for key, value in naming_dict.items():
            setattr(self, key, value)


class HomomerFastaArguments:
    output_directory: str
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
    parent_msa: str
    '''Path to the multiple sequence alignment to be created or a user-generated msa of all parent sequences'''
    collective_fasta: str
    '''Path to the fasta file to be created containing all sequences to be run on Alphafold (parents+chimeras) and their inheritance dicts. Path must have .inh suffix'''
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
    sbatch_group: str
    "group needed to access computing resources through slurm/sbatch"
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
    fasta_dir: Path

    def __init__(self, arg_json,full_use=True):
        self.base_splice = None
        with open(arg_json, 'rb') as jfile:
            self.argument_dict = load(jfile)
        self.has_inheritance = False
        self.operation_toggles = self.argument_dict['operation_toggles']
        self.get_dict_args(HomomerFastaArguments, 'fasta_args', 'fasta_arguments')
        self.get_dict_args(HomomerNamingArguments, 'naming_args', 'naming_arguments')
        self.get_dict_args(HomomerSubmissionArguments, 'submission_args', 'submission_arguments')
        self.get_dict_args(HomomerAnalysisArguments, 'analysis_args', 'analysis_arguments')
        if Path(self.fasta_args.parent_msa).exists():
            self.parent_df = ESM.SequenceDataframe(self.fasta_args.parent_msa)
        elif full_use:
            raise FileNotFoundError(f'Could not find parent msa {self.fasta_args.parent_msa}')
        self.fasta_args.output_directory = Path(self.fasta_args.output_directory)
        # TODO might need to change inh ending
        self.fasta_args.collective_fasta = Path(self.fasta_args.collective_fasta).with_suffix('.inh')
        if self.fasta_args.collective_fasta.exists():
            self.collective_df = ESM.SequenceDataframe(self.fasta_args.collective_fasta)
            self.has_inheritance = True
        self.fasta_dir = self.fasta_args.output_directory.joinpath('Fasta')
        self.pdb_dir = self.fasta_args.output_directory.joinpath('PDB')
        self.alphafold_dir = self.fasta_args.output_directory.joinpath('AlphaFold')

    def get_dict_args(self, dict_class, attr_name, arg_key):
        dict_args = self.argument_dict[arg_key]
        setattr(self, attr_name, dict_class(dict_args))

    def make_fasta_paths(self):
        self.base_splice = AccessiontoAlignment.extract_seq_from_fasta(self.fasta_args.sequence_of_interest)
        self.parent_df['chi_seq'] = self.parent_df.index.map(
            lambda label: AccessiontoAlignment.alignment_finder(self.base_splice, label,
                                                                self.fasta_args.base_identifier,
                                                                self.fasta_args.parent_msa)[0].replace('-', ''))
        self.parent_df['file_stem'] = self.parent_df.index.map(
            lambda label: self.naming_args.WT_convention.replace(self.naming_args.file_stem_placeholder, label))
        self.parent_df['chi_file_stem'] = self.parent_df.index.map(
            lambda label: self.naming_args.chimera_convention.replace(self.naming_args.file_stem_placeholder, label))

        fastas = self.parent_df['file_stem'].map(
            lambda stem: self.fasta_dir.joinpath(stem).with_suffix('.fa').__str__()).to_list()
        chi_fastas = self.parent_df['chi_file_stem'].map(
            lambda stem: self.fasta_dir.joinpath(stem).with_suffix('.fa').__str__()).to_list()
        self.all_fastas = fastas + chi_fastas
        self.make_inheritance()

    def make_inheritance(self):
        for label, data in self.parent_df.iterrows():
            self.collective_df.add_protein(label, data['sequence'], {})
            self.collective_df.add_protein(data['chi_file_stem'],
                                           *AccessiontoAlignment.alignment_finder(self.base_splice, label,
                                                                                  self.fasta_args.base_identifier,
                                                                                  self.fasta_args.parent_msa))
        self.collective_df.dataframe_to_aln(self.fasta_args.collective_fasta)

    def fasta_operations(self):
        self.make_fasta_paths()
        num_of_subunits = self.fasta_args.number_of_subunits
        fasta_dir = self.fasta_args.output_directory.joinpath('Fasta')
        os.makedirs(fasta_dir, exist_ok=True)
        for label, data in self.parent_df.iterrows():
            AccessiontoAlignment.fasta_creation(fasta_dir.joinpath(data['file_stem']).with_suffix('.fa'),
                                                AccessiontoAlignment.create_seq_records(label, data['sequence'],
                                                                                        subunit_count=num_of_subunits))
            AccessiontoAlignment.fasta_creation(fasta_dir.joinpath(data['chi_file_stem']).with_suffix('.fa'),
                                                AccessiontoAlignment.create_seq_records(label, data['chi_seq'],
                                                                                        subunit_count=num_of_subunits))

    def alphafold_submission(self, fastas: list):
        fasta_to_run = ()
        submission_toggles = self.submission_args.submission_toggles
        proteins_per_slurm = self.submission_args.proteins_per_slurm
        # Loops through all fastas created and checks if they are complete by looking for their ranking_debug file
        if submission_toggles['stragglers_or_custom_or_all'] == 'stragglers':
            for fasta in fastas:
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
            fasta_to_run = fastas
        if self.submission_args.submission_toggles['sbatch_slurms']:
            for slurm_index, file_index in enumerate(range(0, len(fasta_to_run), proteins_per_slurm)):
                result = subprocess.run(['sbatch', '-J', str(slurm_index) + self.submission_args.slurm_naming, '-A',
                                         self.submission_args.sbatch_group,
                                         '-e', Path(self.fasta_args.output_directory).joinpath(
                        self.submission_args.slurm_naming + str(slurm_index)).with_suffix(".err"),
                                         self.submission_args.alphafold_shell_script,
                                         ",".join(fasta_to_run[file_index:file_index + proteins_per_slurm]),
                                         self.alphafold_dir], capture_output=True)
                print(result.stdout)

    def analysis_operations(self):
        analysis_toggles = self.analysis_args.analysis_toggles
        if not self.has_inheritance:
            self.make_inheritance()
        if analysis_toggles['analyze_alphafold']:
            Analysis.generate_alphafold_files(self.alphafold_dir, self.fasta_args.output_directory,analysis_toggles['make_plddts'])
            for chi, data in self.collective_df.iterrows():
                chi_plddt = Analysis.get_plddt_dict_from_pdb(
                    self.pdb_dir.joinpath(str(chi)).with_suffix('.pdb'))[self.collective_df.get_sequence(chi)]
                self.collective_df.loc[chi, 'Overall Stability'] = Analysis.overall_confidence(chi_plddt)
                if len(inheritance_dict := self.collective_df.get_description(chi)) == 2:
                    plddt_dict = {}
                    for parent in inheritance_dict.keys():
                        plddt_dict[parent] = Analysis.get_plddt_dict_from_pdb(
                            self.pdb_dir.joinpath(parent).with_suffix('.pdb'))[self.collective_df.get_sequence(parent)]
                    plddt_dict[str(chi)] = chi_plddt
                    relative_stability = Analysis.relative_stabilty(plddt_dict, inheritance_dict)
                    self.collective_df.loc[chi, 'Relative Stability (%)'] = relative_stability
            self.collective_df.get_sequence_identity()
            self.collective_df.to_csv(self.analysis_args.analysis_output_csv)


class ShiftedFastaArguments:
    fasta_toggles: dict
    '''A dictionary that holds optional operations within the fasta operation, they are generally turned on by 
    placing True within the quotes of their value'''
    Make_a_lst_of_created_fasta_files: bool
    Create_an_alignment: bool
    Make_pair_or_combo_heteromers: bool
    '''This toggle either creates every combination between all protein list generated by each json file if "combo" 
    is entered in the json file, if "pair" is entered all chimera lists in the per json file will be paired one to 
    one to one etc. or single file list will be created for homomeric proteins. Switch toggle to pair if you one of 
    your unique proteins will not be a chimera'''
    parent_aln_file: str
    number_of_subunits: int
    output_directory: str
    collective_fasta: str

    def __init__(self, fasta_dict):
        for key, value in fasta_dict.items():
            setattr(self, key, value)


class ShiftedSubmissionArguments:
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
    sbatch_group: str
    "group needed to access computing resources through slurm/sbatch"
    alphafold_shell_script: str
    '''Path to the shell script that runs alphafold, it must take a list of comma-separated fastas 
    and an output as its arguments. It's recommended you use the on provided in github'''
    embedding_script: str

    def __init__(self, fasta_dict):
        for key, value in fasta_dict.items():
            setattr(self, key, value)


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
    submission_args: ShiftedSubmissionArguments
    analysis_args: ShiftedAnalysisArguments

    def __init__(self, shifted_json_file):
        with open(shifted_json_file, 'rb') as jfile:
            self.argument_dict = load(jfile)
        self.get_dict_args(ShiftedFastaArguments, 'fasta_args', 'fasta_arguments')
        self.get_dict_args(ShiftedSubmissionArguments, 'submission_args', 'submission_arguments')
        self.get_dict_args(ShiftedAnalysisArguments, 'analysis_args', 'analysis_arguments')
        self.operation_toggles = self.argument_dict['operation_toggles']

    def get_dict_args(self, dict_class, attr_name, arg_key):
        dict_args = self.argument_dict[arg_key]
        setattr(self, attr_name, dict_class(dict_args))

    def fasta_operations(self):
        # TODO 2 problems 1. theres a issue of the alignment finder not being able to find anything and returning -1 and not throwing error 2. the chimeracombo function will label as 1200-1600 etc even though the protein is only 1300 long
        MIN_SPLICE_PERCENT = .35
        MAX_SPLICE_PERCENT=.65
        seq_dict = AccessiontoAlignment.create_dictionary_from_alignment(self.fasta_args.parent_aln_file)
        # TODO Scanner rate should change somehow
        aln_len=len(list(seq_dict.values())[0])
        aln_dfs=[]
        for scanner in range(round(MIN_SPLICE_PERCENT*aln_len),round(MAX_SPLICE_PERCENT*aln_len)):
            aln_dfs.append(ChimeraGenerator.create_chimera_combinations(self.fasta_args.parent_aln_file,scanner,0,20))
        total=pd.concat(aln_dfs)
        fasta_dir = Path(self.fasta_args.output_directory).joinpath('Fasta')
        fasta_dir.mkdir(exist_ok=True)
        total.__class__=ESM.SequenceDataframe
        total.dataframe_to_aln(self.fasta_args.collective_fasta)
        total.dataframe_to_multi_fa(Path(self.fasta_args.collective_fasta).with_suffix('.mfa'))
        total.make_individual_fasta(fasta_dir, self.fasta_args.number_of_subunits)

    def submission_operations(self):
        submission_toggles = self.submission_args.submission_toggles
        if submission_toggles['get_embeddings']:
            result = subprocess.run(['sbatch',
                                     '-e', Path(self.fasta_args.output_directory).joinpath(
                    self.submission_args.slurm_naming).with_suffix(".err"),
                                     self.submission_args.embedding_script,
                                     self.fasta_args.collective_fasta,
                                     Path(self.fasta_args.collective_fasta).with_suffix('.pkl')], capture_output=True)
            print(result.stdout)

    def analysis_operations(self):
        analysis_toggles = self.analysis_args.analysis_toggles
        if analysis_toggles['analyze_embeddings']:
            embeddings = ESM.SequenceDataframe(Path(self.fasta_args.collective_fasta).__str__(),
                                               Path(self.fasta_args.collective_fasta).with_suffix('.pkl').__str__())
            for func in [ESM.NormType.cosine, ESM.NormType.manhattan, ESM.NormType.euclidean, ESM.NormType.dot_product]:
                embeddings.score_all_embeddings(func, np.mean)
                dist_type = func.name
                embeddings[f'{dist_type}_rank'] = embeddings[dist_type].rank(
                    ascending=func != ESM.NormType.cosine and func != ESM.NormType.dot_product)
            embeddings.get_sequence_identity()
            if self.analysis_args.analysis_output_file:
                embeddings.to_csv(self.analysis_args.analysis_output_file)


#
class NonHomologyFastaArguments:
    fasta_toggles: dict
    '''A dictionary that holds optional operations within the fasta operation, they are generally turned on by 
    placing True within the quotes of their value'''
    Make_a_lst_of_created_fasta_files: bool
    Create_an_alignment: bool
    Make_pair_or_combo_heteromers: bool
    '''This toggle either creates every combination between all protein list generated by each json file if "combo" 
    is entered in the json file, if "pair" is entered all chimera lists in the per json file will be paired one to 
    one to one etc. or single file list will be created for homomeric proteins. Switch toggle to pair if you one of 
    your unique proteins will not be a chimera'''
    parent_mfa_file: str
    number_of_subunits: int
    output_directory: str
    collective_fasta: str

    def __init__(self, fasta_dict):
        for key, value in fasta_dict.items():
            setattr(self, key, value)


class NonHomologySubmissionArguments:
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
    custom_fasta_list: str
    '''Path to file that for list that either be created by toggling on create_file_of_stragglers, or user-generated that should contain 
    line separated fastas of proteins you want predicted'''
    proteins_per_slurm: int
    sbatch_group: str
    "group needed to access computing resources through slurm/sbatch"
    alphafold_shell_script: str
    '''Path to the shell script that runs alphafold, it must take a list of comma-separated fastas 
    and an output as its arguments. It's recommended you use the on provided in github'''
    embedding_script: str

    def __init__(self, fasta_dict):
        for key, value in fasta_dict.items():
            setattr(self, key, value)


class NonHomologyAnalysisArguments:
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


class NonHomologyChimeraArgs:
    argument_dict: dict
    operation_toggles: dict
    run_fasta_operation: str
    alphafold_submission: str
    run_analysis_operation: str
    fasta_args: NonHomologyFastaArguments
    submission_args: NonHomologySubmissionArguments
    analysis_args: NonHomologyAnalysisArguments

    def __init__(self, NonHomology_json_file):
        with open(NonHomology_json_file, 'rb') as jfile:
            self.argument_dict = load(jfile)
        self.get_dict_args(NonHomologyFastaArguments, 'fasta_args', 'fasta_arguments')
        self.get_dict_args(NonHomologySubmissionArguments, 'submission_args', 'submission_arguments')
        self.get_dict_args(NonHomologyAnalysisArguments, 'analysis_args', 'analysis_arguments')
        self.operation_toggles = self.argument_dict['operation_toggles']
        self.fasta_chunk_dir = Path(self.fasta_args.output_directory).joinpath('FastaChunks')
        self.embed_chunk_dir = Path(self.fasta_args.output_directory).joinpath('EmbedChunks')
        self.alphafold_dir=Path(self.fasta_args.output_directory).joinpath('AlphaFold')
        self.fasta_dir =Path(self.fasta_args.output_directory).joinpath('Fasta')
        if Path(self.fasta_args.collective_fasta).exists():
            start=time.time()
            fasta_dfs=[]
            for chunk_fasta in self.fasta_chunk_dir.iterdir():
                fasta_dfs.append(ESM.SequenceDataframe(chunk_fasta))
            self.collective_df = ESM.SequenceDataframe(unconverted_df=pd.concat(fasta_dfs))
            end = time.time()
            self.has_inheritance = True
            print((end-start)/60)
    def get_dict_args(self, dict_class, attr_name, arg_key):
        dict_args = self.argument_dict[arg_key]
        setattr(self, attr_name, dict_class(dict_args))

    def fasta_operations(self):
        MIN_SPLICE_PERCENT = .35
        MAX_SPLICE_PERCENT=.65

        chimera_df=ChimeraGenerator.create_combinations_no_aln(self.fasta_args.parent_mfa_file,(MIN_SPLICE_PERCENT,MAX_SPLICE_PERCENT),self.fasta_args.collective_fasta)
        fasta_dir = self.fasta_chunk_dir
        fasta_dir.mkdir(exist_ok=True)
        num_chunks=chimera_df.shape[0]/2000
        chunks:list[ESM.SequenceDataframe]=[ESM.SequenceDataframe(unconverted_df=chunk) for chunk in np.array_split(chimera_df,num_chunks)]
        aln_df = ESM.SequenceDataframe(self.fasta_args.parent_mfa_file)
        parent1, parent2 = aln_df.index
        aln_df.add_value(parent1, 'description', {parent1: None})
        aln_df.add_value(parent2, 'description', {parent2: None})
        for ind,chunk in enumerate(chunks):
            chunk=ESM.SequenceDataframe(unconverted_df=chunk._append(aln_df))
            chunk.dataframe_to_aln(fasta_dir.joinpath(Path(self.fasta_args.collective_fasta).stem+f'_{str(ind)}.inh'))

    def submission_operations(self):
        submission_toggles = self.submission_args.submission_toggles
        if submission_toggles['get_embeddings']:
            self.embed_chunk_dir.mkdir(exist_ok=True)
            for chunk_fasta in self.fasta_chunk_dir.iterdir():
                result = subprocess.run(['sbatch',
                                         '-e', Path(self.fasta_args.output_directory).joinpath(
                        self.submission_args.slurm_naming).with_suffix(".err"),
                                         self.submission_args.embedding_script,
                                         chunk_fasta,self.embed_chunk_dir.joinpath(chunk_fasta.stem).with_suffix('.pkl')], capture_output=True)
                print(result.stdout,result.stderr)

        if submission_toggles['run_AF']:
            self.alphafold_dir.mkdir(exist_ok=True)
            fasta_dir=self.fasta_dir
            fasta_dir.mkdir(exist_ok=True)
            seq_dict=AccessiontoAlignment.create_dictionary_from_alignment(self.fasta_args.collective_fasta)
            fastas_to_run=[]
            if submission_toggles['run_custom_list']:
                # MAKE MINI INHERITANCE FASTA??
                with open(self.submission_args.custom_fasta_list, 'r') as run_list:
                    chimeras_to_run = [x.replace('\n','') for x in run_list.readlines() if x.replace('\n','')]
                for label in chimeras_to_run:
                    fasta_file=fasta_dir.joinpath(label).with_suffix('.fa')
                    fastas_to_run.append(fasta_file)
                    AccessiontoAlignment.fasta_creation(fasta_file,AccessiontoAlignment.create_seq_records(label,seq_dict[label],subunit_count=self.fasta_args.number_of_subunits))
            else:
                df=ESM.SequenceDataframe(self.fasta_args.collective_fasta)
                df.make_individual_fasta(fasta_dir,self.fasta_args.number_of_subunits)
                for label in df.index:
                    fasta_file = fasta_dir.joinpath(label).with_suffix('.fa')
                    fastas_to_run.append(fasta_file)
            fastas_to_run=[fasta.__str__() for fasta in fastas_to_run]
            for slurm_index, file_index in enumerate(range(0, len(fastas_to_run), self.submission_args.proteins_per_slurm)):
                result = subprocess.run(['sbatch', '-J', str(slurm_index) + self.submission_args.slurm_naming, '-A',
                                         self.submission_args.sbatch_group,
                                         '-e', Path(self.fasta_args.output_directory).joinpath(
                        self.submission_args.slurm_naming + str(slurm_index)).with_suffix(".err"),
                                         self.submission_args.alphafold_shell_script,
                                         ",".join(fastas_to_run[file_index:file_index + self.submission_args.proteins_per_slurm]),
                                         self.alphafold_dir], capture_output=True)
                print(result.stdout)

    def analysis_operations(self):
        analysis_toggles = self.analysis_args.analysis_toggles
        if analysis_toggles['analyze_embeddings']:
            embeddings=[]

            for chunk_fasta in self.fasta_chunk_dir.iterdir():
                embed_file=self.embed_chunk_dir.joinpath(chunk_fasta.stem).with_suffix('.pkl')
                # TODO change
                if not chunk_fasta.exists() and embed_file.exists(): raise FileNotFoundError('One of the files is not like the other')
                embedding_df = ESM.SequenceDataframe(chunk_fasta.__str__(),embed_file.__str__())
                func=ESM.NormType.dot_product
                embedding_df.score_all_per_res(func)
                dist_type = func.name
                embedding_df[f'{dist_type}_rank'] = embedding_df[dist_type].rank(
                    ascending=func != ESM.NormType.cosine and func != ESM.NormType.dot_product)
                embedding_df.drop(columns=['aln_sequence', 'sequence', 'description'])
                del embedding_df.EMBEDDINGS_DICT
                embeddings.append(embedding_df)

            embeddings=pd.concat(embeddings)
            embeddings= embeddings.loc[:, ~embeddings.columns.duplicated()]
            embeddings=ESM.SequenceDataframe(unconverted_df=embeddings)
            embeddings.drop(columns=['aln_sequence'])
            if self.analysis_args.analysis_output_file:
                embeddings.to_csv(self.analysis_args.analysis_output_file)
        if analysis_toggles['analyze_alphafold']:
            for folder in self.alphafold_dir.iterdir():
                chi_label=folder.stem
                pdb=folder.joinpath('ranked_0.pdb')
            # Analysis.generate_alphafold_files(self.alphafold_dir, self.fasta_args.output_directory,
            #                                   analysis_toggles['make_plddts'])
                chi_plddt = Analysis.get_plddt_dict_from_pdb(pdb)[self.collective_df.get_sequence(chi_label)]
                self.collective_df.loc[chi_label, 'Overall Stability'] = Analysis.overall_confidence(chi_plddt)
                if len(inheritance_dict := self.collective_df.get_description(chi_label)) == 2:
                    plddt_dict = {}
                    for parent in inheritance_dict.keys():
                        plddt_dict[parent] = Analysis.get_plddt_dict_from_pdb(
                            self.alphafold_dir.joinpath(parent,'ranked_0').with_suffix('.pdb'))[self.collective_df.get_sequence(parent)]
                    plddt_dict[str(chi_label)] = chi_plddt
                    relative_stability = Analysis.relative_stabilty(plddt_dict, inheritance_dict)
                    self.collective_df.loc[chi_label, 'Relative Stability (%)'] = relative_stability
            self.collective_df=self.collective_df.dropna()
            # self.collective_df.get_sequence_identity()
            self.collective_df.to_csv(self.analysis_args.analysis_output_file)
