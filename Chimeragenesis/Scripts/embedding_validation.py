import subprocess
import numpy as np
import pandas as pd
import ESM,AccessiontoAlignment
from pathlib import Path
no_tag_aln="/scratch/jws6pq/Notebook/ESM/Chimeragenesis/schema_no_tag.aln"
notag_df=ESM.SequenceDataframe(no_tag_aln)

def predict_spired_tm():
    for label in notag_df.index:
        if len(inheritance := notag_df.get_description(label)) == 2:
            for parent in inheritance:
                file_stem=f'{label}w{parent}'
                output=Path(f"/scratch/jws6pq/Notebook/spired_parent_child/{file_stem}")
                output.mkdir(exist_ok=True,parents=True)
                fasta_file=output.joinpath(f'{file_stem}.fa')
                # chimera first
                AccessiontoAlignment.fasta_creation(fasta_file,AccessiontoAlignment.create_seq_records(parent,notag_df.get_aln(parent).replace('-','G'))+AccessiontoAlignment.create_seq_records(label,notag_df.get_aln(label).replace('-','G')))
                result=subprocess.run(['sbatch',"/scratch/jws6pq/Notebook/ESM/Chimeragenesis/spired.slurm",fasta_file,output], text=True, capture_output=True)
                print(result.stdout, result.stderr)

def analyze_tm():
    for label in notag_df.index:
        if len(inheritance := notag_df.get_description(label)) == 2:
            TMs = []
            for parent in inheritance:
                TMs.append(pd.read_csv(f"/scratch/jws6pq/Notebook/spired_parent_child/{f'{label}w{parent}'}/pred.csv").loc[0, 'dTm'])
            notag_df.add_value(label,'Tm',np.mean(TMs))
    notag_df.save_df('/scratch/jws6pq/Notebook/spired_parent_child/tm_data_parent_first.csv')
analyze_tm()
