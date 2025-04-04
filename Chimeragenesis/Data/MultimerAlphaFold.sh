#!/bin/bash

module purge
module load singularity alphafold/2.2.2

export FASTA=$1
export ALPHA_OUT=$2
export SUBUNITS=$3
export DB=$4
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_MEM_FRACTION=20

run --fasta_paths=FASTA \
    --output_dir=ALPHA_OUT \
    --model_preset=SUBUNITS \
    --db_preset=DB \
    --pdb_seqres_database_path=/data/pdb_seqres/pdb_seqres.txt \
    --uniprot_database_path=/data/uniprot/uniprot.fasta \
    --small_bfd_database_path=/data/small_bfd/bfd-first_non_consensus_sequences.fasta \
    --max_template_date=2022-9-30 \
    --num_multimer_predictions_per_model=1 \
    --use_gpu_relax=True
