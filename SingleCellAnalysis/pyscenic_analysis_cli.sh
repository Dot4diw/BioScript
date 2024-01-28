#!/usr/bin/bash

set -e
# Create results directory
OUTS_FOLDER="./ResultsOfPySCENIC"
mkdir -p ${OUTS_FOLDER}

# pySCENIC software path.
PYTHON="/software/miniconda3/bin/python"
PYSCENIC="/software/miniconda3/bin/pyscenic"
RBORETO_WITH_MULTIPROCESSING="/software/miniconda3/lib/python3.9/site-packages/pyscenic/cli/arboreto_with_multiprocessing.py"
ADD_VISUALIZATION="./pySCENIC/add_visualization.py"

# Input Data
IN_DATA="./pySCENIC/placenta_AllCell.loom"
# Output Data
PYSCENIC_GRN_OUT=${OUTS_FOLDER}"/pySCENIC_GRN_adjacencies.csv"
PYSCENIC_CTX_OUT=${OUTS_FOLDER}"/pySCENIC_CTX_regulons.csv"
PYSCENIC_AUCELL_OUT=${OUTS_FOLDER}"/pySCENIC_results_filtered.loom"
VISUALIZATION_LOOM_RESULTS=${OUTS_FOLDER}"/pySCENIC_visualization.loom"

# cisTargetDB
CISTARGETDB="./pySCENIC/cisTargetDB/mc9nr"

# TF_FILE
TF_FILE=${CISTARGETDB}"/allTFs_hg38.txt"
# Feather Database
FEATHER_DATA=${CISTARGETDB}"/*rankings.feather"
# TF Annotations Files
TF_ANNOTATIONS_FILE=${CISTARGETDB}"/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

# Num of Workers
NUM_WORKERS=10

### Step 1
date
# Run pySCENIC get regulatory network from command line interface
${PYSCENIC} grn \
        ${IN_DATA} \
        ${TF_FILE} \
        --num_workers \
        ${NUM_WORKERS} \
        -o ${PYSCENIC_GRN_OUT}

date
echo "regulatory networks defined"

# Using arboreto_with_multiprocessing.py
# Run pySCENIC get regulatory network from command line interface
# This circumvents the use of dask, which seems to cause trouble...
# ${ARBORETO_WITH_MULTIPROCESSING} \
#         ${IN_DATA} \
#         ${TF_FILE} \
#         --num_workers ${NUM_WORKERS} \
#         -o ${PYSCENIC_GRN_OUT}
#         --method grnboost2 \
#         --seed 787878
#
# date
# echo "regulatory networks inferred"


### Step 2
date
# Clean up modules using feather ranking databases
${PYSCENIC} ctx \
        ${PYSCENIC_GRN_OUT} \
        ${FEATHER_DATA} \
        --annotations_fname ${TF_ANNOTATIONS_FILE} \
        --expression_mtx_fname ${IN_DATA} \
        --output ${PYSCENIC_CTX_OUT} \
        --mask_dropouts \
        --num_workers ${NUM_WORKERS}

date
echo "regulon modules defined"

### Step 3
date
# Calculate cell module scores (from filtered loom, to extract acurate regulon incidence matrix later on)
${PYSCENIC} aucell \
        ${IN_DATA} \
        ${PYSCENIC_CTX_OUT} \
         --output ${PYSCENIC_AUCELL_OUT} \
        --num_workers ${NUM_WORKERS}
date
echo "cell scoring completed"

#  ### Step 4
#  date
#  # Add dimensionality reduction (derived from pySCENIC AUCell matrix)
#     ##  always skipped this step, as it can also be performed in R (with proper sample integration)
#  ${PYTHON} ${ADD_VISUALIZATION} \
#       --loom_input ${PYSCENIC_AUCELL_OUT} \
#       --loom_output ${VISUALIZATION_LOOM_RESULTS} \
#       --num_workers ${NUM_WORKERS}
#
#  date
#  echo "dimensionality reduction added"