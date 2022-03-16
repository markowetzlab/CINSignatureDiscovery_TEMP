#### Deriving signatures from input matrix ("3_TCGA_PCAWG_CxS_Matrix.txt").

### 1. Step: Generate signatures

## I use SignatureAnalyzer-GPU which implements the ARD-NMF (Tan and Fevotte, 2012) to derive
## signatures. For this part you need a GPU and Pytorch installed. I have a conda environment
## called bayesnmf that includes all installations.

## I run SignatureAnalyzer-GPU a couple of times (e.g. 50x, 100x, 200x), depending on time constraints
## and conclusiveness of results.


### 2. Step: Decide on best signature solution

## Afterwards I run a R script manually ("2_Decide_signature_solution.R") allowing me first to chose
## the optimal K (rank of NMF, number of signatures) which is the mode of K's over all runs.
## Secondly, it allows me to chose the best solution which is either the most optimal solution
## according to the Kullback-Leibler divergence or the most representative solution according to a
## network graph of all solutions with the same K (number of signatures).


# Standard options
# BASE="/Users/drews01/phd/prjcts/cnsigs2/cnsigs2_revisions"
BASE="/home/drews01/cnsigs2_revisions"
# GPU cluster
# BASE="/Users/drews01/data/phd/cnsigs2_revisions"
OUT="${BASE}/Simulation_2_1_Mutational_processes"

SCRIPTPATH="${BASE}/SignatureAnalyzer-GPU/SignatureAnalyzer-GPU.py"
INPUTMATRIX="${OUT}/Out_4_Step_3_CxS_Matrix_20each_N240_SpecialParam1.txt"
MAXITER=100000
K0=43
OUTPUTDIR="${OUT}/Out_5_NMF_Step_4_CxS_K0-${K0}_20each_SpecialParam1_NoWGD"
REGULARISATIONS="--prior_on_W L1 --prior_on_H L1"

# SignatureAnalyzer-GPU run optionsll
ITERS=15
# Every independent on the full pancancer matrix (6,335 samples x 43 components) run takes up roughly
# 780MB graphics memory.
# ~10500 free ram (not OS) => 13 parallel runs possible
PARALLEL=13

# Thanks to the python script I have to go the directory from where the outputdir is then added.
mkdir -p $OUTPUTDIR
cd $OUTPUTDIR

# SignatureAnalyzer-GPU needs pytorch. I created a conda env with all the right packages. If you
# have conda, then run the following command to install the right env:
# conda create --name bayesnmf --file "${BASE}/scripts/conda_env_requirements.txt"
# If you use a different programme to manage environments or installed pytorch directly, then
# uncomment the following line:
conda activate bayesnmf

# The Python packages needed are found in SignatureAnalyzer-GPU/requirements-py3.txt and can be
# installed with the following command:
# pip3 install -r requirements-py3.txt
# or use the python_requirements.txt file in the scripts folder as pip input.

for ITER in $(seq 1 $ITERS); do

    echo $ITER

    for PAR in $(seq 1 $PARALLEL); do

        UUID=$(cat /dev/urandom | tr -dc 'A-Z0-9' | fold -w 6 | head -n 1)
        OUTPATH="${UUID}"
        mkdir -p $OUTPATH
        python $SCRIPTPATH --data $INPUTMATRIX --max_iter=$MAXITER --output_dir $OUTPATH $REGULARISATIONS --labeled --K0 $K0 > ${OUTPATH}/${UUID}_log.txt &

    done

    sleep 4
    PYRUN=`nvidia-smi | grep -c python`
    while [[ $PYRUN != 0 ]]; do
        sleep 2
        PYRUN=`nvidia-smi | grep -c python`
    done

done
