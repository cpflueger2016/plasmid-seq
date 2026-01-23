#!/bin/bash --login
#SBATCH --job-name=plasmidSeq
#SBATCH --partition=ll
#SBATCH --mem-per-cpu=10G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=04:00:00
#SBATCH --export=NONE
#SBATCH --mail-user=ryan.lister.lab@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

# Start of job
echo $SLURM_JOB_NAME job started at  `date`

# To compile with the GNU toolchain
module load gcc

conda activate /group/llshared/shared_conda_envs/plasmidseq

# Prepend the plasmid-seq environment to the PATH variable to prevent clash with python version from the Module
PATH=/group/llshared/shared_conda_envs/plasmidseq:$PATH
export PATH

# leave in, it lists the environment loaded by the modules
module list


#  Note: SLURM_JOBID is a unique number for every job.
#  These are generic variables
current_date=$(date +%Y-%m-%d)
JOBNAME="${SLURM_JOB_NAME}_${current_date}"
SCRATCH=$MYSCRATCH/$JOBNAME/$SLURM_JOBID
RESULTS=/group/llshared/PlasmidSeq/Results/$JOBNAME/$SLURM_JOBID

# Get the plasmid seq data folder from the first argument
plasmidSeqData=$1

# Divide SLURM_CPUS_PER_TASK by 4
n_plasmids=$((SLURM_CPUS_PER_TASK / 4))


# Check for PL_to_fasta.tsv file and for Reference fasta folder
if [ ! -f "PL_to_fasta.tsv" ] || [ ! -d "Fasta_Reference_Files" ]; then
	echo "<PL_to_fasta.tsv> file AND <Fasta_Reference_Files> folder need to exist in the parent folder. Please provide and try again."
	echo
	exit 1
fi

# Check if sequencing data has been provided

if [ ! -d $plasmidSeqData ]; then
	echo "Please provide the fastq plasmidSeqData folder with all the Projects and the corresponding fastq files in those folders."
    echo
	exit 1
fi

###############################################
# Creates a unique directory in the SCRATCH directory for this job to run in.
if [ ! -d $SCRATCH ]; then 
    mkdir -p $SCRATCH
fi 
echo SCRATCH is $SCRATCH

###############################################
# Creates a unique directory in your GROUP directory for the results of this job
if [ ! -d $RESULTS ]; then 
    mkdir -p $RESULTS
fi
echo the results directory is $RESULTS


################################################
# declare the name of the output file or log file
OUTPUT=${JOBNAME}_${SLURM_JOBID}.log

#############################################


## Copy data to scratch
cp -r PL_to_fasta.tsv Fasta_Reference_Files "$SCRATCH"
rsync -av --exclude='Reports' --exclude='Stats' --exclude="Logs" --exclude="Undetermined*" "$plasmidSeqData"/* "$SCRATCH"


## Navigate to scratch
cd "$SCRATCH"

# Move Fastq files to folders
for i in */PL*R1*; do f=${i%_S*}; mkdir $f; n=${i/_R1_/_R2_}; mv ${i} ${n} ${f}; done

# Match plasmid fasta to samples
/group/llshared/PlasmidSeq/match_plasmid_fasta_refs.bash PL_to_fasta.tsv Fasta_Reference_Files/

# Create jobs file to capture all samples
for i in */PL*/; do ref=`find "$(pwd)/${i}" -name  "*.fa" -o -name "*.fasta" -o -name "na"`; ref=${ref##*/}; one=`ls ${i}*_R1_*gz`; one=${one##*/}; two=`ls ${i}*_R2_*gz`; two=${two##*/}; uID=`echo ${i} | perl -ne '/(PL\d{4,})/ & print $1."\n"'`; echo -ne "${i}\t${ref}\t${one}\t${two}\t${uID}\n" >> jobs.txt; done

# Map all the plasmid files
parallel --verbose -j ${n_plasmids} --colsep "\t" "cd {1}; bash /group/llshared/PlasmidSeq/plasmidseq_mapper_PE_plann_mqc.sh -1 {3} -2 {4} -r {2} -c -m 300 -u {5} -y -s -q 30 ; cd .." :::: jobs.txt

#############################################

# Get an idea of the directory structure and files
ls -lR * >> ${OUTPUT}


#############################################

# Find all folders except "Fasta_Reference_Files" and move them
find . -mindepth 1 -maxdepth 1 -type d ! \( -name "Fasta_Reference_Files" -o -name "Stats" -o -name "Reports" \) -exec mv -t $RESULTS {} +

# Delete SCRATCH
rm -r $MYSCRATCH/$JOBNAME


# Copy results also to Run/Aligned
AlignedData=${plasmidSeqData/fastq*/}"/Aligned"
mkdir -p $AlignedData
cp -r  $RESULTS/* $AlignedData

# Update permissions
USER=${USER:-$(whoami)}
chown -R $USER:llusers $AlignedData $RESULTS 2>/dev/null || true
chmod -R 777 $AlignedData $RESULTS

echo $JOBNAME job finished at  `date`

