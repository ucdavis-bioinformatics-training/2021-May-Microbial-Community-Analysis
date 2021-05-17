#!/bin/bash

#SBATCH --job-name=r16S_amplicon # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=120
#SBATCH --mem=3000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production
#SBATCH --account=workshop
#SBATCH --reservation=workshop
#SBATCH --array=1-8
#SBATCH --output=slurm_out/r16S_amplicon_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=slurm_out/r16S_amplicon_%A_%a.err # File to which STDERR will be written
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail@address.com

start=`date +%s`
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID


export PATH="/share/workshop/mca_htstream/bin/:$PATH"

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" samples.txt`

inpath='00-RawData'
outpath='01-HTS_Preproc'
[[ -d ${outpath} ]] || mkdir -p ${outpath}

echo "SAMPLE: ${sample}"

call=" \
    hts_Stats \
        --stats-file ${outpath}/${sample}.json \
        -1 ${inpath}/${sample}*_R1.fastq.gz \
        -2 ${inpath}/${sample}*_R2.fastq.gz \
        --notes 'Initial Stats' | \
    hts_Overlapper \
        --append-stats-file ${outpath}/${sample}.json \
        --number-of-threads 4 \
        --notes 'Overlap reads' | \
    hts_Primers \
        --append-stats-file ${outpath}/${sample}.json \
        --primers_5p GTGYCAGCMGCCGCGGTAA \
        --primers_3p GGACTACNVGGGTWTCTAAT \
        --min_primer_matches 2 \
        --flip \
        --float 5 \
        --notes 'Single set V3V4 primers' | \
    hts_NTrimmer \
        --append-stats-file ${outpath}/${sample}.json \
        --exclude \
        --notes 'Remove any reads with Ns' | \
    hts_LengthFilter \
        --append-stats-file ${outpath}/${sample}.json \
        --min-length 100 \
        --max-length 400 \
        --notes 'Filter sequences 100 - 400' | \
    hts_Stats \
        --append-stats-file ${outpath}/${sample}.json \
        --force \
        --fastq-output ${outpath}/${sample} \
        --notes 'Final Stats'"
echo $call
eval $call

end=`date +%s`

runtime=$((end-start))

echo $runtime
