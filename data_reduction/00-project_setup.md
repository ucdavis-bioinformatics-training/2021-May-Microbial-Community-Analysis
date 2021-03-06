# Microbial Community Analysis Project Setup

Let's set up a project directory for the week, and talk a bit about project philosophy..

First, create a directory for you and the example project in the workshop share directory:

```bash
cd
mkdir -p /share/workshop/mca_workshop/$USER
```
## DNA Sequence data from a provider

* Fastq files are actually not raw data from the provider, “raw” data is actually bcl files. Sequencing providers will run an application, bcl2fastq, with a sample sheet to produce demultiplexed (by barcode(s)) fastq files.

## Link Raw Fastq files

1. Next, go into that directory, create a raw data directory (we are going to call this 00-RawData) and cd into that directory. Lets then create symbolic links to the sample directories that contains the raw data.

    ```bash
    cd /share/workshop/mca_workshop/$USER
    mkdir 00-RawData
    cd 00-RawData/
    ln -s /share/workshop/mca_htstream/00-RawDataSm/* .
    ```

    This directory now contains a folder for each sample and the fastq files for each sample are in the sample folders.

1. Let's create a sample sheet for the project and store sample names in a file called samples.txt

    ```bash
    ls *R1* > ../samples.txt
    sed -i 's/_R1.fastq.gz//g' ../samples.txt
    cat ../samples.txt
    ```

## Getting To Know Your Raw Data

Lets first spend a minute talking about sequence read files and the [fastq](filetypes) file format.

1. Now, take a look at the raw data directory.

    ```bash
    ls /share/workshop/mca_workshop/$USER/00-RawData
    ```

    Lets get a better look at all the files in all of the directories.

    ```bash
    ls -lah *
    ```

1. Pick a directory and go into it. View the contents of the files using the 'less' command, when gzipped used 'zless' (which is just the 'less' command for gzipped files):

    ```bash
    zless Bs1_2C_A0_R1.fastq.gz
    ```

    Make sure you can identify which lines correspond to a read and which lines are the header, sequence, and quality values. Press 'q' to exit this screen.

1. Then, let's figure out the number of reads in this file. A simple way to do that is to count the number of lines and divide by 4 (because the record of each read uses 4 lines). In order to do this use cat to output the uncompressed file and pipe that to "wc" to count the number of lines:

    ```bash
    zcat Bs1_2C_A0_R1.fastq.gz | wc -l
    ```

    Divide this number by 4 and you have the number of reads in this file.

1. One more thing to try is to figure out the length of the reads without counting each nucleotide. First get the first 4 lines of the file (i.e. the first record):

    ```bash
    zcat Bs1_2C_A0_R1.fastq.gz  | head -2 | tail -1
    ```

    Note the header lines (1st and 3rd line) and sequence and quality lines (2nd and 4th) in each 4-line fastq block.

1. Then, *copy and paste* the DNA sequence line into the following command (replace [sequence] with the line):

    ```bash
    echo -n [sequence] | wc -c
    ```

    This will give you the length of the read. Also can do the bash one liner:

    ```bash
    echo -n $(zcat Bs1_2C_A0_R1.fastq.gz  | head -2 | tail -1) | wc -c
    ```

    See if you can figure out how this command works.

## Prepare our experiment folder for analysis

Now go back to your directory and create two directories called 'slurm_out' and '01-HTS_Preproc':

```bash
cd /share/workshop/mca_workshop/$USER
mkdir slurm_out
mkdir 01-HTS_Preproc
```

The results of all our slurm script will output .out and .err files into the slurm_out folder. The results of our preprocessing steps will be put into the 01-HTS_Preproc directory. The next step after that will go into a "02-..." directory, etc. You can collect scripts that perform each step, and notes and metadata relevant for each step, in the directory for that step. This way anyone looking to replicate your analysis has limited places to search for the commands you used. In addition, you may want to change the permissions on your original 00-RawData directory to "read only", so that you can never accidentally corrupt (or delete) your raw data. We won't worry about this here, because we've linked in sample folders.

Your directory should then look like the below:
```
$ ls
00-RawData  01-HTS_Preproc  samples.txt  slurm_out
```

## Questions you should now be able to answer.

1. How many reads are in the sample you checked?
2. How many basepairs is R1, how many is R2?
3. What is the name of the sequencer this dataset was run on?
4. Which run number is this for that sequencing?
5. What lane was this ran on?
6. Randomly check a few samples, were the all run the same sequencing, run, and lane?
