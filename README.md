# Ensemble Enzyme Prediction Pipeline (E2P2)

The Ensemble Enzyme Prediction Pipeline (E2P2) annotates protein sequences with Enzyme Function classes comprised of full, four-part Enzyme Commission numbers and MetaCyc reaction identifiers. It is the enzyme annotation pipeline used to generate the species-specific metabolic databases at the [Plant Metabolic Network](www.plantcyc.org) since 2013. E2P2 systematically integrates results from two molecular function annotation algorithms using an ensemble classification scheme. For a given genome, all protein sequences are submitted as individual queries against the base-level annotation methods. 

Due to PRIAM's end of development and availability, we've replaced it with [DeepEC](https://bitbucket.org/kaistsystemsbiology/deepec/src/master/) and moved the current E2P2 version to "v5". You can still download the previous version in the "v4" branch.

## Getting Started
The following instuctions are for users to set up the E2P2 pipeline on a Unix machine and start running, developing and/or testing the pipeline.

### Prerequisites
This pipeline is tested on Ubuntu, CentOS, macOS and should theoretically run on all Linux distributions.
* [Python 3](https://www.python.org/downloads/)
* [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* ~~[Java 1.5 or above](https://www.oracle.com/technetwork/java/javase/downloads/index.html)~~
~~**Currently it has been known that PRIAM might not work with Java version 11, we recommend using version 8 instead.~~
* ~~[PRIAM_Search utility V2](http://priam.prabi.fr/REL_JAN18/index_jan18.html)~~
~~**A copy is temperarily included due to the source website being down.~~
* [DeepEC](https://github.com/bxuecarnegie/deepec)
A fork of DeepEC is provided as a submodule, check "environment.yml" for prerequisites.

### Installing

Download E2P2 from [E2P2 at GitHub](https://github.com/carnegie/E2P2)

**Currently the E2P2 pipeline does not support white spaces or other illegal symbols in files paths.**

```
git clone --recurse-submodules https://github.com/carnegie/E2P2.git
```
* If folder "deepec" is empty, run the following command in the 'E2P2' folder
```
git submodule update --init
```


Download Reference Protein Sequence Dataset (RPSD) from https://ftp.dpb.carnegiescience.edu/rpsd/
* Version release_2024_07_31 and up.


Unzip and extract RPSD data (different argument for tar depending on the file extensions)
```
tar -xzf blastdb.tar.gz
tar -xzf deepec.tar.gz
tar -xzf weights.tar.gz
tar -xzf maps.tar.gz

tar -xjf blastdb.tar.bz2
tar -xjf deepec.tar.bz2
tar -xjf weights.tar.bz2
tar -xjf maps.tar.bz2

# If you want to just use the required arguments with a different version of RPSD data, 
# replace the "maps" and "weights" folders under "data" in "E2P2" using the above folders 
# with the same name. And replace the files under "deepec/deepec/data" with the files in deepec.tar.gz.
```

## Usage Example

### config.ini
In the project's data folder is a "config_template.ini".

> **Users need copy the template to the root folder as "config.ini" and edit the environmental variables.**

python3 e2p2.py [-h] --input /PATH/TO/Araport11_genes.201606.pep.repr.fasta -o PATH/TO/output.pf e2p2 --threshold 0.5 

### Required Arguments
    --input INPUT_FILE, -i INPUT_FILE: Path to input protein sequences file

    e2p2: subparser argument, used to separate pipeline arguments and classifier/ensemble arguments

### Optional Arguments Before "e2p2"
    -h, --help            show this help message and exit
    --input INPUT_FILE, -i INPUT_FILE
                        Path to input protein sequences file
    --protein_gene PROTEIN_GENE_PATH, -pg PROTEIN_GENE_PATH
                        Provide a protein to gene map. This can be used to generate a splice variant removed fasta file and output the final version of e2p2.
    --remove_splice_variants, -rm
                        Argument flag to remove splice variants.
    --output OUTPUT_PATH, -o OUTPUT_PATH
                        Path to output file. By Default would be in the same folder of the input.
    --temp_folder TEMP_FOLDER, -tf TEMP_FOLDER
                        Specify the location of the temp folder. By default would be in the same directory of the output.
    --log LOG_PATH, -l LOG_PATH
                        Specify the location of the log file. By default would be "runE2P2.log" in the temp folder.
    --verbose {0,1}, -v {0,1}
                        Verbose level of log output. Default is 0.
                                    0: only step information are logged
                                    1: all information are logged

### Optional Arguments After "e2p2"
    -h, --help: Show help message and exit
    
    --blastp BLASTP, -b BLASTP
                        Command of or path to BLAST+ "blastp".
    --num_threads NUM_THREADS, -n NUM_THREADS
                        Number of threads to run "blastp".
    --blast_db BLAST_DB, -bd BLAST_DB
                        Path to rpsd blast database name. For example, "/PATH/TO/FOLDER/rpsd.fa", where you can find the following files in
                        /PATH/TO/FOLDER:rpsd.fa.phr; rpsd.fa.pin; rpsd.fa.psq
    --blast_e_value BLAST_E_VALUE, -be BLAST_E_VALUE
                        Blastp e-value cutoff
    --blast_weight BLAST_WEIGHT, -bw BLAST_WEIGHT
                        Path to weight file for the blast classifier
    --python_path PYTHON_PATH, -py PYTHON_PATH
                        Command of or path to "java".
    --deepec_path DEEPEC_PATH, -dp DEEPEC_PATH
                        Path to "deepec.py".
    --ec_to_ef_mapping_path EC_TO_EF_MAPPING_PATH, -ee EC_TO_EF_MAPPING_PATH
                        Path to mapping file from ECs to EFs
    --threshold THRESHOLD, -t THRESHOLD
                        Threshold for voting results. Default is 0.5.

### Additional information
- Input protein sequences should be in FASTA format.
- Headers of the FASTA file should begin with the sequence ID followed by a space or '|'.
    For example: >AT1G01010.1 | NAC domain containing protein 1 | Chr1:3760-5630 FORWARD LENGTH=429 | 201606

## Authors

* **Bo Xue** - [bxuecarnegie](https://github.com/bxuecarnegie)

### Previous versions
* **Chuan Wang** - [Chuan Wang](https://github.com/grittyy)
* **Lee Chae**

## Institution
    Plant Metabolic Network
    Department of Plant Biology
    Carnegie Institution for Science
    Stanford, CA 94305


## Acknowledgments

* Special Thanks to
  * Thomas Bernard - *PRIAM*
  * Ludo Cottret - [lipme](https://github.com/lipme) - *Singularity Container For Previous Versions*
