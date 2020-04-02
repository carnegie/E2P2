# Ensemble Enzyme Prediction Pipeline (E2P2)

The Ensemble Enzyme Prediction Pipeline (E2P2) annotates protein sequences with Enzyme Function classes comprised of full, four-part Enzyme Commission numbers and MetaCyc reaction identifiers. It is the enzyme annotation pipeline used to generate the species-specific metabolic databases at the [Plant Metabolic Network](www.plantcyc.org) since 2013. E2P2 systematically integrates results from two molecular function annotation algorithms using an ensemble classification scheme. For a given genome, all protein sequences are submitted as individual queries against the base-level annotation methods.

## Getting Started
The following instuctions are for users to set up the E2P2 pipeline on a Unix machine and start running, developing and/or testing the pipeline.

### Prerequisites
This pipeline is tested on Ubuntu, CentOS, macOS and should theoretically run on all Linux distributions.
* [Python 3](https://www.python.org/downloads/)
* [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [Java 1.5 or above](https://www.oracle.com/technetwork/java/javase/downloads/index.html)
* [PRIAM_Search utility V2](http://priam.prabi.fr/REL_JAN18/index_jan18.html)

### Installing

Download E2P2 from [E2P2 at GitHub](https://github.com/carnegie/E2P2)

```
git clone https://github.com/carnegie/E2P2.git
```

Download Reference Protein Sequence Dataset (RPSD) from https://ftp.dpb.carnegiescience.edu/rpsd/

Unzip and extract RPSD data
```
tar -xzf blastdb.tar.gz
tar -xzf profiles.tar.gz
tar -xzf weights.tar.gz
tar -xzf maps.tar.gz

# If you want to just use the required arguments with a different version of RPSD data, 
# replace the "maps" and "weights" folders under "data" in "E2P2" using the above folders 
# with the same name.
```

## Usage Example

python3 run_pipeline.py [-h] --input /PATH/TO/Araport11_genes.201606.pep.repr.fasta --blastp blastp --java java --priam_search /PATH/TO/PRIAM_search.jar --rpsd /PATH/TO/blastdb/rpsd-4.2.fasta --priam_profile /PATH/TO/profiles

### Required Arguments
    --input INPUT_FILE, -i INPUT_FILE: Path to input protein sequences file
    
    --blastp BLASTP_CMD, -b BLASTP_CMD: Command of or path to NCBI BLAST+ "blastp" executable.
    
    --java JAVA_CMD, -j JAVA_CMD: Command of or path to "java" executable.
    
    --priam_search PRIAM_SEARCH, -ps PRIAM_SEARCH: Path to "PRIAM_search.jar" executable.
    
    --rpsd RPSD_DB, -r RPSD_DB: Path to RPSD BLAST database name.
      For example, "/PATH/TO/FOLDER/rpsd.fa", where you can find the following files in /PATH/TO/FOLDER: rpsd.fa.phr; rpsd.fa.pin; rpsd.fa.psq
      
    --priam_profile PRIAM_PROFILE, -pp PRIAM_PROFILE: Path to PRIAM profiles.
      For example, "/PATH/TO/FOLDER/profiles", where you can find the following in /PATH/TO/FOLDER/profiles:
        files: annotation_rules.xml; genome_rules.xml
        folder: PROFILES: this contains a "LIBRARY" folder and multiple ".chk" files.
### Optional Arguments
    -h, --help: Show help message and exit
    
    --output OUTPUT_PATH, -o OUTPUT_PATH: Path to output file. By Default it would be in the same folder of the input.
    
    --num_threads NUM_THREADS, -n NUM_THREADS: Number of threads to run "blastp". Default is 1
    
    --priam_resume, -pr: Whether or not to resume process if a unfinished PRIAM_search.jar process is found.
    
    --blast_bin BLAST_BIN, -bb BLAST_BIN: Command of or path to NCBI BLAST+ bin folder. By Default, the pipeline would try to retrieve the path from the '--blastp/-b' input path.
    
    --blast_weight BLAST_WEIGHT, -bw BLAST_WEIGHT: Path to weight file for the blast classifier. By default, the path would be defined in the 'definitions' module
    
    --blast_evalue EVALUE, -be EVALUE: blastp evalue cutoff

    --priam_weight PRIAM_WEIGHT, -pw PRIAM_WEIGHT: Path to blast weight for the priam classifier
    
    --efmap EF_MAP, -e EF_MAP: Path to Enzyme function class to Metacyc RXN-ID/EC Number file (efclasses.mapping).
    
    --threshold THRESHOLD, -th THRESHOLD: Threshold for voting results. Default is 0.5.
    
    --temp_folder TEMP_FOLDER, -tf TEMP_FOLDER: Specify the location of the temp folder that contains classifier results. By default it would be in the same directory of the output.
    
    --log LOG_PATH, -l LOG_PATH: Specify the location of the log file. By default it would be "runE2P2.log" in the temp folder.
    
    --protein_gene PROTEIN_GENE_PATH, -pg PROTEIN_GENE_PATH: Provide a protein to gene mapping file. This will be used to generate a splice variant removed fasta file and output our final version of e2p2.
    
    --verbose {0,1}, -v {0,1}: Verbose level of log output. Default is 0.
       0: only step information are logged
       1: all information are logged

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
