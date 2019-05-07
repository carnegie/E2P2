# Ensemble Enzyme Prediction Pipeline (E2P2)

The Ensemble Enzyme Prediction Pipeline (E2P2) annotates protein sequences with Enzyme Function classes comprised of full, four-part Enzyme Commission numbers and MetaCyc reaction identifiers. It is the enzyme annotation pipeline used to generate the species-specific metabolic databases at the [Plant Metabolic Network](www.plantcyc.org) since 2013. E2P2 systematically integrates results from two molecular function annotation algorithms using an ensemble classification scheme. For a given genome, all protein sequences are submitted as individual queries against the base-level annotation methods.

## Getting Started
The following instuctions are for users to set up the E2P2 pipeline on a Unix machine and start running, developing and/or testing the pipeline.

### Prerequisites
* [Python 3](https://www.python.org/downloads/)
* [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [Java](https://www.oracle.com/technetwork/java/javase/downloads/index.html)
* [PRIAM_Search utility](http://priam.prabi.fr/REL_JAN18/index_jan18.html)

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
# The following are not required to run the most current version of E2P2
tar -xzf weights.tar.gz
tar -xzf maps.tar.gz
```

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
