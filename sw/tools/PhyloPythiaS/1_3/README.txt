INSTALLATION
------------

These are the installation instructions for the PhyloPythiaS package.
If you have difficulties please contact us at phylopythias@uni-duesseldorf.de

SOFTWARE REQUIREMENTS
----------------------

Unix-like system [tested on Debian]
glib-2.6.0+
Ruby 1.8+ [Not tested with 1.9]
SQLite3
libsqlite3-ruby/libsqlite3-ruby1.8 (depending upon your ruby version)

On a debian system you can install the above packages using the following command (super user rights required);

	apt-get install package-names

for example all the packages can be installed with;
	
	apt-get install ruby sqlite3 libsqlite3-ruby


HARDWARE REQUIREMENTS
----------------------

~5Gb disk space
~4Gb RAM

You don't need knowledge of databases to use this software, downloading of NCBI taxonomy files and database creation
have been automated. All you need to do is to provide an empty directory where the corresponding files will be stored.

We will assume that commands are being issued from a unix bash shell.

UNPACK BINARIES
-----------------

This package comes with two sets of binary files; bin32.tar.gz and bin64.tar.gz
Please unpack the correct binary file for your system configuration. For example, if you have a 64 bit system;

	tar xvfz bin64.tar.gz

This will create a folder bin with binaries inside. Please make sure that all the binaries have execute permission;

	chmod +x bin/*

THE EASY WAY
-------------

The easiest way to use this package is to change the "sample_config.txt" file in the "examples" directory 
and then run "train.rb" which is the "scripts" directory. For this you will need three empty directories and
a file with either a Newick tree or a list of NCBI taxon identifiers to model;

1. NCBI_DATA_DIR: a directory to store processed NCBI sequence data
2. NCBI_TAX_DIR: a directory to store the NCBI taxonomy and corresponding database
3. PROJECT_DIR: a directory to store a project specific data, like the modeling tree, training data, learned models etc.
4. TREE_FILE: a file containing either a Newick file or a nodes list file with one entry per line.

The first two directories are reusable, i.e. they can be used for different projects. You will need to prepare
a file with list of NCBI Taxonomy identifiers that you want to model (one id per line). You just have to 
give the most specific nodes and rest of the lineage will be extracted automatically. Once you have 
these things ready put this information in the configuration file and run (from this directory);

	scripts/train.rb configuration_file

You will be prompted for some questions, answer them appropriately and you will get the models. The program can run
completely un-attended if you provide empty directories. On our machine with 64-bit linux, 100Gb free disk space and
8Gb RAM, with no competing processes, the default option training takes ~12 hours.

Then run the "predict.rb" script as following to get the predictions using these models (make sure to keep the 
configuration file same);

	scripts/predict.rb fasta_file_to_predict configuration_file

The output will be generated in the same directory where the fasta file is located with an extension ".out" for
raw output and ".PP.out" for PhyloPythia like formatted output. If mentioned in the configuration file (PIE_CHARTS:TRUE)
the program will produce SVG pie charts in the same directory for every taxonomic rank with extension "pie_aRank.svg".
You can use a free software like InkScape (http://www.inkscape.org) to view and manipulate these charts.

It is possible to use sample specific data if available to learn the models. Please see "INCORPORATING SAMPLE SPECIFIC 
DATA" section in this file for more details on this. All of the sample specific data is used in the learning process.

BIT HARDER WAY
--------------

You will want to go the harder way if you want to include sample specific data or change the modeling tree etc.

1. Creating the reference sequence library

The train.rb script can automatically download and process Bacterial and Archaeal genomes from NCBI. But you can 
always do things manually and fine tune it.

A. Download sequences from NCBI

Create a new directory (lets say NCBI_sequences), change to this directory and download sequences from 
NCBI inside this directory. You can use wget command for download. This might take some time as the 
download size is several giga-bytes.

	a. Download the Bacterial and Archaeal genomes and uncompress

	wget ftp://ftp.ncbi.nih.gov/genomes/Bacteria/all.gbk.tar.gz
	tar xvfz all.gbk.tar.gz

	b. Download the WGS data and uncompress [optional but recommended, this takes a lot of time and disk-space]

	wget ftp://ftp.ncbi.nih.gov/genbank/wgs/*.gbff.gz
	gunzip *.gz

B. Process this raw data 

Create another directory (somewhere outside the directory where you downloaded NCBI sequence data) 
to store the processed data (lets say NCBI_processed). Run the script process_ncbi.rb as following 
(located in the scripts directory)

	./process_ncbi.rb path_to_NCBI_sequences path_to_NCBI_processed

You can delete the raw data directory NCBI_sequences if you do not need it yourself.

2. Create the NCBI taxonomy database

Now we use the SQLite3 database for storing and querying the NCBI taxonomy information (before we used MySQL 
which needed some additional user efforts, thanks to Johannes Droege for suggesting this). The program assumes 
that the NCBI taxonomy data is stored in a SQLite3 database file named "ncbitax_sqlite.db" This can be created by
using the "ncbitax2sqlite.rb" script in the scripts directory. It needs a directory name with NCBI taxonomy dump 
files and will download it if the files are not available. If you have NCBI taxonomy dump files (names.dmp and 
nodes.dmp) then you can point the script to the folder containing these. Lets say that the "ncbitax_sqlite.db" file
is in the directory NCBI_TAX_DATA (which also contains the .dmp files). Run the script ncbitax2sqlite.rb as following;

	./ncbitax2sqlite.rb path_to/NCBI_TAX_DATA path_to/NCBI_TAX_DATA/ncbitax_sqlite.db

3. Create a newick tree

If you already have a newick tree (see RESTRICTIONS below) for the model you want to build then skip this step. 
Put all the clades you want to model in a file (one clade per line), lets say "clades.txt". 
You don't have to put higher level clades in the files as they will be automatically determined. 
Use the provided script ncbi2newick as following;  

	./ncbi2newick.rb clades.txt NCBI_TAX_DATA/ncbitax_sqlite.db > tree.newick

4. Generate models

For this you will have to create a configuration file, a sample is provided with this installation.
Please make appropriate changes to the configuration file and run the "train.rb" script with the configuration 
file as a parameter.
	
	./train.rb configuration_file

INCORPORATING SAMPLE SPECIFIC DATA
----------------------------------

Incorporating sample specific data into the models is a two step procedure.

A. Changing the tree

If the clade(s) for the sample specific data is already in the newick tree then you can skip this step. 
Otherwise, you should include these clade(s) in one of the following ways;

	1. By adding them to the clades list and providing it to the configuration file (TREE_FILE). This will work only if the 
	clades are standard NCBI taxon identifiers and have one of the taxonomic ranks provided in the configuration file. 
	In this case please make sure that you give correct taxon id to these clades and not some random NCBI taxon id.

	2. If the previous condition of standard NCBI taxon identifiers is not met, you can manually change a newick tree 
	(which can be generated using the "ncbi2newick.rb" script from a clades list) to include the sample specific clades with 
	some distinct names than already existing nodes in the tree, making sure that they can be traced from the root node. 
	If you do this then PLEASE DO NOT USE PP (PhyloPythia) like output and pie charts, as the NCBI taxon names might get confused 
	if the nodes you added correspond to some valid NCBI taxon id. 

B. Adding the sample specific sequence data

You should create a directory and put the sample specific sequences as fasta files in this directory (see RESTRICTIONS).
Create one fasta file per clade, with a name starting with an id where the sequences should belong to in the tree.
For example 286.fna or 286.1.fna for the sequences that belong to the node with id 286. Provide the name of
this directory via the configuration file (SAMPLE_SPECIFIC_DIR). Please make sure to keep only fasta files inside this directory.

WHAT CLADES TO MODEL
--------------------

If you have some expectations of the contents of your sample (e.g. from 16S rRNA survey) then you should include these clades. 
Also you should include some closely related clades if reference data is available.

It is possible that in some cases there is no reference data available (either in the public domain or sample specific) for some clade, in this case modeling this clade makes no sense. So, there should be either public data or sample specific data available for all the leaf nodes in the modeling tree. In general, we recommend that there should be at-least 100 Kilo-base of data for sample specific clades.


RESTRICTIONS
------------

The binary "fasta2kmers" cannot read a sequence which has more than 50,000 characters in one line. 
So you should fold the sequence into multiple lines before using this binary.

For a FASTA file, please take care that there are no TABs ("\t") in the headers.

The newick tree must have a name for each and every node and that must be a number (an integer). Ideally the number 
should be NCBI taxonomic identifier, so that we can look it up in the database. Further the tree must be 
in a single line and rooted (you can use 1 as the root node, which is also the root of the NCBI taxonomy). 
See the example newick tree provided. 

All of the sample specific data is used when provided via the configuration file. The fragments are created with
a window moving over the sequences (as defined by SAMPLE_SPECIFIC_STEP). 
So if you have too little or too much data then this can lead into
either very little or too much training data. Manual inspection leads to better models in most of the cases so
contact the author if you are looking for collaborative efforts.

Current version reads the complete test file in memory so the number of examples that can be classified at a time is 
restricted by the amount of main memory available on your system. So if you have a very large sample to classify then
please divide it into smaller parts that can fit in your system's memory.


ADDITIONAL INFORMATION
-----------------------

1. Scaffold-contig consistency analysis (this is provided for user convenience, not much testing has been performed)

Here the aim is to analyze the binning consistency of the scaffolds. In this case you know a bunch of contigs and their affiliations to some scaffolds. Create a fasta file with these contigs (not scaffolds) such that each header is in the format

	Scaffold_m_contig_n

Here m is the scaffold identifier and n is the contig identifier. Use PhyloPythiaS to make assignments on this file as described above. Then provide the assignments file (the PP format) the fasta file and the newick tree used to the "scaffold_consistency.rb" with some more needed information. Run the script without any arguments to find out the parameters;

	scripts/scaffold_consistency.rb

Running the script with appropriate parameters will give you the scaffold-contig consistency on the standard output.

If you ask the script to create graph files (parameter "graph") then it will create some more files (actually 4) with extension "_graph_*". Please note that if you want create graph files then the population (parameter "pop") must be provided. These created files then can be passed (only the prefix without "_graph_*") to the R script "scaffold_consistency_barplot.R" to generate the plot. Please check the R script for additional information.

2. The ".out" output file (created by predict.rb) contains sequence ids as the first and assigned taxonomic id as the fifth column. By default taxonomic id "1" means the root node and "-1" means unassigned. Unassignments can occur if the classifiers in an ensemble disagree completely except at the root.
