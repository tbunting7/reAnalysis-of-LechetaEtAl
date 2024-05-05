# reAnalysis-of-LechetaEtAl
This is a small repository acting as a lab notebook during my introduction to RNA-Seq .fastq analysis

[MobaXterm](https://mobaxterm.mobatek.net/) was used on Windows to access a remote server. I can definitely recommend this, it was quite intuitive once I got past the learning curve.

Links I've found useful while learning:

[Intro to Unix and navigation in the command terminal](http://korflab.ucdavis.edu/Unix_and_Perl/current.html)

[Regular expressions for terminal](https://www.regular-expressions.info/quickstart.html)

[FastQ file format info](https://learn.gencore.bio.nyu.edu/ngs-file-formats/fastq-format/)

[FastQC and FastQ format](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon-flipped/lessons/05_qc_running_fastqc_interactively.html)

[STAR software manual](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)

In the repository, there is a command_terminalNotes.txt file that has my notes while I was learning, as well as some commands (with explanations) used for FastQC analysis and STAR indexing and mapping.
There is also a .r file which served as the beginning of differential expression analysis, though not much has been done with it at this point. 
There are also the .html files from FastQC analysis.

# Methods
The steps taken to analyse the date involve: FastQC analysis -> STAR indexing -> STAR mapping -> differential expression analysis. All of this was done to a set of 9 fastq files, the result of RNAseq data, each of which having ~10 millions samples and a read length of 100bp.

FastQC analysis is the first step, and it is done to ensure the quality of the data gathered before further analysis is done. As the name implies, it is a quality control step, meant to identify any problems with the read. These results can be found in the Lecheta_fastqc directory. The results of this were good enough to move onto the next step, the primary factor being that there was a high degree of confidence for each read.

STAR indexing is the next step, done to prepare a genome in order to map the reads to it. This is done because it vastly speeds up the time to map. It uses a known genome to create a new set of files that make it much more efficient to map the read data to, as compared to "manually" looking through the ~180 million basepairs of the genome. For this data, the most recent _Drosophilia Melanogaster_ genome was used.

Once the genome is indexed and prepared with STAR, mapping is ready. At the end of this step, a file is created which has the number of unique reads (from the RNAseq data) that match to any gene in the genome. This file is essentially a table of genes and the number of times that a sequence of RNA was found that matched a gene. 
Additional filters are also put in place that cut out large portions of non-unique reads. These are reads that match to multiple regions of the genome, a primary example being abundant ribosomal RNA.

Once these files are obtained, they can be analyzed. This is what the DESeq2_JAD_counts_TB.r file is. It was modified from a template file to suit the data, and it imports the libraries need to run analysis. After some configuration, parameters were used to compare the different files obtained from mapping. These parameters were temperature, being the cold and heat shock temperatures as well as the room temperature control. This step was not fully completed by me, as I ran into an error from mapping, but I did create a multi-dimentional scaling (MDS) plot. Putting aside the error, this plot was meant to simplify the thousands of genes present into just a single data point from each file (9 files). It essentially creates two new axis to place the data points in a way such that their distance from each other best represents how dissimilar they are. Basically, if the points are closer together on this graph, they are more similar, even if the axes they are placed on tell you nothing about the variables. 
Interestingly, despite the errors in that mapping data, I was able to find clumping of the data points relative to the temperature; the room temperature data was clustered on one side of the graph, with the heat/cold shocked data on the other.
