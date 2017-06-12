# capsim

Capsim is a tool to simulate capture sequencing. It was developed as part
of the Japsa package. 

##Quick installation guide:
On a Linux/Mac machine with `make' and `git' installed, the software can be installed with

    git clone https://github.com/mdcao/japsa
    cd japsa
    make install \
      [INSTALL_DIR=~/.usr/local \] 
      [MXMEM=7000m \] 
      [SERVER=true \]


Details of installation (including for Windows) and usage of Japsa can be found 
in its documentation hosted on [ReadTheDocs](http://japsa.readthedocs.org/en/latest/index.html) 

##Usage

Before running capsim, the probe sequences need to be aligned to the reference sequence (or the genome sequence to simulate sequencing). We recommend using bowtie2 for the alignment.

After alignment, sort the and index the bam file

Capsim takes the bam file as the input:

    jsa.sim.capsim --reference ref.fasta --probe probe.sam --ID ID_OF_THE_DATASET --fmedian 5000  --pacbio output --pblen 3000 --num 20000000
    


