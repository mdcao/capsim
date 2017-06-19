# capsim

Capsim is a tool to simulate capture sequencing. It was developed as part
of the Japsa package. 

## Quick installation guide:
On a Linux/Mac machine with `make' and `git' installed, the software can be installed with

    git clone https://github.com/mdcao/japsa
    cd japsa
    make install \
      [INSTALL_DIR=~/.usr/local \] 
      [MXMEM=7000m \] 
      [SERVER=true \]


Details of installation (including for Windows) and usage of Japsa can be found 
in its documentation hosted on [ReadTheDocs](http://japsa.readthedocs.org/en/latest/index.html) 

## Usage

Before running capsim, the probe sequences need to be aligned to the reference sequence (or the genome sequence to simulate sequencing). We recommend using bowtie2 for the alignment.
   
    #Skip this step if the index has been generated
    bowtie2-build ref.fasta ref
    
    #Align the probes into the reference
    bowtie2 --local --very-sensitive-local --mp 8 --rdg 10,8 --rfg 10,8 -k 10000 -f -x ref -U probes.fa -S probes.sam

Note: for some reason, bowtie2 only accepts the query fasta file (probes.fa) containing one sequence per line.

After alignment, sort the and index the bam file with samtools

     samtools view -bSU probes.sam | samtools sort -o probes.bam -

Capsim takes the bam file as the input:

    jsa.sim.capsim --reference ref.fasta --probe probes.bam --ID someid --fmedian 5000  --pacbio output --pblen 3000 --num 20000000

or 

    jsa.sim.capsim --reference ref.fasta --probe probes.bam --ID someid --fmedian 500  --miseq output --illen 300 --num 20000000



