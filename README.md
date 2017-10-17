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

Options,

<pre>
      --reference=s                    Name of genome to be
                                       (REQUIRED)
      --probe=s                        File containing probes mapped to the reference in bam format
                                       (default='null')
      --logFile=s                      Log file
                                       (default='-')
      --ID=s                           A unique ID for the data set
                                       (default='')
      --miseq=s                        Name of read file if miseq is simulated
                                       (default='null')
      --pacbio=s                       Name of read file if pacbio is simulated
                                       (default='null')
      --fmedian=i                      Median of fragment size at shearing
                                       (default='2000')
      --fshape=d                       Shape parameter of the fragment size distribution
                                       (default='6.0')
      --smedian=i                      Median of fragment size distribution
                                       (default='1300')
      --sshape=d                       Shape parameter of the fragment size distribution
                                       (default='6.0')
      --tmedian=i                      Median of target fragment size (the fragment size of the data).
                                       If specified, will override fmedian and smedian.
                                       Othersise will be estimated
                                       (default='0')
      --tshape=d                       Shape parameter of the effective fragment size distribution
                                       (default='0.0')
      --num=i                          Number of fragments
                                       (default='1000000')
      --pblen=i                        PacBio: Average (polymerase) read length
                                       (default='30000')
      --illen=i                        Illumina: read length
                                       (default='300')
      --ilmode=s                       Illumina: Sequencing mode: pe = paired-end, mp=mate-paired and se=singled-end
                                       (default='pe')
      --seed=i                         Random seed, 0 for a random seed
                                       (default='0')
      --help                           Display this usage and exit
                                       (default='false')
</pre>

CapSim will output sequence reads in fastq format. Users can perform subsequenct analysis by aligning the simulated reads to the reference genome.  

## Off-target analysis

A script file _off_target_probes.sh_ used to identify off target probes is provided in this repository. To run this script file,

    bash off_target_probes.sh -b target_regions.bed -r ref.fasta -a cap_sim.bam -w 1000 -d 10000 -x 500 -q probes.txt -t 4 -p out

where,

<pre>

      -b/--target-bed                  Bed file of the target regions.
      -r/--reference                   Reference genome fasta file.
      -a/--bam                         Bam file of CapSim simulated reads aligned to reference
                                       genome.
      -w/--window-size                 Window size for statistics of the depth of coverage of
                                       the off target regions (default 1000).
      -d/--min-depth                   Minimum depth of base coverage of the off target regions
                                       to analyse (default 10000).
      -x/--padding-size                Extension/padding size to the up and downstream of the
                                       target regions (default 500).
      -q/--probe-seq                   Text file containing the probe ID and sequence.
      -t/--threads                     Number of threads for alignment (default 1).
      -p/--prefix                      Prefix of the output files (default ./out).
</pre>
    
The following tools should be installed and added to system path.

* [samtools](http://samtools.sourceforge.net/)
* [bwa](http://bio-bwa.sourceforge.net/)
* [bedtools](https://github.com/arq5x/bedtools2)

**_The script file consists of 12 main commands which could be run step by step._**

1\. generate reference index file,
    
    samtools faidx ref.fasta

2\. generate genome file for bedtools,

    awk '{print $1"\t"$2}' ref.fasta.fai > ref.genome

3\. generate fasta file of probe sequences,

    tail -n +2 probes.txt | awk '{printf(">%s:%s\n%s\n",$1,$2,$3)}' > probes.fa

4\. generate target bed file with extra regions upstream and downstream,

    bedtools slop -i target_regions.bed -g ref.genome -b 500 > target_regions_500bp.bed
    
5\. generate a bed file of regions in a genome that are not covered by the target bed file,

    bedtools complement -i target_regions_500bp.bed -g ref.genome | awk -v w=100 '{num=($3-$2)/w; for(i=0;i<num-1;i++) print $1"\t"($2+w*i)"\t"($2+w*(i+1)-1); if(w*int(num)!=$3-$2) print $1"\t"($2+w*int(num))"\t"$3;}' > target_regions_complement.bed
    
This file will be used to calculate the depth of coverage of the bam file across the off target regions. Window size other than 1Kb could be used here. A smaller window size will generally result in more precise statistics but will be more time-consuming. The window size could be changed by the awk parameter _w_.
 
6\. generate a bam file of alignments which overlap with the bed file,

    bedtools intersect -abam cap_sim.bam -b target_regions_complement.bed > out_complement.bam

7\. generate an alignment index,
    
    samtools index out_complement.bam out_complement.bai
    
8\. generate a base coverage,

    samtools bedcov target_regions_complement.bed out_complement.bam | awk -v d="10000" '{if ($4>d) {print $1"\t"$2"\t"$3"\t"}}' > out_off_target_regions.bed
    
The threshold of the base coverage filtering could be specified by the awk parameter _d_.

9\. generate a fasta sequence file for off_target regions,

    bedtools getfasta -fi ref.fasta -bed out_off_target_regions.bed -fo out_off_target_regions.fas
    
10\. generate bwa index files for the off target region fasta file,

    bwa index -a bwtsw out_off_target_regions.fas

11\. align the probe sequences to the off target region fasta file,

    bwa mem -t 4 -Y -c 1000 out_off_target_regions.fas probes.fa | samtools view -@4 -bS | samtools sort -@4 -o out_off_target_regions.bam && samtools index out_off_target_regions.bam out_off_target_regions.bai
    
12\. extract the off target alignment probe sequences,

    samtools view -F4 out_off_target_regions.bam | cut -f1 | sort -u > out_off_target_probes.txt
    
