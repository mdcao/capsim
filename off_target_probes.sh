#!/bin/bash

printUsage() {
    echo ""
    echo ""
    echo "  Usage: bash off_target_probes.sh [options] -b <bed> -r <fasta> -a <bam> -q <fastq>" 
    echo "  Options:"
    echo "      -b/--target-bed                  Bed file of the target regions. "
    echo "      -r/--reference                   Reference genome fasta file. "
    echo "      -a/--bam                         Bam file of CapSim simulated reads aligned to reference"
    echo "                                       genome. "
    echo "      -w/--window-size                 Window size for statistics of the depth of coverage of"
    echo "                                       the off target regions (default 1000)."
    echo "      -d/--min-depth                   Minimum depth of base coverage of the off target regions"
    echo "                                       to analyse (default 10000). "
    echo "      -x/--padding-size                Extension/padding size to the up and downstream of the "
    echo "                                       target regions (default 500). "
    echo "      -q/--probe-seq                   Text file containing the probe ID and sequence."
    echo "      -t/--threads                     Number of threads for alignment (default 1)."
    echo "      -p/--prefix                      Prefix of the output files (default ./out)."
    echo ""
}

TARGET_BED=""
REFERENCE=""
BAM=""
WINDOW_SIZE=1000
MIN_DEPTH=10000
PADDING_SIZE=500
PROBE_SEQ=""
THREADS=1
PREFIX_OUT="./out"

while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
        -b|--target-bed)
            TARGET_BED="$2"
            shift
            shift
            ;;
        -r|--reference)
            REFERENCE="$2"
            shift
            shift
            ;;
        -a|--bam)
            BAM="$2"
            shift
            shift
            ;;
        -w|--window-size)
            WINDOW_SIZE="$2"
            shift
            shift
            ;;
        -d|--min-depth)
            MIN_DEPTH="$2"
            shift
            shift
            ;;
        -x|--padding-size)
            PADDING_SIZE="$2"
            shift
            shift
            ;;
        -q|--probe_seq)
            PROBE_SEQ="$2"
            shift
            shift
            ;;
        -t|--threads)
            THREADS="$2"
            shift
            shift
            ;;
        -p|--prefix)
            PREFIX_OUT="$2"
            shift
            shift
            ;;
        *) # unknown option
            echo ""
            echo "!!!ERROR: unknown options $1"
            printUsage
            exit
            ;;
    esac
done

if [[ $TARGET_BED == "" ]]
then
    echo ""
    echo "!!!ERROR: A bed file for targe regions is required: -b/--target-bed"
    printUsage
    exit
fi

if [[ $REFERENCE == "" ]]
then
    echo ""
    echo "!!!ERROR: A reference genome fasta file is required: -r/--reference"
    printUsage
    exit
fi

if [[ $BAM == "" ]]
then
    echo ""
    echo "!!!ERROR: A bam file XXX is required: -a/--bam"
    printUsage
    exit
fi

if [[ $PROBE_SEQ == "" ]]
then
    echo ""
    echo "!!!ERROR: A text file for probe sequences is required: -q/--probe-seq"
    printUsage
    exit
fi

echo ""
echo "  Options provided:"
echo "      TARGET BED     = $TARGET_BED"
echo "      REFERENCE      = $REFERENCE"
echo "      BAM            = $BAM"
echo "      WINDOW SIZE    = $WINDOW_SIZE"
echo "      MIN DEPTH      = $MIN_DEPTH"
echo "      PADDING SIZE   = $PADDING_SIZE"
echo "      PROBE SEQ      = $PROBE_SEQ"
echo "      THREADS        = $THREADS"
echo "      PREFIX OUT     = $PREFIX_OUT"

STEP=1

echo "["`date`"]" STEP $((STEP++)): Generate reference index
samtools faidx $REFERENCE 

echo "["`date`"]" STEP $((STEP++)): Generate genome file for bedtools
awk '{print $1"\t"$2}' "$REFERENCE".fai > `echo $(basename $REFERENCE) | sed 's/.fasta$//g'`.genome

echo "["`date`"]" STEP $((STEP++)): Generate fasta file of probe sequences
tail -n +2 $PROBE_SEQ | awk '{printf(">%s:%s\n%s\n",$1,$2,$3)}' > `echo $(basename $PROBE_SEQ) | sed 's/.txt$//g'`.fa

echo "["`date`"]" STEP $((STEP++)): Generate target bed file with extra regions upstream and downstream
bedtools slop -i $TARGET_BED -g `echo $(basename $REFERENCE) | sed 's/.fasta$//g'`.genome -b $PADDING_SIZE > `echo $(basename $TARGET_BED) | sed 's/.bed$//g'`_"$PADDING_SIZE"bp.bed

echo "["`date`"]" STEP $((STEP++)): Generate a bed file of regions in a genome that are not covered by the target bed file
bedtools complement -i `echo $(basename $TARGET_BED) | sed 's/.bed$//g'`_"$PADDING_SIZE"bp.bed -g `echo $(basename $REFERENCE) | sed 's/.fasta$//g'`.genome | awk -v w="$WINDOW_SIZE" '{num=($3-$2)/w; for(i=0;i<num-1;i++) print $1"\t"($2+w*i)"\t"($2+w*(i+1)-1); if(w*int(num)!=$3-$2) print $1"\t"($2+w*int(num))"\t"$3;}' > `echo $(basename $TARGET_BED) | sed 's/.bed$//g'`_complement.bed

echo "["`date`"]" STEP $((STEP++)): Generate a bam file of alignments which overlap with the bed file
bedtools intersect -abam $BAM -b `echo $(basename $TARGET_BED) | sed 's/.bed$//g'`_complement.bed > "$PREFIX_OUT"_complement.bam

echo "["`date`"]" STEP $((STEP++)): Generate an alignment index
samtools index "$PREFIX_OUT"_complement.bam "$PREFIX_OUT"_complement.bai

echo "["`date`"]" STEP $((STEP++)): Generate a base coverage
samtools bedcov `echo $(basename $TARGET_BED) | sed 's/.bed$//g'`_complement.bed "$PREFIX_OUT"_complement.bam | awk -v d="$MIN_DEPTH" '{if ($4>d) {print $1"\t"$2"\t"$3"\t"}}' > "$PREFIX_OUT"_off_target_regions.bed

echo "["`date`"]" STEP $((STEP++)): Generate a fasta sequence file for off_target regions
bedtools getfasta -fi $REFERENCE -bed "$PREFIX_OUT"_off_target_regions.bed -fo "$PREFIX_OUT"_off_target_regions.fas

echo "["`date`"]" STEP $((STEP++)): Generate index sequences for the fasta sequence
bwa index -a bwtsw "$PREFIX_OUT"_off_target_regions.fas 

echo "["`date`"]" STEP $((STEP++)): Generate an alignement of off_target regions
bwa mem -t $THREADS -Y -c 1000 "$PREFIX_OUT"_off_target_regions.fas `echo $(basename $PROBE_SEQ) | sed 's/.txt$//g'`.fa | samtools view -@"$THREADS" -bS | samtools sort -@"$THREADS" -o "$PREFIX_OUT"_off_target_regions.bam && samtools index "$PREFIX_OUT"_off_target_regions.bam "$PREFIX_OUT"_off_target_regions.bai

echo "["`date`"]" STEP $((STEP++)): Generate a list of probes contributing to off_target alignment
samtools view -F4 "$PREFIX_OUT"_off_target_regions.bam | cut -f1 | sort -u > "$PREFIX_OUT"_off_target_probes.txt


