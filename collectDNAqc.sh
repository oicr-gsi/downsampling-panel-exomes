!#/bin/bash
module load samtools;

OUTDIR=190606_NB551056_0089_AHGVK3BGXB_precapture_pooling_round_1_Kayla
downsampling=`cat downsampling_bam_2.txt | head -1 | tr "\t" "\n" | tail -n +8`

for d in $downsampling; do
	echo $d
	mkdir -p $OUTDIR/$d
	less downsampling_bam_2.txt | tr -s "\t" "," | csvcut -c $d,Sample_ID | tr "," "\t" | column -t | awk '{print $1,$NF}' | tr " " "\t" > $OUTDIR/$d/$d.txt
done


CURRDIR=`pwd`
for DS in `ls $OUTDIR/Downsample*/Downsample*.txt`; do
samples=`less $DS | tail -n +2 | cut -f2`;
dirt=`echo $DS | cut -d/ -f2`
for s in $samples; do
frac=`cat $DS | grep $s | cut -f1`;
bam=`ls $CURRDIR/$OUTDIR/*${s}*.bam`;
echo "module load samtools/1.9; samtools view -b -s $frac $bam > $CURRDIR/$OUTDIR/$dirt/${s}.bam; samtools index $CURRDIR/$OUTDIR/$dirt/${s}.bam" > scripts/02_${s}_${dirt}_downsample.sh;
chmod +x $CURRDIR/scripts/02_${s}_${dirt}_downsample.sh;
qsub -V -l h_vmem=32G -q production -o $CURRDIR/logs -e $CURRDIR/logs -N ${s}_${dirt} $CURRDIR/scripts/02_${s}_${dirt}_downsample.sh
done
done

### cd into each downsampling directory ###

# run BamQC
mkdir -p BamQC
module load picard-tools/1.89
REF=/scratch2/groups/tgl/GBS-1313_old/hg19.fa
PERL=/.mounts/labs/seqprodbio/private/seqware/hsqwprod/provisioned-bundles/Workflow_Bundle_BamQC_2.5_SeqWare_1.1.0/Workflow_Bundle_BamQC/2.5/bin/perl-5.14.1/perl
BAMQC=/.mounts/labs/seqprodbio/private/seqware/hsqwprod/provisioned-bundles/Workflow_Bundle_BamQC_2.5_SeqWare_1.1.0/Workflow_Bundle_BamQC/2.5/bin/bamqc/bamqc.pl
TARG=/.mounts/labs/PDE/data/interval-files/MA38-pad10.bed
for BWAMEMBAM in `ls *.bam`; do
        SNAME=`echo $BWAMEMBAM | cut -d\. -f1`
        SORTEDBAM=$SNAME.sorted.bam
        java -Xmx8g -jar $PICARD_TOOLS_ROOT/SortSam.jar I=$BWAMEMBAM O=$SORTEDBAM SO=coordinate CREATE_INDEX=true TMP_DIR=picardTmp
        samtools view $SORTEDBAM | $PERL $BAMQC -s 1001 -i 1500 -q 30 -r $TARG -j BamQC/$SNAME.meta.json > BamQC/$SNAME.json;
done

# collect Picard metrics
mkdir -p PicardMetrics
TARG=/.mounts/labs/PDE/data/interval-files/MA38-pad10.bed
REF=/.mounts/labs/PDE/data/gatkAnnotationResources/hg19_random.fa
for BWAMEMBAM in `ls *.bam`; do
        SNAME=`echo $BWAMEMBAM | cut -d\. -f1`
        SORTEDBAM=$SNAME.sorted.bam
       java -jar $TOOLS/PicardTools/picard-tools-2.5.0/picard.jar CollectHsMetrics \
                I=$SORTEDBAM \
                O=PicardMetrics/$SNAME.hsmetrics \
                R=$REF \
                BAIT_INTERVALS=$TARG.interval_list \
                TARGET_INTERVALS=$TARG.interval_list
done
