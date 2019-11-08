module load samtools;

SAMTOOLS=/.mounts/labs/PDE/modulator/sw/Debian8.11/samtools-1.9/bin/samtools
# bam csv identified as project_name,sample_id,path_to_bam_file
BAMCSV=
# DOWNSAMPLING RATIO
DOWNSAMPLE=
# OUTPUT DIR
WORKDIR= # path to working dir
# ARGUMENT TO SIMULATE EXPECTED RESULTS
SIM= # yes/no
# interval bed file
INTERVAL_BED=
SEED= # set seed
NUM_ITERS= # number of iterations

# set up the working directory
mkdir -p $WORKDIR
TMP=$WORKDIR/tmp
mkdir -p $TMP
EXPECTED=$WORKDIR/expected_ouput
mkdir -p $WORKDIR
if [ $SIM == "no" ]; then
  DOWNSAMPLED_TMP=$TMP/downsampled_bam # this is where the downsampled bam files sit
  mkdir -p $DOWNSAMPLED_TMP
  OBSERVED=$WORKDIR/observed_output # this is where all the final QC metrics are collected
  mkdir -p $OBSERVED
fi

# first get expected results
for record in `cat $BAMCSV`; do
  project_name=`echo $record | cut -d, -f1`
  sample_id=`echo $record | cut -d, -f2`
  original_bam=`echo $record | cut -d, -f3`
  # compute original QC metrics
  # Read these metrics in R script; divide by the dsFraction
  # write this to file $EXPECTED/downsample
done

# read downampling (set seed to downsampling index)
seed=0
for dsFrac in `cat $DOWNSAMPLE`; do
  seed=$(( $seed + 1))
  downsample_dir_name=`echo $dsFrac | tr "." "_"`
  downsample_bam_loc=$DOWNSAMPLED_TMP/"downsample_to_frac_"${downsample_dir_name}
  for record in `cat $BAMCSV`; do
    project_name=`echo $record | cut -d, -f1`
    sample_id=`echo $record | cut -d, -f2`
    original_bam=`echo $record | cut -d, -f3`

    # launch downsampling jobs









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
done
