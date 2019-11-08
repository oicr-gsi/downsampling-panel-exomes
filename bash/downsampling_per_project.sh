project=$1
downsample_fraction=$2
frac_string=`echo ${downsample_fraction} | tr "." "_"`
workdir=/scratch2/groups/tgl/GBS-1572/$project; mkdir -p $workdir
outdir=$workdir/output; mkdir -p $outdir
prov=/.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz
# get all bwa mem bam files
input=$workdir/input_logs_bmpp.txt
if [[ ! -f $input ]]; then
  zgrep $project $prov | grep "BamMerge" | cut -f14,19,47 | grep "bam" | sort | uniq | tr "\t" "," | grep -v "_M" | grep -v "bai" | grep -E "EX|TS" > $input
fi
# compute base-metrics
job_string=""
for record in `cat ${input}`; do
  sample_name=`echo ${record} | cut -d, -f1`
  run_id=`echo ${record} | cut -d, -f2`
  bam_loc=`echo ${record} | cut -d, -f3`
  BASEDIR=${workdir}
  echo ${sample_name}
  output=$outdir/$project/${sample_name}; mkdir -p $output
  scriptDir=${output}/script; mkdir -p $scriptDir
  logDir=${output}/logs; mkdir -p $logDir
  tmp=${output}/tmp; mkdir -p $tmp
  # set up various downsampling fractions
  # make ini
  INI=${BASEDIR}/${sample_name}_downsample_${frac_string}.ini
  echo "#!/bin/bash" > $INI
  echo "project_name=\"$project\"" >> $INI
  echo "sample_name=\"${sample_name}\"" >> $INI
  echo "input_bam=\"${bam_loc}\"" >> $INI
  echo "downsample_fraction=${downsample_fraction}" >> $INI
  echo "seed=1" >> $INI
  echo "temp_dir=$tmp" >> $INI
  echo "output_dir=$output" >> $INI
  echo "frac_string=${frac_string}" >> $INI
  echo "output_bam=${output}/${sample_name}_downsampled_to_${frac_string}.bam" >> $INI
  echo "sorted_bam=${output}/${sample_name}_downsampled_to_${frac_string}.sorted.bam" >> $INI
  echo "output_picard_metrics=${output}/${sample_name}_downsampled_to_${frac_string}.hsmetrics" >> $INI
  echo "bamqc_json=${output}/${sample_name}_downsampled_to_${frac_string}.json" >> $INI
  # echo "SETTINGS=${BASEDIR}/downsample_settings.ini" >> $INI
  script=$scriptDir/${sample_name}_${frac_string}.sh
  echo "#!/bin/bash" > $script
  echo "SETTINGS=/scratch2/groups/tgl/GBS-1572/Downsampling/bash/downsample_setings.ini; source /scratch2/groups/tgl/GBS-1572/Downsampling/bash/design_functions.sh \$SETTINGS $INI"  >> $script
  echo "# start stringing" >> $script
  echo "00_checkINI" >> $script
  echo "01_downsample_bam" >> $script
  echo "02_runBamQC" >> $script
  echo "03_collectPicardMetrics" >> $script
  chmod +x $script
  jname=${sample_name}_${frac_string}
  qsub -V -l h_vmem=32G -N ${jname} -e $logDir -o $logDir $script
  job_string="${job_string}${job_string:+,}${jname}"
done

echo "#!/bin/bash" > $workdir/summarize_init_${frac_string}.sh
echo "for hsmetrics in \`ls /scratch2/groups/tgl/GBS-1572/${project}/output/${project}/*/*${frac_string}.hsmetrics\`; do sample_name=\`basename \${hsmetrics} \".hsmetrics\"\`; cat \$hsmetrics | tail -n +7 | head -1 | awk '{print \"sample_name\",\$0}' | tr \" \" \"\t\"; cat \$hsmetrics | tail -n +8 | head -1 | awk -v var=\${sample_name} '{print var,\$0}' | tr \" \" \"\t\"; done | sort | uniq | sort -r > $workdir/initQC_${frac_string}.txt" >> $workdir/summarize_init_${frac_string}.sh
echo "sort -r $workdir/initQC_${frac_string}.txt > .tmp_${frac_string}; mv .tmp_${frac_string} $workdir/initQC_${frac_string}.txt" >> $workdir/summarize_init_${frac_string}.sh
# echo "cp $workdir/initQC.txt /.m" >> $workdir/summarize_init.sh
chmod +x $workdir/summarize_init_${frac_string}.sh
echo "#!/bin/bash" > run_${frac_string}.sh
echo "qsub -hold_jid ${job_string} -cwd $workdir/summarize_init_${frac_string}.sh" >> run_${frac_string}.sh
chmod +x run_${frac_string}.sh
./run_${frac_string}.sh
