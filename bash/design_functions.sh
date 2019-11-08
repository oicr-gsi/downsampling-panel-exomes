#!/bin/bash
SETTINGS=$1
source ${SETTINGS}
SAMPLE_DOWNSAMPLE_INI=$2
source ${SAMPLE_DOWNSAMPLE_INI}

echo "Printing settings ini params .... \n \n"
for entry in `cat $SETTINGS`; do
  echo $entry
done

echo "Printing downsampling ini params .... \n \n"
for rec in `cat ${SAMPLE_DOWNSAMPLE_INI}`; do
  echo $rec
done

# utils

function 0r_generateIntervals {
  source ${SETTINGS}
  source ${SAMPLE_DOWNSAMPLE_INI}

  TARG=$1

  $JAVA -jar ${PICARD} BedToIntervalList \
  I=${TARG} \
  O=${TARG}.interval_list \
  SD=${REF_DICT}
}

function 0r_generateRefDict {
  source ${SETTINGS}
  source ${SAMPLE_DOWNSAMPLE_INI}
  REF=$1
  $JAVA -jar ${PICARD} CreateSequenceDictionary \
  R=${REF} \
  O=${REF_DICT}
}

# check if all ini files are exisitent
function 00_checkINI {
  source ${SETTINGS}
  source ${SAMPLE_DOWNSAMPLE_INI}
  # function to check setting up of all ini params
  if [[ ! -f ${input_bam} ]];then
    echo -e "Missing input bam file"
  fi

  if [[ -z ${sample_name} ]];then
    echo -e "enter sample name or setting sample name to the name of input bam file "
    sample_name=`basename ${input_bam}`
  fi

  if [[ -z ${project_name} ]];then
    echo -e "enter project name or setting sample name to \"test\" "
    project_name="test"
  fi

  if [[ -z ${downsample_fraction} ]];then
    echo -e "No downsampling fraction provided. Setting it to 1.00"
    downsample_fraction=1
    frac_string="none"
  fi

  if [[ -z ${output_bam} ]];then
    echo -e "No path for output bam file is provided. Creating a new output bame file name"
    output_bam=${output_dir}/${sample_name}_downsampled_to_${frac_string}.bam
  fi

  if [[ -z ${output_picard_metrics} ]];then
    echo -e "No path for output_picard_metrics file is provided. Creating a new output_picard_metrics file name"
    output_picard_metrics=${output_dir}/${sample_name}_downsampled_to_${frac_string}.hsmetrics
  fi

  if [[ -z ${bamqc_json} ]];then
    echo -e "No path for bamqc_json file is provided. Creating a new bamqc_json file name"
    bamqc_json=${output_dir}/${sample_name}_downsampled_to_${frac_string}.json
  fi

  if [[ -z ${seed} ]];then
    echo -e "No seed provided. Setting seed to 1"
    seed=1
  fi

  if [ -z $REF ] || [ ! -f ${REF} ]; then
    echo -e "No reference file invalid"
  fi

  if [ -z ${REF_DICT} ]; then
    if [ ! -f ${REF_DICT} ] ; then
    echo -e "Generating reference dict for $REF \n This process will take a while ... "
    0r_generateRefDict $REF
    fi
    echo -e "Provide path to REF_DICT "
  fi

  if [ ! -f ${TARG}.interval_list ]; then
    echo -e "Generating target intervals for the target bed file \n This process will take a while ... "
    0r_generateIntervals $TARG
  fi
}

function 01_downsample_bam {
  module load samtools/1.9
  source ${SETTINGS}
  source ${SAMPLE_DOWNSAMPLE_INI}
  # input_bam=$1
  # downsample_fraction=$2
  #statements
  if [[ ! -f ${input_bam} ]]; then
    echo "no input bam provided to downsample"
  elif [ -z ${downsample_fraction} ] || [ ${downsample_fraction} == 1.0 ]; then
    echo "No downsampling will be performed"
    echo "Launching samtools flagstat for this input ..."
    $SAMTOOLS flagstat ${input_bam} > ${output_dir}/${sample_name}_flagstat.txt;
    $JAVA -jar $PICARD CollectHsMetrics \
             I=${input_bam} \
             O=${output_picard_metrics} \
             R=$REF \
             BAIT_INTERVALS=${TARG}.interval_list \
             TARGET_INTERVALS=${TARG}.interval_list;
    # exit -1
  else
    module load samtools/1.9
    $SAMTOOLS view -b -s ${downsample_fraction} ${input_bam} > ${output_bam};
    samtools index ${output_bam};
  fi
}

function 02_runBamQC {
  source ${SETTINGS}
  source ${SAMPLE_DOWNSAMPLE_INI}
  if [[ ! -f ${output_bam} ]]; then
    echo "No input bam provided to compute bamqc metrics"
    # exit -1
  else
    mkdir -p ${temp_dir}/BamQC;
    mkdir -p ${temp_dir}/picardTmp;
    $JAVA -Xmx16g -jar ${PICARD_TOOLS_ROOT}/SortSam.jar I=${output_bam} O=${sorted_bam} SO=coordinate CREATE_INDEX=true TMP_DIR=${temp_dir}/picardTmp;
    $SAMTOOLS view ${sorted_bam} | ${PERL} ${BAMQC} -s 1001 -i 1500 -q 30 -r ${TARG} -j ${temp_dir}/BamQC/${sample_name}.meta.json > ${temp_dir}/BamQC/${sample_name}.json;
    cp ${temp_dir}/BamQC/${sample_name}.json ${bamqc_json}
  fi
}


function 03_collectPicardMetrics {
  source ${SETTINGS}
  source ${SAMPLE_DOWNSAMPLE_INI}
  if [[ ! -f ${sorted_bam} ]]; then
    echo "No sorted bam provided to compute picard metrics"
    # exit -1
  else
    $JAVA -jar $PICARD CollectHsMetrics \
             I=${sorted_bam} \
             O=${output_picard_metrics} \
             R=$REF \
             BAIT_INTERVALS=${TARG}.interval_list \
             TARGET_INTERVALS=${TARG}.interval_list
  fi
}
