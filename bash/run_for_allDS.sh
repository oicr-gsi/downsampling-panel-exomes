#!/bin/bash
project=TGL40
downsample_txt=/scratch2/groups/tgl/GBS-1572/Downsampling/bash/downsampling_fractions.txt
SCRP=/scratch2/groups/tgl/GBS-1572/Downsampling/bash
for downsample_fraction in `cat ${downsample_txt}`; do
  $SCRP/downsampling_per_project.sh $project ${downsample_fraction}
done
