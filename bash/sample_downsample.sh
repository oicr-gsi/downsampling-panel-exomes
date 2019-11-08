!#/bin/bash
BASEDIR=/Users/prath/Documents/TGL/TGL_support_project/WDL/Downsampling/bash/
INI=${BASEDIR}/example.ini
SETTINGS=${BASEDIR}/downsample_settings.ini
source ${BASEDIR}/bash/design_functions.sh $SETTINGS $INI
# start stringing

00_checkINI()
01_downsample_bam()
02_runBamQC()
03_collectPicardMetrics()
