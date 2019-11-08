for x in `ls */*.txt`; do
  total_reads=`cat $x | grep total | cut -d+ -f1 | tr -d " "`;

done
