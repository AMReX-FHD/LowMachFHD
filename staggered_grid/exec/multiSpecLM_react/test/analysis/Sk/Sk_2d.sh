#!/bin/bash

# this script gathers S_k data for all pairs into res.Sk
# pair=1 (which is usually for rhotot) is excluded.
# (but you can easily change the script)
# S(0) values are excluded.
# 
# first col = kx
# i-th col = S(kx) for pair=i

# get filePrefix
if [ $# -eq 1 ]
then
  filePrefix=$1
else
  echo "ERROR: (usage) $0 filePrefix"
  exit
fi

# get info
firstfile=${filePrefix}.S_k.pair=1.Re.dat
if [ -f $firstfile ]
then
  # number of lines including the first two comment lines
  nline=$(wc -l < $firstfile)

  # number of files
  nfile=$(ls ${filePrefix}.S_k.pair=*.Re.dat | wc -w)

  # which column to get (corresponding to ky=0)
  col=$(((nline-2)/2+1))

  # which line to skip (corresponding to kx=0)
  skipline=$(((nline-2)/2+2))
else
  echo "ERROR: $firstfile does not exist"
  exit
fi

# get data
for ((i=1;i<=nfile;i++))
do
  infile=${filePrefix}.S_k.pair=${i}.Re.dat
  outfile=tmp$i

  if [ $i -eq 1 ]
  then
    awk -v var1=$skipline -v var2=1 '(NR>2 && NR!=var1){printf "%e\n", $var2}' $infile > $outfile
    list="$outfile" 
  else
    sed -i -e 's/ 0.00000000/#0.00000000/g' $infile
    awk -v var1=$skipline -v var2=$col '(NR>2 && NR!=var1){printf "%e\n", $var2}' $infile > $outfile 
    list="$list $outfile"
  fi 
done

output=res.Sk
paste $list > $output 
rm $list
echo "$output generated"
