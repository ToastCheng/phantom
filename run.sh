#!/bin/bash
##################################
#              RUN!              #
##################################
#
#run the MonteCarlo and parameter calculator automatic.
dir=`pwd`
if [ $# = 1 ]; then
  echo "created by SHI-CHENG TU 2017.09.13"
  echo ""
  echo "This is a bash script which will call two program while working."
  echo ""
  echo "input the parameters for calculate using the form:"
  echo "[min size(um)] [max size(um)] [step] [concentration(g/g)]"
  echo "e.g. bash run.sh 0.9 1.1 0.01 0.025"
  echo "means"
  echo "min size of sphere [um]: 0.9"
  echo "max size of sphere [um]: 1.1"
  echo "step between min & max:  0.01"
  echo "concentration [#/um^3]:     0.025"
  exit 0
fi
min=$1
max=$2
step=$3
con[0]=$4
con[1]=$5
con[2]=$6
con[3]=$7
con[4]=$8
echo $#
declare -i count=5
#while [ "$count" -le "$#+1" ]
while (("$count"<="$#+1"))
do
declare -i x=${count}-5
declare -i index=0
check=1
while [ $check = 1 ]
do
  size=$(bc <<< "$min+($index*$step)")
  check=$(bc <<< "$min+($index*$step)<$max")

  cd $dir/phase_function
  ./mie $size ${con[$x]}
  cd $dir/MCML_Phantom
  ./mc
  index=index+1
done
count=count+1
done
echo "================================================"
echo "Monte Carlo simulation complete."
#echo "press any key to contine for python plotting program."
#read -n 1 -s -r
cd $dir
#python plot.py $1 $2 $3
