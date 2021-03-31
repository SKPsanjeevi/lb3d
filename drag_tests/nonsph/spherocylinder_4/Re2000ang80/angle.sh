#!/bin/bash
# deg=45
# while read f1 f8 f9 rest; do
deg=`echo $1`
INPUT=input-default
  nx=`grep "^nx" $INPUT | awk '{print $3}' | sed "s/'//g"`
  ny=`grep "^ny" $INPUT | awk '{print $3}' | sed "s/'//g"`
  nz=`grep "^nz" $INPUT | awk '{print $3}' | sed "s/'//g"`
  f1=$((nx/2))
  f2=$((ny/2))
#  f3=$((nz/4))
  f3=$((nz/3))
  f4="0.0 0.0 0.0"
  f7=0  
  f10="0 0 0"
  pi=`echo "4*a(1)" | bc -l`
  rad=`echo "$deg*($pi/180)" | bc -l` 
  f8=`echo "-1*s($rad)" | bc -l`
  f9=`echo "1*c($rad)" | bc -l`
  printf '%s\n' "$f1 $f2 $f3    $f4    $f7 $f8 $f9    $f10" > init.cfg
# done < test.cfg


# while read f1 f2 f3 f4 f5 f6 f7 f8 f9 rest; do
#  printf '%s\n' "$f1 $f2 $f3 $f4 $f5 $f6 $f7 $f8 $f9 $rest" &> test.cfg
