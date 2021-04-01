#!/bin/bash
# rm -rf Re*
# declare -a array=(0 10 30 60 80 90)
declare -a array=(0 10 30 45 60 80 90)
# declare -a array=(30)
# declare -a array=(45)
# declare -a array=(0 10 20 30 40 50 60 70 80 90)
# declare -a array=(30)
# declare -a array=(0 90)
for deg in "${array[@]}"; do
  target=Re2000ang$deg
  mkdir -p $target
  cp in* $target/
  cp myjob.cmd $target/
  cp angle.sh $target/
  cd $target

  sh angle.sh $deg
  sbatch myjob.cmd
  cd ..
done
