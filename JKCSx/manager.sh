#!/bin/bash
#sh manager.sh INPUT
start=`date +%s`

# locate TOOLS path
scriptpath="$( cd "$(dirname "$0")" ; pwd -P )"
toolspath="$scriptpath/../TOOLS"
source ~/.JKCSusersetup.txt

input=$1
if [ -z "$input" ];then echo Missing input;exit;fi
if [ ! -e "$input" ];then echo Input file does not exist;exit;fi
start_line=$2 
if [ -z "$start_line" ];then start_line=1;fi
line_num=$start_line

function submitME {
  #sbatch -p q24,q36,q40,q48 --time 30:00 -n 1 -J "manager" JKsend 
  submit_MANAGER sh $scriptpath/manager.sh $input $line_num "SUB"
}
if [ -z "$3" ]; then submitME; exit; fi

max_walltime=`echo 24*60*60 | bc` #in seconds
max_SingleCommandExecutionTime=`echo 60*60 | bc` #in seconds
walltime_cutoff=`echo $max_walltime-$max_SingleCommandExecutionTime | bc`

function check_if_all_finished {
  test=0
  while [ $test -eq 0 ]
  do
    ### Resubmit myself if I am running too long
    end=`date +%s`;runtime=$((end-start))
    if [ $runtime -gt $walltime_cutoff ]
    then
      submitME
      exit
    else
      sleep 120
    fi
    ###
    runningcpus=`squeue -u $USER --array -o "%C" | awk 'BEGIN{j=-1;c=0}{j+=1;if (j>0) c+=$1}END{print c}'`
    if [ $runningcpus -gt 3600 ]
    then
      sleep 120
      continue
    fi
    test=1
    ### TESTING VIA JKcheck
    JKcheck -num -this > .test 2>>output
    testlines=`wc -l .test | awk '{print $1}'`
    for l in `seq 1 $testlines`
    do
      cl=`head -n $l .test | tail -n 1`
      if [ "$cl" != "100.00" ]
      then
        test=0
      fi
    done
    ### TESTING WHETHER FOLDER IS IN QUEUE
    thisfolder=$(basename "$PWD")
    check2=$(squeue --format="%.100j %.1T" -u $USER | grep "$thisfolder ")
    if  [ ! -z "$check2" ]
    then
      test=0
    fi
    ###
  done
  return
}

lines=`wc -l $input | awk '{print $1}'`
original_directory=$PWD
for line_num in `seq $start_line $lines`
do
  if [ $line_num -ne $start_line ] && [ ! -z "$last_ID" ]
  then
    SBATCH_PREFIX+=" --dependency=afterok:$last_ID "
    #echo $SBATCH_PREFIX
    submit_MANAGER sh $scriptpath/manager.sh $input $line_num "SUB"
    exit
  fi
  cd $original_directory
  #check_if_all_finished
  line=`head -n $line_num $input | tail -n 1`
  echo $line > commandsNOW.txt
  echo $line >> commandsSENT.txt
  job=$(/bin/bash commandsNOW.txt)
  #if [ ! -z "$job" ] ;then echo $job;fi
  #echo $job | awk '{print $4}'  >> .jobs.txt
  echo -e "$job" >> output
  last_ID=`echo -e "$job" | grep "Submitted batch job" | awk '{print $4}' | xargs | sed "s/ /,/g"`
  rm commandsNOW.txt
done

