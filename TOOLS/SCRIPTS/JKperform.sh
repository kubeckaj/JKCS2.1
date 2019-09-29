#!/bin/bash

# FUNCTIONS:
# waitcheck = .wait
# getcommand = DOING/DONE/WAITING
# fincommand = DOING -> DONE
# timming

starttime=`date +%s`

# WAITING QUESTION. MAYBE SOME OTHER GUY IS WORKING WITH FILE
#####################
function waitcheck {
  if [ $1 -eq 0 ] 
  then
    echo 0 > .wait
    return
  fi
  local waittest=0
  while [ $waittest -eq 0 ]
  do
    if [ -e .wait ]
    then
      local value=`cat .wait`
      if [ "$value" == "0" ]
      then
        echo "2 $MY_output" > .wait
        sleep 1
        local value=`cat .wait`
        if [ "$value" == "2 $MY_output" ]
        then
          echo 1 > .wait
          waittest=1
        fi
      else
        sleep 3
      fi
    else
      echo 1 > .wait
      waittest=1
    fi
  done
}
######################
function getcommand {
  waitcheck 1
  # SOLVING all commands in commands_TODO.txt
  if [ -e commands_TODO.txt ]
  then
    local jobscount=`wc -l commands_TODO.txt | awk '{print $1}'`
    #echo jobscount=$jobscount
    if [ $jobscount -eq 0 ]
    then
      echo "No jobs in commands_TODO.txt. [EXITING]" >> $MY_output
    fi
    local doingcount=0
    local linenumber=0
    local test=0
    local linktest=0 #testing how many calc. are already linked 
    for i in `seq 1 $jobscount`
    do
      local linenumber=`echo $linenumber+1 | bc`
      local line=`head -n $i commands_TODO.txt | tail -n 1`
      local firstword=`echo $line | awk '{print $1}'`
      if [ "$firstword" == "DOING" ]
      then
        doingcount=`echo $doingcount+1 | bc`
        continue
      elif [ "$firstword" == "DONE" ]
      then
        continue
      elif [ "$firstword" == "WAITING" ]
      then
        if [ $doingcount -eq 0 ]
        then
          test=1
          line=${line:8}
          sed "${linenumber}s/WAITING//" commands_TODO.txt > .commands_TODO.txt_help; mv .commands_TODO.txt_help commands_TODO.txt
          sed "${linenumber}s/^/DOING $MY_ID /" commands_TODO.txt > .commands_TODO.txt_help; mv .commands_TODO.txt_help commands_TODO.txt
          break
        else
          waitcheck 0
          exit
        fi
      elif [ "$firstword" == "LINK" ]
      then
        test=1
        linkentered=`echo $line | awk '{print $2}'`
        if [ "$linkentered" == "0" ]
        then
          line=${line:7}
          sed "${linenumber}s/LINK 0//" commands_TODO.txt > .commands_TODO.txt_help; mv .commands_TODO.txt_help commands_TODO.txt
          sed "${linenumber}s/^/LINK $MY_ID /" commands_TODO.txt > .commands_TODO.txt_help; mv .commands_TODO.txt_help commands_TODO.txt
          break
        elif [ "$linkentered" == "$MY_ID" ]
        then
          sed "${linenumber}s/LINK/DONE/" commands_TODO.txt > .commands_TODO.txt_help; mv .commands_TODO.txt_help commands_TODO.txt           
          waitcheck 0
          exit
        else
          doingcount=`echo $doingcount+1 | bc`
          continue
        fi 
      else
        test=1
        sed "${linenumber}s/^/DOING $MY_ID /" commands_TODO.txt > .commands_TODO.txt_help; mv .commands_TODO.txt_help commands_TODO.txt
        break
      fi
    done
    if [ $test -eq 0 ]
    then
      waitcheck 0
      echo "Everything is finished or already calculating. [EXITING]" >> $MY_output
      exit
    fi
#  # I SHOULD NOT CARE ABOUT commands.txt
#  elif [ -e commands.txt ]
#  then
#    wc -l commands.txt 
  else
    echo "No commands in $PWD. [EXITING]" >> $MY_output
    waitcheck 0
    exit
  fi
  waitcheck 0
  MY_LINE=$line
  MY_LINENUMBER=$linenumber 
}
####################
function fincommand { 
  waitcheck 1
  sed "${MY_LINENUMBER}s/DOING //" commands_TODO.txt > .commands_TODO.txt_help; mv .commands_TODO.txt_help commands_TODO.txt
  sed "${MY_LINENUMBER}s/LINK //" commands_TODO.txt > .commands_TODO.txt_help; mv .commands_TODO.txt_help commands_TODO.txt
  sed "${MY_LINENUMBER}s/^/DONE /" commands_TODO.txt > .commands_TODO.txt_help; mv .commands_TODO.txt_help commands_TODO.txt
  waitcheck 0 
}
####################
function timing {  
  # current runtime
  endtime=`date +%s`
  runtime=$((endtime-jobtime))
  totruntime=$((endtime-starttime))
  echo runtime = $runtime >> $MY_output
  echo totruntime = $totruntime >> $MY_output
  
  # time need for one job
  if [ ! -e .maxtime ]
  then
    echo $runtime > .maxtime
    maxtime=$runtime
  else
    maxtime=`cat .maxtime`
    maxtimetest=`echo "$maxtime < $runtime" | bc -l`
    if [ $maxtimetest -eq 1 ]
    then
      echo $runtime > .maxtime
      maxtime=$runtime
    fi
  fi

  # runtime / req.time
  percentage=`echo "100*($totruntime+$maxtime)/($MY_reqtime*3600)" | bc -l`
  percentagetest=`echo "$percentage < 90" | bc -l`
  #echo $percentage $percentagetest
  if [ $percentagetest -eq 0 ]
  then
    echo "UPGRADE THIS PLS. OVER 90% OF REQUESTED TIME USED. [EXITING]"
    exit
  fi 
}
#############################################################
#############################################################
# WHO AM I
MY_ID=$SLURM_JOBID
if [ -z "$MY_ID" ]; then MY_ID="LOC"; echo 1 local task is running; fi
MY_output=output.$MY_ID
MY_motherdir=$PWD
#TODO link hours to requesting time
MY_reqtime=24 #in hours

echo JOBID = $MY_ID >> $MY_output
echo mother dir = $MY_motherdir >> $MY_output
# LOAD USER SETUP
source ~/.JKCSusersetup.txt

MY_jobnumber=0
while [ 1 -eq 1 ]
do
  echo "################" >> $MY_output 
  jobtime=`date +%s`
  MY_jobnumber=`echo $MY_jobnumber+1 | bc`
  echo "Job number: $MY_jobnumber" >> $MY_output
  # LOADING COMMAND
  getcommand #gives: MY_LINE, MY_LINENUMBER
  echo "Line to be done: $MY_LINE" >> $MY_output
  echo "Linenumber to be done: $MY_LINENUMBER" >> $MY_output
  # PERFORMING THE LINE
  eval $MY_LINE
  # RETURN TO MOTHER DIR
  cd $MY_motherdir
  # UPDATING COMMANDS
  fincommand 
  # CHECK IF I AM NOT OUT OF TIME
  timing
  echo "################" >> $MY_output 
done


