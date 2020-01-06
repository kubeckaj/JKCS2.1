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
    local test=0             #did I found some job todo
    local Qwaiting=0         #did I pass some waiting job?
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
        Qwaiting=1
        if [ $doingcount -eq 0 ]
        then
          test=1
          line=${line:8}
          sed "${linenumber}s/WAITING//" commands_TODO.txt > .commands_TODO.txt_help; mv .commands_TODO.txt_help commands_TODO.txt
          sed "${linenumber}s/^/DOING $MY_ID /" commands_TODO.txt > .commands_TODO.txt_help; mv .commands_TODO.txt_help commands_TODO.txt
          ### OPEN LINKS ###
          for j in `seq $i $jobscount`    
          do
            local WAITINGline=`head -n $j commands_TODO.txt | tail -n 1`
            local WAITINGfirstword=`echo $WAITINGline | awk '{print $1}'`
            if [ $WAITINGfirstword == "LINK" ]
            then
              WAITINGnewcurrentdir=`echo $WAITINGline | awk '{print $3}'`
              cd $WAITINGnewcurrentdir
              WAITINGlinkteststatus=`grep -c "LINK 0 $MY_motherdir" commands_TODO.txt`
              if [ $WAITINGlinkteststatus -eq 1 ]
              then
                waitcheck 1
                linkline=`grep -n "LINK 0 $MY_motherdir" commands_TODO.txt | sed 's/:/ /' | awk '{print $1}'`
                sed "${linkline}s/LINK 0/LINK 1/" commands_TODO.txt > .commands_TODO.txt_help; mv .commands_TODO.txt_help commands_TODO.txt
                waitcheck 0
              fi
              cd $MY_motherdir
            fi
          done
          ###     
          break
        else
          continue
        fi
      elif [ "$firstword" == "LINK" ]
      then
        # LINK STATUS DIR
        #STATUS: 0 closed, 1 open
        linkstatus=`echo $line | awk '{print $2}'`
        newcurrentdir=`echo $line | awk '{print $3}'`
        if [ "$linkstatus" == "0" ]
        then
          continue
        else
          test=1
          if [ $doingcount -eq 0 ]
          then
            sed "${linenumber}s/LINK 1/DONE/" commands_TODO.txt > .commands_TODO.txt_help; mv .commands_TODO.txt_help commands_TODO.txt
            cd $newcurrentdir
            linkteststatus=`grep -c "LINK 1 $MY_motherdir" commands_TODO.txt`
            if [ $linkteststatus -eq 1 ]
            then
              waitcheck 1
              linkline=`grep -n "LINK 1 $MY_motherdir" commands_TODO.txt | sed 's/:/ /' | awk '{print $1}'`
              sed "${linkline}s/LINK 1/DONE/" commands_TODO.txt > .commands_TODO.txt_help; mv .commands_TODO.txt_help commands_TODO.txt
              waitcheck 0
            fi
            linkteststatus=`grep -c "LINK 0 $MY_motherdir" commands_TODO.txt`
            if [ $linkteststatus -eq 1 ]
            then
              waitcheck 1
              linkline=`grep -n "LINK 1 $MY_motherdir" commands_TODO.txt | sed 's/:/ /' | awk '{print $1}'`
              sed "${linkline}s/LINK 0/DONE/" commands_TODO.txt > .commands_TODO.txt_help; mv .commands_TODO.txt_help commands_TODO.txt
              waitcheck 0
            fi
            cd $MY_motherdir
            echo "MY OLD MY_motherdirNEW = $MY_motherdirNEW" >> $MY_output
            MY_motherdirNEW=$newcurrentdir
          else
            cd $newcurrentdir
            linkteststatus=`grep -c "LINK 1 $MY_motherdir" commands_TODO.txt`
            if [ $linkteststatus -eq 1 ]
            then
              waitcheck 1 
              linkline=`grep -n "LINK 1 $MY_motherdir" commands_TODO.txt | sed 's/:/ /' | awk '{print $1}'`
              sed "${linkline}s/LINK 1/LINK 0/" commands_TODO.txt > .commands_TODO.txt_help; mv .commands_TODO.txt_help commands_TODO.txt
              waitcheck 0
            fi
            cd $MY_motherdir
            MY_motherdirNEW=$newcurrentdir         
          fi
          line="echo Linking to $MY_motherdirNEW" >> $MY_output
          test=1
          break
        fi 
      elif [ $Qwaiting -eq 0 ]
      then
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
}

function timingTEST {
  endtime=`date +%s`
  totruntime=$((endtime-starttime))
  if [ -e .maxtime ]
  then 
    maxtime=`cat .maxtime`
  else
    maxtime=0
    echo "maxtime=0 (hopefully it is ok to assume then next job will not take so much time)" >> $MY_output
  fi

  # runtime / req.time
  percentage=`echo "100*($totruntime+$maxtime)/($MY_reqtime*3600)" | bc -l`
  percentagetest=`echo "$percentage < 90" | bc -l`
  #echo $percentage $percentagetest
  if [ $percentagetest -eq 0 ]
  then
    echo "SORRY I HAVE TO TURN OF :-( OVER 90% OF REQUESTED TIME USED. [EXITING]"
    exit
  fi 
}

function check_jobs_availability {
  if [ -e .crealTODO.txt ]
  then
    if [ ! -e commands_TODO.txt ]
    then
      cp .crealTODO.txt commands_TODO.txt
    fi
    N1=`wc -l .crealTODO.txt | awk '{print $1}'`
    #
    if [ -e .crealJOBS.txt ]
    then
      N2=`wc -l .crealJOBS.txt | awk '{print $1}'`
    else
      N2=0
    fi
    #
    if [ $N1 -gt $N2 ]
    then
      echo 1
    else
      echo 0
    fi
    # 
  else
    echo 0
  fi
}

function perform_job {
  #sign in
  echo $MY_ID >> .crealJOBS.txt
  jobline=`grep -n $MY_ID .crealJOBS.txt | tail -n 1 | sed 's/:/ /' | awk '{print $1}'`
  N1=`wc -l .crealTODO.txt | awk '{print $1}'`
  if [ $N1 -ge $jobline ]
  then
    command=`head -n $jobline .crealTODO.txt | tail -n 1`
    dir=$PWD
    eval $command
    cd $dir
    echo $jobline >> .crealDONE.txt
    change_commandsTODO $jobline
    N2=`wc -l .crealDONE.txt | awk '{print $1}'`
    if [ $N1 -eq $N2 ]
    then
      finish_all
    fi 
  fi
}

function Isleep {
  sleep $sleeptime
  if [ $sleeptime -lt 5 ]
  then
    sleeptime=`echo $sleeptime+1 |bc`
  fi
}

function perform_command {
  if [ -e ".waitcommand" ]
  then
    Isleep
  else
    testjobs=`check_finished_jobs`
    # Check if jobs are all done
    if [ $testjobs -eq 0 ] 
    then     
      Isleep 
    else
      echo $MY_ID >> .waitcommand
      sleep 1
      testw=`grep -n "$MY_ID" .waitcommand | head -n 1 | sed 's/:/ /' | awk '{print $1}'`
      if [ $testw -eq 1 ]
      then
        rm .creal*.txt .link.txt 2> /dev/null
        #sign in
        echo $MY_ID >> .cbigJOBS.txt
        jobline=`grep -n $MY_ID .cbigJOBS.txt | tail -n 1 | sed 's/:/ /' | awk '{print $1}'`
        N1=`wc -l .cbigTODO.txt | awk '{print $1}'`
        if [ $N1 -ge $jobline ]
        then
          dir=$PWD
          command=`head -n $jobline .cbigTODO.txt | tail -n 1`
          if [[ "$command" == *"JKCS"* ]]; 
          then 
            eval "$command -maxtasks 0"
          else
            eval "$command"
          fi
          cd $dir
          echo $jobline >> .cbigDONE.txt
          rm .waitcommand
          #change_commandsTODO $jobline
          #N2=`wc -l .cDONE.txt | awk '{print $1}'`
          #if [ $N1 -eq $N2 ]
          #then
          #  finish_all
          #fi
        fi
        sleeptime=1
      else
        Isleep
      fi
    fi
  fi
}

function change_commandsTODO {
  update=$1
  if [ -e commands_TODO.txt ] 
  then
    test=`wc -l commands_TODO.txt | awk '{print $1}'`
    t1n=`grep -c "DONE" commands_TODO.txt`
    t2n=`wc -l .crealDONE.txt | awk '{print $1}'`
    t1n=`echo $t1n+1| bc`
  else
    test=0
  fi
  if [ $test -gt 0 ] && [ $test -ge $update ] && [ $t1n -eq $t2n ]
  then
    finish_line $update commands_TODO.txt
  else
    finish_all
  fi 
}

function finish_all {
  comfile=`newfile .commands_TODO.txt_redo`
  cp .crealTODO.txt $comfile
  list=`cat .crealDONE.txt | xargs`
  for i in $list
  do
    finish_line $i $comfile
  done
  n1test=`echo $list | xargs -n 1 | wc -l`
  n2test=`wc -l .crealDONE.txt | awk '{print $1}'`
  if [ $n1test -eq $n2test ]
  then
    mv $comfile commands_TODO.txt
  else
    rm $comfile
  fi
}

function finish_line {
  update=$1
  file=$2
  mvfile=`newfile .commands_TODO.txt_help`
  sed "${update}s/^/DONE $MY_ID /" $file > $mvfile ;
  mv $mvfile $file 2> /dev/null
}

function check_finished_jobs {
  if [ -e .crealDONE.txt ]
  then
    N3=`wc -l .crealDONE.txt | awk '{print $1}'`
  else
    N3=0
  fi
  if [ -e .crealTODO.txt ]
  then
    N4=`wc -l .crealTODO.txt | awk '{print $1}'`
  else
    N4=0
  fi

  ##
  if [ $N3 -eq $N4 ]
  then
    rm .creal*.txt 2> /dev/null
    echo 1
  else
    echo 0
  fi
}

function check_finished_commands {
  if [ -e .cbigDONE.txt ]
  then
    N3=`wc -l .cbigDONE.txt | awk '{print $1}'`
  else
    N3=0
  fi
  if [ -e .cbigTODO.txt ]
  then
    N4=`wc -l .cbigTODO.txt | awk '{print $1}'`
  else
    N4=0
  fi

  ##
  if [ $N3 -eq $N4 ] 
  then 
    rm .cbig*.txt 2> /dev/null
    echo 1
  else
    echo 0
  fi
}

function check_link {
  if [ -e .link.txt ]
  then
    lines=`wc -l .link.txt | awk '{print $1}'`
    for i in `seq 1 $lines`
    do
      line=`head -n $i .link.txt | tail -n 1  2> /dev/null`
      opentest=`echo $line | awk '{print $1}' 2> /dev/null`
      if [ ! -z "$opentest" ]
      then
        if [ $opentest -eq 1 ]
        then
          newdir=`echo $line | awk '{print $2}'`
          if [ ! -d $newdir ]
          then
            sed "${i}s/1 ${newdir%/*}\//0 ${newdir%/*}\//" .link.txt > .linkbcp
            mv .linkbcp .link.txt
          else
            echo $newdir
            break
          fi
        fi
      fi
    done
    echo 0
  else
    echo 0
  fi 
}

function close_above_link {
  if [ -e ../.link.txt ]
  then 
    thisdir=${PWD##*/}
    count=`grep -c " ${thisdir}/" ../.link.txt`
    if [ $count -gt 0 ]
    then
      bcpfile=`newfile .linkbcp`
      mvfile=`newfile .linkmv`
      cp ../.link.txt $bcpfile
      testjobs=`check_finished_jobs`
      if [ $testjobs -eq 1 ]
      then
        sleep 1
        sed "s/1 ${thisdir}\//0 ${thisdir}\//g" ../.link.txt > ${mvfile}B
        sed "s/2 ${thisdir}\//0 ${thisdir}\//g" ${mvfile}B > $mvfile
        rm ${mvfile}B
        rm .creal* .link.txt 2>/dev/null 
      else
        sed "s/1 ${thisdir}\//2 ${thisdir}\//g" ../.link.txt > $mvfile
      fi
      mv $mvfile ../.link.txt 2> /dev/null
      check=`wc -l ../.link.txt | awk '{print $1}'`
      if [ $check -eq 0 ]
      then
        cp $bcpfile ../.link.txt
      fi
      rm $mvfile $bcpfile 2> /dev/null
    fi
  fi
}

function newfile {
  base=$1
  exist=1
  while [ $exist -eq 1 ]
  do
    file=$base$RANDOM
    if [ ! -e $file ]
    then
      exist=0
    fi
  done
  echo $file
} 


#############################################################
#############################################################
# FILES TO KNOW: commands_TODO.txt, commands.txt
# HIDDEN FILES: .crealTODO.txt, .crealJOBS.txt, .crealDONE.txt 
#               .cbigTODO.txt, .cbigJOBS.txt, cbigDONE.txt
#               .link.txt
#               .wait
# WHO AM I
sleeptime=1
MY_ID=$SLURM_JOBID
if [ -z "$MY_ID" ]; then MY_ID="LOC"; echo "1 local task (LOC) is running"; fi
#MY_output=.output.$MY_ID
MY_dir=$PWD
#MY_motherdirNEW=$PWD
#TODO link hours to requesting time
#MY_reqtime=24 #in hours

#echo JOBID = $MY_ID >> $MY_output
#echo mother dir = $MY_motherdir >> $MY_output
# LOAD USER SETUP
source ~/.JKCSusersetup.txt

#MY_jobnumber=0

##
#RUN
while [ 1 -eq 1 ]
do
  cd $MY_dir

  ### IS HERE A LINK? ###
  test=`check_link`
  if [ "$test" != "0" ]
  then
    cd $test
    MY_dir=$PWD
    sleeptime=1
    continue
  fi

  ### IS HERE A JOB? ###
  test=`check_jobs_availability`
  if [ $test -gt 0 ]
  then
    perform_job
    sleeptime=1
    continue
  fi

  ### IS HERE BIG COMMAND? ###  
  test=`check_finished_commands`
  if [ $test -eq 0 ] 
  then 
    perform_command
    continue
  fi
 
  #close above link
  `close_above_link`
 
  ### IS EVERYTHING DONE HERE? ###
  if [ ! -e input.txt ]
  then
    cd ..
    MY_dir=$PWD
    sleeptime=1
    continue
  fi

  ### I CAN GO      
  echo $MY_ID is done.
  exit
done


#COMMANDNUMBER="0"
#while [ 1 -eq 1 ]
#do
#  NT1=`grep -c "DONE" commands_TODO.txt`
#  NT2=`wc -l commands_TODO.txt | awk '{print $1}'`
#  if [ $NT1 -ne $NT2 ] 
#  then
#    echo "################" >> $MY_output 
#    jobtime=`date +%s`
#    MY_jobnumber=`echo $MY_jobnumber+1 | bc`
#    echo "Job number: $MY_jobnumber" >> $MY_output
#    timingTEST #check if I have still energy to calculate
#    # LOADING COMMAND
#    getcommand #gives: MY_LINE, MY_LINENUMBER
#    echo "Line to be done: $MY_LINE" >> $MY_output
#    echo "Linenumber to be done: $MY_LINENUMBER" >> $MY_output
#    # PERFORMING THE LINE
#    eval $MY_LINE
#    if [ "$MY_motherdir" == "$MY_motherdirNEW" ]
#    then
#      # RETURN TO MOTHER DIR
#      cd $MY_motherdir
#      # UPDATING COMMANDS
#      fincommand 
#      timing
#      echo "################" >> $MY_output 
#    else
#      cd $MY_motherdirNEW
#      MY_motherdir=$MY_motherdirNEW
#    fi
#  else
#    exit
#  fi
#done


