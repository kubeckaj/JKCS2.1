#!/bin/bash

#printing COMMAND also with comments
function echoCOMMAND() {
whitespace="[[:space:]]"
  for i in "$@"
  do
    if [[ $i =~ $whitespace ]]
    then
        i=\"$i\"
    fi
    printf "%s" "$i "
  done
  printf "\n"
}

#checking if path is full or absolute
function is_absolute() {
  case "$1" in
    ///* | //) echo 1;;
          //*) echo 0;; # on some systems, //foo is special and is
                       # not an absolute path. // alone is /
           /*) echo 1;;
            *) echo 0
  esac
}

#return lowercased letters
function lowercase() {
	echo "$1" | sed "y/ABCDEFGHIJKLMNOPQRSTUVWXYZ/abcdefghijklmnopqrstuvwxyz/"
}

function JKos {
  OS=`lowercase \`uname\``

  if [ "$OS" == "darwin" ]; then
    OS="osx"
  elif [ "$OS" == "windowsnt" ] ; then
    OS="vs"
  elif [ "${OS:0:5}" == "mingw" -o "$OS" == "msys_nt-6.3" ]; then
    OS="msys2"
  elif [ "$OS" == "linux" ]; then
    ARCH=`uname -m`
    if [ "$ARCH" == "i386" -o "$ARCH" == "i686" ] ; then
      OS="linux"
    elif [ "$ARCH" == "x86_64" ] ; then
      OS="linux64"
    elif [ "$ARCH" == "armv6l" ] ; then
      OS="linuxarmv6l"
    elif [ "$ARCH" == "armv7l" ] ; then
      OS="linuxarmv7l"
    else
    # We don't know this one, but we will try to make a reasonable guess.
      OS="linux"$ARCH
    fi
  fi
  
  echo "$OS"
}

#need Qprint; evaluate as JKecho 1 "string"
function JKecho {
  if [ $Qprint -ge $1 ]
  then
    if [ $Qsymbol -eq 1 ]
    then
      printf "😘  $scriptfile: ${@:2}\n"
    else
      printf    "$scriptfile: ${@:2}\n"
    fi
    OSTYPE=`JKos`
    if [ "$OSTYPE" == "osx" ]
    then
      printf "$scriptfile: ${@:2}\n" | sed -E "s/\x1B\[([0-9]{1,3}((;[0-9]{1,3})*)?)?[m|K]//g" >> output
    else
      printf "$scriptfile: ${@:2}\n" | sed -r "s/\x1B\[([0-9]{1,3}((;[0-9]{1,3})*)?)?[m|K]//g" >> output
    fi
  fi
}

#loading -help
function JKloadhelp {
  arguments_help=()
  # first search for -print
  ADDITION=$1
  for i in "${!arguments[@]}"
  do
    if [ "${arguments[i]}"  == "-help$ADDITION" ] || [ "${arguments[i]}"  == "--help$ADDITION" ];
    then
      help$ADDITION
    else
      arguments_help+=( "${arguments[i]}" )
    fi
  done
  arguments=("${arguments_help[@]}")
}

#load charge and multiplicity 
function JKloadCHARM {
  if [ -e $inputfile ]
  then
    CHARGE=`grep "TotalCharge" $inputfile | awk '{print $2}'`
    MULTIPLICITY=`grep "TotalMultiplicity" $inputfile | awk '{print $2}'`
  else
    JKecho 1 "Using charge 0 and multiplicity 1"
    CHARGE=0
    MULTIPLICITY=1
  fi
  JKecho 2 "Loading charge ($CHARGE) and multiplicity ($MULTIPLICITY)"
}

#load parameters for supercomputer
function JKloadsupercomputer {
  JKecho 2 "Loading default parameters for supercomputer."
  METHODsupercomputer=$program
  if [ -z "$METHODsupercomputer" ]
  then
    METHODsupercomputer="XTB"
  fi
  arguments_help=()
  # first search for -loc
  for i in "${!arguments[@]}"
  do
    if [ "${arguments[i]}"  == "-loc" ]
    then
      METHODsupercomputer="loc"
      arguments_help_pass+=( "${arguments[i]}" )
      continue
    fi
    arguments_help+=( "${arguments[i]}" )
  done
  arguments=("${arguments_help[@]}")
  JKecho 2 "Type of method used: $METHODsupercomputer"
  if [ -e $inputfile ]
  then
    # extracting info from userinputfile/inputfile
    lines=`grep -n "======================================================" $inputfile | sed "s/:/ /" | awk '{print $1}' | xargs`
    line1=`echo $lines | awk '{print $1}'`
    line2=`echo $lines | awk '{print $2}'`
    line2m1=`echo $line2-$line1-1 | bc`
    line2=`echo $line2-1 | bc`
    head -n $line2 $inputfile | tail -n $line2m1 > .${inputfile}_supercomputer
    # extracting SC=supercomputer variables
    methodtest=`grep -c "$METHODsupercomputer" .${inputfile}_supercomputer`
    if [ $methodtest -gt 0 ]
    then
      JKecho 2 "Extracting values from ${cfYELLOW}${inputfile}${cfDEF}"
      supercomputerline=`grep "$METHODsupercomputer" .${inputfile}_supercomputer`
      if [ -z "$supercomputerline" ]
      then
        JKecho 0 "Method $METHODsupercomputer does not exist in file ${cfYELLOW}${inputfile}${cfDEF}. ${cfRED}[EXITING]${cfDEF}"
        exit
      else
        JKecho 2 "Supercomputer parameters: PROGRAM SCtasks SCcpu SCnodes SCtime SCpar SCmem"
        JKecho 2 "Default supercomputer parameters: `echo $supercomputerline | column -t`"
      fi
      NoC=`grep "## Number of Combinations" ${inputfile} | awk '{print $6}'`
    else
      supercomputerline="$METHODsupercomputer 1 1 1 24:00:00 small 4000mb"
    fi
  else
    supercomputerline="$METHODsupercomputer 1 1 1 24:00:00 small 4000mb"
  fi
  ###
  SCtasks=`echo $supercomputerline | awk '{print $2}'`
  SCcpu=`echo $supercomputerline | awk '{print $3}'`
  SCnodes=`echo $supercomputerline | awk '{print $4}'`
  SCtime=`echo $supercomputerline | awk '{print $5}'`
  SCpar=`echo $supercomputerline | awk '{print $6}'`
  SCmem=`echo $supercomputerline | awk '{print $7}'`
  # checking script arguments for change
  JKecho 2 "Loading parameters given by user"
  arguments_help=()
  last=''
  for i in "${!arguments[@]}"
  do
    # -maxtasks
    if [ "${arguments[i]}"  == "-tasks" ] || [ "${arguments[i]}"  == "-maxtasks" ]; then
      last="-maxtasks"
      continue
    fi
    if [ "$last" == "-maxtasks" ]; then
      SCtasks=`echo "${arguments[i]}" | sed "s/M/$M/" | sed "s/NoC/$NoC/" | bc`
      last=''
      continue
    fi
    # -cpu
    if [ "${arguments[i]}"  == "-cpu" ] || [ "${arguments[i]}"  == "-CPU" ]; then
      last="-cpu"
      continue
    fi
    if [ "$last" == "-cpu" ]; then
      SCcpu=`echo "${arguments[i]}" | sed "s/M/$M/" | sed "s/NoC/$NoC/" | bc`
      last=''
      continue
    fi
    # -nodes
    if [ "${arguments[i]}"  == "-nodes" ] || [ "${arguments[i]}"  == "-N" ]; then
      last="-nodes"
      continue
    fi
    if [ "$last" == "-nodes" ]; then
      SCnodes=`echo "${arguments[i]}" | sed "s/M/$M/" | sed "s/NoC/$NoC/" | bc`
      last=''
      continue
    fi
    # -time
    if [ "${arguments[i]}"  == "-time" ] || [ "${arguments[i]}"  == "-time" ]; then
      last="-time"
      continue
    fi
    if [ "$last" == "-time" ]; then
      SCtime="${arguments[i]}"
      last=''
      continue
    fi
    # -partition
    if [ "${arguments[i]}"  == "-par" ] || [ "${arguments[i]}"  == "-partition" ]; then
      last="-par"
      continue
    fi
    if [ "$last" == "-par" ]; then
      SCpar="${arguments[i]}"
      last=''
      continue
    fi
    # -mem
    if [ "${arguments[i]}"  == "-mem" ] || [ "${arguments[i]}"  == "-memory" ]; then
      last="-mem"
      continue
    fi
    if [ "$last" == "-mem" ]; then
      SCmem="${arguments[i]}"
      if [[ ${SCmem,,} == *"mb"* ]] || [[ ${SCmem,,} == *"gb"* ]] || [[ ${SCmem,,} == *"mw"* ]] || [[ ${SCmem,,} == *"w"* ]] || [[ ${SCmem,,} == *"MB"* ]] || [[ ${SCmem,,} == *"GB"* ]] || [[ ${SCmem,,} == *"MW"* ]] || [[ ${SCmem,,} == *"W"* ]]
      then
        if [[ ${SCmem,,} == *"*"* ]]
        then
          JKecho 0 "Cannot procces memory requirement: ${cfYELLOW}$SCmem${cfDEF}. ${cfRED}[EXITING]${cfDEF}"
          JKecho 1 "It is easy to programmier, I was just lazy. JK"
          exit
        fi
      else
        SCmem=`echo $SCmem | sed "s/M/$M/" | sed "s/NoC/$NoC/" | bc`       
      fi
      last=''
      continue
    fi
    arguments_help+=( "${arguments[i]}" )
  done
  arguments=("${arguments_help[@]}")
  JKecho 2 "Final supercomputer parameters: `echo $METHODsupercomputer $SCtasks $SCcpu $SCnodes $SCtime $SCpar $SCmem | column -t`"
  if [ "$METHODsupercomputer" == "loc" ]
  then
    SC_command="sh $toolspath/SCRIPTS/JKsend "
  else
    SC_command=`program_SBATCH "$currentdir" $SCpar $SCtime $SCnodes $SCmem $SCcpu`" $toolspath/SCRIPTS/JKsend "
    #SC_command="sbatch -J "$currentdir" -p $SCpar --time $SCtime -N $SCnodes --mem-per-cpu $SCmem -n $SCcpu $SBATCHuseradd $toolspath/SCRIPTS/JKsend "
  fi
}

#loading directories from argument and entering them
function JKloaddirs {
  arguments_help=()
  # first search for -print
  for i in "${!arguments[@]}"
  do
    L4=`echo "${arguments[i]}" | cut -c1-4`
    if [ "$L4"  == "${folderbasename}_" ] 
    then
      motherdir=$PWD
      uploadfolder=`Cbrackets ${arguments[i]}`
      uploadfolder=`Cdash $uploadfolder`
      for j in $uploadfolder
      do
        if [ -d "$j" ]
        then
          folders+=" $j"
        else
          JKecho 0 "Folder ${cfGREEN}${j}${cfDEF} does not exist. ${cfRED}[EXITING]${cfDEF}"
          exit
        fi
      done
    else
      arguments_help+=( "${arguments[i]}" )
    fi
  done
  arguments=("${arguments_help[@]}")
  # checking variable folders
  if [ -z "$folders" ]
  then
    if [ ! -e $inputfile ]
    then
      motherdir=$PWD
      folders=`ls -d ${folderbasename}_* 2>/dev/null`
      folders=`echo $folders | xargs`
      if [ -z "$folders" ]
      then
        JKecho 0 "Neither ${cfYELLOW}$inputfile${cfDEF} nor ${cfBLUE}${folderbasename}${cfDEF}_ exist! ${cfRED}[NO - EXITING]${cfDEF}"
        JKecho 2 "Continuing anyway."
        #exit
      else
        JKecho 0 "  All subfolders = ${cfBLUE}$folders${cfDEF}"
        JKecho 1 "I, this script, will enter to all subfolders. :-D"
      fi
    fi
  fi
  # folders
  for i in $folders
  do
    JKecho 1 "Entering directory ${cfGREEN}${i}${cfDEF}"
    cd $motherdir/$i
    commandarguments1=`printf '"%s" ' "${arguments[@]}"` 
    commandarguments2=`printf '"%s" ' "${arguments_help_pass[@]}"`
    #echo $commandarguments1
    #echo $commandarguments2
    #echo $command
    command=`echo $scriptpath/$scriptfilecommand $commandarguments1 $commandarguments2`
    JKecho 2 "Evaluating command ${cfYELLOW}$command${cfDEF}"
    #echo "LINK 1 $motherdir/$i" >> ../commands_TODO.txt
    #echo "LINK 1 $motherdir" >> commands_TODO.txt
    eval $command
    cd ..
    JKecho 1 "Leaving directory ${cfRED}${i}${cfDEF}"
  done
  if [ ! -z "$folders" ]
  then
    JKecho 1 "I, this script, did all my job :-D"
    exit 
  fi
}

#loading program from arguments, checking if the function exist
function JKloadprogram {
  program=XTB
  arguments_help=()
  # first search for -print
  saveQprogram="no"
  for i in "${!arguments[@]}"
  do
    # this is just when -print is argument
    if [ "$saveQprogram" == "yes" ]
    then
      saveQprogram="no"
      arguments_help_pass+=( "${arguments[i]}" )
      program=${arguments[i]};
      continue
    fi
    if [ "${arguments[i]}"  == "-p" ] || [ "${arguments[i]}"  == "-program" ];
    then
      arguments_help_pass+=( "${arguments[i]}" )
      saveQprogram="yes"
    else
      arguments_help+=( "${arguments[i]}" )
    fi
  done
  #now check if function program_${program} exist?
  JKecho 2 "Program set: $program"
  if [ `type -t program_$program`"" == 'function' ]; 
  then 
    JKecho 2 "Function ${cfGREEN}program_$program${cfDEF} from ${cfYELLOW}~/.JKCSusersetup.txt${cfDEF} is utilized."
  else
    JKecho 0 "Sorry, I did not find function ${cfGREEN}program_$program${cfDEF} in ${cfYELLOW}~/.JKCSusersetup.txt${cfDEF}. [${cfRED}EXITING${cfDEF}]"
    exit 
  fi
  arguments=("${arguments_help[@]}")
}


#loading Qprint from arguments
function JKloadprint {
  Qprint=1
  arguments_help=()
  # first search for -print
  saveQprint="no"
  for i in "${!arguments[@]}"
  do
    # this is just when -print is argument
    if [ "$saveQprint" == "yes" ]
    then
      saveQprint="no"
      arguments_help_pass+=( "${arguments[i]}" )
      Qprint=${arguments[i]};
      continue
    fi
    if [ "${arguments[i]}"  == "-print" ];
    then
      arguments_help_pass+=( "${arguments[i]}" )
      saveQprint="yes"
    elif [ "${arguments[i]}"  == "-nosymbol" ]
    then
      arguments_help_pass+=( "${arguments[i]}" )
      Qsymbol=0
    else
      arguments_help+=( "${arguments[i]}" )
    fi
  done
  arguments=("${arguments_help[@]}")
  JKecho 2 "Output printing set to level $Qprint"
}

#loading -nocolours from argument
function JKloadcolours {
  arguments_help=()
  # first search for -nocolours
  for i in "${!arguments[@]}"
  do
    if [ "${arguments[i]}"  == "-nocolours" ]; 
    then 
      Qcolours=0; 
      arguments_help_pass+=( "${arguments[i]}" )
    else
      arguments_help+=( "${arguments[i]}" )
    fi
  done
  arguments=("${arguments_help[@]}")
  # load colours
  source $toolspath/LOADING/colours.txt
  # re-set scriptfile
  scriptfilecommand=$scriptfile
  scriptfile="${cfBLUE}J${cfRED}K${cfGREEN}C${cfYELLOW}S${cfDEF}"`echo $scriptfile | cut -c5-`
  JKecho 2 "Colours loaded."
}

#CALCULATE NUMBER OF ELEMENTS
function Felements {
  N=0
  for i in $*
  do
     N=`echo $N+1 |bc`
  done
  echo $N
}

# JKCS1 - replace brackets by number in composition 
function Cbrackets {
  outputTOT=""
  for input in $*
  do
    if [[ "$input" == *"("*")"* ]];
    then
      n1=`echo "$input" | sed -n "s/[\(].*//p" | wc -c | xargs`
      n2=`echo "$input" | sed -n "s/[\)].*//p" | wc -c | xargs`
      #n1=$(echo $(expr index "$input" "\("))
      #n2=$(echo $(expr index "$input" "\)"))
      n1p=`echo $n1+1 |bc`
      n2p=`echo $n2+1 |bc`
      diff=`echo $n2-$n1p | bc`
      if [ $n1 -eq 0 ]
      then
        before=''
      else
        n1m=`echo $n1-1 |bc`
        before=`echo ${input:0:$n1m}`
      fi
      in=`echo ${input:$n1:$diff} | sed 's/,/ /g'`
      after=`echo ${input:$n2}`
      output=''
      for i in $in
      do
        new=$before$i$after
        output+=" $new"
      done
      output2=''
      for j in $output
      do
        output20=`Cbrackets "$j"`
        output2+=" $output20"
      done
      outputTOT+=" $output2"
    else
      outputTOT+=" $input"
    fi
  done
  echo $outputTOT
}

# JKCS1 - replace dash by serie
function Cdash {
  outputTOT=""
  for input in $*
  do
    if [[ "$input" == *"-"* ]];
    then
      c="$input"
      n1=`echo "$input" | sed -n "s/[-].*//p" | wc -c | xargs`
      #n1=$(echo $(expr index "$input" "-"))
      n1m=$(echo $n1-1|bc)
      t1=`echo ${c:0:$n1m}`
      if [[ "$t1" == *"_"* ]];
      then
        F1=`echo ${t1%_*}`;
        F1="${F1}_"
        l=${#F1};
        F2=`echo ${t1:$l} `;
      else
        F1=""
        F2="$t1"
      fi
      t2=`echo ${c#*-}`;
      if [[ "$t2" == *"_"* ]];
      then
        F4=`echo ${t2#*_}`;
        n1=`echo "$t2" | sed -n "s/[_].*//p" | wc -c | xargs`
        #n1=$(echo $(expr index "$t2" "_"))
        n1m=`echo $n1-1 | bc`
        F3=`echo ${t2:0:$n1m}`;
        F4="_${F4}"
      else
        F3="$t2"
        F4=""
      fi
      output=''
      F3str=`grep -Eo '[[:alpha:]]+' <<< $F3`
      F3=`grep -Eo '[0-9]+' <<< $F3`
      for j in `seq $F2 $F3`
      do
        output+=" $F1$j$F4$F3str"
      done
      output2=''
      for i in $output
      do
        output2h=`Cdash $i`
        output2+=" $output2h"
      done
      outputTOT+=" $output2"
    else
      outputTOT+=" $input"
    fi
  done
  echo $outputTOT
}

function strindex { 
    local str=$1
    local search=$2
    let local n=1
    local retval=1 # here, 1 is failure, 0 success
    Qfound=0
    for col in $str; do # $str unquoted -> whitespace tokenization!
      if [ $col = $search ]; then
        echo $n
        Qfound=1
        break
      else
        ((n++))
      fi
    done
    if [ $Qfound -eq 0 ]; then echo -1; fi
}
