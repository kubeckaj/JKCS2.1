#!/bin/bash
CURRDIR=$PWD

cat TOOLS/label_small

echo ___________________
echo TESTING JKCS SETUP:

### ~/.JKCSusersetup.txt
if [ -e .log ]
then
  rm .log
fi

printf "== testing existence of ~/.JKCSusersetup.txt:"
printf "== testing existence of ~/.JKCSusersetup.txt:\n" >> .log
if [ ! -e ~/.JKCSusersetup.txt ]
then
  printf "\nDid you run already setup.sh? [The ~/..JKCSusersetup.txt is missing] [EXITING]\n"
  printf "\nDid you run already setup.sh? [The ~/..JKCSusersetup.txt is missing] [EXITING]\n" >> .log
  exit
else
  printf " SUCCESFULL\n"
  printf " SUCCESFULL\n" >> .log
fi
source ~/.JKCSusersetup.txt
source TOOLS/LOADING/colours.txt

### WRKDIR
printf "== testing existence of WRKDIR:"
printf "== testing existence of WRKDIR:\n" >> .log
if [ ! -d "$WRKDIR" ]
then 
  printf "Working directory does not exist. ($WRKDIR) <== THIS MIGHT BE PROBLEM"
else
  printf " ${cfGREEN}SUCCESFULL${cfDEF}\n"
  printf "   testing permission to write:"
  printf "   testing permission to write:\n" >> .log
  touch $WRKDIR/.test 2>/dev/null
  if [ ! -e "$WRKDIR/.test" ]
  then
    printf "${cfRED}UNSUCCESFULL${cfDEF} //There is no permission for creating file in the WRKDIR <== THIS MIGHT BE SOMETIMES OK\n"
    printf "UNSUCCESFULL //There is no permission for creating file in the WRKDIR <== THIS MIGHT BE SOMETIMES OK\n" >> .log
  else
    printf " ${cfGREEN}SUCCESFULL${cfDEF}\n"
    printf " SUCCESFULL\n" >> .log
    rm $WRKDIR/.test
  fi
fi
##########

### PYTHON
if [ ! -d TEST ]; then mkdir TEST; fi
cd TEST
testdir=$PWD
printf "== testing Python:"
printf "== testing Python:\n" >> ../.log

if [ -e .help ]; then rm .help; fi
program_PYTHON --version > .help 2>&1
result=`cat .help | tail -n 1 | cut -c-9`
if [ "$result" == "Python 3." ]
then
  printf " ${cfGREEN}SUCCESFULL${cfDEF}\n"
  printf " SUCCESFULL\n" >> ../.log
 
  #Testing numpy
  echo "import numpy" > .test.py
  echo "a=numpy.matrix('1 2;3 4')" >> .test.py
  echo "print(a*a)" >> .test.py
  program_PYTHON .test.py > .test.out 2>&1
  result=`grep -c "[15 22]]" .test.out`
  if [ $result -eq 1 ]
  then
    echo -e "   -- numpy: ${cfGREEN}SUCCESFULL${cfDEF}"
    echo "   -- numpy: SUCCESFULL" >> ../.log
  else
    echo -e "   -- numpy: ${cfRED}UNSUCCESFULL${cfGREEN}"
    echo "   -- numpy: UNSUCCESFULL" >> ../.log
    cat .test.out >> ../.log
    echo "   :: see .log for the error"
    echo "   :: your python version probably does not have numpy libraries"
  fi
  rm .test.py .test.out

  #Testing numpy
  echo "import pandas" > .test.py
  program_PYTHON .test.py > .test.out 2>&1
  result=`wc -l .test.out | awk '{print $1}' `
  if [ $result -eq 0 ]
  then
    echo -e "   -- pandas: ${cfGREEN}SUCCESFULL${cfDEF}"
    echo "   -- pandas: SUCCESFULL" >> ../.log 
  else
    echo -e "   -- pandas: ${cfRED}UNSUCCESFULL${cfDEF}" 
    echo "   -- pandas: UNSUCCESFULL" >> ../.log
    cat .test.out >> ../.log
    echo "   :: see .log for the error"
    echo "   :: your python version probably does not have numpy libraries"
  fi
  rm .test.py .test.out

  #Testing ase
  echo "import ase" > .test.py
  program_PYTHON .test.py > .test.out 2>&1
  result=`wc -l .test.out | awk '{print $1}' `
  if [ $result -eq 0 ]
  then
    echo -e "   -- ase: ${cfGREEN}SUCCESFULL${cfDEF}"
    echo "   -- ase: SUCCESFULL" >> ../.log
  else
    echo -e "   -- ase: ${cfRED}UNSUCCESFULL${cfDEF}" 
    echo "   -- ase: UNSUCCESFULL" >> ../.log
    cat .test.out >> ../.log
    echo "   :: see .log for the error"
    echo "   :: your python version probably does not have numpy libraries"
  fi
  rm .test.py .test.out
else
  printf " ${cfRED}UNSUCCESFULL${cfDEF}\n" 
  cat .help >> ../.log
  echo "####################################" >> ../.log
  echo "   :: see .log for the error"
  echo "   :: open ~/.JKCSusersetup -> check program_PYTHON -> check the function/setup paths"
fi
rm .help
##############

### JKQCpickle
cd $testdir
printf "== testing JKQCpickle:"
printf "== testing JKQCpickle:\n" >> ../.log
program_JKQCpickle > .test.out 2>&1
cd $testdir
result=`grep -c "No inputs. No" .test.out`
if [ $result -eq 1 ]
then
  printf " ${cfGREEN}SUCCESFULL${cfDEF}\n"
  printf " SUCCESFULL\n" >> ../.log
else
  printf " ${cfRED}UNSUCCESFULL${cfDEF}\n"
  printf " UNSUCCESFULL\n" >> ../.log
  cat .test.out >> ../.log
  echo "####################################" >> ../.log
  echo "   :: see .log for the error"
  echo "   :: open ~/.JKCSusersetup -> check program_JKQCpickle"
fi
rm .test.out
##########

### ABCluster
cd $testdir
printf "== testing ABCluster:"
printf "== testing ABCluster:\n" >> ../.log
touch .calc.inp
program_ABC .calc.inp 2> .calc.out
cd $testdir
result=`grep -c "Cannot read the cluster file name." .calc.out`
if [ $result -eq 1 ]
then 
  printf " ${cfGREEN}SUCCESFULL${cfDEF}\n"
  printf " SUCCESFULL\n" >> ../.log
else
  printf " ${cfRED}UNSUCCESFULL${cfDEF}\n"
  printf " UNSUCCESFULL\n" >> ../.log
  cat .calc.out >> ../.log
  echo "####################################" >> ../.log
  echo "   :: see .log for the error"
  echo "   :: open ~/.JKCSusersetup -> check program_ABC or setup properly path PATH_ABCluster"
fi
rm .calc.inp .calc.out
##########

### XTB
cd $testdir
printf "== testing XTB:"
printf "== testing XTB:\n" >> ../.log
if [ -e .test.xyz ]
then 
  rm .test.xyz
fi
touch .test.xyz
program_XTB .test.xyz > .test.log 2> .test.log
cd $testdir
if [ ! -e .test.log ]; then touch .test.log; fi
result=`grep -c "abnormal termination of xtb" .test.log`
if [ $result -eq 1 ]
then
  printf " ${cfGREEN}SUCCESFULL${cfDEF}\n"
  printf " SUCCESFULL\n" >> ../.log
else
  printf " ${cfRED}UNSUCCESFULL${cfDEF}\n" 
  printf " UNSUCCESFULL\n" >> ../.log
  cat .test.log >> ../.log
  echo "####################################" >> ../.log
  echo "   :: see .log for the error"
  echo "   :: open ~/.JKCSusersetup -> check program_XTB or setup properly path PATH_XTB"
fi
rm .test.xyz .test.log
##########

### GAUSSIAN
cd $testdir
printf "== testing Gaussian(G16):"
printf "== testing Gaussian(G16):\n" >> ../.log
touch .test.com
program_G16 .test.com > .test.log 2> .test.log
cd $testdir
if [ ! -e .test.log ]; then touch .test.log; fi
result=`grep -c "Route card not found." .test.log`
if [ $result -eq 1 ]
then
  printf " ${cfGREEN}SUCCESFULL${cfDEF}\n"
  printf " SUCCESFULL\n" >> ../.log
else
  printf " ${cfRED}UNSUCCESFULL${cfDEF}\n"
  printf " UNSUCCESFULL\n" >> ../.log
  cat .test.log >> ../.log
  echo "####################################" >> ../.log
  echo "   :: see .log for the error"
  echo "   :: open ~/.JKCSusersetup -> check program_G16 or setup properly path PATH_G16"
fi
rm .test.log .test.com
##########

### ORCA
cd $testdir
printf "== testing ORCA:"
printf "== testing ORCA:\n" >> ../.log
touch .test.inp 
program_ORCA .test.inp > .test.out  2> .test.out
cd $testdir
if [ ! -e .test.out ]; then touch .test.out; fi
result=`grep -c "You must have a \[COORDS\] ... \[END\] block in your input" .test.out`
if [ $result -ne 0 ]
then
  printf " ${cfGREEN}SUCCESFULL${cfDEF}\n"
  printf " SUCCESFULL\n" >> ../.log
else
  printf " ${cfRED}UNSUCCESFULL${cfDEF}\n"
  printf " UNSUCCESFULL\n" >> ../.log
  cat .test.out >> ../.log
  echo "####################################" >> ../.log
  echo "   :: see .log for the error"
  echo "   :: open ~/.JKCSusersetup -> check program_ORCA or setup properly path PATH_ORCA"
fi
rm .test.inp .test.inp.out .test.out .test.xyz  2>/dev/null

cd $testdir
cd ..
if [ -e TEST ]; then rm -r TEST; fi
echo ___________________
