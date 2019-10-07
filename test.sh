#!/bin/bash
echo ___________________
echo TESTING JKCS SETUP:

if [ ! -e ~/.JKCSusersetup.txt ]
then
  echo "Did you run already setup.sh? [The ~/..JKCSusersetup.txt is missing]"
fi
if [ -e .log ]
then
  rm .log
fi
source ~/.JKCSusersetup.txt

##########
if [ ! -d "$WRKDIR" ]
then 
  echo "!!!!!!!!!!!! Working directory does not exist. ($WRKDIR)"
else
  touch $WRKDIR/.test
  if [ ! -e "$WRKDIR/.test" ]
  then
    echo "There is no permission for creating file in the WRKDIR"
  else
    rm $WRKDIR/.test
  fi
fi
##########
printf "== testing Python2.x:"
if [ -e .help ]; then rm .help; fi
program_PYTHON2 --version > .help 2>&1
result=`cat .help | tail -n 1 | cut -c-9`
if [ "$result" == "Python 2." ]
then
  printf " SUCCESFULL\n"
  echo "import numpy" > .test.py
  echo "a=numpy.matrix('1 2;3 4')" >> .test.py
  echo "print(a*a)" >> .test.py
  program_PYTHON2 .test.py > .test.out 2> .test.out
  result=`grep -c "[15 22]]" .test.out`
  if [ $result -eq 1 ]
  then
    echo "   -- numpy: SUCCESFULL"
  else
    echo "   -- numpy: UNSUCCESFULL"
    echo "   :: see .log for the error"
    echo "   :: your python version probably does not have numpy libraries"
  fi
  rm .test.py .test.out
else
  printf " UNSUCCESFULL\n" 
  cat .help >> .log
  echo "####################################" >> .log
  echo "   :: see .log for the error"
  echo "   :: open ~/.JKCSusersetup -> check program_PYTHON2 -> check the function/setup paths"
fi
rm .help
##########
##########
printf "== testing Python3.x:"
if [ -e .help ]; then rm .help; fi
program_PYTHON3 --version > .help 2>&1
result=`cat .help| tail -n 1 |cut -c-9`
if [ "$result" == "Python 3." ]
then
  printf " SUCCESFULL\n"
  echo "import numpy" > .test.py
  echo "a=numpy.matrix('1 2;3 4')" >> .test.py
  echo "print(a*a)" >> .test.py
  program_PYTHON3 .test.py > .test.out 2> .test.out
  result=`grep -c "[15 22]]" .test.out`
  if [ $result -eq 1 ]
  then
    echo "   -- numpy: SUCCESFULL"
  else
    echo "   -- numpy: UNSUCCESFULL"
    echo "   :: see .log for the error"
    echo "   :: your python version probably does not have numpy libraries"
  fi
  rm .test.py .test.out
else
  printf " UNSUCCESFULL\n"
  cat .help >> .log
  echo "####################################" >> .log
  echo "   :: see .log for the error"
  echo "   :: open ~/.JKCSusersetup -> check program_PYTHON3 -> check the function/setup paths"
fi
rm .help
##########
printf "== testing ABCluster:"
touch .calc.inp
program_ABC .calc.inp 2> .calc.out
result=`grep -c "Cannot read the cluster file name." .calc.out`
if [ $result -eq 1 ]
then 
  printf " SUCCESFULL\n"
else
  printf " UNSUCCESFULL\n"
  cat .calc.out >> .log
  echo "####################################" >> .log
  echo "   :: see .log for the error"
  echo "   :: open ~/.JKCSusersetup -> check program_ABC or setup properly path PATH_ABCluster"
fi
rm .calc.inp .calc.out
##########
printf "== testing XTB:"
if [ -e .test.xyz ]
then 
  rm .test.xyz
fi
touch .test.xyz
program_XTB .test.xyz > .test.log 2> .test.log
result=`grep -c "#ERROR! no atoms!" .test.log`
if [ $result -eq 1 ]
then
  printf " SUCCESFULL\n"
else
  printf " UNSUCCESFULL\n"
  cat .test.log >> .log
  echo "####################################" >> .log
  echo "   :: see .log for the error"
  echo "   :: open ~/.JKCSusersetup -> check program_XTB or setup properly path PATH_XTB"
fi
rm .test.xyz .test.log
##########
printf "== testing Gaussian(G16):"
touch .test.com
program_G16 .test.com > .test.log 2> .test.log
result=`grep -c "Route card not found." .test.log`
if [ $result -eq 1 ]
then
  printf " SUCCESFULL\n"
else
  printf " UNSUCCESFULL\n"
  cat .test.log >> .log
  echo "####################################" >> .log
  echo "   :: see .log for the error"
  echo "   :: open ~/.JKCSusersetup -> check program_G16 or setup properly path PATH_G16"
fi
rm .test.log .test.com
##########
printf "== testing GoodVibes:"
touch .test.log
program_GoodVibes .test.log > .test.out 2> /dev/null 
result=`grep -c "Warning! Couldn't " .test.out`
if [ $result -eq 1 ]
then
  printf " SUCCESFULL\n"
else
  printf " UNSUCCESFULL\n"
  cat .test.out >> .log
  echo "####################################" >> .log
  echo "   :: see .log for the error"
  echo "   :: open ~/.JKCSusersetup -> check program_GoodVibes or setup properly path PATH_GoodVibes"
fi
rm .test.out .test.log
if [ -e Goodvibes_output.dat ]; then rm Goodvibes_output.dat; fi
##########
printf "== testing ORCA:"
touch .test.inp
program_ORCA .test.inp > .test.out 2> .test.out
result=`grep -c "You must have a \[COORDS\] ... \[END\] block in your input" .test.out`
if [ $result -ne 0 ]
then
  printf " SUCCESFULL\n"
else
  printf " UNSUCCESFULL\n"
  cat .test.out >> .log
  echo "####################################" >> .log
  echo "   :: see .log for the error"
  echo "   :: open ~/.JKCSusersetup -> check program_ORCA or setup properly path PATH_ORCA"
fi
rm .test.inp .test.out .test.xyz
echo ___________________
