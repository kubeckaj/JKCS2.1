
echo ___________________
echo TESTING JKCS SETUP:

if [ ! -e ~/.JKCSusersetup.txt ]
then
  echo "Did you run already setup.sh? [The ~/..JKCSusersetup.txt is missing]"
fi
source ~/.JKCSusersetup.txt

##########
printf "== testing Python2.x:"
if [ -e .help ]; then rm .help; fi
program_PYTHON2 --version > .help 2>&1
result=`cat .help| cut -c-9`
if [ "$result" == "Python 2." ]
then
  printf " SUCCESSFULL\n"
else
  printf " UNSUCCESSFULL\n" 
  cat .help
  echo ":: open ~/.JKCSusersetup -> find program_PYTHON2 -> check the function/setup paths"
  exit
fi
rm .help
##########
##########
printf "== testing Python3.x:"
if [ -e .help ]; then rm .help; fi
program_PYTHON3 --version > .help 2>&1
result=`cat .help| cut -c-9`
if [ "$result" == "Python 3." ]
then
  printf " SUCCESSFULL\n"
else
  printf " UNSUCCESSFULL\n"
  cat .help
  echo ":: open ~/.JKCSusersetup -> find program_PYTHON3 -> check the function/setup paths"
fi
rm .help
##########
printf "== testing ABCluster:"
touch .calc.inp
program_ABC .calc.inp
result=`grep -c "Cannot read the cluster file name." .calc.out`
if [ $result -eq 1 ]
then 
  printf " SUCCESSFULL\n"
else
  printf " UNSUCCESSFULL\n"
  cat .calc.out
  echo ":: open ~/.JKCSusersetup -> find program_ABC -> check the function/or re-setup path PATH_ABCluster"
fi
rm .calc.inp .calc.out
##########
printf "== testing XTB:"
touch .test.xyz
program_XTB .test.xyz
result=`grep -c "#ERROR! no atoms!" .test.log`
if [ $result -eq 1 ]
then
  printf " SUCCESSFULL\n"
else
  printf " UNSUCCESSFULL\n"
  cat .test.log
  echo ":: open ~/.JKCSusersetup -> find program_XTB -> check the function/or re-setup path PATH_XTB"
fi
rm .test.xyz .test.log
##########
printf "== testing Gaussian(G16):"
touch .test.com
program_G16 .test.com
result=`grep -c "Route card not found." .test.log`
if [ $result -eq 1 ]
then
  printf " SUCCESSFULL\n"
else
  printf " UNSUCCESSFULL\n"
  cat .test.log
  echo ":: open ~/.JKCSusersetup -> find program_G16 -> check the function/or re-setup path PATH_G16"
fi
rm .test.xyz .test.log .test.com
##########
printf "== testing GoodVibes:"
touch .test.log
program_GoodVibes .test.log > .test.out 2> /dev/null 
result=`grep -c "Warning! Couldn't find frequency information ..." .test.out`
if [ $result -eq 1 ]
then
  printf " SUCCESSFULL\n"
else
  printf " UNSUCCESSFULL\n"
  cat .test.out
  echo ":: open ~/.JKCSusersetup -> find program_GoodVibes -> check the function/or re-setup path PATH_GoodVibes"
fi
rm .test.out .test.log
if [ -e Goodvibes_output.dat ]; then rm Goodvibes_output.dat; fi
#PATH_GoodVibes
echo ___________________
