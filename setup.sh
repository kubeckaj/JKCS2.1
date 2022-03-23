# USE
# ./setup.sh
# ./setup.sh -r  //rewrite ~/.JKCSusersetup.txt + run JKQC install
# ./setup.sh -r2 //rewrite ~/.JKCSusersetup.txt

### MODIFY ###

PYTHON="python3.9"    #Please modify this and use python version >3.8.0 but <4.0.0
MODULE_PYTHON="module load python-env/3.9-1"  #Is there some module required to load python?

###########################################################################################################
## DO NOT MODIFY

cat TOOLS/label_big

# JKCS python environment
if [ ! -e JKQC/JKCS ] || [ "$1" == "-r" ]
then
  cd JKQC
  sh .install.sh "$PYTHON" "$MODULE_PYTHON"
  if [ ! -e JKCS ]
  then
    exit
  fi
  cd ..
fi

# START MODIFYING ~/.bashrc
echo "-----------------------"
echo 'Now, I will just test what is written in your ~/.bashrc file and modify it'

# PATH DECLARATION
path=$PWD
path1="${path}/JKCSx"
path2="${path}/TOOLS/MANIPULATE"

# DOES ~/.bashrc exist?
if [ ! -e ~/.bashrc ]
then
  touch ~/.bashrc
fi

function writetobashrc {
  command="$1"
  test=`grep -c "$command" ~/.bashrc`
  if [ $test -ne 0 ]
  then 
    sed  "\#$command#d" ~/.bashrc > ~/.bashrc_help
    mv ~/.bashrc_help ~/.bashrc
  fi
  echo "$command" >> ~/.bashrc
}

command1="export PATH=$path1:\$PATH"
writetobashrc "$command1"

command2="export PATH=$path2:\$PATH"
writetobashrc "$command2"

Qcolours=1
source TOOLS/LOADING/colours.txt

echo -e "${cfRED}Write following command:${cfDEF} "
echo -e "          ${cfYELLOW}source ~/.bashrc${cfDEF}"
echo "-----------------------"

if [ ! -e ~/.JKCSusersetup.txt ] || [ "$1" == "-r" ] || [ "$1" == "-r2" ]
then
  if [ -e ~/.JKCSusersetup.txt ]; then cp ~/.JKCSusersetup.txt ~/.oldJKCSusersetup.txt; fi
  cp TOOLS/.JKCSusersetup.txt .help1
  sed 's,REPLACE_jkcs_path,'"$PWD"',g' .help1 > .help2
  sed 's,REPLACE_python,'"$PYTHON"',g' .help2 > .help3
  sed 's,REPLACE_module_python,'"$MODULE_PYTHON"',g' .help3 > ~/.JKCSusersetup.txt
  rm .help1 .help2 .help3
  echo -e "${cfRED}Please, change all required user settings (e.g. paths) in file ~/.JKCSusersetup.txt${cfDEF}"
else
  echo -e "File ~/.JKCSusersetup.txt already exists. However, check, if all paths in this file are correct."
  echo "...or use the following command the rewrite the old ~/.JKCSusersetput.txt:"
  echo -e "          ${cfYELLOW}sh setup.sh -r${cfDEF}"
fi
echo "-----------------------"
echo "Anyway, you can check if everything is working by running:"
echo -e "          ${cfYELLOW}sh test.sh${cfDEF}"
