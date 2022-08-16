# USE
# ./setup.sh
# ./setup.sh -r  //rewrite ~/.JKCSusersetup.txt + run JKQC install
# ./setup.sh -r2 //rewrite ~/.JKCSusersetup.txt

### MODIFY ###
#default setup

PYTHON="python3.9"    #Please modify this and use python version >3.8.0 but <4.0.0
MODULE_PYTHON=""  #Is there some module required to load python?
PATH_ABC="/users/kubeckaj/ABCluster-2.0-Linux/"
MODULE_ABC="module load gcc"
PATH_XTB="/users/kubeckaj/XTB6.0/"
MODULE_XTB=""
PATH_G16="/appl/soft/chem/gaussian/G16RevC.01/"
MODULE_G16="module load gaussian/G16RevC.01"
PATH_ORCA="/users/kubeckaj/ORCA/orca_4_2_0_linux_x86-64_shared_openmpi314/"
MODULE_ORCA="module load intel/19.0.4 hpcx-mpi/2.4.0 intel-mkl/2019.0.4"
SBATCH_PREFIX=""
WRKDIR="./"

###########################################################################################################
## DO NOT MODIFY

cat TOOLS/label_big

Qr=0
for i in "$@"
do 
  if [ "$i" == "-help" ]
  then
    echo "see https://jkcs.readthedocs.io/en/latest/JKCSSetupAndInstallation.html"
    exit
  fi
  if [ "$i" == "grendel" ]
  then
    PYTHON="python3.9"
    MODULE_PYTHON="module load python/3.9.4"
    PATH_ABC="/home/kubeckaj/Applications/ABCluster-2.0-Linux/"
    MODULE_ABC="module load gcc"
    PATH_XTB="/home/kubeckaj/Applications/XTB6.4/"
    MODULE_XTB=""
    PATH_G16="/comm/groupstacks/gaussian/gaussian/gaussian16/Rev.B.01/"
    MODULE_G16="source /comm/groupstacks/gaussian/bin/modules.sh --silent; module load gaussian16/Rev.B.01;"
    PATH_ORCA="/comm/groupstacks/chemistry/apps/orca/5.0.3/"
    MODULE_ORCA="source /comm/groupstacks/gaussian/bin/modules.sh --silent; module load orca/5.0.3"
    SBATCH_PREFIX=""
    WRKDIR="/scratch/\\\$SLURM_JOB_ID/"
    continue
  fi
  if [ "$i" == "mahti" ]
  then
    PYTHON="python3.9"                             #Please modify this and use python version >3.8.0 but <4.0.0
    MODULE_PYTHON="module load python-data/3.9-3"    #Is there some module required to load python?
    PATH_ABC="/users/kubeckaj/ABCluster-2.0-Linux/"
    MODULE_ABC="module load gcc"
    PATH_XTB="/users/kubeckaj/XTB6.4/"
    MODULE_XTB=""
    PATH_G16="/appl/soft/chem/gaussian/G16RevC.01/"
    MODULE_G16="module load gaussian/G16RevC.01"
    PATH_ORCA="/users/kubeckaj/ORCA/orca_4_2_0_linux_x86-64_shared_openmpi314/"
    MODULE_ORCA="module load intel/19.0.4 hpcx-mpi/2.4.0 intel-mkl/2019.0.4"
    project=`csc-projects | grep Owner | awk '{print $2}' | grep -v $USER | grep -v gaussian`
    SBATCH_PREFIX="--account=$project "
    WRKDIR="./"
    continue
  fi
  if [ "$i" == "puhti" ]
  then
    PYTHON="python3.9"    #Please modify this and use python version >3.8.0 but <4.0.0
    MODULE_PYTHON="module load python-data/3.9-1"  #Is there some module required to load python?
    PATH_ABC="/users/kubeckaj/ABCluster-2.0-Linux/"
    MODULE_ABC="module load gcc"
    PATH_XTB="/users/kubeckaj/XTB6.0/"
    MODULE_XTB=""
    PATH_G16="/appl/soft/chem/gaussian/G16RevC.01/"
    MODULE_G16="module load gaussian/G16RevC.01"
    PATH_ORCA="/users/kubeckaj/ORCA/orca_4_2_0_linux_x86-64_shared_openmpi314/"
    MODULE_ORCA="module load intel/19.0.4 hpcx-mpi/2.4.0 intel-mkl/2019.0.4"
    project=`csc-projects | grep Owner | awk '{print $2}' | grep -v $USER | grep -v gaussian`
    SBATCH_PREFIX="--account=$project "
    WRKDIR="./"
    continue
  fi
  if [ "$i" == "-qml" ] || [ "$i" == "qml" ]
  then
    ADD="-qml"
    continue
  fi
  if [ "$i" == "-r" ]
  then
    Qr=1
    continue
  fi
  if [ "$i" == "-r2" ]
  then
    Qr=2
    continue
  fi 
  echo "I do not understand argument: $i [EXITING]"
  exit
done

# JKCS python environment
if [ ! -e JKQC/JKCS ] || [ "$Qr" == "1" ] 
then
  cd JKQC
  rm -r JKCS 2>/dev/null 
  sh .install.sh "$PYTHON" "$MODULE_PYTHON" "$ADD"
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

printf "${cfRED}Write following command:${cfDEF}\n"
printf "          ${cfYELLOW}source ~/.bashrc${cfDEF}\n"
echo "-----------------------"

if [ ! -e ~/.JKCSusersetup.txt ] || [ $Qr -gt 0 ] 
then
  if [ -e ~/.JKCSusersetup.txt ]; then cp ~/.JKCSusersetup.txt ~/.oldJKCSusersetup.txt; fi
  cp TOOLS/.JKCSusersetup.txt .help1
  sed 's,REPLACE_jkcs_path,'"$PWD"',g' .help1 > .help2
  sed 's,REPLACE_python,'"$PYTHON"',g' .help2 > .help3
  sed 's,REPLACE_module_python,'"$MODULE_PYTHON"',g' .help3 > .help4
  sed 's,REPLACE_abc,'"$PATH_ABC"',g' .help4 > .help5
  sed 's,REPLACE_module_abc,'"$MODULE_ABC"',g' .help5 > .help6
  sed 's,REPLACE_xtb,'"$PATH_XTB"',g' .help6 > .help7
  sed 's,REPLACE_module_xtb,'"$MODULE_XTB"',g' .help7 > .help8
  sed 's,REPLACE_g16,'"$PATH_G16"',g' .help8 > .help9
  sed 's,REPLACE_module_g16,'"$MODULE_G16"',g' .help9 > .help10
  sed 's,REPLACE_orca,'"$PATH_ORCA"',g' .help10 > .help11
  sed 's,REPLACE_module_orca,'"$MODULE_ORCA"',g' .help11 > .help12
  sed 's,REPLACE_sbatch_prefix,'"$SBATCH_PREFIX"',g' .help12 > .help13
  sed 's,REPLACE_wrkdir,'"$WRKDIR"',g' .help13 > .help14
  mv .help14 ~/.JKCSusersetup.txt
  rm .help*
  printf "${cfRED}Please, change all required user settings (e.g. paths) in file ~/.JKCSusersetup.txt${cfDEF}\n"
else
  echo "File ~/.JKCSusersetup.txt already exists. However, check, if all paths in this file are correct."
  echo "...or use the following command the rewrite the old ~/.JKCSusersetput.txt:"
  printf "          ${cfYELLOW}sh setup.sh -r${cfDEF}\n"
fi
echo "-----------------------"
echo "Anyway, you can check if everything is working by running:"
printf "          ${cfYELLOW}sh test.sh${cfDEF}\n"
