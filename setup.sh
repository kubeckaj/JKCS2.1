# USE
# ./setup.sh
# ./setup.sh -r  //rewrite ~/.JKCSusersetup.txt + run JKQC install
# ./setup.sh -r2 //rewrite ~/.JKCSusersetup.txt

### MODIFY ###
#default setup

PYTHON="python3.9"    #Please modify this and use python version >3.8.0 but <4.0.0
MODULE_PYTHON=""  #Is there some module required to load python?
PATH_ABC="/users/kubeckaj/ABCluster-2.0-Linux/"
PATH_ABC3="/users/kubeckaj/ABCluster-3.1-Linux/"
MODULE_ABC="module load gcc"
PATH_XTB="/users/kubeckaj/XTB6.0/"
MODULE_XTB=""
PATH_CREST="/users/kubeckaj/crest/"
PATH_G16="/appl/soft/chem/gaussian/G16RevC.01/"
MODULE_G16="module load gaussian/G16RevC.01"
PATH_ORCA="/users/kubeckaj/ORCA/orca_4_2_0_linux_x86-64_shared_openmpi314/"
MODULE_ORCA="module load intel/19.0.4 hpcx-mpi/2.4.0 intel-mkl/2019.0.4"
EXTRA_ORCA_LINES=""
MODULE_GCC="module load gcc"
SBATCH_PREFIX=""
WRKDIR="./"

time1="72:00:00"
time2="330:00:00"
queue1="small"
queue2="longrun"

###########################################################################################################
## DO NOT MODIFY

QUpython=""
QUmodulepython=""
Qr=0
ADD=""
last=""
for i in "$@"
do 
  if [ "$i" == "-help" ]
  then
    echo """    sh setup.sh [:cluster:] [:arguments:] [:extra packages:]

OPTIONS (cluster):
  <empty> ...... an example setting
  grendel ...... (J.Elm AU [Grendel])
  molas ........ (N.Myllys UH [Puhti])
  puhti/mahti .. (H.Vehkamaki UH [Puhti/Mahti])
 
OPTIONS (arguments):
  -r .............. rewrite the old installation
  -r2 ............. rewrite the old installation but skip Python installation
  -up ............. (only) updates python libraries
  -python <str> ... define how do you call python (e.g. python3.9)
  -module \"<str>\" . define python module (e.g. \"module load python/3.9.4\")

OPTIONS (extra packages):
  -qml ............ quantum machine learning program
  -nn ............. schnetpack [!needs: module load gcc]
  -descriptors .... dscribe library
  -calculators .... python TBlite, XTB, ORCA
  -all ............ all above = qml,nn,descriptors,calculators
  OTHERS:
  -aimnet .......... AIMNet
  -physnet ........ physnet
  -mbdf ........... MBDF for categorization trick
  experimental: -qml-lightning

EXAMPLE:
    sh setup.sh -python python3.9 -module \"module load python\"
    sh setup.sh grendel -r -qml -descriptors
    sh setup.sh -up grendel -nn
    """
    exit
  fi
  if [ "$i" == "grendel" ]
  then
    PYTHON="python3.11"
    MODULE_PYTHON="source /com/bin/modules.sh;module load python/3.11.1"
    #PATH_ABC="/home/kubeckaj/Applications/ABCluster-2.0-Linux/"
    #PATH_ABC3="/home/kubeckaj/Applications/ABCluster-3.1-Linux/"
    PATH_ABC="/home/kubeckaj/Applications/ABCluster-3.3pre-Linux/"
    PATH_ABC3="/home/kubeckaj/Applications/ABCluster-3.3pre-Linux/"
    MODULE_ABC="module load gcc"
    PATH_XTB="/home/kubeckaj/Applications/XTB6.7.0/"
    MODULE_XTB=""
    PATH_CREST="/home/kubeckaj/Applications/crest/"
    PATH_G16="/comm/groupstacks/gaussian/gaussian/gaussian16/Rev.B.01/"
    MODULE_G16="source /com/bin/modules.sh;source /comm/groupstacks/gaussian/bin/modules.sh --silent; module load gaussian16/Rev.B.01;"
    #PATH_ORCA="/comm/groupstacks/chemistry/apps/orca/5.0.3/"
    #MODULE_ORCA="source /comm/groupstacks/gaussian/bin/modules.sh --silent; module load orca/5.0.3"
    #PATH_ORCA="/home/kubeckaj/Applications/orca_5_0_4_linux_x86-64_shared_openmpi411/"
    #MODULE_ORCA="module load gcc; module load openmpi"
    PATH_ORCA="/comm/groupstacks/chemistry/apps/orca/5.0.4/"
    #MODULE_ORCA="ulimit -c 0;source /com/bin/modules.sh;source /comm/groupstacks/chemistry/bin/modules.sh;ml orca/5.0.4;ml gcc/9.1.0;ml openmpi/3.1.3;"
    MODULE_ORCA="ulimit -c 0;source /comm/groupstacks/chemistry/bin/modules.sh; source /com/bin/modules.sh;ml orca/5.0.4;ml gcc/9.1.0;ml openmpi;export OMPI_MCA_btl_openib_allow_ib=1;export OMPI_MCA_btl=^openib;export OMPI_MCA_mca_base_component_show_load_errors=0"
    MODULE_GCC="module load gcc/11.1.0"
    SBATCH_PREFIX=""
    WRKDIR="/scratch/\\\$SLURM_JOB_ID/"
    time1="10-00:00:00"
    time2="10-00:00:00"
    queue1="q64,q48,q40,q36"
    queue2="q64,q48,q40,q36"
    continue
  fi
  if [ "$i" == "mahti" ]
  then
    PYTHON="python3.9"                             #Please modify this and use python version >3.8.0 but <4.0.0
    MODULE_PYTHON="module load python-data/3.9-3"    #Is there some module required to load python?
    PATH_ABC="/users/kubeckaj/ABCluster-2.0-Linux/"
    PATH_ABC3=$PATH_ABC"../ABCluster-3.1-Linux/"
    MODULE_ABC="module load gcc"
    PATH_XTB="/users/kubeckaj/XTB6.4/"
    MODULE_XTB=""
    PATH_CREST=""
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
    MODULE_PYTHON="module load python-data/3.9-22.04"  #Is there some module required to load python?
    PATH_ABC="/users/kubeckaj/ABCluster-3.0-Linux/"
    PATH_ABC3="/users/ineefjes/Applications/ABCluster-3.0-Linux/"
    MODULE_ABC="module load gcc"
    PATH_XTB="/users/ineefjes/Applications/xtb-6.5.1/"
    MODULE_XTB=""
    PATH_CREST=""
    PATH_G16="/appl/soft/chem/gaussian/G16RevC.02/"
    MODULE_G16="module load gaussian"
    PATH_ORCA="/users/ineefjes/Applications/orca_5_0_4_linux_x86-64_shared_openmpi411"
    MODULE_ORCA="module purge; module load gcc/11.3.0 openmpi/4.1.4 intel-oneapi-mkl/2022.1.0"
    EXTRA_ORCA_LINES="ORTERUN=\\\\\`which orterun\\\\\`\nln -sf \\\\\${ORTERUN}  \\\\\${SLURM_SUBMIT_DIR}/mpirun\nexport PATH=\\\\\${SLURM_SUBMIT_DIR}:\\\\\${PATH}\n"
    project="hvehkama"
    #project=`csc-projects | grep Owner | awk '{print $2}' | grep -v $USER | grep -v gaussian`
    SBATCH_PREFIX="--account=$project "
    WRKDIR="./"
    continue
  fi
  if [ "$i" == "molas" ]
  then
    PYTHON="python3.9"    #Please modify this and use python version >3.8.0 but <4.0.0
    MODULE_PYTHON=""  #Is there some module required to load python?
    PATH_ABC="/projappl/project_2006166/APP/ABCluster3.1"
    PATH_ABC3="/projappl/project_2006166/APP/ABCluster3.1"
    MODULE_ABC="module load gcc"
    PATH_XTB="/projappl/project_2006166/APP/XTB6.5.1"
    MODULE_XTB=""
    PATH_CREST="/projappl/project_2006166/APP/crest/"
    PATH_G16="/appl/soft/chem/gaussian/G16RevC.02/"
    MODULE_G16="module load gaussian"
    PATH_ORCA="/projappl/project_2006166/APP/ORCA5.0.3"
    #MODULE_ORCA="module load intel/19.0.4 hpcx-mpi/2.4.0 intel-mkl/2019.0.4"
    #MODULE_ORCA="module purge; module load gcc/11.3.0 openmpi/4.1.4 intel-oneapi-mkl/2022.1.0"
    #MODULE_ORCA="module load intel/19.0.4 hpcx-mpi/2.4.0 intel-mkl/2019.0.4"
    MODULE_ORCA="module purge; module load gcc/11.3.0 openmpi/4.1.4 intel-oneapi-mkl/2022.1.0"
    EXTRA_ORCA_LINES="ORTERUN=\\\\\`which orterun\\\\\`\nln -sf \\\\\${ORTERUN}  \\\\\${SLURM_SUBMIT_DIR}/mpirun\nexport PATH=\\\\\${SLURM_SUBMIT_DIR}:\\\\\${PATH}\n"
    project="project_2001079"  #project_2006166
    #project=`csc-projects | grep Owner | awk '{print $2}' | grep -v $USER | grep -v gaussian`
    SBATCH_PREFIX="--account=$project "
    WRKDIR="./"
    continue
  fi
  if [ "$i" == "-all" ]
  then
    ADD+=" -qml -xtb -descriptors -nn"
    continue
  fi
  if [ "$i" == "-aimnet" ]
  then
    ADD+=" -aimnet "
    continue
  fi
  if [ "$i" == "-physnet" ]
  then
    ADD+=" -physnet "
    continue
  fi  
  if [ "$i" == "-xtb" ] || [ "$i" == "-calculators" ]
  then
    ADD+=" -xtb "
    continue
  fi
  if [ "$i" == "-qml" ] || [ "$i" == "qml" ]
  then
    ADD+=" -qml "
    continue
  fi
  if [ "$i" == "-mbdf" ] || [ "$i" == "mbdf" ]
  then
    ADD+=" -mbdf "
    continue
  fi
  if [ "$i" == "-descriptors" ] || [ "$i" == "descriptors" ]
  then
    ADD+=" -descriptors "
    continue
  fi
  if [ "$i" == "-nn" ] || [ "$i" == "nn" ]
  then
    ADD+=" -nn "
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
  if [ "$i" == "-up" ] || [ "$i" == "-update" ] || [ "$i" == "update" ] || [ "$i" == "up" ]
  then
    Qr=3
    continue
  fi
  if [ "$i" == "-python" ]
  then
    last="-python"
    continue
  fi
  if [ "$last" == "-python" ]
  then
    QUpython="$i"
    last=""
    continue
  fi
  if [ "$i" == "-module" ]
  then
    last="-module"
    continue
  fi
  if [ "$last" == "-module" ]
  then
    QUmodulepython="$i"
    last=""
    continue
  fi
  echo "I do not understand argument: $i [EXITING]"
  exit
done
if [ ! -z "$QUpython" ]; then PYTHON="$QUpython"; fi
if [ ! -z "$QUmodulepython" ]; then MODULE_PYTHON="$QUmodulepython"; fi

cat TOOLS/label_big

# JKCS python environment
if [ ! -e JKQC/JKCS ] || [ "$Qr" == "1" ] || [ "$Qr" == "3" ]
then
  cd JKQC
  if  [ "$Qr" != "3" ]
  then
    rm -rf JKCS
    if [ -e ../JKML/src/PhysNet_DER ]
    then
      rm -rf ../JKML/src/PhysNet_DER
    fi 
  fi
  sh .install.sh "$PYTHON" "$MODULE_PYTHON" "$ADD"
  if [ ! -e JKCS ]
  then
    exit
  fi
  cd ..
fi

if  [ "$Qr" != "3" ]
then
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
  
  if [ ! -e TOOLS/LOADING/JK_ACDC/ACDC_input/acdc_2020_07_20.pl ]
  then 
    cd TOOLS/LOADING/JK_ACDC/ACDC_input
    tar -xzf acdc_2020_07_20.pl.tar.gz
    cd -
  fi
  
  printf "${cfGREEN}Write following command:${cfDEF}\n"
  printf "          ${cfYELLOW}source ~/.bashrc${cfDEF}\n"
  echo "-----------------------"
  
  if [ ! -e ~/.JKCSusersetup.txt ] || [ $Qr -gt 0 ] 
  then
    if [ -e ~/.JKCSusersetup.txt ]; then cp ~/.JKCSusersetup.txt ~/.oldJKCSusersetup.txt; fi
    cp TOOLS/.JKCSusersetup.txt .help1
    sed 's,REPLACE_jkcs_path,'"$PWD"',g' .help1 > .help2
    sed 's,REPLACE_python,'"$PYTHON"',g' .help2 > .help3
    sed 's,REPLACE_module_python,'"$MODULE_PYTHON"',g' .help3 > .help4
    sed 's,REPLACE_abc3,'"$PATH_ABC3"',g' .help4 > .help4b
    sed 's,REPLACE_abc,'"$PATH_ABC"',g' .help4b > .help5
    sed 's,REPLACE_module_abc,'"$MODULE_ABC"',g' .help5 > .help6
    sed 's,REPLACE_xtb,'"$PATH_XTB"',g' .help6 > .help6b
    sed 's,REPLACE_crest,'"$PATH_CREST"',g' .help6b > .help7
    sed 's,REPLACE_module_xtb,'"$MODULE_XTB"',g' .help7 > .help8
    sed 's,REPLACE_g16,'"$PATH_G16"',g' .help8 > .help9
    sed 's,REPLACE_module_g16,'"$MODULE_G16"',g' .help9 > .help10
    sed 's,REPLACE_orca,'"$PATH_ORCA"',g' .help10 > .help11
    sed 's,REPLACE_module_orca,'"$MODULE_ORCA"',g' .help11 > .help11b
    sed 's,REPLACE_extra_orca_lines,'"$EXTRA_ORCA_LINES"',g' .help11b > .help12
    sed 's,REPLACE_sbatch_prefix,'"$SBATCH_PREFIX"',g' .help12 > .help13
    sed 's,REPLACE_wrkdir,'"$WRKDIR"',g' .help13 > .help14
    sed 's,REPLACE_time1,'"$time1"',g' .help14 > .help15
    sed 's,REPLACE_time2,'"$time2"',g' .help15 > .help16
    sed 's/REPLACE_queue1/'"$queue1"'/g' .help16 > .help17
    sed 's/REPLACE_queue2/'"$queue2"'/g' .help17 > .help18
    sed 's,REPLACE_gcc,'"$MODULE_GCC"',g' .help18 > .help19
    mv .help19 ~/.JKCSusersetup.txt
    rm .help*
    printf "${cfGREEN}Please, change all required user settings (e.g. paths) in file ~/.JKCSusersetup.txt${cfDEF}\n"
  else
    echo "File ~/.JKCSusersetup.txt already exists. However, check, if all paths in this file are correct."
    echo "...or use the following command the rewrite the old ~/.JKCSusersetput.txt:"
    printf "          ${cfYELLOW}sh setup.sh -r${cfDEF}\n"
  fi
  echo "-----------------------"
  echo "Anyway, you can check if everything is working by running:"
  printf "          ${cfYELLOW}sh test.sh${cfDEF}\n"
else
  printf "  Update finished\n"
fi
