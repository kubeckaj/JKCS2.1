#!/bin/bash
PYTHON=$1         #python3.9
MODULE_PYTHON=$2  #"module load python-data/3.9-1"

eval $MODULE_PYTHON
### MAKING OWN ENVIRONMENT
#echo """name: jkcs
#
#channels:
#  - conda-forge
#  - defaults
#""" > jkcs.yaml

#####################################################################
### Just checking if python exists ##################################
version=$($PYTHON -V 2>&1 | awk '{print $2}')                       #
if [[ -z "$version" ]]                                              #
  then echo "No Python!"; exit; fi                                  #
#####################################################################
### Just checking if appropriate version is used ####################
parsedVersion=$(echo "${version//./}")                              #
if [[ "$parsedVersion" -gt "400" || "$parsedVersion" -lt "380" ]]   #
  then echo "Invalid Python version"; exit;fi                       #
#####################################################################

echo "- creating environment"
#conda env create --name JKCS --file jkcs.yaml
$PYTHON -m venv --system-site-packages JKCS

#echo "- loading anaconda"
#module load anaconda3/5.0.1

#echo "- checking path"
#path=`conda info | grep "base environment" | awk '{print $4}'`

#echo "- starting conda"
#. ${path}/etc/profile.d/conda.sh

#echo "- listing environmnets"
#conda env list

echo "- activating environment"
#source /home/kubeckaj/.conda/envs/JKCS/bin/activate
#source my-venv/bin/activate
#conda activate JKCS
source JKCS/bin/activate

#echo "- installing python 3.8"
#echo y | conda install -c anaconda python=3.8

#echo y | conda install xtb-python
#echo y | conda install numba
PIP="$PYTHON -m pip"
#$PIP --version
$PIP install pip --force-reinstall --upgrade
$PIP install pathlib #Perhaps this one is not necessary
$PIP install numexpr==2.7.0 --force-reinstall --upgrade
$PIP install numpy==1.21.4 --force-reinstall --upgrade
$PIP install pandas==1.3.4 --force-reinstall --upgrade
$PIP install ase
#$PIP install sklearn
#$PIP install cffi
#$PIP install dscribe
#$PIP install qml

#echo "- exporting environment"
#conda env export > environment.yml

#echo "- removing environment"
#conda remove --name JKCS --all

# update your jupyter:
#   python -m ipykernel install --user --name=jkcs

### ADDING JKpython to bashrc
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

echo " - saving JKpython alias to ~/.bashrc"
PyPATH=$(which $PYTHON)

command="alias JKpython='$Module; source $PWD/JKCS/bin/activate; alias python=$PyPATH'"
writetobashrc "$command"

#search=$(grep -c JKpython ~/.bashrc)
#if [ $search -eq 0 ]
#then
#  echo "alias JKpython='$Module; source $PWD/JKCS/bin/activate; alias python=$PyPATH'" >> ~/.bashrc
#  #echo "export PATH=$PWD/xtb-python/:\$PATH" >> ~/.bashrc
#fi
