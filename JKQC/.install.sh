#!/bin/bash
PYTHON=$1         #python3.9
if [ ! -z "$2" ]
then
  last=`echo $2 | rev | cut -c-1`
  if [ "$last" == ";" ]
  then
    MODULE_PYTHON="$2"  #"module load python-data/3.9-1"
  else
    MODULE_PYTHON="$2;"  #"module load python-data/3.9-1"
  fi
fi

echo Let me install some python packages
eval $MODULE_PYTHON
### MAKING OWN ENVIRONMENT
#echo """name: jkcs
#
#channels:
#  - conda-forge
#  - defaults
#""" > jkcs.yaml

$PYTHON -V

#####################################################################
### Just checking if python exists ##################################
version=$($PYTHON -V 2>&1 | awk '{print $2}')                       #
if [[ -z "$version" ]]                                              #
  then echo "No Python!"; exit; fi                                  #
#####################################################################
### Just checking if appropriate version is used ####################
parsedVersion=$(echo "${version//./}" | cut -c-3)                   #
if [[ "$parsedVersion" -gt "400" || "$parsedVersion" -lt "380" ]]   #
  then echo "Invalid Python version ($version)"; exit;fi            #
#####################################################################

echo "- creating environment"
#conda env create --name JKCS --file jkcs.yaml
$PYTHON -m venv JKCS
#$PYTHON -m venv --no-site-packages JKCS
#$PYTHON -m venv --system-site-packages JKCS

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
#TMPDIR='$PWD/JKCS/'

#echo "- installing python 3.8"
#echo y | conda install -c anaconda python=3.8

#ADD="--force-reinstall --upgrade"

#echo y | conda install xtb-python
#echo y | conda install numba
PIP="$PYTHON -m pip --no-cache-dir" #--cache-dir=$PWD/JKCS/"

if [[ "$*" == *"-mbdf"* ]]
then
  echo "======================"
  echo "MBDF requires some old numpy, here you go:"
  $PIP install numpy==1.21.4
  echo "======================"
  currdir=$PWD
  cd JKCS/lib64/py*/site-packages/
  git clone https://github.com/dkhan42/MBDF.git
  cp MBDF/MBDF.py .
  git clone https://github.com/dkhan42/qml2.git
  cd qml2
  cp readme.rst README.rst
  $PYTHON setup.py install
  cd $currdir
fi

#pyarrow==15.0.0 ????

echo "======================"
$PIP install --upgrade pip
echo "======================"
$PIP install numpy==1.26.3
echo "======================"
$PIP install scipy==1.9.3   #I need this version for ArbAlign 
echo "======================"
$PIP install joblib==1.3.2  
echo "======================"
$PIP install ase==3.22.1
echo "======================"
$PIP install pathlib==1.0.1 #Perhaps this one is not necessary
echo "======================"
$PIP install numexpr==2.8.4 #2.7.0
echo "======================"
$PIP install pyarrow==15.0.0 #this bullshit is required by new pandas
echo "======================"
$PIP install pandas==2.2.0  #
echo "======================"
$PIP install lapjv==1.3.1
echo "======================"
$PIP install xlsxwriter     #important only for IMoS output but cheap to install
if [[ "$*" == *"-descriptors"* ]]
then
  echo "======================"
  $PIP install scikit-learn
  echo "======================"
  $PIP install cffi
  echo "======================"
  $PIP install dscribe
fi
if [[ "$*" == *"-qml "* ]]
then
  if [[ "$*" != *"-descriptors"* ]]
  then
    echo "======================"
    $PIP install scikit-learn
  fi
  echo "======================"
  $PIP install qml          #$PIP install git+https://github.com/qmlcode/qml@develop $ADD
fi
if [[ "$*" == *"-mbdf"* ]]
then
  echo "======================"
  currdir=$PWD
  cd JKCS/lib64/py*/site-packages/
  git clone https://github.com/dkhan42/MBDF.git
  cp MBDF/MBDF.py .
  git clone https://github.com/dkhan42/qml2.git
  cd qml2
  cp readme.rst README.rst
  $PYTHON setup.py install 
  cd $currdir
fi
if [[ "$*" == *"-qml-lightning"* ]]
then
  echo "======================"
  currdir=$PWD
  cd JKCS/lib64/py*/site-packages/
  git clone https://github.com/nickjbrowning/qml-lightning.git
  #LDFLAGS=-L/usr/local/cuda-11.4
  $PYTHON setup.py build
  cd $currdir
fi
if [[ "$*" == *"-xtb"* ]]
then
  echo "======================"
  export CC=gcc
  $PIP install tblite
  ## Use following in Python script:
  #from tblite.ase import TBLite 
  #atoms.calc = TBLite(method="GFN1-xTB")
fi
if [[ "$*" == *"-nn"* ]]
then
  echo "======================"
  $PIP install torch
  echo "======================"
  $PIP install lightning
  echo "======================"
  $PIP install hydra-core
  echo "======================"
  $PIP install schnetpack
  #AGOX was not able to use schnet calculator, this will resolve it:
  sed -i "s/elif type(atoms) == ase.Atoms:/elif type(atoms) == ase.Atoms or issubclass(type(atoms), ase.Atoms):/" JKCS/lib64/pyth*/site-packages/schnetpack/interfaces/ase_interface.py
  sed -i "s/elif type(atoms) == ase.Atoms:/elif type(atoms) == ase.Atoms or issubclass(type(atoms), ase.Atoms):/" JKCS/lib/pyth*/site-packages/schnetpack/interfaces/ase_interface.py
  echo "======================"
  $PIP install pytorch-lightning==2.0.6
  echo "======================"
  $PIP install tensorboard
fi

echo "======================"
#ArbAlign stuff:
cp ../TOOLS/SCRIPTS/modifiedArbAlign.py JKCS/lib/$PYTHON/site-packages/ArbAlign.py

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

command="alias JKpython='$MODULE_PYTHON source $PWD/JKCS/bin/activate; alias python=$PyPATH'"
writetobashrc "$command"

#search=$(grep -c JKpython ~/.bashrc)
#if [ $search -eq 0 ]
#then
#  echo "alias JKpython='$Module; source $PWD/JKCS/bin/activate; alias python=$PyPATH'" >> ~/.bashrc
#  #echo "export PATH=$PWD/xtb-python/:\$PATH" >> ~/.bashrc
#fi
