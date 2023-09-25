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
#TMPDIR='$PWD/JKCS/'

#echo "- installing python 3.8"
#echo y | conda install -c anaconda python=3.8

#ADD="--force-reinstall --upgrade"

#echo y | conda install xtb-python
#echo y | conda install numba
PIP="$PYTHON -m pip " #--cache-dir=$PWD/JKCS/"
#$PIP --version
$PIP install --upgrade pip
$PIP install scipy==1.9.3 $ADD 
$PIP install pathlib #Perhaps this one is not necessary
$PIP install numexpr==2.7.0 $ADD
#NUMPY: the numpy is reinstalled later in the script because lapjv would not run with this version
#this version is however needed for QML. When QML is installed with the newest version, it works ...I guess
$PIP install numpy==1.21.4 $ADD #1.23.0 had some issues for qml np.distutils or something like that
$PIP install pandas==1.3.4 $ADD
$PIP install ase 
if [[ "$*" == *"-descriptors"* ]]
then
 $PIP install sklearn
 $PIP install cffi
 $PIP install dscribe
 #ADD="--force-reinstall --upgrade"
 #$PIP install install git+https://github.com/qmlcode/qml@develop $ADD
fi
if [[ "$*" == *"-qml "* ]]
then
  #$PIP install sklearn
  #$PIP install cffi
  #$PIP install dscribe
  #ADD="--force-reinstall --upgrade"
  $PIP install git+https://github.com/qmlcode/qml@develop $ADD
fi
if [[ "$*" == *"-mbdf"* ]]
then
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
  currdir=$PWD
  cd JKCS/lib64/py*/site-packages/
  git clone https://github.com/nickjbrowning/qml-lightning.git
  #LDFLAGS=-L/usr/local/cuda-11.4
  $PYTHON setup.py build
  cd $currdir
fi
if [[ "$*" == *"-xtb"* ]]
then
  export CC=gcc
  $PIP install tblite
  ## Use following in Python script:
  #from tblite.ase import TBLite 
  #atoms.calc = TBLite(method="GFN1-xTB")
fi
$PIP install numpy==1.23.0
if [[ "$*" == *"-nn"* ]]
then
  $PIP install torch
  $PIP install lightning
  $PIP install hydra-core
  $PIP install schnetpack
  #AGOX was not able to use schnet calculator, this will resolve it:
  sed -i "s/elif type(atoms) == ase.Atoms:/elif type(atoms) == ase.Atoms or issubclass(type(atoms), ase.Atoms):/" JKCS/lib64/pyth*/site-packages/schnetpack/interfaces/ase_interface.py
  sed -i "s/elif type(atoms) == ase.Atoms:/elif type(atoms) == ase.Atoms or issubclass(type(atoms), ase.Atoms):/" JKCS/lib/pyth*/site-packages/schnetpack/interfaces/ase_interface.py
  $PIP install pytorch-lightning==2.0.6
  $PIP install tensorboard
fi

#NOT IMPORTANT FOR ALL BUT CHEAP TO INSTALL
$PIP install xlsxwriter

#ArbAlign stuff:
$PIP install lapjv
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
