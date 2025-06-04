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
then                                                                #
  if [[ "$parsedVersion" -lt "310" || "$parsedVersion" -gt "320" ]] #
  then                                                              #
    echo "Invalid Python version ($version/$parsedVersion)";        #
    exit;                                                           #
  fi                                                                #
fi                                                                  #
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
PIP="$PYTHON -m pip --no-cache-dir --disable-pip-version-check " #--cache-dir=$PWD/JKCS/"

echo "======================"
$PIP install --upgrade pip
#$PIP install pip==22.0.4
#$PIP install ipython
rm -r JKCS/lib/python*/site-packages/~* 2>/dev/null
if [[ "$*" == *"-mbdf"* ]]
then
  echo "======================"
  echo "MBDF requires some old numpy, here you go:"
  $PIP install numpy==1.21.4 --force-reinstall
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
echo "======================"
$PIP install scipy==1.9.3   #I need this version for ArbAlign 
#QML
if [[ "$*" == *"-qml "* ]]
then
  echo "======================"
  
  ####TODO might require pip downgrade
  ####$PIP install numpy==1.21.4 --force-reinstall
  ####$PIP uninstall numpy
  ###$PIP install numpy==1.25.2 
  ####$PIP install numpy==1.21.1
  ###echo "======================"
  ###$PIP install scikit-learn
  ###echo "======================"
  ####$PYTHON -m pip install qml          #
  ####$PIP install git+https://github.com/qmlcode/qml@develop $ADD
  ####$PIP install git+https://github.com/qmlcode/qml@develop
  ####$PYTHON -m pip install qml  --global-option="build" --global-option="--compiler=intelem" --global-option="--fcompiler=intelem"
  ###cd JKCS
  ###if [ ! -e qml ]
  ###then
  ###  git clone https://github.com/qmlcode/qml.git
  ###fi
  ###cd qml
  echo "======================"
  $PIP install scikit-learn
  if ! command -v gcc 2>&1 >/dev/null
  then
    echo "I cannot see gcc so I am trying to: module load gcc"
    module load gcc
  fi
  ###sed -i 's/ -lpthread/-lpthread/' setup.py   # remove the stray space
  ###$PIP install .
  ###cd ../..
  if ! command -v lapack 2>&1 >/dev/null
  then
    echo "I cannot see lapack so I am trying to: module load lapack"
    module load lapack
  fi
  $PIP install --no-deps qmllib
fi
echo "======================"
$PIP install lapjv
echo "======================"
$PIP install numpy==1.25.2
echo "======================"
$PIP install joblib==1.3.2  
echo "======================"
$PIP install ase==3.24.0
#FOR UMBRELLA SAMPLING
test=`grep -c "CS = atoms.constraints" JKCS/lib/python*/site-packages/ase/calculators/calculator.py`
if [ $test -eq 0 ]
then
  #792 for xxx
  #881 for 3.24.0
  sed -i '881s/.*/            CS = atoms.constraints\n            del atoms.constraints\n            self.atoms = atoms.copy()\n            atoms.set_constraint(CS)/' JKCS/lib/python*/site-packages/ase/calculators/calculator.py
fi
echo "======================"
$PIP install pathlib==1.0.1 #Perhaps this one is not necessary
echo "======================"
$PIP install numexpr==2.8.4 #2.7.0
#echo "======================"
#$PIP install pyarrow==15.0.0 #this bullshit is required by new pandas
echo "======================"
$PIP install pandas==1.3.4  #2.2.0  #for the newer one, I have to get rid of append
echo "======================"
$PIP install psutil
echo "======================"
$PIP install xlsxwriter     #important only for IMoS output but cheap to install
echo "======================"
$PIP install rdkit     #important only for IMoS output but cheap to install
#LAPJV can be installed directlty:
#echo "======================"
#module load gcc 2>/dev/null
#$PIP install git+https://github.com/src-d/lapjv
#$PIP install lapjv==1.3.1  #this bitch has some issues to see numpy


if [[ "$*" == *"-descriptors"* ]]
then
  echo "======================"
  $PIP install scikit-learn
  echo "======================"
  $PIP install cffi
  echo "======================"
  $PIP install dscribe
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


if [[ "$*" == *"-xtb"* ]] || [[ "$*" == *"-calcs"* ]]
then
  echo "======================"
  export CC=gcc
  module load gcc 2>/dev/null
  $PIP install tblite
  $PIP install XTB
  $PIP install ORCA
  cp ../TOOLS/LOADING/JKase.py JKCS/lib/python*/site-packages/tblite/ase.py
  sed -i "s/print(_nentries)//" JKCS/lib64/pyth*/site-packages/tblite/library.py
  sed -i "s/print(_dict_py)//" JKCS/lib64/pyth*/site-packages/tblite/library.py
  ## Use following in Python script:
  #from tblite.ase import TBLite 
  #atoms.calc = TBLite(method="GFN1-xTB")
fi


if [[ "$*" == *"-nn"* ]]
then
  echo "======================"
  $PIP install sympy==1.13.2
  echo "======================"
  $PIP install torch==2.4.0
  echo "======================"
  $PIP install torchmetrics==1.0.1
  echo "======================"
  $PIP install torch_ema==0.3
  echo "======================"
  $PIP install pytorch-lightning==2.0.6
  echo "======================"
  $PIP install lightning
  echo "======================"
  $PIP install hydra-core
  echo "======================"
  $PIP install tensorboard==2.17.1
  echo "======================"
  $PIP install tensorboardX==2.6.2.2
  echo "======================"
  $PIP install distlib pygments platformdirs pathspec nodeenv mypy-extensions mdurl identify fasteners dirsync colorlog click cfgv virtualenv scipy markdown-it-py h5py black rich pre-commit hydra-colorlog matscipy 
  echo "======================"
  $PIP install --no-deps schnetpack==2.0.4
  #AGOX was not able to use schnet calculator, this will resolve it:
  sed -i "s/elif type(atoms) == ase.Atoms:/elif type(atoms) == ase.Atoms or issubclass(type(atoms), ase.Atoms):/" JKCS/lib64/pyth*/site-packages/schnetpack/interfaces/ase_interface.py
  sed -i "s/elif type(atoms) == ase.Atoms:/elif type(atoms) == ase.Atoms or issubclass(type(atoms), ase.Atoms):/" JKCS/lib/pyth*/site-packages/schnetpack/interfaces/ase_interface.py
  #THIS IS NEEDED FOR SCHNETPACK + PHYSNET TRICK
  $PIP install tad-dftd4 tad-mctc tad-multicharge tad-dftd3
fi

if [[ "$*" == *"-physnet"* ]]
then
  echo "======================"
  $PIP install torch==2.4.0
  echo "======================"
  $PIP install tensorboardX==2.6.2.2
  echo "======================"
  $PIP install numpy==1.25.2
  echo "======================"
  if [ ! -e ../JKML/src/PhysNet_DER ];
  then
    cd ../JKML/src
    git clone https://github.com/LIVazquezS/PhysNet_DER.git
    sed -i "s/np.float/float/g" PhysNet_DER/DataContainer.py
    sed -i "s/return tensor/segment_ids = segment_ids.to(torch.int64)\n    tensor = torch.zeros(*shape,device=device).scatter_add(0, segment_ids, data)\n    return tensor/" PhysNet_DER/layers/utils.py 
    sed -i "s/tensor = torch.zeros(/data = data.to(torch.float32)\n    segment_ids = segment_ids.to(torch.int64)\n    tensor = torch.zeros(/g" PhysNet_DER/layers/utils.py 
    #sed -i "s/np.Inf/np.inf/g" PhysNet_DER/run_train.py
    #sed -i "s/filename='train.log', //g" PhysNet_DER/run_train.py
    cp ../../TOOLS/LOADING/PN/run_train.py PhysNet_DER/
    cp ../../TOOLS/LOADING/PN/PNcalculator.py PhysNet_DER/
    cp ../../TOOLS/LOADING/PN/Neural_Net_evid.py PhysNet_DER/
    rm PhysNet_DER/calculator.py	
    cd -
  fi
fi
if [[ "$*" == *"-aimnet"* ]]
then
  #Installation of aimnet:
  #  conda create -n aimnet2 python=3.11
  #  conda activate aimnet
  #Install PyTorch with a proper CUDA version according to instructions at pytorch.org. E.g.
  #  conda install pytorch pytorch-cuda=12.4 -c pytorch -c nvidia
  #Install other dependencies.
  #  conda install -c conda-forge -c pytorch -c nvidia -f requirements.txt
  #Finally, install using setuptools.
  #  python setup.py install
  echo "======================"
  currdir=$PWD
  cd JKCS/lib64/py*/site-packages/
  git clone https://github.com/zubatyuk/aimnet2.git
  cd aimnet2
  $PYTHON setup.py build
  cd $currdir	
  echo "======================"
fi

#ArbAlign stuff:
cp ../TOOLS/SCRIPTS/modifiedArbAlign.py JKCS/lib/$(basename $PYTHON)/site-packages/ArbAlign.py

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
