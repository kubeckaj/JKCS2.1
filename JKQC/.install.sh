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
$PYTHON -m venv JKCS
#$PYTHON -m venv --no-site-packages JKCS
#$PYTHON -m venv --system-site-packages JKCS

echo "- activating environment"
source JKCS/bin/activate

#ADD="--force-reinstall --upgrade"
PIP="$PYTHON -m pip --no-cache-dir --disable-pip-version-check " #--cache-dir=$PWD/JKCS/"

#Remove some non functional libraries
rm -r JKCS/lib/python*/site-packages/~* 2>/dev/null

####################################
echo "======================"
$PIP install --upgrade pip
#$PIP install pip==22.0.4
####################################

#echo "======================"
#$PIP install ipython
#echo "======================"
#$PIP install pyarrow==15.0.0 #this bullshit is required by new pandas
#echo "======================"
#$PYTHON -m ipykernel install --user --name=jkcs #update your jupyter

####################################
echo "======================"
$PIP install numpy==1.25.2
echo "======================"
$PIP install scipy==1.9.3   #I need this version for ArbAlign 
echo "======================"
$PIP install joblib==1.3.2  
echo "======================"
$PIP install pathlib==1.0.1 #Perhaps this one is not necessary
echo "======================"
$PIP install numexpr==2.8.4 #2.7.0
echo "======================"
$PIP install psutil
echo "======================"
$PIP install xlsxwriter     #important only for IMoS output but cheap to install
echo "======================"
$PIP install rdkit     #important only for IMoS output but cheap to install
####################################
echo "======================"
#module load gcc 2>/dev/null; $PIP install git+https://github.com/src-d/lapjv
#$PIP install lapjv==1.3.1  #this bitch has some issues to see numpy
$PIP install lapjv
####################################
echo "======================"
#$PIP install pandas==1.3.4  #2.2.0  #for the newer one, I have to get rid of append
$PIP install pandas==1.5.3  #2.2.0  #for the newer one, I have to get rid of append
####################################
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
####################################
#ArbAlign stuff:
cp ../TOOLS/SCRIPTS/modifiedArbAlign.py JKCS/lib/$(basename $PYTHON)/site-packages/ArbAlign.py
####################################

if [[ "$*" == *"-descriptors"* ]]
then
  echo "======================"
  $PIP install scikit-learn
  echo "======================"
  $PIP install cffi
  echo "======================"
  $PIP install dscribe
fi

####################################

if [[ "$*" == *"-xtb"* ]] || [[ "$*" == *"-calcs"* ]]
then
  #############
  echo "======================"
  export CC=gcc
  module load gcc 2>/dev/null
  $PIP install tblite
  cp ../TOOLS/LOADING/JKase.py JKCS/lib/python*/site-packages/tblite/ase.py
  sed -i "s/print(_nentries)//" JKCS/lib64/pyth*/site-packages/tblite/library.py
  sed -i "s/print(_dict_py)//" JKCS/lib64/pyth*/site-packages/tblite/library.py
  #############
  echo "======================"
  $PIP install XTB
  #############
  echo "======================"
  $PIP install py-cpuinfo urllib3 typing-extensions toolz platformdirs ndindex msgpack idna charset-normalizer certifi requests 
  #############
  echo "======================"
  $PIP install --no-deps ORCA blosc2 tables
  #############
fi

####################################

if [[ "$*" == *"-qml "* ]]
then
  echo "======================"
  $PIP install scikit-learn
  if ! command -v gcc 2>&1 >/dev/null
  then
    echo "======================"
    echo "I cannot see gcc so I am trying to: module load gcc"
    module load gcc
  fi
  ###sed -i 's/ -lpthread/-lpthread/' setup.py   # remove the stray space
  if ! command -v lapack 2>&1 >/dev/null
  then
    echo "======================"
    echo "I cannot see lapack so I am trying to: module load lapack"
    module load lapack
  fi
  echo "======================"
  $PIP install --no-deps qmllib
fi

####################################
## SCHNETPACK ##
# it has its own environmnet

if [[ "$*" == *"-nn"* ]]
then
  currdir=$PWD
  cd JKCS
  mkdir SCHNETPACK 2>/dev/null
  cd SCHNETPACK
  mkdir lib 2>/dev/null; mkdir lib/$PYTHON 2>/dev/null; mkdir lib/$PYTHON/site-packages 2>/dev/null
  export PYTHONPATH=$PWD/lib/$PYTHON/site-packages/
  echo "======================"
  $PIP install --prefix ./ sympy==1.13.2 torch==2.4.0 torchmetrics==1.0.1 torch_ema==0.3 pytorch-lightning==2.0.6 lightning hydra-core tensorboard==2.17.1 tensorboardX==2.6.2.2 distlib pygments platformdirs pathspec nodeenv mypy-extensions mdurl identify fasteners dirsync colorlog click cfgv virtualenv scipy markdown-it-py h5py black rich pre-commit hydra-colorlog matscipy scikit-learn
  echo "======================"
  $PIP install --prefix ./ --no-deps schnetpack==2.0.4
  #AGOX was not able to use schnet calculator, this will resolve it:
  sed -i "s/elif type(atoms) == ase.Atoms:/elif type(atoms) == ase.Atoms or issubclass(type(atoms), ase.Atoms):/" lib/$PYTHON/site-packages/schnetpack/interfaces/ase_interface.py
  #THIS IS NEEDED FOR SCHNETPACK + PHYSNET TRICK
  #$PIP install --target JKCS/lib/SCHNETPACK tad-dftd4 tad-mctc tad-multicharge tad-dftd3
  cd $currdir
  export PYTHONPATH=""
fi

####################################
## AIMNET ##
# it has its own environmnet

if [[ "$*" == *"-aimnet"* ]]
then
  currdir=$PWD
  cd JKCS 
  mkdir AIMNET 2>/dev/null
  cd AIMNET
  mkdir lib 2>/dev/null; mkdir lib/$PYTHON 2>/dev/null; mkdir lib/$PYTHON/site-packages 2>/dev/null
  export PYTHONPATH=$PWD/lib/$PYTHON/site-packages/
  ######
  $PIP install --prefix ./ torch==2.7.0
  echo "======================"
  $PIP install --no-deps --prefix ./ pytorch-ignite wandb==0.19.11 torchvision torchaudio h5py wheel scikit-learn threadpoolctl pyyaml click omegaconf tqdm requests urllib3 idna certifi appdirs protobuf==5.29.2 charset_normalizer platformdirs pydantic==2.10.6 pydantic_core==2.27.2 typing_inspection annotated_types docker-pycreds sentry-sdk numba llvmlite
  echo "======================"
  $PIP install --prefix ./ antlr4-python3-runtime==4.9.3 --upgrade --use-pep517 --no-build-isolation
  echo "======================"
  $PIP install --no-deps --prefix ./ torch-cluster -f https://data.pyg.org/whl/torch-2.7.0+cpu.html
  $PIP install --no-deps --prefix ./ torch-cluster -f https://data.pyg.org/whl/torch-2.7.0+cu118.html
  $PIP install --no-deps --prefix ./ torch-cluster -f https://data.pyg.org/whl/torch-2.7.0+cu126.html
  $PIP install --no-deps --prefix ./ torch-cluster -f https://data.pyg.org/whl/torch-2.7.0+cu128.html
  echo "======================"
  if [ -e ./lib/$PYTHON/site-packages/aimnet ]
  then
    echo "Requirement already satisfied: aimnet in ../lib/python3.11/site-packages"
  else
    $PIP install --no-deps --prefix ./ git+https://github.com/isayevlab/aimnetcentral.git
    sed -i 's/if loss_fn\.weights is not None:/if hasattr(loss_fn, "weights") and loss_fn.weights is not None:/g' lib/$PYTHON/site-packages/aimnet/train/utils.py
    if [ ! -e ../bin/aimnet ]
    then
      ln -s $PWD/bin/aimnet ../bin/aimnet
    fi
    awk '
    /if cfg.wandb.watch_model:/ {
        print "    best_val = {\"loss\": float(\"inf\")}"
        print "    "
        print "    @validator.on(Events.EPOCH_COMPLETED)"
        print "    def update_best_val_loss(engine):"
        print "        current = engine.state.metrics.get(\"loss\")"
        print "        if current is not None and current < best_val[\"loss\"]:"
        print "            best_val[\"loss\"] = current"
        print "            wandb.run.summary[\"best_val/loss\"] = current"
        print "    "
    }
    { print }
    ' lib/$PYTHON/site-packages/aimnet/train/utils.py > tmp
    mv tmp lib/$PYTHON/site-packages/aimnet/train/utils.py
  fi
  ######################
  cd $currdir
  export PYTHONPATH=""
fi

####################################
####################################
####################################

if [[ "$*" == *"-physnet"* ]]
then
  echo " AJAJAJAJ --------------- NOT TESTED FOR AGES"
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

####################################

if [[ "$*" == *"-mbdf"* ]]
then
  echo " AJAJAJAJ --------------- NOT TESTED FOR AGES"
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

if [[ "$*" == *"-mbdf"* ]]
then
  echo " AJAJAJAJ --------------- NOT TESTED FOR AGES"
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

####################################

if [[ "$*" == *"-qml-lightning"* ]]
then
  echo " AJAJAJAJ --------------- NOT TESTED FOR AGES"
  echo "======================"
  currdir=$PWD
  cd JKCS/lib64/py*/site-packages/
  git clone https://github.com/nickjbrowning/qml-lightning.git
  #LDFLAGS=-L/usr/local/cuda-11.4
  $PYTHON setup.py build
  cd $currdir
fi

####################################
####################################
####################################

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

echo "=================================================================="
echo " - saving JKpython alias to ~/.bashrc"
PyPATH=$(which $PYTHON)

command="alias JKpython='$MODULE_PYTHON source $PWD/JKCS/bin/activate; alias python=$PyPATH'"
writetobashrc "$command"

