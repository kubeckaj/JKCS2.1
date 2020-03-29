#!/bin/bash
export OMP_NUM_THREADS=1 
module load intel/19.0.4 hpcx-mpi/2.4.0 intel-mkl/2019.0.4 2>/dev/null
export LD_LIBRARY_PATH=/users/kubeckaj/ORCA/orca_4_2_1_linux_x86-64_shared_openmpi314/:/appl/spack/install-tree/gcc-9.1.0/hpcx-mpi-2.4.0-dnpuei/lib:/appl/opt/cluster_studio_xe2019/compilers_and_libraries_2019.4.243/linux/tbb/lib/intel64_lin/gcc4.7:/appl/opt/cluster_studio_xe2019/compilers_and_libraries_2019.4.243/linux/compiler/lib/intel64_lin:/appl/opt/cluster_studio_xe2019/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64_lin:/appl/spack/install-tree/gcc-4.8.5/gcc-9.1.0-vpjht2/lib64:/appl/spack/install-tree/gcc-4.8.5/gcc-9.1.0-vpjht2/lib

#CREATING WORKING DIRECTORY
if [ ! -d .//TMP ]; then mkdir .//TMP; fi
ADD=""
test=0
while [ $test -eq 0 ]
do
  CALC_NAME=.//TMP/ORCA${SLURM_JOBID}${ADD}
  if [ -d $CALC_NAME ]; then ADD="_${RANDOM}"
  else test=1;fi
done

#dirs
mkdir $CALC_NAME

if [ -e .test.xyz ]; then cp .test.xyz $CALC_NAME/; fi
cp .test.inp $CALC_NAME/
cd $CALC_NAME

/users/kubeckaj/ORCA/orca_4_2_1_linux_x86-64_shared_openmpi314//orca .test.inp > .test.out 2> .test.out

#COPYING RESULTS BACK
if [ -e /projappl/hvehkama/kubeckaj/Apps/JKCS2.1/.test.out ]
then
  mv /projappl/hvehkama/kubeckaj/Apps/JKCS2.1/.test.out /projappl/hvehkama/kubeckaj/Apps/JKCS2.1/.testO.out
fi
cp .test.out /projappl/hvehkama/kubeckaj/Apps/JKCS2.1/
cp *.xyz /projappl/hvehkama/kubeckaj/Apps/JKCS2.1/ 2>/dev/null
cp *.err /projappl/hvehkama/kubeckaj/Apps/JKCS2.1/ 2>/dev/null
cd /projappl/hvehkama/kubeckaj/Apps/JKCS2.1

#CLEANING
rm -rf $CALC_NAME
