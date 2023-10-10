#test=0;ADD=0;
#while [ $test -eq 0 ]
#do
#  random_folder=XTB_${ADD}
#  if [ -d $random_folder ]; then ADD="_${RANDOM}"; else test=1;fi
#done
#mkdir $random_folder
#cd $random_folder
cp ../param_gfn1-xtb.txt .

NCORES=$1

test=$(cat parameter.txt)
TEST=($test)

NUMBA=$(grep "AAA" param_gfn1-xtb.txt -o | wc -l)
for ((i = 1; i <= NUMBA; i++)); do
    loop_counter=$((i - 1))
    sed -i "s/AAA$i /${TEST[$loop_counter]} /g" param_gfn1-xtb.txt
done

JKQC ../input.pkl -xyz

LW=$(ls *.xyz | xargs)

# initialize a semaphore with a given number of tokens
open_sem(){
    mkfifo pipe-$$
    exec 3<>pipe-$$
    rm pipe-$$
    local i=$1
    for((;i>0;i--)); do
        printf %s 000 >&3
    done
}

# run the given command asynchronously and pop/push tokens
run_with_lock(){
    local x
    # this read waits until there is something to read
    read -u 3 -n 3 x && ((0==x)) || exit $x
    (
     ( "$@"; )
    # push the return code of the command to the semaphore
    printf '%.3d' $? >&3
    )&
}

# Setting XTB 
export XTBHOME=$PWD
export XTBPATH=$PWD
source ~/.JKCSusersetup.txt AAAlocAAA
$MODULE_XTB 2>/dev/null
ulimit -s unlimited
export OMP_STACKSIZE=1G
export OMP_NUM_THREADS=1
module_xtb_loaded="loaded"
open_sem $NCORES
for i in $LW
do
  run_with_lock program_XTB $i --gfn 1
  #program_XTB $i --gfn 1 &
  # allow to execute up to $N jobs in parallel
  #if [[ $(jobs -r -p | wc -l) -ge $NCORES ]]; then
    # now there are $N jobs already running, so wait here for any job
    # to be finished so there is a place to start next one.
  #  wait -n
  #fi
done
wait

#array_i=($LW)

#AS=${#array_i[@]}  # Size of the array
#AB=$NCORES   # Size of each split

## Calculate the number of splits required
#num_splits=$((AS / AB + (AS % AB > 0)))

# Loop through the splits
#for ((i = 0; i < num_splits; i++)); do
#    start=$((i * AB))
#    end=$((start + AB - 1))
#    # Adjust the end index for the last split
#    if [ $end -ge $AS ]; then
#        end=$((AS - 1))
#    fi
#    # Extract the split from the array
#    split=("${array_i[@]:$start:$((end - start + 1))}")
#
#    # Print the split
#    # echo "Split $i: ${split[@]}"
#    for x in ${split[@]}
#    do
#	program_XTB $x --gfn 1 &
#    done
#    wait
#done

#Get results
JKQC *.log -add HIGH ../B_El.txt -ct -el -extra HIGH -natoms -formation -noex -unit > FORM.txt
#ddt_MAE=`cat FORM.txt | awk '{printf("%.9f\n",($2-$3)^2/$4)}' | xargs | sed 's/ /+1.*/g' | bc -l | awk '{printf("%.9f\n",$1^0.5)}'`
#t_MAE=`cat FORM.txt | awk 'BEGIN{n=0;sum=0}{sum+=(($2-$3)/$4)^2;n++}END{print (sum/n)^0.5}'`
#t_MAE=`cat FORM.txt | awk -F'[-,]' 'function abs(v) {return v < 0 ? -v : v}BEGIN{n=0;sum=0}{sum+=(abs($2-$3)/$4;n++}END{print (sum/n)^0.5}'`
t_MAE=`cat FORM.txt | awk 'BEGIN{n=0;sum=0}{sum+=sqrt(($2-$3)^2)/$4;n++}END{print (sum/n)}'`
f_MAE=`grep "GRADIENT NORM " *.log | awk '{print $5}' | awk 'BEGIN{n=0;sum=0}{sum+=($1);n++}END{print (sum/n)}'`
tf_MAE=`echo $t_MAE+$f_MAE | bc -l`
#echo "$t_MAE + $f_MAE = $tf_MAE"
t_Failed=`grep -c "nan" FORM.txt`
t_Failed2=`grep "GRADIENT NORM " *.log | awk '{print $5}' | grep -c "NaN"`
rm *.log *.xyz

if [ "$t_Failed" -gt 0 ] || [ "$t_Failed2" -gt 0 ]
then
	echo "FAILED" > result
        exit
else
	echo $tf_MAE > result
        printf "." >> ../output
        #echo "$test"
        #echo MAE = $t_MAE kcal/mol >> ../output
fi

#cd ..
#rm -r $random_folder
