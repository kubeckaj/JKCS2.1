if [ "$1" == "-help" ]
then
  echo "ARGUMENTS: <file.pkl> [<int=gfn version> <0/1=print modified input>]"
  echo "e.g. sh /home/kubeckaj/Applications/JKCS2.1/TOOLS/SCRIPTS/XTB1_prepare_parameter_file.sh 0 1"
  exit
fi

input=$1
version=$2
modify=$3
if [ -z "$input" ]; then echo Missing input; exit; fi 
if [ ! -e "$input" ]; then echo File does not exist; exit; fi 
if [ -z "$version" ]; then version=1; fi
if [ -z "$modify" ]; then modify=0; fi
eval `grep PATH_XTB ~/.JKCSusersetup.txt | head -n 1`;
if [ ! -e $PATH_XTB/share/xtb/param_gfn$version-xtb.txt ]; then echo Parameter file not found; exit; fi 
cp $PATH_XTB/share/xtb/param_gfn$version-xtb.txt .
search=`JKQC $input -atoms`
lines=""
lines+=`grep -ne '^\$globpar' param_gfn$version-xtb.txt | sed "s/:/ /" | awk '{print $1}'`
lines+=" "
lines+=`grep -ne '^\$pairpar' param_gfn$version-xtb.txt | sed "s/:/ /" | awk '{print $1}'`
for i in $search
do
  lines+=" "`grep -ne "^\\\$Z=[[:space:]]*${i}[[:space:]]" param_gfn$version-xtb.txt | sed "s/:/ /" | awk '{print $1}'`
done
#echo $lines
any_overlap() {
    for e1 in $search
    do
        for e2 in 9 17 35 53 85 117
        do
            if [[ "$e1" = "$e2" ]]
            then
                echo "1"
                exit
            fi
        done
    done
    echo "0"
}
both_present() {
    for e1 in $*
    do
        t=0
        for e2 in $search
        do
            if [[ "$e1" = "$e2" ]]
            then
                t=1
            fi
        done
        if [ $t -eq 0 ]; then echo "0"; exit; fi
    done
    echo "1"
}
Qhalogens=`any_overlap`
#echo $Qhalogens

vars=1
search_num=0
initial=""
for i in $lines
do
  search_num=`echo $search_num+1 | bc`
  test=0
  linenum=$i
  while [ $test -eq 0 ]
  do
    linenum=`echo $linenum+1 | bc`
    line=`head -n ${linenum} param_gfn$version-xtb.txt | tail -n 1`
    #echo $line
    if [ "$line" == '$end' ]; then test=1; continue; fi
    columns=`echo $line | wc -w`
    if [ $search_num -eq 2 ]
    then
      newline=`echo $line | awk '{print $1 " " $2}'`
      Qbothpresent=`both_present $newline`
      if [ $Qbothpresent -eq 0 ]; then continue; fi
      startcol=3
    else
      newline=`echo $line | awk '{print $1}'`
      startcol=2
    fi
    #XTB1: perhaps useless ipeashift|s9|cnd2
    #XTB2: perhaps useless kdiff|ipeashift|gam3d2
    #if [[ "${newline}" == @(ipeashift|s9|cnd2) ]]; then continue; fi
    if [ $Qhalogens -eq 0 ] && [[ "${newline}" == @(xbdamp|xbrad) ]]; then continue; fi
    for j in `seq $startcol $columns`
    do
      initial+=" `echo $line | awk -v var=$j '{print $var}'`"
      newline+=" AAA$vars "
      vars=`echo $vars+1 | bc`
    done
    sed -i "${linenum}s/.*/$newline/" param_gfn$version-xtb.txt
  done
done
echo $initial > initial.txt
if [ $modify -eq 1 ]
then
  if [ -e mod_initial.txt ]; then rm mod_initial.txt; fi
  lines=`cat initial.txt | wc -w`
  for i in `seq 1 $lines`
  do
    number=`cat initial.txt | awk -v el=$i '{print $el}'`
    echo "AAA$i $number" >> mod_initial.txt
    modulo=$((i % 5))
    if [ $modulo -eq 0 ]
    then
      echo "" >> mod_initial.txt
    fi
  done
fi 
