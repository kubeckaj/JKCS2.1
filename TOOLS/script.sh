savefile=collectionXTB.txt
resultfile=resultsXTB.dat
movie=movieXTB.xyz

GRID=50  # grid points on each axis
Qprint=1 # print comments?
Qs=3     # sort column

###########################
#############################
S=`grep -c "info" $0`
S=`echo $S-2 | bc`
S0=0
function info {
  if [ $Qprint -eq 1 ]
  then
    S0=`echo $S0+1 |bc`
    echo $S0/$S:$1 
  fi
}

# 3 data columns?
info "Checking amount of columns"
c4=`head -n 1 $savefile | awk '{print $4}'`
if [ -z "$c4" ]
then 
  t4=0
else
  t4=1
fi
# remove if for whatever reason exist these files
info "Removing old files .help1_$resultfile & .help2_$resultfile (if they exist)"
if [ -e .help1_$resultfile ]; then rm .help1_$resultfile; fi
if [ -e .help2_$resultfile ]; then rm .help2_$resultfile; fi
# MIN MAX VALUES
info "Searchin for MIN and MAX value in each column"
cat $savefile | awk '{print $2}' > .help0_$resultfile
value=`head -n 1 .help0_$resultfile`
max1=`awk -v max=$value '$1 > max { max = $1 } END { print max }' .help0_$resultfile` 
min1=`awk -v min=$value '$1 < min { min = $1 } END { print min }' .help0_$resultfile` 
cat $savefile | awk '{print $3}' > .help0_$resultfile
value=`head -n 1 .help0_$resultfile`
max2=`awk -v max=$value '$1 > max { max = $1 } END { print max }' .help0_$resultfile` 
min2=`awk -v min=$value '$1 < min { min = $1 } END { print min }' .help0_$resultfile` 
if [ $t4 -eq 1 ]
then
  cat $savefile | awk '{print $4}' > .help0_$resultfile
  value=`head -n 1 .help0_$resultfile`
  max3=`awk -v max=$value '$1 > max { max = $1 } END { print max }' .help0_$resultfile` 
  min3=`awk -v min=$value '$1 < min { min = $1 } END { print min }' .help0_$resultfile` 
else
  min3="";max3="";
fi 
#echo $min1 $max1 $min2 $max2 $min3 $max3
info "Sorting ... this might take some time"
if [ $t4 -eq 1 ]
then
  sed '/^[ \t]*$/d' "$savefile" | awk -v var=$GRID -v var1a=$min1 -v var1b=$max1 -v var2a=$min2 -v var2b=$max2 -v var3a=$min3 -v var3b=$max3 '{printf("%s %.0f %.0f %.0f\n"),$1,($2-var1a)/(var1b-var1a)*var,($3-var2a)/(var2b-var2a)*var,($4-var3a)/(var3b-var3a)*var}' | sort -u -n -k 2 -k 3 -k 4 > .help1_$resultfile
else
  sed '/^[ \t]*$/d' "$savefile" | awk -v var=$GRID -v var1a=$min1 -v var1b=$max1 -v var2a=$min2 -v var2b=$max2 '{printf("%s %.0f %.0f\n"),$1,($2-var1a)/(var1b-var1a)*var,($3-var2a)/(var2b-var2a)*var}' | sort -u -n -k 2 -k 3 > .help1_$resultfile
fi

# GREPING WHOLE LINES
info "Searching selected lines for $resultfile ... this might take also some time"
lines=`wc -l .help1_$resultfile | awk '{print $1}'`
for i in `seq 1 $lines`
do
  if [ $Qprint -eq 1 ]
  then
    echo -en "\r line $i/$lines"
  fi
  text=`head -n $i .help1_$resultfile | tail -n 1 | awk '{print $1}'`
  grep $text $savefile >> .help2_$resultfile
done
if [ $Qprint -eq 1 ]
then
  echo -en "\n"
fi
#echo "UNIQUENESS: $resultfile - $lines files ( u1= $u1 ,u2= $u2 ,u3= $u3 ,sort= $Qs )"
#echo "UNIQUENESS: $resultfile - $lines files ( u1= $u1 ,u2= $u2 ,u3= $u3 ,sort= $Qs )" >> FILTER.txt

# SORTING
info "Sorting + formation of $resultfile."
cat .help2_$resultfile | sort -nk $Qs > $resultfile
#RED='\033[0;31m';NC='\033[0m';printf "I ${RED}love${NC} Stack Overflow\n" > o

# MOVIE
info "Creating file $movie."
if [ -e $movie ]; then rm $movie; fi
o=`cat $resultfile | awk '{print $1}'`;
for i in $o;
do
  cat ${i%.*}.xyz >> $movie
done

# REMOVE
info "Removing unnecessary files."
rm .help1_$resultfile .help2_$resultfile .help0_$resultfile

info "Done"

