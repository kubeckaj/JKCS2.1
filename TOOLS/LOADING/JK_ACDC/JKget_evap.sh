if [ ! -e RUN/get_evap.m ]
then
  echo "The RUN/get_evap.m file must exist. Run a short ACDC simulation first."
else
  cat RUN/get_evap.m | grep % | awk '{gsub(";", "", $3); print $5,$3}' | awk '{sum[$1] += $2} END {for (name in sum) print name, sum[name]}' | column -t | grep -v coefficients
fi

