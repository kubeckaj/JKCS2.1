function [] = ACDCinit(Temp,rh)
%ACDCINIT Summary of this function goes here
%   Detailed explanation goes here

  command="perl ../ACDC/create_acdc_2016_09_30.pl";
  command=command+" --temperature "+num2str(Temp);
  command=command+" --e ../ACDC_input/dG.txt";
  command=command+" --i ../ACDC_input/input2.inp";
  command=command+" --rh "+num2str(rh);
  %command=command+" --dip ../ACDC_input/dipoles.txt";
  %command=command+" --no_generic_ions";
  ok=system(command);
  if ok ~= 0,cd ..;error('Trouble with perl script!');end

end

