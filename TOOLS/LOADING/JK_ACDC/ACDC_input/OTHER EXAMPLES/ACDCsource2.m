function [] = ACDCsource(c1)
%ACDCSOURCE Summary of this function goes here
%   Detailed explanation goes here
  fileID = fopen('sources.txt','w');
  fprintf(fileID,"constant 1A %3.3e \n",c1);
  if fclose(fileID)==1
    error
  end
end

