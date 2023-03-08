function [] = ACDCsource(c1,c2)
  fileID = fopen('sources.txt','w');

  fprintf(fileID,"constant 1A %3.3e -1A1N-1A2N-1A3N-1A4N\n",c1);
  fprintf(fileID,"constant 1N %3.3e \n",c2);
  fprintf(fileID,"source neg 3.000000 \n");
  fprintf(fileID,"source pos 3.000000 \n");
  fprintf(fileID,"wall enhancement 1.000000e+00 \n");
  
  
  if fclose(fileID)==1
    error
  end
end

%fprintf(fileID,"constant 1A %3.3e -1A1MSA-1A1MSA1N-1A1MSA2N-1A2MSA-1A2MSA1N-1A3MSA-1A1N-1A2N-1A3N\n",c1);
%constant 1A 1.000000e+08 -1A1N-1A2N-1A3N-1A4N
%constant 1N 2.621047e+10
%source neg 3.000000e+00
%source pos 3.000000e+00
%wall enhancement 1.000000e+00 

