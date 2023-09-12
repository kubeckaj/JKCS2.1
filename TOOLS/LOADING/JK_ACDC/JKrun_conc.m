clc
clear all
if not(isfolder("RUN"))
  error("Folder RUN missing. Have you run e.g. JKrun_single.m?")
  %exit
end
init
%%

%%%%%%%%%%%%%%
%%% LOADING THE RESULTS
%%%%%%%%%%%%%%

load("variables.mat")
fprintf("JK: Loading ACDC output.\n");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% READ ACDC CLUSTERS
fp=fopen('RUN/driver_acdc.m','r');
line=0;
while line ~= -1
    line=[fgetl(fp),' '];
    if ~isempty(strfind(line, 'clust = {'))
        line=strrep(line,'clust','clust_acdc');
        eval(line);
        break
    end
end
fclose(fp);
% -> clust_acdc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:size(C,2)
  fprintf('{\"%s\",',clust_acdc{i})
  fprintf("%f},",C(end,i))
end
fprintf("\n")
