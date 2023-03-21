clc
clear all
cd ../../
init

cd 'RUN/'
load("variables.mat")
fprintf("JK: Loading ACDC output.\n");

fp=fopen('driver_acdc.m','r');
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
get_evap; get_coll; get_wl; get_cs;
source=zeros([size(C(end,:),2),1]);
[flux_2d, coll_evap_2d, sources, flux_3d, flux_out, clust_flux]=dofluxes(C(end,:)*10^6,K,E,WL,CS,source,0,0,0);

out=[];
max=size(flux_out,1);
total=0;
for i=1:max
 for j=1:max   
   f = flux_out(i,j,1);
   total=total+f;
 end
end
for i=1:max
 for j=1:max   
   f = flux_out(i,j,1);
   if f > 0.05*total %This means that you do not care about less than 5 % growth paths
       out=[out;[clust_acdc{i}," + ",clust_acdc{j}," --> ",JKcombine_clusters(clust_acdc{i}+"",clust_acdc{j}+""),num2str(round(1000*f/total)/10)+"%"]];
   end
 end
end
out
"Total flux is "+num2str(total/10^6)+" cm-3s-1"
