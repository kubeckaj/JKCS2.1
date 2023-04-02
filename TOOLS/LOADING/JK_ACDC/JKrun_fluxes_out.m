clc
clear all
if not(isfolder("RUN")) 
  error("Folder RUN missing. Have you run e.g. JKrun_single.m?")
  %exit
end
init
%%
PERCENTAGE = 0.02; %0.05 = 5 %

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
[flux_2d, coll_evap_2d, sources, flux_3d, flux_out, clust_flux]=dofluxes(C(end,:)*10^6,K,E,WL,CS,source,0);
%%
out=[];
out2=[];
maxi=size(flux_out,1);
total=0;
for i=1:maxi
 for j=1:maxi   
   f = flux_out(i,j,1);
   total=total+f;
 end
end
for i=1:maxi
 for j=1:maxi   
   f = flux_out(i,j,1);
   if f > PERCENTAGE*total %This means that you do not care about less than 5 % growth paths
     if i>j
       rem=i;
     else
       rem=j;
     end
     out=[out;[clust_acdc{i},clust_acdc{j},JKcombine_clusters(clust_acdc{i}+"",clust_acdc{j}+""),f/total*100,rem]];     
   end
 end
 f2 = sum(flux_out(i,:,1))+sum(flux_out(:,i,1));
 if f2 > PERCENTAGE*total
    out2=[out2;[clust_acdc{i},JKcombine_clusters(clust_acdc{i}+"",clust_acdc{i}+""),f2/total*100,i]];
 end
end
if size(out,1) > 1
%  out=sort(out,1,"descend");
   [~,order] = sort(str2double(out(:,4)),"descend");
   out = out(order,:);
end
if size(out2,1) > 1
%  out=sort(out,1,"descend");
   [~,order2] = sort(str2double(out2(:,3)),"descend");
   out2 = out2(order2,:);
end
fprintf("\nSingle fluxes:\n")
for i = 1:size(out,1)
  fprintf("%20s + %20s --> %20s (%.2f %s)\n",out(i,1),out(i,2),out(i,3),out(i,4),"%")
end
fprintf("\nNext collision outgrowing fluxes (distinguish monomers and clusters):\n")
for i = 1:size(out2,1)
  fprintf("%20s (%.2f %s)\n",out2(i,1),out2(i,3),"%")
end
fprintf("\nTotal flux is %f cm-3s-1\n\n",total/10^6)

%%
num_monomers=0;
max_cluster=0;
for i=1:size(clust_acdc,2)
   cluster=clust_acdc{i};
   if sum(str2double(regexp(cluster, '\d+', 'match')))==1
      num_monomers=num_monomers+1;
   end
   if sum(str2double(regexp(cluster, '\d+', 'match')))>max_cluster
      max_cluster=sum(str2double(regexp(cluster, '\d+', 'match')));
   end
end
test = transpose(str2double(unique(out(:,5))));
for idx0 = test
  idx = idx0;
  myflux = flux_2d;
  output=[clust_acdc(idx)];
  for i = 1:max_cluster-2
    max_flux_in_coll=0;
    this_cluster=clust_acdc{idx};
    this_cluster_num=sum(str2double(regexp(this_cluster, '\d+', 'match')));
    for j = 1:size(clust_acdc,2)
       test_cluster=clust_acdc{j};
       test_cluster_num=sum(str2double(regexp(test_cluster, '\d+', 'match')));
       if test_cluster_num<this_cluster_num && test_cluster_num>1
          if myflux(j,idx)>max_flux_in_coll
              idx1=j;
              max_flux_in_coll=myflux(j,idx);
          end
       end
    end
    %[~, idx1] = max(myflux(num_monomers+1:idx,idx));
    if idx1 ~= idx
      %idx=idx1+num_monomers;
      idx=idx1;
      %display(idx); 
      output=[output," -> ",clust_acdc(idx)];
    end
  end
  fprintf(join([join(output(end:-1:1)),"\n"]));
end

%%
%cellArr = output(end:-1:1);
%newCell = {'->'};
%n = numel(cellArr) + numel(newCell) * (numel(cellArr)-1);
%newArr = reshape([cellArr; repmat(newCell, [1 numel(cellArr)-1])], [1 n]);
%disp(newArr);

