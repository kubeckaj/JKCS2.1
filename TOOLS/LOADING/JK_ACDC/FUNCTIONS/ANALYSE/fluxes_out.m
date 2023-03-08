
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
"Total flux is "+num2str(total)+" cm-3s-1"


            