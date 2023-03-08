out=[];
for i=[1 5 17]
   newcluster=JKcombine_clusters("1A",clust_acdc{i}+"");
   test=0;
   for ii=1:75
    if compare_clusters(newcluster,clust_acdc{ii})
      test=ii;
      break
    end
   end
   if test==0
    test=81
   end
   if flux_2d(1,test) > 0.001
     newcluster=JKcombine_clusters("1A",clust_acdc{i}+"");
     out=[out;["1A",newcluster,test,flux_2d(1,test)]];
   end
end

for repeat=1:8
out2=[];
s=size(out,2);
for j=1:size(out,1)
 if str2num(out(j,s-1)) > 0 && str2num(out(j,s-1)) < 81
  for i=[1 5 17]
   newcluster=JKcombine_clusters(clust_acdc{str2num(out(j,s-1))},clust_acdc{i}+"");   
   test=0;
   for ii=1:75
    if compare_clusters(newcluster,clust_acdc{ii})
     test=ii;
     break
    end
   end 
   if test==0
    test=81;
   end
   if flux_2d(str2num(out(j,s-1)),ii) > 0.1
     out2=[out2;[out(j,1:s),newcluster,test,flux_2d(str2num(out(j,s-1)),test)]];
   %else
   %  out2=[out2;[out(j,1:s),"","",""]];
   end
  end
 else
  out2=[out2;[out(j,1:s),0,0,0]];
 end  
end
out=out2
end
out2