function [largest,smallest,maxmols,minmols] = max_mols(clus1,clus2)
% Takes in names of clusters, calculates which one has the largest
% number of molecules, and returns the name of the largest cluster

clusters={clus1 clus2};
totmols=zeros(1,2);

for i=1:2
    totmols(i)=calc_mols(clusters{i});
end

[maxmols,nmax]=max(totmols);
if (maxmols==0)
    error(['Maximum number of total molecules is zero? ',clusters])
end
largest=clusters{nmax};

if nmax==1
    nmin=2;
else
    nmin=1;
end

smallest=clusters{nmin};
minmols=totmols(nmin);

end