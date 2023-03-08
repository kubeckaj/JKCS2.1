function [clus_ch] = solve_charge(clus)
% Takes in the name of a cluster and solves the charge
% ch -> charge: 1=neutral, 2=neg., 3=pos.

% Names of the generic ions
genneg={'neg','nion'};
genpos={'pos','pion'};
% Names of the ions in a cluster
neg={'B','R','RP'};
pos={'P'};
clus_ch=NaN;

if ismember(1,strcmpi(clus,genneg))
    clus_ch=2;
elseif ismember(1,strcmpi(clus,genpos))
    clus_ch=3;
else
    [molnames,~]=parse_cluster(clus);
    negidx=ismember(molnames,neg);
    posidx=ismember(molnames,pos);
    if ismember(1,negidx)
        clus_ch=2;
    elseif ismember(1,posidx)
        clus_ch=3;
    else
        clus_ch=1;
    end
end

end