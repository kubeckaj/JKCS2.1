function [lcluster]=check_cluster(namelabel)
% Takes in a namelabel and checks if it is a cluster,
% i.e. of the form [numbers][letters][numbers][letters]...
% Returns a logical lcluster, which is 1 if the thing is a cluster, and 0 if not

lcluster=0;

matchstr=regexp(namelabel,'([0-9]+[A-Za-z]+)+','match');
if strcmp(matchstr,namelabel)
    lcluster=1;
end

end