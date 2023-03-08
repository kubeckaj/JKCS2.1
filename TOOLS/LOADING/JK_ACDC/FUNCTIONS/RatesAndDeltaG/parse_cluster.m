function [molnames,nmols]=parse_cluster(namelabel)
% takes in a namelabel of a cluster and returns vectors containing the
% molecule names and corresponding molecule numbers

molnames=regexp(namelabel,'[A-Za-z]+','match');
nmols=regexp(namelabel,'[0-9]+','match');
nmols=cellfun(@str2num,nmols);

end