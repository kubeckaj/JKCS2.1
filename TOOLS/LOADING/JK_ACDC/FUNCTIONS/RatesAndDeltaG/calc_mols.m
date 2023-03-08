function [totmols]=calc_mols(namelabel)
% Takes in a name of a cluster and calculates the total number of molecules
% in it (a proton/removed proton is NOT considered to be a separate molecule)

posion='P'; % this is assumed to be a PROTON
negion='RP'; % this is assumed to be a missing proton

[molnames,nmols]=parse_cluster(namelabel);
totmols=sum(nmols);
% don't take the proton to the sum
npos=strmatch(posion,molnames,'exact');
nneg=strmatch(negion,molnames,'exact');

if ~isempty(npos)
    totmols=totmols-nmols(npos);
end

if ~isempty(nneg)
    totmols=totmols-nmols(nneg);
end

end