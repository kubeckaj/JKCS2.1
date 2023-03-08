function [combined]=combine_clusters(clus1,clus2)
% Combines two clusters and returns the name label of the result cluster
% Input: name labels of 2 clusters, e.g. 2A1B1N
% NB!!! name labels of generic negative and positive ions are assumed to be neg and pos, or nion and pion
% NB!!! The order of the molecules is not standard, remember to check it by e.g. using the function compare_clusters
% NB!!! If a sulfuric acid -containing cluster gets negatively ionized, the default is to replace one acid (A)
% with a bisulfate (B); however, the code works also for recombinations of a cluster containing a missing proton (RP) (the missing proton is removed)

% Names of acids and bases, used when checking if clusters can get ionized/recombined
% Remember to update!!!
acid='A';
bases={'N' 'D' 'M' 'T'};
negacid='B'; % this is the ionized acid

posion='P'; % this is assumed to be a proton
negion='RP'; % this is assumed to be a missing proton

clusters={clus1 clus2};

% Logicals: 1 when true, 0 when not
% check if the cluster is a generic negative ion
lneg=zeros(1,2);
% generic positive ion
lpos=zeros(1,2);
% cluster containing a negative ion
lnegclus=zeros(1,2);
% cluster containing a positive ion
lposclus=zeros(1,2);
% cluster containing acid
lhasacid=zeros(1,2);
% cluster containing base
lhasbase=zeros(1,2);

for i=1:2
    
    if(strcmpi(clusters{i},'neg') || strcmpi(clusters{i},'nion'))
       lneg(i)=1;
    end
    
    if(strcmpi(clusters{i},'pos') || strcmpi(clusters{i},'pion'))
       lpos(i)=1;
    end
    
end

% create vectors containing the names of different molecules and the corresponding numbers of molecules for both clusters
nmols=cell(1,2);
molnames=cell(1,2);

for i=1:2
    if (~lneg(i) && ~lpos(i))
        [molnames{i},nmols{i}]=parse_cluster(clusters{i});
        
        % no names or numbers
        if (isempty(molnames{i}) || isempty(nmols{i}))
            combined='undef';
            return;
        end
        
        % find out if this is negative or positive, or has acid or base
        if ismember(negion,molnames{i}) || ismember(negacid,molnames{i})
            lnegclus(i)=1;
        end
        
        if ismember(posion,molnames{i})
            lposclus(i)=1;
        end
        
        if ismember(acid,molnames{i})
            lhasacid(i)=1;
        end
        
        for j=1:length(bases)
            if ismember(bases{j},molnames{i})
                lhasbase(i)=1;
                break;
            end
        end
    end
end

% check all possible collisions, ionizations and recombinations

% if one of the clusters is a generic negative ion
if (lneg(1) || lneg(2))

    % check which one is the ion
    if lneg(1)
        nclus=2;
    else
        nclus=1;
    end
    
    % check all the possibilities, depending on what kind of cluster the other party is
    % these don't produce a result cluster
    if (lneg(nclus) || lnegclus(nclus) || lpos(nclus))
        
        combined='undef';
        return;
    
    % if this is a positive cluster, it gets deprotonated
    elseif (lposclus(nclus))
        
        % find the proton and remove it
        nposion=strcmp(molnames{nclus},posion);
        nmols{nclus}(nposion)=nmols{nclus}(nposion)-1;
    
        combined=create_label(nmols{nclus},molnames{nclus});
        return;
        
    else
        % if this is a neutral cluster, it gets ionized if it has acid
        if (lhasacid(nclus))
            % remove one acid
            nacid=strcmp(molnames{nclus},acid);
            nmols{nclus}(nacid)=nmols{nclus}(nacid)-1;

            % add a corresponding negative ion
            molnames{nclus}{length(molnames{nclus})+1}=negacid;
            nmols{nclus}(length(nmols{nclus})+1)=1;
            
            combined=create_label(nmols{nclus},molnames{nclus});
            return;
            
        else
            
            combined='undef';
            return;

        end
    
    end

% if one of the clusters is a generic positive ion
elseif (lpos(1) || lpos(2))
    
    if lpos(1)
        nclus=2;
    else
        nclus=1;
    end    
    
    if (lpos(nclus) || lposclus(nclus))
        
        combined='undef';
        return;
    
    % if this is a negative cluster, it gets neutralized
    elseif (lnegclus(nclus))
        
        if ismember(negion,molnames{nclus})
            
            % if there is a missing proton, just remove it
            nnegion=strcmp(molnames{nclus},negion);
            nmols{nclus}(nnegion)=nmols{nclus}(nnegion)-1;
            
        else
            
            % add one acid
            nacid=strcmp(molnames{nclus},acid);
            if ismember(1,nacid)
                nmols{nclus}(nacid)=nmols{nclus}(nacid)+1;
            else
                molnames{nclus}{length(molnames{nclus})+1}=acid;
                nmols{nclus}(length(nmols{nclus})+1)=1;
            end
        
            % remove the negative ionic acid
            nnegacid=strcmp(molnames{nclus},negacid);
            nmols{nclus}(nnegacid)=nmols{nclus}(nnegacid)-1;
        
        end
    
        combined=create_label(nmols{nclus},molnames{nclus});
        return;
        
    else

        % if this is a neutral cluster, it gets protonated if it has base
        if (lhasbase(nclus))
   
            % add the proton
            molnames{nclus}{length(molnames{nclus})+1}=posion;
            nmols{nclus}(length(nmols{nclus})+1)=1;
            
            combined=create_label(nmols{nclus},molnames{nclus});
            return;
            
        else
            
            combined='undef';
            return;

        end
    
    end

% if both parties are molecular clusters
elseif (lnegclus(1) || lnegclus(2))
    
    if lnegclus(1)
        nnegclus=1;
        nclus=2;
    else
        nnegclus=2;
        nclus=1;
    end
    
    if (lnegclus(nclus))
        
        combined='undef';
        return;   
        
    elseif (lposclus(nclus))
        % if this is a recombination of charged clusters of opposite signs, the clusters get neutralized
        
        if ismember(negion,molnames{nnegclus})
        
            % if there is a missing proton, just remove it
            nnegion=strcmp(molnames{nnegclus},negion);
            nmols{nnegclus}(nnegion)=nmols{nnegclus}(nnegion)-1;
            
        else
        
            % add one acid
            nacid=strcmp(molnames{nnegclus},acid);
            if ismember(1,nacid)
                nmols{nnegclus}(nacid)=nmols{nnegclus}(nacid)+1;
            else
                molnames{nnegclus}{length(molnames{nnegclus})+1}=acid;
                nmols{nnegclus}(length(nmols{nnegclus})+1)=1;
            end
        
            % remove the negative ionic acid
            nnegacid=strcmp(molnames{nnegclus},negacid);
            nmols{nnegclus}(nnegacid)=nmols{nnegclus}(nnegacid)-1;
        
        end
        
        % remove the proton
        nposion=strcmp(molnames{nclus},posion);
        nmols{nclus}(nposion)=nmols{nclus}(nposion)-1;
        
        % calculate the molecules together
        [molsum,sumnames]=add_mols(nmols,molnames);
        combined=create_label(molsum,sumnames);
        return;
        
    else
        % if the other party is neutral
        [molsum,sumnames]=add_mols(nmols,molnames);
        combined=create_label(molsum,sumnames);
        return;
        
    end
    
elseif (lposclus(1) || lposclus(2))    
    
    if lposclus(1)
        nclus=2;
    else
        nclus=1;
    end    

    if (lposclus(nclus))
        
        combined='undef';
        return;
        
    else
        % the only option left here is a neutral cluster (all the others should have been checked already)
        [molsum,sumnames]=add_mols(nmols,molnames);
        combined=create_label(molsum,sumnames);
        return;
        
    end

% if both clusters are neutral
else
    
    [molsum,sumnames]=add_mols(nmols,molnames);
    combined=create_label(molsum,sumnames);
    return;  
    
end


end

function [label]=create_label(nums,names)
% takes in 2 vectors containing the numbers and names of the molecules and creates a name label
% e.g. [2 1 1], {'A' 'B' 'N'} -> 2A1B1N

if (length(nums)~=length(names))
    error(['The amount of molecule numbers is not equal to the amount of molecule names: ',nums,', ',names])
end

label='';

for i=1:length(nums)
    if (nums(i)>0)
        label=[label,num2str(nums(i)),names{i}];
    end
end


end

function [molsum,sumnames]=add_mols(nmols,molnames)
% takes in 2 cells containg vectors for the numbers and names of the molecules of different molecular clusters
% and combines the clusters by adding the molecules together

if (~iscell(nmols) || ~iscell(molnames))
    error(['Input should be cells: ',nmols,', ',molnames])
end

for i=1:(length(nmols))
    if (length(nmols{i})~=length(molnames{i}))
        error(['The amount of molecule numbers is not equal to the amount of molecule names: ',nmols{i},', ',molnames{i}])
    end
end

% create the vector for the names e.g. from the names in the 1st cluster
sumnames=molnames{1};
% check if the other clusters have more molecule types
for i=2:length(molnames)
    for j=1:length(molnames{i})
        % check if this is already in the name vector, if not, add it
        lfound=strmatch(molnames{i}{j},sumnames,'exact');
        if (isempty(lfound))
            sumnames=[sumnames,molnames{i}{j}];
        end
    end
end

molsum=zeros(1,length(sumnames));

% calculate together all the molecules that are of the same type
for i=1:length(nmols)
    for j=1:length(nmols{i})
        if (nmols{i}(j)>0)
            nmol=strmatch(molnames{i}{j},sumnames,'exact');
            molsum(nmol)=molsum(nmol)+nmols{i}(j);
        end
    end
end

end


