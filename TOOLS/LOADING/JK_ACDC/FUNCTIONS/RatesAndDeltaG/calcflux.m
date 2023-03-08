function [varargout] = calcflux(tar,direction,clust_flux,flux,varargin)
% Checks the fluxes into or out from a specific cluster (tar) and draws a pie chart of them
% Input: cluster name, direction (i.e. if we want to track the flux into or out from the cluster), names of the fluxes (clust_flux) and flux matrix (from ACDC),
% additional input for plotting and the criterion for if a flux is taken into account individually or put under the label 'others' in the plot  
% Output arguments (if specified) are the fluxes and the corresponding cluster names, and the fluxes and names
% where the fluxes that don't meet the criterion are put together under the label 'others'

% Names of the generic ions
genions={'neg','pos','nion','pion'};

% logicals that tell the direction of the flux
lfrom=0;
lto=0;

if strcmpi(direction,'to')
    lto=1;
elseif strcmpi(direction,'from')
    lfrom=1;
    % flip the flux matrix
    flux=flux.';
else
    error('Input argument ''direction'' not well defined; should be ''to'' or ''from''.')
end

% possible other input arguments
% default values
lpie=0;
crit=NaN;
lsame_ch=0;
lpairs=0;
lgrowth=0; % when the direction is "from", growth to sizes that are larger than
% 2 times the number of molecules in the growing ("target") cluster, i.e. if A+B->C, C_mols<=2*A_mols, are
% classified as self-coagulation

if ~isempty(varargin)
    for i=1:length(varargin)
        if strcmpi(varargin{i},'crit')
            crit=varargin{i+1};
        elseif strcmpi(varargin{i},'pie')
            lpie=1;
        elseif strcmpi(varargin{i},'lsame_ch')
            if length(varargin)>i && isnumeric(varargin{i+1})
                lsame_ch=varargin{i+1};
            else
                lsame_ch=1;
            end
            if lsame_ch
                tar_ch=solve_charge(tar);
            end
        elseif strcmpi(varargin{i},'clust')
            clust=varargin{i+1}; % only the "real" clusters, i.e. no source, out, wall etc.
        elseif strcmpi(varargin{i},'lpairs')
            if length(varargin)>i && isnumeric(varargin{i+1})
                lpairs=varargin{i+1};
            else
                lpairs=1;
            end
        elseif strcmpi(varargin{i},'lgrowth')
            if length(varargin)>i && isnumeric(varargin{i+1})
                lgrowth=varargin{i+1};
            else
                lgrowth=1;
            end
        end
    end
end

% find the column of the target cluster
col=strmatch(tar,clust_flux,'exact');
if isempty(col)
    error(['Cluster not found: ',tar])
end

% Size of the target cluster
tar_mols=calc_mols(tar);

% See if the target cluster is a boundary cluster outside the true system
lbound_tar=0;
if exist('clust','var')
    if check_cluster(tar) && ~ismember(tar,clust)
        lbound_tar=1;
    end
end

% ltaken is a logical that tells if the cluster in question has already been taken into account in the flux
% (because fluxes from collisions are otherwise counted twice)
ltaken=zeros(1,length(clust_flux));

if lgrowth
    if ~lfrom
        error('lgrowth can be set to 1 only if you are tracking the fluxes FROM a cluster!')
    end
    % if we're studying cluster growth, save the self-coagulation fluxes to variable self_coag_flux
    self_coag_flux=0.0;
end

% index in the flux vector
k=0;

for i=1:size(flux,1)
    
    if ((flux(i,col)~=0) && (ltaken(i)==0))
        
        % first check if this is a flux from a generic ion; then it should not be taken into account
        % (since the corresponding recombination/ionization flux will be missing, i.e. this is not a net flux)
        lgenion=strcmpi(clust_flux{i},genions);
        if ismember(1,lgenion)
            ltaken(i)=1;
            %disp(['Discarded flux from ',clust_flux{i}])
            continue
        end
        
        % If we're tracking the flux from a boundary cluster outside the system, we are not interested in the evaporating monomers
        % (this has to be excluded here since combining the product clusters doesn't work because there may be more than one monomer)
        if lbound_tar && lfrom
            [~,nmols]=parse_cluster(clust_flux{i});
            if length(nmols)==1
                if nmols(1)==1
                    ltaken(i)=1;
                    disp([tar,': Discarded flux to ',clust_flux{i}])
                    continue
                end
            end
        end
        
        lfound=0;
        
        if check_cluster(clust_flux{i})
            clust_mols=calc_mols(clust_flux{i});
            if clust_mols > tar_mols
                % count this as self-coagulation if the result cluster is "too" large
                if lgrowth && clust_mols >= 2*tar_mols
                    %disp(['Self-coagulation flux to ',clust_flux{i},' from ',tar])
                    if clust_mols > 2*tar_mols
                        self_coag_flux=self_coag_flux+flux(i,col);
                    else
                        % if the cluster collides with another similar cluster,
                        % one of them is considered to "grow", and the other one goes to self-coagulation
                        self_coag_flux=self_coag_flux+flux(i,col)./2;
                        
                        k=k+1;
                        fluxname{k,1}=clust_flux{i};
                        fluxval(k)=flux(i,col)./2;
                        if lpairs
                            fluxname{k,2}='';
                        end
                    end
                else
                    k=k+1;
                    fluxname{k,1}=clust_flux{i};
                    fluxval(k)=flux(i,col);
                    if lpairs
                        fluxname{k,2}='';
                    end
                end
                % no need to try to find the other collider, it just makes the script slow (and can't find it anyway)
                continue
            end
        end
        
        % if the flux is from a collision in the case of flux in, or from an evaporation in the case of flux out
        %(or ionization), find the other party
        for j=(i+1):size(flux,1)
            if ltaken(j)==0
                % combine the clusters
                combined=combine_clusters(clust_flux{i},clust_flux{j});
                % check if the result is the target cluster
                lsame=compare_clusters(combined,tar);
                if lsame
                    lfound=1;
                    % this cluster has now been taken into account
                    ltaken(j)=1;
                    
                    k=k+1;
                    % name of the cluster
                    [largest,smallest]=max_mols(clust_flux{i},clust_flux{j});
                    if lsame_ch==1
                        % choose the one with the same charge as the target cluster
                        clust_ch=solve_charge(largest);
                        if clust_ch==tar_ch
                            fluxname{k,1}=largest;
                            if lpairs
                                fluxname{k,2}=smallest;
                            end
                        else
                            clust_ch=solve_charge(smallest);
                            if clust_ch==tar_ch
                                fluxname{k,1}=smallest;
                                if lpairs
                                    fluxname{k,2}=largest;
                                end
                            else
                                % no match, choose the bigger one
                                fluxname{k,1}=largest;
                                if lpairs
                                    fluxname{k,2}=smallest;
                                end
                            end
                        end
                    else
                        % choose the bigger one
                        fluxname{k,1}=largest;
                        if lpairs
                            fluxname{k,2}=smallest;
                        end
                    end
					
                    % numerical value of the flux
                    fluxval(k)=flux(i,col);
                    
                    break
                end
            end
        end
        if (lfound==0)
            % this flux is not from a collision (for flux in)/evaporation (for flux out)...
            k=k+1;
            fluxname{k,1}=clust_flux{i};
            
            % ... except if it's two identical clusters colliding
            % then the flux needs to be divided by 2, since it's from the
            % point of view of the colliding parties
            combined=combine_clusters(clust_flux{i},clust_flux{i});
            lsame=compare_clusters(combined,tar);
            if lsame
                fluxval(k)=flux(i,col)./2;
                if lpairs
                    fluxname{k,2}=clust_flux{i};
                end
                %disp(['Divided: ',clust_flux{i}])
            else
                fluxval(k)=flux(i,col);
                if lpairs
                    fluxname{k,2}='';
                end
            end
        end
    end

end

if ~exist('fluxval','var')
    % no fluxes found??
    disp(['No fluxes ',direction,' ',tar,'!'])
    fluxval=[];
    fluxname={};
    roundfluxval=[];
    roundfluxname={};
else
    if lgrowth && self_coag_flux > 0
        k=k+1;
        fluxname{k,1}='self_coag';
        fluxval(k)=self_coag_flux;
        if lpairs
            fluxname{k,2}='';
        end
    end
    
    % sort the fluxes in descending order
    [fluxval,ind]=sort(fluxval,'descend');
    fluxname=fluxname(ind,:);

    % to make the plot clearer, just take the largest fluxes and put the rest
    % under the label 'others'
    if isnan(crit)
        [roundfluxname,roundfluxval] = get_significant(fluxname,fluxval);
    else
        [roundfluxname,roundfluxval] = get_significant(fluxname,fluxval,'crit',crit);
    end
end

% Return the output parameters, if there are any
if nargout>0
    output={fluxval,fluxname,roundfluxval,roundfluxname};
    varargout=output(1:nargout);
end

if lpie
    % do a pie plot
    figure(12)
    set(gca,'LineWidth',2,'FontWeight','bold','FontSize',12);
    pie(roundfluxval/sum(fluxval))
    hold on
    legend(roundfluxname{:,1},'Location','BestOutside')
    title(['Total flux ',direction,' ',tar,' ',sprintf('%0.3g',sum(fluxval)),' cm^{-3}s^{-1}'])
    drawnow
    hold off
end

end