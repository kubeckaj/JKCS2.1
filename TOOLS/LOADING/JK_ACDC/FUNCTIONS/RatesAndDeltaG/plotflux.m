function [allfluxnames,allfluxes,mainfluxnames,mainfluxes,mainroute] = plotflux(name_mon,C_mon,tar,direction,ltrack,varargin)
% Find the fluxes to or from the wanted cluster; also the complete main flux chain to or from the cluster can be tracked
% Input:
%   name_mon: a cell string containing the names of the monomers, e.g. {'A','N'}
%   C_mon: a vector containing the concentrations of the monomers (see also C_mon_notation below), e.g. [1e6,100]
%   tar: the name of the cluster for which the fluxes are calculated, e.g. '2A1N'
%   direction: 'to' corresponding to fluxes into the cluster, or 'from', corresponding to fluxes out of the cluster
%   ltrack: 1 if all the flux chain is tracked, 0 if not
% Output:
%   allfluxnames: cell string containing the names of the clusters (or other creatures) where all the fluxes go to/come from
%   allfluxes: vector containing the absolute values for all the fluxes
%   mainfluxnames: cell string containing the names of the clusters (or other creatures) where all the fluxes that form at least
%                   a certain fraction of the total flux, determined by the variable crit, go to/come from
%   mainfluxes: vector containing the absolute values for the fluxes that form a certain fraction of the total flux
%   mainroute: cell string containing the names of the clusters along the main flux chain (if ltrack=1)


% Default values for the input variables

C_mon_notation=cell(1,length(name_mon)); C_mon_notation(:)={'e'}; % Notation for the concentration (fixed-point: 'f', exponential: 'e')
suffix=''; % suffix for the files

crit=0.01; % criteria for the most significant fluxes (fraction of the total flux)

lshowall=0; % print all the fluxes to the screen (instead of just the main fluxes)?
lsame_ch=0; % track the fluxes based on the charge and not on the absolute number of molecules
lpairs=0; % return both parties of a collision/evaporation (only the larger one is returned by default)


% See what is given as the function input

if ~isempty(varargin)
    for i=1:length(varargin)
        if strcmpi(varargin{i},'C_mon_notation') | strcmpi(varargin{i},'notation')
            C_mon_notation=varargin{i+1};
        elseif strcmpi(varargin{i},'suffix')
            suffix=varargin{i+1};
        elseif strcmpi(varargin{i},'dirpath') | strcmpi(varargin{i},'dp')
            dirpath=varargin{i+1};
        elseif strcmpi(varargin{i},'crit')
            crit=varargin{i+1};
        elseif strcmpi(varargin{i},'lshowall')
            if length(varargin)>i && isnumeric(varargin{i+1})
                lshowall=varargin{i+1};
            else
                lshowall=1;
            end
        elseif strcmpi(varargin{i},'lsame_ch')
            if length(varargin)>i && isnumeric(varargin{i+1})
                lsame_ch=varargin{i+1};
            else
                lsame_ch=1;
            end
        elseif strcmpi(varargin{i},'lpairs')
            if length(varargin)>i && isnumeric(varargin{i+1})
                lpairs=varargin{i+1};
            else
                lpairs=1;
            end
        elseif strcmpi(varargin{i},'temp')
            temp=varargin{i+1};
        elseif strcmpi(varargin{i},'IPR')
            IPR=varargin{i+1};
        end
    end
end

% Define the notation for the concentrations in the filenames
C_mon_str=cell(1,length(name_mon));
for i=1:length(name_mon)
    if strcmp(C_mon_notation{i},'f')
        C_mon_str{i}=num2str(C_mon(i));
    elseif strcmp(C_mon_notation{i},'e')
        C_mon_str{i}=sprintf('%.2e',C_mon(i));
    end
end

% Default directory path, if one is not given
if ~exist('dirpath','var') && exist('temp','var')
    dirpath=['.\',cell2mat(name_mon),'_',num2str(temp),'\'];
end

% Filename
ffile=['fluxes',cell2mat(strcat('_',name_mon,C_mon_str))];
% Possible IPR in the filename (this has been used sometimes)
if exist('IPR','var')
    ffile=[ffile,'_I',num2str(round(IPR))];
end
% Rest of the filename
ffile=[ffile,suffix,'.txt'];

% Complete path
fn=[dirpath,ffile];
if ~exist(fn,'file')
    % this might be the old typo
    ffile=regexprep(ffile,['_',name_mon{1}],['_',name_mon{1},'_']);
    fn=[dirpath,ffile];
end
% if ~exist(fn,'file')
%     % this might be some sprintf conversion problem
%     fn=regexprep(fn,'e\+0','e\+00');
% end

%%%% ~o~ %%%%

% Read the flux data

% read the header containing the cluster names + other flux names
fid=fopen(fn,'r');
s=textscan(fid,'%s',2,'delimiter','\n');
fclose(fid);
s=strtrim(s);
s2=regexp(s{1,1}{1},'\s+','split');
% ignore the percent mark in the beginning of the header
if strcmp(s2{1},'%')
    clust_flux=s2(2:end);
else
    error([s2{1},'????'])
end
% if there is another header line, assume it's the names of the explicitly simulated clusters
s3=regexp(s{1,1}{2},'\s+','split');
if strcmp(s3{1},'%')
    clust=s3(2:end);
    disp('Found info on the second line, assuming it''s the names of the simulated clusters')
end

flux=load(fn);
if (size(flux,2) ~= length(clust_flux))
    error(['Number of columns in ',fn,' doesn''t match with the header of the data set!'])
end

%%%% ~o~ %%%%

% Find the fluxes

info_str='';
for i=1:length(name_mon)
    info_str=[info_str,'[',name_mon{i},'] = ',C_mon_str{i},' (cm^-3 or ppt)    '];
end
disp(info_str)

% mainroute contains the names of the clusters in the main flux chain;
% corresponding elements in mainfluxes are the largest fluxes to/from this
% cluster and the names corresponding to these fluxes
mainroute={tar};
mainfluxes=cell(1);
mainfluxnames=cell(1);
allfluxes=cell(1);
allfluxnames=cell(1);

count=0;

% use the calcflux function
if ltrack

    lcluster=1;
    
    % track the mainfluxes back to a non-cluster, e.g. out, bound, wall...
    while lcluster
        
        if exist('clust','var')
            [fluxval,fluxname,roundfluxval,roundfluxname]=calcflux(tar,direction,clust_flux,flux,'pie','crit',crit,'lsame_ch',lsame_ch,'lpairs',lpairs,'clust',clust);
        else
            [fluxval,fluxname,roundfluxval,roundfluxname]=calcflux(tar,direction,clust_flux,flux,'pie','crit',crit,'lsame_ch',lsame_ch,'lpairs',lpairs);
        end
        
        % find the main flux and the corresponding cluster name
        [~,nmain]=max(fluxval);
        tar=fluxname{nmain,1};
        
        % check if the cluster is already in the route (i.e. the flux is circulating around the same clusters)
        lfound=strmatch(tar,mainroute,'exact');
        if isempty(lfound)
            count=count+1;
            mainfluxes{count}=roundfluxval;
            mainfluxnames{count}=roundfluxname;
            allfluxes{count}=fluxval;
            allfluxnames{count}=fluxname;
            % this will now be the next cluster in the main route
            mainroute{count+1}=tar;
            [lcluster]=check_cluster(tar);
        else
            disp(['Came back to ',tar,', exiting the loop'])
            break
        end
        
    end
    
else
    
    if exist('clust','var')
        [fluxval,fluxname,roundfluxval,roundfluxname]=calcflux(tar,direction,clust_flux,flux,'pie','crit',crit,'lsame_ch',lsame_ch,'lpairs',lpairs,'clust',clust);
    else
        [fluxval,fluxname,roundfluxval,roundfluxname]=calcflux(tar,direction,clust_flux,flux,'pie','crit',crit,'lsame_ch',lsame_ch,'lpairs',lpairs);
    end
    
    count=count+1;
    mainfluxes{count}=roundfluxval;
    mainfluxnames{count}=roundfluxname;
    allfluxes{count}=fluxval;
    allfluxnames{count}=fluxname;
    
end

for i=1:max(length(mainroute)-1,1)
        
    str=['Main fluxes ',direction,' ',mainroute{i},': '];        
    for j=1:length(mainfluxes{i})
        str=[str,'  ',mainfluxnames{i}{j,1},' ',num2str(round(mainfluxes{i}(j)/sum(allfluxes{i})*100)),'% '];
    end
    str=[str,'  Total flux ',sprintf('%0.3g',sum(allfluxes{i})),' cm^-3s^-1'];
    disp(str)
        
    if lshowall
        disp(['All fluxes ',direction,' ',mainroute{i},': ']);        
        for j=1:length(allfluxes{i})
            disp([allfluxnames{i}{j,1},' ',sprintf('%0.3g',allfluxes{i}(j)/sum(allfluxes{i})*100),'% '])
        end
    end
        
end

end
