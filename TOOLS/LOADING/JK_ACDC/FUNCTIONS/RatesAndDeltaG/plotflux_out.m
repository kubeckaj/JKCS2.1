function [alloutname,allout,roundoutname,roundout] = plotflux_out(name_mon,C_mon,varargin)
% Find the fluxes out from the system
% Input:
%   name_mon: a cell string containing the names of the monomers, e.g. {'A','N'}
%   C_mon: a vector containing the concentrations of the monomers (see also C_mon_notation below), e.g. [1e6,100]
% Output:
%   alloutname: cell string containing the names of the cluster pairs contributing to the outgoing flux, e.g. '5A5N+1A'
%   allout: vector containing the absolute values for all the outgoing fluxes
%   roundoutname: cell string containing the names of the outgrowing cluster pairs for all the fluxes that form at least
%                   a certain fraction of the total flux, determined by the variable crit
%   roundout: vector containing the absolute values for the fluxes that form a certain fraction of the total flux


% Default values for the input variables

C_mon_notation=cell(1,length(name_mon)); C_mon_notation(:)={'e'}; % Notation for the concentration (fixed-point: 'f', exponential: 'e')
suffix=''; % suffix for the files

crit=0.01; % criteria for the most significant fluxes (fraction of the total flux)

lshowall=0; % print all the outgoing fluxes to the screen (instead of just the main fluxes)?
lclusters_out=0; % list also the clusters that got out?

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
            lshowall=1;
        elseif strcmpi(varargin{i},'lclusters_out')
            lclusters_out=1;
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
ffile=['outmat',cell2mat(strcat('_',name_mon,C_mon_str))];
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

% read the header containing the cluster names
fid=fopen(fn,'r');
s=textscan(fid,'%s',1,'delimiter','\n');
fclose(fid);
s=strtrim(s);
s2=regexp(s{1,1},'\s+','split');
% ignore the percent mark in the beginning of the header
if strcmp(s2{1,1}{1},'%')
    clust=s2{1,1}(2:end);
else
    error([s2{1,1}{1},'????'])
end

outflux=load(fn);
if size(outflux,1) ~= length(clust) || size(outflux,2) ~= length(clust)
    error(['Number of rows or columns in ',fn,' doesn''t match with the header of the data set!'])
end

%%%% ~o~ %%%%

% Find the fluxes

info_str='';
for i=1:length(name_mon)
    info_str=[info_str,'[',name_mon{i},'] = ',C_mon_str{i},' (cm^-3 or ppt)    '];
end
disp(info_str)

outflux_clean=removedoubles(outflux,clust);
outfluxsum=sum(sum(outflux_clean));

allout=[];
alloutname={};
alloutclust={};
l=0;
for i=1:size(outflux_clean,1)
    for j=1:size(outflux_clean,2)
        if outflux_clean(i,j)~=0
            l=l+1;
            allout(l)=outflux_clean(i,j);
            alloutname{l}=[clust{i},'+',clust{j}];
            alloutclust{l}=combine_clusters(clust{i},clust{j});
        end
    end
end

% sort them in descending order
[allout,ind]=sort(allout,'descend');
alloutname=alloutname(ind);
alloutclust=alloutclust(ind);

% to make the plot clearer, just take the largest fluxes and put the rest under the label 'Others'
if isnan(crit)
    [roundoutname,roundout] = get_significant(alloutname,allout);
else
    [roundoutname,roundout] = get_significant(alloutname,allout,'crit',crit);
end

disp(['Total flux out ',sprintf('%0.3g',outfluxsum),' cm^-3s^-1'])
disp('Main collisions leading out of the system: ');
for i=1:length(roundout)
        disp([roundoutname{i},' ',num2str(round(roundout(i)/outfluxsum*100)),'% '])
end
if lshowall
    disp('All collisions out: ');
    for i=1:length(allout)
        disp([alloutname{i},' ',sprintf('%0.3g',allout(i)/outfluxsum*100),'% '])
    end
end

figure(11)
set(gca,'LineWidth',2,'FontWeight','bold','FontSize',12);
pie(roundout/outfluxsum)
hold on;
legend(roundoutname,'Location','BestOutside')
title(['Total flux out ',sprintf('%0.3g',outfluxsum),' cm^{-3}s^{-1}'])
hold off;

% list the clusters getting out?
if lclusters_out
    % first see if the same outgoing cluster forms from different collisions
    allout_forclust=allout;
    for i=1:length(alloutclust)
        if ~strcmpi(alloutclust{i},'')
            for j=(i+1):length(alloutclust)
                lsame=compare_clusters(alloutclust{i},alloutclust{j});
                if lsame
                    % add the value to the current element and reset the value of the other element to zero
                    allout_forclust(i)=allout_forclust(i)+allout_forclust(j);
                    allout_forclust(j)=0;
                    alloutclust{j}='';
                end
            end
        end
    end
    % remove empty elements
    ind=find(allout_forclust==0);
    allout_forclust(ind)=[];
    alloutclust(ind)=[];
    
    % sort them in descending order
    [allout_forclust,ind]=sort(allout_forclust,'descend');
    alloutclust=alloutclust(ind);
    
    [roundoutclust,roundout_forclust] = get_significant(alloutclust,allout_forclust,'crit',crit);

    disp('Clusters getting out of the system: ');
    for i=1:length(roundout_forclust)
        disp([roundoutclust{i},' ',num2str(round(roundout_forclust(i)/outfluxsum*100)),'% '])
    end
    
    figure(2)
    set(gca,'LineWidth',2,'FontWeight','bold','FontSize',12);
    pie(roundout_forclust/outfluxsum)
    hold on;
    legend(roundoutclust)
    title(['Total flux out ',sprintf('%0.3g',outfluxsum),' cm^{-3}s^{-1}'])
    hold off;
end

end

function [outflux_clean] = removedoubles(outflux,clust)
% OLDER FORMAT: element (j,i) should be the same as (i,j); set the one corresponding to the smaller cluster
% to zero; in this case also fluxes due to collisions of two similar clusters need to be divided by two
% NEWER FORMAT: just see that the flux is in the element corresponding to the larger cluster
% (and the element corresponding to the smaller cluster is zero); change this, if needed
outflux_clean=zeros(size(outflux));
lnew_format=0;
for i=1:length(clust)
    for j=i:length(clust)
        if outflux(i,j)~=0
            if i~=j
                if outflux(j,i)==0
                    lnew_format=1;
                elseif outflux(i,j)~=outflux(j,i)
                    error(['Fluxes don''t match: ',clust{i},' and ',clust{j},' fluxes out ',sprintf('%0.3g',outflux(i,j)),' and ',sprintf('%0.3g',outflux(j,i))])
                end
                
                largest=max_mols(clust{i},clust{j});
                if strcmpi(largest,clust{i})
                    outflux_clean(i,j)=outflux(i,j);
                    outflux_clean(j,i)=0;
                else
                    outflux_clean(j,i)=outflux(i,j);
                    outflux_clean(i,j)=0;
                end
            else
                outflux_clean(i,j)=outflux(i,j);
            end
        end
    end
end
if lnew_format == 0
    disp('Using the older outflux_matrix format! i.e. dividing fluxes from similar clusters by two')
    for i=1:length(clust)
        % division by 2 since the fluxes are always from the point of view of the smaller cluster
        outflux_clean(i,i)=outflux_clean(i,i)/2;
    end
end

end
