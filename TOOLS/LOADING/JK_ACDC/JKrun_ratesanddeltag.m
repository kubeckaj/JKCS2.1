%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INTIALIZATION %%
init
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MANUAL %%
%1) in ACDC_input folder, create: input.inp, dG.txt (or dHdS.txt)
%2) in JK folder, adjust: ACDCinit, ACDCsource
%3) run single job: JKrun_ratesanddeltag
%4) plot some properties: plotJT, plotCT, plotCd
%   plotSG, plotAG, plotKX, plotKY, plotKXY, plotE, plotKXYdivE, plotFLUXES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER AREA %%

temp=280;  
temp_conv=280;
%define temp_conv in case of dH and dS, otherwise based on the 2.input line
    
sa_conc=10^7; %cmmc=cm^3   10^6-10^8
am_conc=100;  %ppt         20-150
%msa_conc=10^7;  %ppt      20-150

Xname='A';
Yname='N';
%Zname='MSA';
monomers_concentrations={{'A',sa_conc},{'N',ppt2cmmc(am_conc,temp)}};
%monomers_concentrations={{'A',sa_conc},{'N',ppt2cmmc(dma_conc,temp)},{'MSA',msa_conc}}; % put in cm^3 units (cmmc)

% create your map and check it first!
max_X=5;  %e.g. max_X=2 -> max 2SA
max_Y=5; %e.g. max_Y=3 -> max 3AM
map=cell(max_Y+1,max_X+1);
for y=0:max_Y
    for x=0:max_X
        name='';
        %name=[name,[num2str(1),Zname]]; %ADD EXTRA MOLECULE
        if y>0; name=[name,[num2str(y),Yname]];end
        if x>0; name=[name,[num2str(x),Xname]];end       
        map{y+1,x+1}=name;
    end
end
map=flip(map) %visualize the final map


% WHAT TO CALCULATE  1=yes,0=no
CALCstandarddG=1;   %this will be calculated anyway
CALCactualdG=1;     %required for plotSG
CALCevaporations=1;
CALCcollisions=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER AREA RUN %%
fprintf("JK: ACDC run.\n");

ACDCinit(temp);                    %adjust ACDCinit
ACDCsource(sa_conc,ppt2cmmc(am_conc,temp)) %adjust ACDCsource
[T,C,J_out]=ACDCrun(); 

get_evap; get_coll; get_wl; get_cs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD DATA %% --DO NOT TOUCH
fprintf("JK: Loading free energies files.\n");
fn_dg="../ACDC_input/HS298.15K.txt";

% Read the free energy file
fid=fopen(fn_dg);
% Read the pressure and temperature (1st and 2nd lines)
s=textscan(fid,'%f',2,'Delimiter','\n');
if ~isnumeric(s{1})
    error(['The 1st and 2nd lines in ',fn_dg,' should be the pressure and the temperature!'])
else
    pres=s{1}(1);
    temp=s{1}(2);
end
% Read in the DeltaGs
% Find the number of columns
testline=fgetl(fid);
while isempty(testline) || ismember(1,regexp(testline,'\s*#'))
    testline=fgetl(fid);
end
ndatacols=numel(regexp(testline,'\S+'));

frewind(fid)
if ndatacols == 3
    % DeltaHs and DeltaSs
    fprintf('\nAssuming that DeltaHs and DeltaSs are given instead of DeltaGs\n')
    if exist('temp_conv','var')
        temp=temp_conv;
    end
    M=textscan(fid,'%s %f %f','HeaderLines',2,'CommentStyle','#');
    clust=M{1,1};
    deltah=M{1,2};
    deltas=1e-3.*M{1,3}; % cal/molK to kcal/molK
    deltag=deltah-temp*deltas;
elseif ndatacols == 2
    % DeltaGs
    M=textscan(fid, '%s %f','HeaderLines',2,'CommentStyle','#');
    clust=M{1,1};
    deltag=M{1,2};
end
fclose(fid);
% -> press,temp(=temp_conv or 2nd line),
%    M,clust,deltag,(deltah,deltas)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CHECK MONOMERS %% --DO NOT TOUCH
% Reference DeltaGs should be zero for monomers (they are probably not
% listed in the DeltaG file since they are assumed to be zero by default)
fprintf("JK: Searching for monomers.\n");

%find all types of units
monomers={};
for i=1:length(clust)
    [i1,i2]=parse_cluster(clust{i});
    for j=1:length(i1)
        test=0;
        for k=1:length(monomers)
            if strcmpi(i1{j},monomers{k})
                test=1;
            end
        end
        if test==0
           monomers_size=length(monomers)+1;
           monomers{monomers_size}=i1{j}; 
        end
    end
end
%if monomer is not defined in the input file, let us add it to the end
for i=1:length(monomers)
    test=0;
    for j=1:length(clust)
        if strcmpi(clust{j},['1',monomers{i}])
            test=1;
        end
    end
    if test==0
        nclust=length(clust)+1;
        clust{nclust}=['1',monomers{i}];
        deltag(nclust)=0;
    end
end
% -> adjust: clust, deltag 
%    monomers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% READ ACDC CLUSTERS %% --DO NOT TOUCH
fprintf("JK: Loading ACDC output.\n");

fp=fopen('../RUN/driver_acdc.m','r');
line=0;
while line ~= -1
    line=[fgetl(fp),' '];
    if ~isempty(strfind(line, 'clust = {'))
        line=strrep(line,'clust','clust_acdc');
        eval(line);
        break
    end
end
fclose(fp);
% -> clust_acdc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RECALCULATE dG %% --DO NOT TOUCH
fprintf("JK: Recalculating dG.\n");

% DeltaG surface in the given concentrations
deltag_new=NaN(size(deltag));
for i=1:length(deltag)
    deltag_new(i)=deltag(i);
    [i1,i2]=parse_cluster(clust{i});
    for j=1:length(monomers_concentrations)
        [j1,j2]=parse_cluster(monomers_concentrations{j}{1});
        j3=monomers_concentrations{j}{2};
        for k=1:length(i1)
           if strcmpi(i1{k},j1)
              deltag_new(i)=deltag_new(i)-kB*temp*i2(k)*log(cmmc2pa(j3,temp)/pres)*J2kcal; 
           end
        end
    end
end
%take XY monomer conc.
for j=1:length(monomers_concentrations)
	[j1,j2]=parse_cluster(monomers_concentrations{j}{1});
	j3=monomers_concentrations{j}{2};
    if compare_clusters(['1',Xname],['1',j1{1}])
       Xconc=j3;
    end
    if compare_clusters(['1',Yname],['1',j1{1}])
       Yconc=j3;
    end
end
% -> deltag_new,Xconc,Yconc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAP ALL VARIABLES %% --DO NOT TOUCH
fprintf("JK: Mapping all variables.\n");

deltag_mat=NaN(max_Y+1,max_X+1);     % Reference and actual DeltaGs
deltag_mat_new=NaN(max_Y+1,max_X+1); % (+1 because of the elements with 0 molec.)    
evap_tot_mat=NaN(max_Y+1,max_X+1);   % Total evaporation frequency (s^-1)
monX_coll_mat=NaN(max_Y+1,max_X+1);  % Monomer collision frequencies (s^-1)
monY_coll_mat=NaN(max_Y+1,max_X+1);  % Monomer collision frequencies (s^-1)
clust_mat=cell(max_Y+1,max_X+1);     % Names of the clusters

if CALCactualdG
fprintf("JK: - deltaG\n");
for x=0:max_X
    for y=0:max_Y
        for i=1:length(clust)
            if compare_clusters(map{y+1,x+1},clust{i})
                deltag_mat(y+1,x+1)=deltag(i);
                deltag_mat_new(y+1,x+1)=deltag_new(i);              
            end
        end
    end
end
end
if CALCevaporations
fprintf("JK: - evaporation\n");
for x=0:max_X
    for y=0:max_Y
        for i=1:length(clust_acdc)
            if compare_clusters(map{y+1,x+1},clust_acdc{i})
                for j=1:length(clust_acdc)
                    for k=j:length(clust_acdc)
                        %[clust_acdc{j}," ",clust_acdc{k}]
                        if compare_clusters(combine_clusters(clust_acdc{j},clust_acdc{k}),clust_acdc{i})
                            if isnan(evap_tot_mat(y+1,x+1)); evap_tot_mat(y+1,x+1)=0;end
                            %x
                            %y
                            if not(isnan(E(j,k)))
                              evap_tot_mat(y+1,x+1)=evap_tot_mat(y+1,x+1)+E(j,k);
                            end
                            %if compare_clusters(clust_acdc{j},['1',Xname]) || compare_clusters(clust_acdc{k},['1',Xname])
                            %    if isnan(monX_coll_mat(y+1,x+1)); monX_coll_mat(y+1,x+1)=0;end
                            %    monX_coll_mat(y+1,x+1)=monX_coll_mat(y+1,x+1)+Xconc*K(j,k)/kB/temp;
                            %end
                            %if compare_clusters(clust_acdc{j},['1',Yname]) || compare_clusters(clust_acdc{k},['1',Yname])
                            %    if isnan(monY_coll_mat(y+1,x+1)); monY_coll_mat(y+1,x+1)=0;end
                            %    monY_coll_mat(y+1,x+1)=monY_coll_mat(y+1,x+1)+Yconc*K(j,k)/kB/temp;
                            %end
                        end
                    end
                end
            end
        end
    end
end
end
if CALCcollisions
fprintf("JK: - collision\n");
for x=0:max_X
    for y=0:max_Y
        for i=1:length(clust_acdc)
            if compare_clusters(map{y+1,x+1},clust_acdc{i})
                for j=1:length(clust_acdc)
                    %for k=j:length(clust_acdc)
                        %if compare_clusters(combine_clusters(clust_acdc{j},clust_acdc{i}),clust_acdc{k})
                            %if isnan(evap_tot_mat(y+1,x+1)); evap_tot_mat(y+1,x+1)=0;end
                            %evap_tot_mat(y+1,x+1)=evap_tot_mat(y+1,x+1)+E(j,k);
                            if compare_clusters(clust_acdc{j},['1',Xname])
                                if isnan(monX_coll_mat(y+1,x+1)); monX_coll_mat(y+1,x+1)=0;end
                                monX_coll_mat(y+1,x+1)=monX_coll_mat(y+1,x+1)+cmmc2pa(Xconc,temp)*K(j,i)/kB/temp;
                            end
                            if compare_clusters(clust_acdc{j},['1',Yname]) 
                                if isnan(monY_coll_mat(y+1,x+1)); monY_coll_mat(y+1,x+1)=0;end
                                monY_coll_mat(y+1,x+1)=monY_coll_mat(y+1,x+1)+cmmc2pa(Yconc,temp)*K(j,i)/kB/temp;
                            end
                        %end
                    %end
                end
            end
        end
    end
end
end
fprintf("JK: Done\n");

cd ..