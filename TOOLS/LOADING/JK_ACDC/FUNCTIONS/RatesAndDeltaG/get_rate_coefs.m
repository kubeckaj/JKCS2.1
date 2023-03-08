function [varargout]=get_rate_coefs(clus,ratetype,varargin)
% Get the collision or evaporation rates from the ACDC Matlab files
% Input:
% clus = the name of the cluster for which the rate is wanted;
% if clus is a cell string, it contains the colliding/evaporating parties,
% if clus is a string, it is the evaporating cluster
% ratetype = either 'coll' or 'evap'
% keyword 'ldisp' = print the results to the Matlab workspace
% keyword 'lrun' = run the ACDC Perl script to create the rate files (see the settings below)
% keyword 'dp' and the full directory path = directory for the rate files, if it's not the current one
% keyword 'suffix' and the suffix string = suffix in the name of the rate files,
% if  they are not simply 'get_coll.m' and 'get_evap.m'
% Output (optional):
% rate = the rate(s) in SI units (m^3 s^-1 for collisions, s^-1 for evaporations)
% cluslabels = names of the colliding/evaporating parties


ldisp=0;
lrun_perl=0;
dp='';
suffix='';

if ~isempty(varargin)
    for i=1:length(varargin)
        if strcmpi(varargin{i},'ldisp')
            ldisp=1;
        elseif strcmpi(varargin{i},'lrun_perl') | strcmpi(varargin{i},'lrun')
            lrun_perl=1;
        elseif strcmpi(varargin{i},'dp')
            dp=varargin{i+1};
        elseif strcmpi(varargin{i},'suffix')
            suffix=varargin{i+1};
        end
    end
end

if lrun_perl
    disp('Running the Perl script: check that the Perl call is correct!')
    % DeltaG and dip_pol files
    Gfile='\\ATKK\home\n\nmyllys\Desktop\RateConstants\Guanidine-RC\G298.15K_AG.txt';
    DPfile='\\ATKK\home\n\nmyllys\Desktop\RateConstants\Guanidine-RC\G298.15K_AG.txt';
    inputfile='\\ATKK\home\n\nmyllys\Desktop\RateConstants\Guanidine-RC\ACDC_input\input_AGN.inp';
    rh=0;
  
    commandstr=['perl create_acdc_2016_09_30.pl --i ',inputfile,' --e ',Gfile,' --dip ',DPfile,' --rh ',num2str(rh)];
    lrun=dos(commandstr);
    if lrun ~= 0
        error('Running Perl failed!')
    end
    rehash pathreset
end

% "Type" of the process
lcoll=0;
levap=0;
if strcmp(ratetype,'coll')
    lcoll=1;
    ratefile=[dp,'get_coll',suffix,'.m'];
elseif strcmp(ratetype,'evap')
    levap=1;
    ratefile=[dp,'get_evap',suffix,'.m'];
end

% Charged cluster: lch=1, neutral cluster: lch=0
lch=0;

% Is this the rate for specific daughter clusters?
% If not, then it's all possible evaporations of the cluster
lspecific=0;
if levap && ischar(clus)
    cluslabels={};
    % The comments in the files should be in this format
    str0=['% ',clus,' -> '];
    % Is this an ion evaporation?
    if solve_charge(clus) ~= 1
       lch=1;
    end
elseif iscellstr(clus) && length(clus)==2
    lspecific=1;
    cluslabels={[clus{1},' + ',clus{2}]};
    % The comments in the files should be in this format
    str1=[' ',clus{1},' + ',clus{2}];
    str2=[' ',clus{2},' + ',clus{1}];
    % Is this an ion collision/evaporation?
    for i=1:2
        if solve_charge(clus{i}) ~= 1
            lch=1;
        end
    end
end

if lch==1
    % Enhancement factor file for charged clusters
    fcrfile=['get_fcr',suffix,'.m'];
    if exist('fcrfile','file')
        disp('You seem to have an FCR file (i.e. a really antique version of ACDC??);')
        error('if you want to use it, comment out the row where lch forced to 0 by default')
    end
end
% lch is now always forced to 0 by default, as there is no longer a separate FCR file in the Perl output
lch=0;

% imax is the number of files that are read (1 for neutral clusters, 2 for ions if there is an FCR file)
if lch==0
    imax=1;
else
    imax=2;
end
for i=1:imax
    if lspecific
        value=NaN;
    else
        if i==1
            j=0;
            cluslabels=cell(1,1);
            no_enh_rate=NaN(1,1);
        else
            enh=NaN(size(no_enh_rate));
        end
    end
    if i==1
        filename=ratefile;
    else
        filename=fcrfile;
    end
    fid=fopen(filename,'r');
    line=fgetl(fid);
    while ischar(line)
        if lspecific
            found1=strfind(line,str1);
            found2=strfind(line,str2);
            if (~isempty(found1) && (length(str1)==length(line(found1:end))...
                    || isempty(regexp(line(found1+length(str1)),'[A-Za-z0-9]','match'))))...
                    || (~isempty(found2) && (length(str2)==length(line(found2:end))...
                    || isempty(regexp(line(found2+length(str2)),'[A-Za-z0-9]','match'))))
                %disp(line)
                split1=regexp(line,';','split');
                split2=regexp(split1{1},'=','split');
                split3=regexp(split1{end},'%','split');
                cluslabels_info=split3{end};
                eval(['value=',split2{end},';']);
                break
            end
        else
            if i==1
                if ~isempty(strfind(line,str0))
                    %disp(line)
                    split1=regexp(line,';','split');
                    split2=regexp(split1{1},'=','split');
                    split3=regexp(split1{end},'->','split');
                    j=j+1;
                    cluslabels{j}=split3{end};
                    eval(['no_enh_rate(j)=',split2{end},';']);
                end
            else
                for j=1:length(cluslabels)
                    found0=strfind(line,cluslabels{j});
                    if ~isempty(found0) && length(cluslabels{j})==length(line(found0:end))
                        %disp(line)
                        split1=regexp(line,';','split');
                        split2=regexp(split1{1},'=','split');
                        eval(['enh(j)=',split2{end},';']);
                        break
                    end
                end
            end
        end
        line = fgetl(fid);
    end
    fclose(fid);
    if lspecific
        if i==1
            if isnan(value)
                disp(['Couldn''t find the rate from ',filename,' for ',clus])
            end
            no_enh_rate=value;
        else
            if isnan(value)
                value=1;
                disp(['Couldn''t find FCR for ',clus])
            end
            enh=value;
        end
    else
        if i==1
            if isequaln(no_enh_rate,NaN(1,1))
                disp(['Couldn''t find any rates for ',clus])
            end
        else
            if isequaln(enh,NaN(size(enh)))
                enh=ones(size(enh));
                disp(['Couldn''t find any FCRs for ',clus])
            else
                for j=1:length(enh)
                    if isnan(enh(j))
                        enh(j)=1;
                        disp(['Couldn''t find FCR for ',cluslabels{j}])
                    end
                end
            end
        end
    end
end

if lch==0
    rate=no_enh_rate;
else
    rate=enh.*no_enh_rate;
end

if lspecific && ldisp && exist('cluslabels_info','var')
    if lcoll
        disp([cluslabels_info,'    ',sprintf('%.6e',rate),' m^3 s^-1'])
    elseif levap
        disp([cluslabels_info,'    ',sprintf('%.6e',rate),' s^-1'])
    end
elseif ~lspecific
    [rate,ind]=sort(rate,'descend');
    cluslabels=cluslabels(ind);
    if lch==1
        no_enh_rate=no_enh_rate(ind);
        enh=enh(ind);
    end
    for j=1:length(cluslabels)
        if ldisp
            disp([cluslabels{j},'    ',sprintf('%.6e',rate(j)),' s^-1'])
        end
        if ~isempty(strfind(cluslabels{j},','))
            split1=regexp(cluslabels{j},',','split');
            cluslabels{j}=split1{1};
        end
    end
    if ldisp
        disp(['Total rate    ',sprintf('%.6e',sum(rate,'omitnan')),' s^-1'])
    end
end

if nargout>0
    output={rate,cluslabels};
    varargout=output(1:nargout);
end

end