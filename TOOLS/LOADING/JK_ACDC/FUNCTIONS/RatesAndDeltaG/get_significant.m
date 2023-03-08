function [roundnames,roundvalues] = get_significant(names,values,varargin)
% Reads in two vectors: first has the names of some quantities and second has numerical values corresponding to the names
% Takes all the values that are more than the set limit value of the sum of the number vector
% and puts the rest under the label 'Others', and returns these new vectors
% Optional 3rd input variable is the limit value (in the format 'crit',value); if it's not given, the default is 3% of the sum

% default value for the limit (1% of the sum)
crit=0.01;

if ~isempty(varargin)
    for i=1:length(varargin)
        if strcmpi(varargin{i},'crit')
            crit=varargin{i+1};
        end
    end
end

if size(names,1)==1
    names=names';
end

roundvalues=[];
roundnames=cell(1,size(names,2));
k=0;
for i=1:length(values)
    % take just the significant values
    if (values(i)/sum(values,'omitnan') >= crit)
        k=k+1;
        roundvalues(k)=values(i);
        roundnames(k,:)=names(i,:);
    end
end
k=k+1;
roundnames{k,1}='Others';
roundvalues(k)=sum(values,'omitnan')-sum(roundvalues(1:(k-1)),'omitnan');


end