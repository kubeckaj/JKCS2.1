% Function for creating a neat output string for scientific notation
function [es_str] = get_es_str(num,prec)
    es_exp=floor(log10(num));
    es_coef=num/10^es_exp;
    if strcmp(sprintf(['%0.',num2str(prec),'f'],es_coef),'10')
        es_str=['10^{',num2str(es_exp+1),'}'];
    elseif strcmp(sprintf(['%0.',num2str(prec),'f'],es_coef),'1')
        es_str=['10^{',num2str(es_exp),'}'];
    else
        es_str=[sprintf(['%0.',num2str(prec),'f'],es_coef),'\cdot','10^{',num2str(es_exp),'}'];
    end
end