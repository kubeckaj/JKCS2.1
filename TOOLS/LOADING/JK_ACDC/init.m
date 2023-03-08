clc
clear all

addpath ACDC_input
addpath FUNCTIONS 
           %pa2cmmc = Pa -> cm-3
           %ppt2cmmc = ppt -> cm-3
addpath FUNCTIONS/RatesAndDeltaG -end
JKconstants %Na,R,kB,pres_atm
            %J2kcal,amu2kg, Torr2Pa
            
%RUN folder
%if not(isfolder("RUN")); mkdir('RUN'); end; cd RUN; 
pwdf=pwd;if pwdf(end-2:end)~='RUN';if not(isfolder("RUN")); mkdir('RUN'); end; cd RUN;end
