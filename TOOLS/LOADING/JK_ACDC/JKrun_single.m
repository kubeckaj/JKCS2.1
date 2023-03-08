%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INTIALIZATION %%
init
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MANUAL %%
%1) in ACDC_input folder, create: input.inp, dG.txt (or dHdS.txt)
%2) in ACDC folder, adjust: ACDCinit, ACDCsource
%3) run single job: JKrun_single
%4) plot some properties: plotJT, plotCT, plotCd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER AREA %%
temp=280; %Temperature

sa_conc=10^7; %cmmc=cm^-3 
am_conc=10; %ppt

ACDCinit(temp);                    %adjust ACDCinit
ACDCsource(sa_conc,ppt2cmmc(am_conc,temp)) %adjust ACDCsource
[T,C,J_out]=ACDCrun(); 
J_out(end)

%get_evap; get_coll; get_wl; get_cs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%leave back from RUN folder
cd '../'