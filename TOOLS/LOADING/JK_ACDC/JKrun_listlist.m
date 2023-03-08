%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INTIALIZATION %%
init
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MANUAL %%
%1) in ACDC_input folder, create: input.inp, dG.txt (or dHdS.txt)
%2) in ACDC folder, adjust: ACDCinit, ACDCsource
%3) run 
%4) plot VAR1(different lines) and VAR2(x-axis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER AREA %%
temp=298.15; %Temperature

resultss=[];
VAR1range=[10 100 1000]; %different lines
for VAR1=VAR1range
    results=[];
    VAR2range=6:0.5:8; %x-axis
    VAR2range=10.^VAR2range;
    for VAR2=VAR2range
        fprintf("->"+num2str(VAR1)+"/"+num2str(VAR2));
        
        sa_conc=VAR2; %cmmc=cm^-3 
        am_conc=VAR1; %ppt
        
        ACDCinit(temp);                     %adjust ACDCinit
        ACDCsource(VAR2,ppt2cmmc(VAR1,temp)) %adjust ACDCsource
        [T,C,J_out]=ACDCrun();  
        
        results=[results J_out(end)];
    end
    resultss=[resultss; results];
end

figure()
loglog(VAR2range,resultss)
%loglog(cmmc2ppb(VAR2range,temp),resultss)
%title("SA-W nucleation")
xlabel("Concentration of H_2SO_4 (cm^{-3})")
ylabel("Nucleation rate (cm^{-3}s^{-1})")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%leave back from RUN folder
cd '../'
            
