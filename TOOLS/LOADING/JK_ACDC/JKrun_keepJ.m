%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INTIALIZATION %%
init
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MANUAL %%
%1) in ACDC_input folder, create: input.inp, dG.txt (or dHdS.txt)
%2) in ACDC folder, adjust: ACDCinit, ACDCsource
%3) run 
%4) plot VAR1(y-axis = varry) and VAR2(x-axis) keeping constant J
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER AREA %%
temp=298.15; %Temperature

resultss=[];
Jrange=[10^15]; %different lines
for J=Jrange
    results=[];
    VAR2range=[10^8 10^9 10^10 10^11 10^12]; %x-axis
    for VAR2=VAR2range
        fprintf("->"+num2str(J)+"/"+num2str(VAR2));
        
        converged=0;factor=10;
        VAR1=0.000001; %has to be lower than resltant value %y-axis
        while converged==0
    
            ACDCinit(temp);                     %adjust ACDCinit
            ACDCsource(VAR2,rh2cmmc(VAR1,temp)) %adjust ACDCsource
            [T,C,J_out]=ACDCrun();
      
            %J test and adjustment
            if J_out(end)<J
                VAR1=VAR1*factor;
            else
                VAR1=VAR1/factor;
                factor=factor^.5;
                VAR1=VAR1*factor;
                if factor < 1.001
                    converged=1;
                end
            end
        end
        
        results=[results VAR1];
    end
    resultss=[resultss; results];
end

figure()
loglog(VAR2range,resultss)
title("SA-W nucleation")
xlabel("SA  conc [cm^{-3}]")
ylabel("W  conc [cm^{-3}]")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%leave back from RUN folder
cd '../'
            
