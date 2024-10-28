function D12=getD12(molname1,molname2)
%A function that return the diffusion coefficient D12 of teh couple of
%molecule molname1, molname2(string).
%VAlues are from Joseph Oakland Hirschfelder, Charles Francis Curtiss et
% Robert Byron Bird, Molecular Theory of Gases and Liquids,
% John Wiley and Sons, 1966 page 579
%D12 is in m2 s-1


name=strcat(molname1,molname2);

if name=="ArN2"|name=="N2Ar"
    D12=0.2;
elseif name=="ArCO2"|name=="CO2Ar"
    D12=0.14;
elseif name=="ArO2"|name=="O2Ar"
    D12=0.20;
elseif name=="N2O2"|name=="O2N2"
    D12=0.181;
elseif name=="N2CO2"|name=="CO2N2"
    D12=0.144;
elseif name=="CO2O2"|name=="O2CO2"
    D12=0.14;

else
    D12=0;
    disp(strcat("Not diffusion coefficient found for ",name))

end

D12=D12*1e-4;

end