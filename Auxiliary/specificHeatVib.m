function c=specificHeatVib(nmole,nmode,T,model)
%c=specificHeatVib(mol1,T)
%A function that return the specific heat of the vibrationnal mode "nmode"
%of the molecule "nmole" under the temperature T (K). c is in J/mol/K?????
%% Loading parameters
load_physic_constant

if isstring(nmole)
    gas_p=model.gasT;
    vib_p=model.vibT;

    mode=vib_p(vib_p.MolName==nmole & vib_p.VibrationMode==nmode,:);
    mole=gas_p(gas_p.Name==nmole,:);
else
    mode=nmode;
    mole=nmole;
end
%%
c=mode.g*R*(mode.theta/T)^2*exp(mode.theta/T)/(exp(mode.theta/T)-1)^2;
end