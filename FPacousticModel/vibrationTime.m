function tau_vib=vibrationTime(nmole1,nmode1,nmole2,nmode2,T)
%tau_trans=vibrationTime(nmol1,nmode1,nmol2,nmode2,T)
%A function that compute the vibrationnal relaxation time tau_vib(s)
%between the nmode1 vibrationnal mode of the molecule nmole1 and the
%nmode2 mode of the molecule nmole2.
%% Loading parameters
load_physic_constant

if isstring(nmole2)
    disp(strcat("vibrationTime: Reading table for",nmole2))
    gas_p=readtable("gas_parametersSIModel.csv");
    vib_p=readtable("vibration_mode_parametersSIModel.csv");

    mode2=vib_p(string(vib_p.MolName)==nmole2 & string(vib_p.VibrationMode)==nmode2,:);
    mol2=gas_p(gas_p.Name==nmole2,:);
else
    mode2=nmode2;
    mol2=nmole2;
end

%%
inv_tau=mol2.MolarFract*mode2.g*collisionRate(nmole1,nmole2,T)*transition_probability(nmole1,nmode1,nmole2,nmode2,T);
%inv_tau=mode2.g*collisionRate(nmole1,nmole2,T)*transition_probability(nmole1,nmode1,nmole2,nmode2,T);

tau_vib=1/(inv_tau);

end