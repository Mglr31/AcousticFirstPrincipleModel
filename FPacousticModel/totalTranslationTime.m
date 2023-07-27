function tau=totalTranslationTime(nmole,nmode,model)
%tau=totalTranslationTime(nmole,nmode,T)
%A function that compute the total translation time for the vibration mode nmode of the molecule nmole

T=model.T;
%% Loading parameters
load_physic_constant

if isstring(nmole)
    disp(strcat("totalTranslationTime: Reading table for",nmole))

    gas_p=readtable("gas_parametersSIModel.csv");
    vib_p=readtable("vibration_mode_parametersSIModel.csv");

    mode=vib_p(vib_p.MolName==nmole & vib_p.VibrationMode==nmode,:);
    mole=gas_p(gas_p.Name==nmole,:);
else
    mode=nmode;
    mole=nmole;
end
%%


%gas_p=readtable("gas_parametersSIModel.csv");
gas_p=model.gasT;
n=height(gas_p);
list_of_tau=zeros(1,n);
for i=1:n
    list_of_tau(i)= translationTime(mole,mode,gas_p(i,:),T);
end
somme=sum(gas_p.MolarFract./list_of_tau');
%somme=sum(1./list_of_tau');

tau=1/somme;
end