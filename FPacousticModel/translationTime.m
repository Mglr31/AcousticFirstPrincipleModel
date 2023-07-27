function tau_trans=translationTime(nmole1,nmode1,nmole2,T)
%tau_trans=translationTime(nmol1,nmode1,nmol2,nmode2,T)
%A function that compute the translationnal relaxation time tau_trans(s)
%between the nmode1 vibrationnal mode of the molecule nmole1 and the
%molecule nmole2.
%% Loading parameters
load_physic_constant

if isstring(nmole1)
vib_p=readtable("vibration_mode_parametersSI.csv");

mode1=vib_p(vib_p.MolName==nmole1 & vib_p.VibrationMode==nmode1,:);

else
    mode1=nmode1;
    
end
%%
inv_tau=collisionRate(nmole1,nmole2,T)*transition_probability(nmole1,nmode1,nmole2,"",T,'TransitionType',"V-T")...
    *(1-exp(-h*mode1.n/(kb*T)));

if nmole1.Name=="CO2"
    if nmole1.VibrationMode=="nu2" & nmole2.Name=="H2O"
        inv_tau=collisionRate(nmole1,nmole2,T)*3.78e-2...
    *(1-exp(-h*mode1.n/(kb*T)));
        transition_probability(nmole1,nmode1,nmole2,"",T,'TransitionType',"V-T")
    end
end
tau_trans=1/(inv_tau);

end