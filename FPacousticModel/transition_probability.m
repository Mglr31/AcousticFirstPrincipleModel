function p=transition_probability(nmol1,nmode1,nmol2,nmode2,T,varargin)
%p=transition_probability(nmol1,nmode1,nmol2,nmode2,T,varargin)
%A function that compute the probability of a "V-V"(or a "V-T" or "T-V" if specified in the parameter)
% transition between two species at the temperature T.
%
%Inputs
%   nmol1: molecule name for the first species
%   nmode1: vibration mode name for the first species
%   nmol1: molecule name for the second species
%   nmode1: vibration mode name for the second species
%   T: temperature in K
%   Optionnal
%       'TransitionType',"V-V" "T-V" or "V-T": transition type.
%Outputs
%    p: The transition probability

%% Input Parser
p = inputParser;
%Required input
addRequired(p,'nmol1');
addRequired(p,'nmode1');
addRequired(p,'nmol2');
addRequired(p,'nmode2');
addRequired(p,'T');
%Optionnal input
default_transition_type="V-V";
addParameter(p,'TransitionType',default_transition_type);
parse(p,nmol1,nmode1,nmol2,nmode2,T,varargin{:});
%% Loading parameters
load_physic_constant

if isstring(nmol1)
    disp(strcat("transition_probability: Reading table for",nmol1))

    gas_p=readtable("gas_parametersSI.csv");
    vib_p=readtable("vibration_mode_parametersSI.csv");

    mode1=vib_p(vib_p.MolName==nmol1 & vib_p.VibrationMode==nmode1,:);
    mol1=gas_p(gas_p.Name==nmol1,:);


else
    mode1=nmode1;
    mol1=nmol1;

end

if isstring(nmol2)
    disp(strcat("transition_probability: Reading table for",nmol2))

    gas_p=readtable("gas_parametersSI.csv");
    vib_p=readtable("vibration_mode_parametersSI.csv");



    mode2=vib_p(string(vib_p.MolName)==nmol2 & string(vib_p.VibrationMode)==nmode2,:);
    mol2=gas_p(gas_p.Name==nmol2,:);
else

    mode2=nmode2;
    mol2=nmol2;
end
%% Computing LJ fit
[rc ,alpha ,Estar]=fitLJ(nmol1,nmode1,nmol2,nmode2,T,'TransitionType',p.Results.TransitionType);




%% Computing pairwise parameters
mu=mol1.mass*mol2.mass/(mol1.mass+mol2.mass);% reduced mass
epsilon=sqrt(mol1.epsilon*mol2.epsilon);% pairwise potential depth
coeff=getPolarCoeff(mol1.Name,mol2.Name);
epsilon=epsilon*coeff;
sigma=(mol1.sigma+mol2.sigma)/2;% pairwise collision diameter




%% Getting the species parameters
C=(mol1.SutherlandC+mol2.SutherlandC)/2;%Sutherland constant Averaging???: +mol2.SutherlandC)/2;
ksi=Estar/(kb*T);


if p.Results.TransitionType=="V-V"
    [rc1 ,alpha1 ,Estar1]=fitLJ(nmol1,nmode1,nmol1,nmode1,T,'TransitionType',p.Results.TransitionType);
    [rc2 ,alpha2 ,Estar2]=fitLJ(nmol2,nmode2,nmol2,nmode2,T,'TransitionType',p.Results.TransitionType);

    Vsquare1=(alpha1)^2*mode1.VibAmplCoeff*h_barre/2/mode1.omega;
    Vsquare2=(alpha2)^2*mode2.VibAmplCoeff*h_barre/2/mode2.omega;
    DeltaE=h_barre*(mode1.omega-mode2.omega);
elseif p.Results.TransitionType=="V-T"
    [rc1 ,alpha1 ,Estar1]=fitLJ(nmol1,nmode1,nmol1,nmode1,T,'TransitionType',p.Results.TransitionType);

    Vsquare1=alpha1^2*mode1.VibAmplCoeff*h_barre/2/mode1.omega;
    Vsquare2=1;
    DeltaE=h_barre*(mode1.omega);

elseif p.Results.TransitionType=="T-V"
    [rc2 ,alpha2 ,Estar2]=fitLJ(nmol2,nmode2,nmol2,nmode2,T,'TransitionType',p.Results.TransitionType);

    Vsquare1=1;
    Vsquare2=alpha2^2*mode2.VibAmplCoeff*h_barre/2/mode2.omega;
    DeltaE=h_barre*(-mode2.omega);

else
    disp('Error: possible transition type are "V-V" "T-V" and "V-T"' )
    return
end

ksi_alpha=(DeltaE^2*mu*pi*pi/(2*alpha^2*h_barre^2*kb*T))^(1/3);

if isstring(mode2)%V-T transition, a P0 is needed even if this dos not make sense because P0 is characteristic of a vibration mode
    clear mode2
    mode2.P0=2/3;%Default value?
end



p=mode1.P0*mode2.P0*1.364/(1+(C/T))*(rc/sigma)^2*Vsquare1*Vsquare2*8*sqrt(pi/3)*(2*pi*mu*DeltaE/(alpha^2*h_barre^2))^2 ...
    *sqrt(ksi)*exp(-3*ksi+DeltaE/(2*kb*T)+epsilon/(kb*T));

end