function[rc ,alpha ,Estar,percentage_diff]=fitLJ(nmol1,nmode1,nmol2,nmode2,T,varargin)
%[rc alpha Estar]=fitLJ(sigma1,epsilon1,sigma2,epsilon2)
%A function that fit an exponential to the Lennard Jone potential between
%species 1 and 2
%Inputs
%   nmol1: molecule name for the first species(or the corresponding model table
%   line)
%   nmode1: vibration mode name for the first species(or the corresponding model table
%   line)
%   nmol1: molecule name for the second species(or the corresponding model table
%   line)
%   nmode1: vibration mode name for the second species(or the corresponding model table
%   line)
%   T: temperature in K
%   Optionnal
%       'TransitionType',"V-V" "T-V" or "V-T": transition type.
%Outputs
%    rc: classical turning point [m]
%    alpha: decay parameter [m-1]
%    Estar: incident kinetic energy [J]
%The fit potential is given by Vexp=(epsilon+Estar)exp(alpha(rc-r))-epsilon

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
%Loading parameters
load_physic_constant

if isstring(nmol1)
    gas_p=readtable("gas_parametersSI.csv");
    vib_p=readtable("vibration_mode_parametersSI.csv");

    mode1=vib_p(vib_p.MolName==nmol1 & vib_p.VibrationMode==nmode1,:);
    mol1=gas_p(gas_p.Name==nmol1,:);

    
else
    mode1=nmode1;
    mol1=nmol1;
    
end

if isstring(nmol2)
    gas_p=readtable("gas_parametersSI.csv");
    vib_p=readtable("vibration_mode_parametersSI.csv");

   

    mode2=vib_p(string(vib_p.MolName)==nmol2 & string(vib_p.VibrationMode)==nmode2,:);
    mol2=gas_p(gas_p.Name==nmol2,:);
else
    
    mode2=nmode2;
    mol2=nmol2;
end
%% Computing pairwise parameters
mu=mol1.mass*mol2.mass/(mol1.mass+mol2.mass);% reduced mass
epsilon=sqrt(mol1.epsilon*mol2.epsilon);% pairwise potential depth
sigma=(mol1.sigma+mol2.sigma)/2;% pairwise collision diameter

%% Finding a good (alpha,rc,E) triplet
%First try
RC=linspace(0.1*sigma,1.5*sigma,100000);%array of all possible rc
ESTAR=4*epsilon*((sigma*(RC.^-1)).^12-(sigma*(RC.^-1)).^6);%Corresponding  kinetic energies
ALPHA=log(epsilon*(epsilon+ESTAR).^-1)./(RC-sigma);%Corresponding decay parameters

% if p.Results.TransitionType=="V-V"
%     BB=kb*T*(2*pi^4*mu*(mode1.omega-mode2.omega)^2/(kb*T))^(1/3);
% elseif p.Results.TransitionType=="V-T"
%     BB=kb*T*(2*pi^4*mu*(mode1.omega)^2/(kb*T))^(1/3);
% elseif p.Results.TransitionType=="T-V"
%     BB=kb*T*(2*pi^4*mu*(-mode2.omega)^2/(kb*T))^(1/3);
% else
%     disp('Error: possible transition type are "V-V" "T-V" and "V-T"' )
%     return
% end
if p.Results.TransitionType=="V-V"
    BB=kb*T*(pi^2*mu*(mode1.omega-mode2.omega)^2/(2*kb*T))^(1/3);
elseif p.Results.TransitionType=="V-T"
    BB=kb*T*(pi^2*mu*(mode1.omega)^2/(2*kb*T))^(1/3);
elseif p.Results.TransitionType=="T-V"
    BB=kb*T*(pi^2*mu*(mode2.omega)^2/(2*kb*T))^(1/3);
else
    disp('Error: possible transition type are "V-V" "T-V" and "V-T"' )
    return
end
ESTARalpha=BB*ALPHA.^(-2/3);%ESTAR computed from the decay parameter alpha;

ESTARdiff=abs((ESTAR-ESTARalpha));%./ESTAR).^2;

[minimum ,iapprox]=min(ESTARdiff);
Estar_approx=ESTAR(iapprox);

%Second search around this first appoximative value of rc
RCsmall=linspace(RC(iapprox-1),RC(iapprox+1),100000);%array of all possible rc
ESTAR=4*epsilon*((sigma*(RCsmall.^-1)).^12-(sigma*(RCsmall.^-1)).^6);%Corresponding  kinetic energies
ALPHA=log(epsilon*(epsilon+ESTAR).^-1)./(RCsmall-sigma);%Corresponding decay parameters
ESTARalpha=BB*ALPHA.^(-2/3);%ESTAR computed from the decay parameter alpha;

ESTARdiff=abs((ESTAR-ESTARalpha));%./ESTAR).^2;

[minimum ,i]=min(ESTARdiff);

rc=RCsmall(i);
alpha=ALPHA(i);
Estar=ESTARalpha(i);
percentage_diff=ESTARdiff(i)/Estar*100;
end