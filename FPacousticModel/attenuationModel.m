function result=attenuationModel(model,varargin)
%result=attenuationModel(model,varargin)
%A function that compute the attenuation in a mixture of gases with a
%first principle model based on the work of Richard Lueptow and Andi
%Petculescu(DOI: 10.1121/1.1828547, 10.1121/1.1352087).
%
% Inputs:
%   model: a structure containing the parameter of the simulation:
%    Required parameters:
%       model.listOfMole:       the list of molecule present in the model. Each
%                               molecule name must be repeated for each associated vibration mode.
%                               Example:["CH4"    "CH4"    "N2"    "H2O"]
%
%       model.listOfMode:       the list of vibration modes used to compute the
%                               relaxational attenuation. Example: ["nu2"    "nu4"    "nu"    "nu2"]
%
%       model.MolarFrac:        A list of the molar fraction of each
%                               molecule in the mixture. The molar fraction must be in the same
%                               order than model.listOfMole and each molar fraction appears once
%                               for the corresponding species.Example:[0.940.03 0.03]
%
%       model.T:                The temperature of the gas mixture in K.
%                               Ex:297
%
%       model.P:                The pressure of the gas mixture in Pa. Ex:
%                               101325
%
%      Optionnal parameters:
%       model.model_csv:        The name of a csv file containing a
%                               reference model. Example: "Tri_CH4_model_297.csv".
%
%       model.experiment_csv:   The name of a csv file containing
%                               experimental data for the model. Example: "Tri_CH4_exp_297.csv".
%
%   Other inputs(optionnal Name-Value pairs):
%       'FreqArray':            An array of the frequency for which the
%                               attenuation will be computed. Default value: logspace(-2 ,10,3000);
%
%
%   Outputs:
%       result: a structure containing the result of the model.
%       result's fields:
%           result.f:           An array of the frequency for which the
%                               attenuation was computed.
%
%           result.alpha:       The attenuation coefficient for each
%                               frequency in result.f.
%
%           result.alpha_c:     The classical attenuation coefficient for each
%                               frequency in result.f.
%
%           result.alpha_r:     The relaxational attenuation coefficient for each
%                               frequency in result.f.
%
%           result.alpha_rot:   The rotationnal attenuation coefficient for each
%                               frequency in result.f.
%
%           result.fmax:        The frequency for which the relaxational
%                               attenuation coeff is maximum.
%
%           result.c:           The sound of speed for each
%                               frequency in result.f.
%
%

%% Input parser
p = inputParser;
%Required input
addRequired(p,'model');
%Optionnal input
default_freqarray=logspace(-2 ,10,3000);
addParameter(p,'FreqArray',default_freqarray);

parse(p,model,varargin{:});


%% Complete the model and tables
[gas_t,vib_t]=makeModelTables(model);%create and save tables corresponding to the model molecule and modes parameters
model.gasT=gas_t;
model.vibT=vib_t;
model=addDerivativePArameter(model);%Compute secondary parameters(density,viscosity,cv,cp...) and add them to the model

%% Compute relaxational attenuation and sound of speed
s1=tic;
[A ,q]=getA(model);
t1=toc(s1);
f=p.Results.FreqArray;
omegas=2*pi*f;
s2=tic;
nombredonde=arrayfun(@(x) getK(x,A,q,model),omegas);
t2=toc(s2);
alpha_r=(imag(nombredonde));

c=omegas./real(nombredonde);
lambda=c./f;

[~, i]=max(alpha_r.*lambda);
fmax=f(i);
[TF,P] = islocalmax(alpha_r.*lambda);
fmaxs.freq=f(TF);
fmaxs.prominence=P(TF);


%% Compute the classical attenuation
s3=tic;
alpha_c=computeClassicalAtt(model,c,f);
t3=toc(s3);


%% Make result
alpha=alpha_r+alpha_c;
result.f=f;
result.alpha=alpha;
result.alpha_c=alpha_c;
result.alpha_r=alpha_r;
result.fmax=fmax;
result.fmaxs=fmaxs;
result.c=c;

result.MatrixComputationTime=t1;
result.FrequencyComputationTime=t2;
result.ClassicalAndRotComputationTime=t3;




end
