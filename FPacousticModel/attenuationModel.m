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
%       'PlotOption':           Default:'none', nothing will be plotted.
%                               Possible values ["alpha_r","alpha_c","alpha_rot","c","alpha","a","all"] for plotting
%                               respectively the relaxation attenuation,
%                               the classical attenuation,the rotationnal
%                               attenuation
%                               the speed of sound, the total attenuation and all attenuation. Example:
%                               'PlotOption',["alpha_r","c"].
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
default_plotoption="none";
addParameter(p,'PlotOption',default_plotoption);

parse(p,model,varargin{:});

%% Check the current directory
% list=ls;
% list=deblank(string(list));
% flag=sum(list=="Tables");
% if flag==0
%     %disp("Error: the function must be executed in a directory that contains a folder named Tables with all the gas and vibration mode tables in it")
%     %return
%     dir=pwd;
%     cd("D:\m.gillier\Documents\MATLAB\Model")
% end

%% Complete the model and tables
[gas_t,vib_t]=makeModelTables(model);%create and save tables corresponding to the model molecule and modes parameters
model.gasT=gas_t;
model.vibT=vib_t;
model=addDerivativePArameter(model);%Compute secondary parameters(density,viscosity,cv,cp...) and add them to the model
%% Debug water vapor
l=gas_t(gas_t.Name=="H2O",:);
result.Nwater=l.N;




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


%% Compute the classical and rotationnal attenuation(and diffusion)
s3=tic;
alpha_c=computeClassicalAtt(model,c,f);
%ccoeff=alpha_c./f.^2.*c.^3
alpha_rot=computeRotationnalAtt(model,c,f);
%diff=computeDiffusionAttenuation(model,c,f);
%alpha_diff=diff.alpha;
%diff.name

t3=toc(s3);

%% Compute the dust attenuation

if isfield(model,'dust')
    a2pilambda=getDustAtt(model.dust,f,model);
    a_dust=a2pilambda*2*pi./lambda;
else
    a_dust=[];
end
%% Make result
alpha=alpha_r+alpha_c+alpha_rot;
result.f=f;
result.alpha=alpha;
result.alpha_c=alpha_c;
result.alpha_r=alpha_r;
result.alpha_rot=alpha_rot;
result.alpha_dust=a_dust;
result.fmax=fmax;
result.fmaxs=fmaxs;
result.c=c;

result.MatrixComputationTime=t1;
result.FrequencyComputationTime=t2;
result.ClassicalAndRotComputationTime=t3;


%% Getting back to original directory
if flag==0

    cd(dir)
end
%% Plotting(if required)
if p.Results.PlotOption ~="none"% if plotting was required
    %% Check for reference model and experimental data
    model_flag=0;
    exp_flag=0;
    if isfield(model,'model_csv')
        model_flag=1;
        tab=readtable(model.model_csv,'Delimiter',';','DecimalSeparator',',');
        f_ref=tab.Var1;
        alambda_ref=tab.Var2;
    end

    if isfield(model,'experiment_csv')
        exp_flag=1;
        tab2=readtable(model.experiment_csv,'Delimiter',';','DecimalSeparator',',');
        f_exp=tab2.Var1;
        al_exp=tab2.Var2;
    end
    n=length(p.Results.PlotOption);
    titre="";
    u=unique(model.listOfMole,"stable");
    for k=1:length(u)
        titre=strcat(titre,u(k),"(",string(model.MolarFrac(k)),"), ");

    end
    titre=strcat(titre," T=",string(model.T),", P=",string(model.P));
    fp=f/(model.P/101325);
    for i=1:n
        if p.Results.PlotOption(i)=="none"
            disp("No plotting")
        elseif p.Results.PlotOption(i)=="alpha_r"
            figure
            semilogx(fp,alpha_r.*lambda,'DisplayName',"Model",'LineWidth',2)
            if exp_flag
                hold on
                semilogx(f_exp,al_exp,'+','DisplayName',"Experiment",'LineWidth',2,'MarkerSize',15)
            end
            if model_flag
                hold on
                semilogx(f_ref,alambda_ref,'--','DisplayName',"Reference Model(DL)",'LineWidth',2)
            end
            legend show
            grid on
            xlabel("Frequency/Pressure(Hz/atm)");
            ylabel("Alpha\_relax*Lambda")
            set(gca,'FontSize',22)
            title(titre)
        elseif p.Results.PlotOption(i)=="alpha_c"
            figure
            semilogx(f,alpha_c,'DisplayName',"Model",'LineWidth',2)
            legend show
            grid on
            xlabel("Frequency(Hz)");
            ylabel("Alpha classical(m-1)")
            set(gca,'FontSize',22)
            title(titre)

        elseif p.Results.PlotOption(i)=="alpha"
            figure
            semilogx(fp,alpha.*lambda,'DisplayName',"Model",'LineWidth',2)
            if exp_flag
                hold on
                semilogx(f_exp,al_exp,'+','DisplayName',"Experiment",'LineWidth',2,'MarkerSize',15)
            end
            if model_flag
                hold on
                semilogx(f_ref,alambda_ref,'--','DisplayName',"Reference Model(DL)",'LineWidth',2)
            end
            legend show
            grid on
            xlabel("Frequency/Pressure(Hz/atm)");
            ylabel("Alpha*Lambda")
            set(gca,'FontSize',22)
            title(titre)

        elseif p.Results.PlotOption(i)=="c"
            figure
            semilogx(f,c,'DisplayName',"Model",'LineWidth',2)
            legend show
            grid on
            xlabel("Frequency(Hz)");
            ylabel("Speed of Sound(m/s)")
            set(gca,'FontSize',22)
            title(titre)


        elseif p.Results.PlotOption(i)=="a"
            figure
            loglog(f,alpha,'DisplayName',"Model",'LineWidth',2)
            if exp_flag
                hold on
                loglog(f_exp,al_exp,'+','DisplayName',"Experiment",'LineWidth',2,'MarkerSize',15)
            end
            if model_flag
                hold on
                loglog(f_ref,alambda_ref,'--','DisplayName',"Reference Model(DL)",'LineWidth',2)
            end
            legend show
            grid on
            xlabel("Frequency(Hz)");
            ylabel("\alpha(m^{-1})")
            set(gca,'FontSize',22)
            title(titre)

        elseif p.Results.PlotOption(i)=="alpha_rot"
            figure
            loglog(f,alpha_rot,'DisplayName',"Model",'LineWidth',2)
            if exp_flag
                hold on
                loglog(f_exp,al_exp,'+','DisplayName',"Experiment",'LineWidth',2,'MarkerSize',15)
            end
            if model_flag
                hold on
                loglog(f_ref,alambda_ref,'--','DisplayName',"Reference Model(DL)",'LineWidth',2)
            end
            legend show
            grid on
            xlabel("Frequency(Hz)");
            ylabel("Alpha(m-1)")
            set(gca,'FontSize',22)
            title(titre)

        elseif p.Results.PlotOption(i)=="all"
            figure
            loglog(f,alpha_rot,'DisplayName',"alpha\_rot",'LineWidth',2)
            hold on
            loglog(f,alpha_r,'DisplayName',"alpha\_relax",'LineWidth',2)
            hold on
            loglog(f,alpha_c,'DisplayName',"alpha\_classical",'LineWidth',2)
            hold on
       %     loglog(f,alpha_diff,'DisplayName',"alpha\_diff",'LineWidth',2)
            hold on
            loglog(f,alpha,'DisplayName',"alpha\_total",'LineWidth',2)
            if exp_flag
                hold on
                %loglog(f_exp,al_exp,'+','DisplayName',"Experiment",'LineWidth',2,'MarkerSize',15)
            end
            if model_flag
                hold on
               % loglog(f_ref,alambda_ref,'--','DisplayName',"Reference Model(DL)",'LineWidth',2)
            end
            legend show
            grid on
            xlabel("Frequency(Hz)");
            ylabel("Alpha(m^{-1})")
            set(gca,'FontSize',22)
            title(titre)
        end
    end
end
end
