%% Environment definition

load_physic_constant%loading a bunch of physics constants(kb,R,h...)

clear model

model.listOfMole=["N2","CO2","CO2","CO2"];
model.listOfMode=["nu","nu1","nu2","nu3"];
model.MolarFrac=[0.02 0.98] ;
model.T=297;
model.P=101325; %1 atm

%%
model=modelTri;
makeModelTables(model);%create and save tables corresponding to the model molecule and modes parameters

model=addDerivativePArameter(model)
%% Computation of the relaxation equation parameter
tic
[A ,q]=getA(model);
toc
%% Wave number computation
f=logspace(-2 ,10,3000);

omegas=2*pi*f;
tic
nombredonde=arrayfun(@(x) getK(x,A,q,model),omegas);
gammas=arrayfun(@(x) sum(getGamma(A,q,x))/4,omegas);
toc
r=real(gammas);
alpha=(imag(nombredonde));
c=omegas./real(nombredonde);
lambda=c./f;
[m i]=max(alpha.*lambda);
fmax=f(i)
%% Plotting
tab=readtable("CO2_model_air.csv",'Delimiter',';','DecimalSeparator',',');
tab2=readtable("CO2_exp_air.csv",'Delimiter',';','DecimalSeparator',',');
tab2=readtable("N2_exp_T_297.csv",'Delimiter',';','DecimalSeparator',',');

f_ref=tab.Var1;
alambda_ref=tab.Var2;
f_exp=tab2.Var1;
al_exp=tab2.Var2;

figure
semilogx(f,alpha.*lambda,'DisplayName',"Model",'LineWidth',2)
hold on

semilogx(f_exp,al_exp,'+','DisplayName',"Experiment",'LineWidth',2,'MarkerSize',15)
hold on
semilogx(f_ref,alambda_ref,'--','DisplayName',"Reference Model(DL)",'LineWidth',2)

legend show
grid on
xlabel("Frequency/Pressure(Hz/atm)");
ylabel("Alpha*Lambda")
set(gca,'FontSize',22)
%xlim([1e4 1e6])
%%
figure
semilogx(f,c)