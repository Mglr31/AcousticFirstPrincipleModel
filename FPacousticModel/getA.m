function [A,q]=getA(model)
%A function that compute the A matrix for the molecule, vib mode pressure and temp of model
molList=model.listOfMole;
modList=model.listOfMode;

T=model.T;

load_physic_constant
vib_p=model.vibT;%table with all the vibration mode ppties
gasT=model.gasT;%table with all the gas properties


%Making the table of all the vibration mode with all the ppties(mol and
%mode related)
vib_p.Properties.VariableNames{1} = 'Name';
vib_p = removevars(vib_p, 'MolarFract');
mol_tab=join(vib_p,gasT,'Keys','Name');
mol_tab(mol_tab.VibrationMode=="none",:)=[];%removing the molecules that do not vibrate
n=height(mol_tab);
A=zeros(n,n);
q=zeros(n,1);
tvibs=45*ones(n,n);
 for i=1:n
     for j=1:n
         tvibs(i,j)=vibrationTime(mol_tab(i,:),mol_tab(i,:),mol_tab(j,:),mol_tab(j,:),T);
     end
 end
ttransls=45*ones(n,1);
for j=1:n
    ttransls(j)=totalTranslationTime(mol_tab(j,:),mol_tab(j,:),model);
end
ns=mol_tab.n;

%filling A and q
for j=1:n
    tau_vib=tvibs(j,:);
    q(j)=1/ttransls(j)+sum(tau_vib.^-1.*(1-exp(-h*ns(j)/(kb*T)))./(1-exp(-h*ns'/kb/T)).*(1-ns/ns(j))');
    for k=1:n
        if j==k
            A(j,k)=1/ttransls(j)+sum(tau_vib.^-1.*(1-exp(-h*ns(j)/(kb*T)))./(1-exp(-h*ns'/kb/T)));
        else
            A(j,k)=-1/tau_vib(k)*(1-exp(-h*ns(j)/(kb*T)))/(1-exp(-h*ns(k)/kb/T))*ns(k)/ns(j);

        end
    end
end

 
 
