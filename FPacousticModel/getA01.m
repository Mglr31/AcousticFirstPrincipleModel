function [A,q]=getA01(model)
%A function that compute the A matrix for the molecule, vib mode pressure and temp of model
molList=model.listOfMole;
modList=model.listOfMode;
do_not_vib=(modList=="none");%get the index of the species that are not vibrating
molList(do_not_vib)=[];%getting ride of this molecule and mode
modList(do_not_vib)=[];

T=model.T;
n=length(molList);
A=zeros(n,n);
q=zeros(n,1);
load_physic_constant
%vib_p=readtable("vibration_mode_parametersSIModel.csv");
vib_p=model.vibT;
vib_p(do_not_vib,:)=[];
gasT=model.gasT;%table with all the gas properties
list_n=vib_p.n;
for j=1:n
    molj=molList(j);
    modj=modList(j);
    


    tau_vib_j=arrayfun(@(x,y) vibrationTime(molj,modj,x,y,T),molList,modList);
    nj=vib_p(vib_p.MolName==molj & vib_p.VibrationMode==modj,:).n;

    other_tau_vib=tau_vib_j;
    other_tau_vib(j)=[];
    list_nj=list_n;
    list_nj(j)=[];
    q(j)=1/totalTranslationTime(molj,modj,T)+sum(other_tau_vib.^-1.*(1-exp(-h*nj/(kb*T)))./(1-exp(-h*list_nj'/kb/T)).*(1-list_nj/nj)');


    for k=1:n

        molk=molList(k);
        modk=modList(k);
        nk=vib_p(vib_p.MolName==molk & vib_p.VibrationMode==modk,:).n;


        if j==k
            other_tau_vib=tau_vib_j;
            other_tau_vib(j)=[];
            list_nj=list_n;
            list_nj(j)=[];
            A(j,k)=1/totalTranslationTime(molj,modj,T)+sum(other_tau_vib.^-1.*(1-exp(-h*nj/(kb*T)))./(1-exp(-h*list_nj'/kb/T)));
        else
            A(j,k)=-1/tau_vib_j(k)*(1-exp(-h*nj/(kb*T)))/(1-exp(-h*list_n(k)/kb/T))*list_n(k)/nj;
        end
    end

end