function k=getK(omega,A,q,model)
%A function that compute the dispersion relation k=f(omega) for the
%parametrs of model. A and q must have benn computed using [A,q]=getA(model)

Gamma=getGamma(A,q,omega);

x=model.MolarFractMode;
cvib=model.cvib;

modList=model.listOfMode;
do_not_vib=(modList=="none");%get the index of the species that are not vibrating
x(do_not_vib)=[];%getting ride of the associated molar frac
cvib(do_not_vib)=[];%getting ride of the associated cvib

squareroot=sqrt(model.rho/model.P*(model.cv+sum(x.*cvib.*(Gamma-1)'))/(model.cp+sum(x.*cvib.*(Gamma-1)')));
if imag(squareroot)<0
    squareroot=-squareroot;
end
k=omega*squareroot;
end