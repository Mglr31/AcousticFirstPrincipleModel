function Amixt=wilkeExpression(x,A,phi)
%w=wilkeExpression(x,A,phi)
%A function that compute the desired transport quantity A(thermal cond or
%viscosity) for a mixture of gas. x are the molar frac,phi a parameter of
%Wilke expression(Wilke,1950)
n=length(x);
Amixt=0;
for a=1:n
    Amixt=Amixt+x(a)*A(a)/sum(x.*phi(a,:));
end

end