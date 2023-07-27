function Gamma=getGamma(A,q,omega)
%gamma=getGAmma(A,q,omega)
%A function that return the vector gamma from the matrix a and the vector q
%and the pulsation omega
%A=R(n,n), q=R(n),Gamma=C(n)
n=length(q);
Id=eye(n);
Gamma=(A+1i*omega*Id)\(q);

end