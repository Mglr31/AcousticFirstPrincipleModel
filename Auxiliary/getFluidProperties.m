function [cv,cp,eta,kappa]=getFluidProperties(P,T,molname)
%A function that get some physical properties of fluids from the webbokk
%website(https://webbook.nist.gov)
% Inputs:
%   P:          the pressure at which the fluid properties are valid(can be a scalar or a
%       vector) [PA]
%   T:          the temperature at which the fluid properties are valid(can be a scalar or a
%       vector)[K]
%   molname:       the name of the considered fluid("O2","CH4"...)
%
% Outputs:
%   cv:         the isochoric heat capacity of the fluid at the requested P
%               and T(can be a scalar, a vector or a matrix) [J/mol*K]
%   cp:         the isobaric heat capacity of the fluid at the requested P
%               and T(can be a scalar, a vector or a matrix)[J/mol*K]
%   eta:        the viscosity of the fluid at the requested P
%               and T(can be a scalar, a vector or a matrix)[Pa*s]
%   kappa:      the thermal conductivity of the fluid at the requested P
%               and T(can be a scalar, a vector or a matrix)[W/m*K]
%
%
%   If P and T are of size m and n, the result will be a matrix of size
%   (m+1,n+1) with the first column being the pressure points in Pa and thz
%   first line the temperature point in K.

%Unit conversion
P=1e-6*P;
%URL chunks
requestIsoT1="https://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=";
requestIsoT2="&Type=IsoTherm&Digits=5&PLow=";
requestIsoT3="&RefState=DEF&TUnit=K&PUnit=MPa&DUnit=mol%2Fm3&HUnit=kJ%2Fmol&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm";

requestIsoP2="&Type=IsoBar&Digits=5&P=";

%Some flags
scalar=0;
IsoT=0;
IsoP=0;
matrix=0;

% Get molecule ID
WebbookID=readtable('WebbookID.csv');
line=WebbookID(WebbookID.MolName==molname,:);
ID=line.ID;
if isempty(line)
    disp("Error: this molecule is not available")
    return
end
if max(size(P))==1&max(size(T))==1 %Request at a single (P,T) point, outputs are scalar
    scalar=1;
    request=strcat(requestIsoT1,ID,requestIsoT2,string(P),"&PHigh=",string(P),"&PInc=0&T=",string(T),requestIsoT3);
    try
        out_table = getTableFromWeb_mod(request,1);
    catch
        if molname=="CO2"
            disp("CO2: not available in the webbok, using extrapolations")
            [cv,cp,eta,kappa]=getCO2Properties(T);
        elseif molname=="H2O"
            disp("H2O: not available in the webbok, using extrapolations")

            [cv,cp,eta,kappa]=getH2OProperties(T);
        else
            disp(strcat(molname," : the physical properties for the required P and T were not in webbook. The default values were used"))
            gas_p=readtable("gas_parametersSI.csv");
            line=gas_p(gas_p.Name==molname,:);
            cv=line.cv;
            cp=line.cp;
            eta=line.etha;
            kappa=line.kappa;
        end
        return
    end
    cv=double(string(out_table{2,8}));
    cp=double(string(out_table{2,9}));
    eta=double(string(out_table{2,12}));
    kappa=double(string(out_table{2,13}));

end

if max(size(P))>1&max(size(T))==1 %Request at a range of P with one T, outputs are colum vectors
    IsoT=1;
    request=strcat(requestIsoT1,ID,requestIsoT2,string(P(1)),"&PHigh=",string(P(end)),"&PInc=",string(P(2)-P(1)),"&T=",string(T),requestIsoT3);
    out_table = getTableFromWeb_mod(request,1);
    table=string(out_table);
    cv=double(table(2:end,[2,8]));
    cp=double(table(2:end,[2,9]));
    eta=double(table(2:end,[2,12]));
    kappa=double(table(2:end,[2,13]));
    %Conversion from MPa to Pa
    cv(:,1)=cv(:,1)*1e6;
    cp(:,1)=cp(:,1)*1e6;
    eta(:,1)=eta(:,1)*1e6;
    kappa(:,1)=kappa(:,1)*1e6;
end

if max(size(P))==1&max(size(T))>1 %Request at a range of T with one P, outputs are line vectors
    IsoP=1;
    request=strcat(requestIsoT1,ID,requestIsoP2,string(P),"&THigh=",string(T(end)),"&TLow=",string(T(1)),"&TInc=",string(T(2)-T(1)),requestIsoT3);
    out_table = getTableFromWeb_mod(request,1);
    table=string(out_table);
    cv=double(table(2:end,[1,8]))';
    cp=double(table(2:end,[1,9]))';
    eta=double(table(2:end,[1,12]))';
    kappa=double(table(2:end,[1,13]))';

end

if max(size(P))>1&max(size(T))>1 %Request at a range of T and P, outputs are matrices
    %Get the first column(P vector)
    request=strcat(requestIsoT1,ID,requestIsoT2,string(P(1)),"&PHigh=",string(P(end)),"&PInc=",string(P(2)-P(1)),"&T=",string(T(1)),requestIsoT3);
    out_table = getTableFromWeb_mod(request,1);
    table=string(out_table);
    cv=double(table(2:end,2));
    cp=double(table(2:end,2));
    eta=double(table(2:end,2));
    kappa=double(table(2:end,2));
    cv=[NaN;cv];
    cp=[NaN;cp];
    eta=[NaN;eta];
    kappa=[NaN;kappa];
    n=length(T);
    for i=1:n
        request=strcat(requestIsoT1,ID,requestIsoT2,string(P(1)),"&PHigh=",string(P(end)),"&PInc=",string(P(2)-P(1)),"&T=",string(T(i)),requestIsoT3);
        out_table = getTableFromWeb_mod(request,1);
        table=string(out_table);
        cv=[cv [T(i);double(table(2:end,8))]];
        cp=[cp [T(i);double(table(2:end,9))]];
        eta=[eta [T(i);double(table(2:end,12))]];
        kappa=[kappa [T(i);double(table(2:end,13))]];
    end

end





end