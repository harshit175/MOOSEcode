%this file calculated input parameters for MOOSE 

% DCPD PDMS Wire  -- order of entries

clc; clear

rho = [980 1030 8920 8050];

k = [0.152 .27 401 50.2];

Cp = [1440 1100 386 490];

Hr = 460e3;

%% 

L = 8e-3; %Characteristic length

Avec = 10.^[4.409 5.281 6.649];

Evec = [48.76 51.1 58.16]*1e3;

nvec = [1.616 1.927 2.192];

Kcatvec = [0.436 0.365 0.370];

i=2;  %i=2 is medium

A =  Avec(i);  %from kessler's paper

E = Evec(i);

n= nvec(i);

Kcat = Kcatvec(i);

Ttrig =  339.44;  %Ignition temperatue in degree C

Tinit = 20;  %Intial Temperature in degree C

%% 

%DCPD-PDMS-WIRE

Tconductivity = k./rho(1)/Cp(1)/A/L^2;

Hrinput(1) = Hr/Cp(1)/(Ttrig-Tinit);

Hrinput(2:4)=0;

rhocp = rho.*Cp./(rho(1)*Cp(1)); %goes in specific heat

%% Print
i=1;
fid = fopen('matprop.txt','wt');

fprintf(fid, 'Order is DCPD, PDMS, Copper Wire, Steel wire \n \n');

fprintf(fid, 'Entries for DCPDnonDgeneral \n');
fprintf(fid,'Ttrig = %1.2f \n Tintl = %1.2f \n',Ttrig,Tinit);
fprintf(fid,'_E = %1.4d \n _n = %1.4f \n _Kcat = %1.4f \n\n',E,n,Kcat);

for i=1:4
    fprintf(fid,'prop_names = specific_heat density Hr TConductivity \n');
    fprintf(fid,'prop_values = %1.1f %1.4f %1.4f %1.4d  # %1.0f \n \n', 1, rhocp(i),Hrinput(i), Tconductivity(i),i);
end




