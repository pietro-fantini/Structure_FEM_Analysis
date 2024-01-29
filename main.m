%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mechanical System Dynamics
% FEM script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; 
close all; 
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load input file and assemble structure
[file_i,xy,nnod,sizee,idb,ndof,incid,l,gamma,m,EA,EJ,T,posit,nbeam,pr]=loadstructure;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Draw structure
dis_stru(posit,l,gamma,xy,pr,idb,ndof);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assemble mass and stiffness matricies
[M,K] = assem(incid,l,m,EA,EJ,T,gamma,idb);

%kx1 contribution
i_ndof1(1) = idb(1,1);
kx1 = 5e4; %N/m
K(i_ndof1(1),i_ndof1(1))=K(i_ndof1(1),i_ndof1(1)) + kx1;

%ky1 contribution
i_ndof1(2) = idb(1,2);
ky1 = 5e4; %N/m
K(i_ndof1(2),i_ndof1(2))=K(i_ndof1(2),i_ndof1(2)) + ky1;

%kx8 contribution
i_ndof8(1) = idb(8,1);
kx8 = 5e4; %N/m
K(i_ndof8(1),i_ndof8(1))=K(i_ndof8(1),i_ndof8(1)) + kx8;

%ky8 contribution
i_ndof8(2) = idb(8,2);
ky8 = 5e4; %N/m
K(i_ndof8(2),i_ndof8(2))=K(i_ndof8(2),i_ndof8(2)) + ky8;

%contribution of lumped mass
i_ndof12=idb(12,2);
Mc=5; %kg
M(i_ndof12,i_ndof12)=M(i_ndof12,i_ndof12) + Mc;

%contribution of internal spring kc
i_ndof11=idb(11,2);
i_ndof12=idb(12,2);
i_dofkc= [i_ndof11 i_ndof12];
kc=1e5; %N/m
K_k2=[kc -kc; -kc kc];
K(i_dofkc,i_dofkc)=K(i_dofkc,i_dofkc) + K_k2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute natural frequencies and mode shapes
MFF = M(1:ndof,1:ndof);
MCF = M(ndof+1:end,1:ndof);
MFC = M(1:ndof,ndof+1:end);
MCC = M(ndof+1:end,ndof+1:end);

KFF = K(1:ndof,1:ndof);
KCF = K(ndof+1:end,1:ndof);
KFC = K(1:ndof,ndof+1:end);
KCC = K(ndof+1:end,ndof+1:end);

[modes, omega] = eig(MFF\KFF);
omega = sqrt(diag(omega));
[omega,i_omega] = sort(omega);
freq0 = omega/2/pi;
modes = modes(:,i_omega);

nmodes = 4;
scale_factor = 1;
for ii=1:nmodes
    mode = modes(:,ii);
    figure();
    diseg2(mode,scale_factor,incid,l,gamma,posit,idb,xy)
    title(['mode ',num2str(ii),' freq ',num2str(freq0(ii)),' Hz']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add contribution of internal damping 
alpha = 1.5; %s^-1
beta = 9*10^-5; %s
R = alpha*M + beta*K;
%add concetrated damper 
r = 50; %Ns/m
R_r = [r -r; -r r]; 
i_ndof11 = idb(11,2);
i_ndof12 = idb(12,2);
i_dofr = [i_ndof11 i_ndof12];
R(i_dofr,i_dofr) = R(i_dofr,i_dofr) + R_r;
%partitioning of R
RFF = R(1:ndof,1:ndof);
RCF = R(ndof+1:end,1:ndof);
RFC = R(1:ndof,ndof+1:end);
RCC = R(ndof+1:end,ndof+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%FRF displacement
F0 = zeros(ndof,1);
F0(idb(11,2)) = 1;
om = (0:0.01:20)*2*pi; %rad/s
for ii=1:length(om)
    A = -om(ii)^2*MFF + 1i*om(ii)*RFF +KFF; %1/G
    X(:,ii) = A\F0; %(=1/A=G*F0)
    Xpp(:,ii)=-om(ii)^2*X(:,ii);
end

figure
%plot FRF horizontal displ node 1
nexttile
semilogy(om/(2*pi),abs(X(idb(1,1),:))) %forse meglio solo om in rad/s
title("FRF magnitude of horizontal displacement of A")
xlabel('Frequency [Hz]')
ylabel('magnitude(FRF) [m/N]')
grid on
nexttile
plot(om/(2*pi),angle(X(idb(1,1),:)))
title("FRF phase of horizontal displacement of A")
xlabel('Frequency [Hz]')
ylabel('phase(FRF) [rad]')
grid on
hold off

figure
%plot FRF vertical displ node 4
nexttile
semilogy(om/(2*pi),abs(X(idb(4,2),:)))
title("FRF magnitude of vertical displacement of D")
xlabel('Frequency [Hz]')
ylabel('magnitude(FRF) [m/N]')
grid on
nexttile
plot(om/(2*pi),angle(X(idb(4,2),:)))
title("FRF magnitude of vertical displacement of D")
xlabel('Frequency [Hz]')
ylabel('phase(FRF) [rad]')
grid on
hold off

figure
%plot FRF horizontal acceleration node 1
nexttile
semilogy(om/(2*pi),abs(Xpp(idb(1,1),:)))
title("FRF magnitude of horizontal acceleration of A")
xlabel('Frequency [Hz]')
ylabel('magnitude(FRF) [m/(s^2*N)]')
grid on
nexttile
plot(om/(2*pi),angle(Xpp(idb(1,1),:)))
title("FRF phase of horizontal acceleration of A")
xlabel('Frequency [Hz]')
ylabel('phase(FRF) [rad]')
grid on
hold off

figure
%plot FRF vertical acceleration node 4
nexttile
semilogy(om/(2*pi),abs(Xpp(idb(4,2),:)))
title("FRF magnitude of vertical acceleration of D")
xlabel('Frequency [Hz]')
ylabel('magnitude(FRF) [m/(s^2*N)]')
grid on
nexttile
plot(om/(2*pi),angle(Xpp(idb(4,2),:)))
title("FRF phase of vertical acceleration of D")
xlabel('Frequency [Hz]')
ylabel('phase(FRF) [rad]')
grid on
hold off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Internal actions
nel=4;
idof_A = idb(2,:); % : take all elements along second row 
idof_Ar = idb(5,:);

lambda= [cos(gamma(nel)) sin(gamma(nel)) 0;
         -sin(gamma(nel)) cos(gamma(nel)) 0;
         0 0 1];

XiG = X(idof_A,:); %node2 global
XjG = X(idof_Ar,:); %node5 global
XiL = lambda*XiG;
XjL = lambda*XjG;

L_el = l(nel);
csi = L_el; %internal action in the specific csi position
c = -3/L_el^2*XiL(2,:)+3/L_el^2*XjL(2,:)-2/L_el*XiL(3,:)-1/L_el*XjL(3,:);
d = 2/L_el^3*XiL(2,:)-2/L_el^3*XjL(2,:)+1/L_el^2*XiL(3,:)+1/L_el^2*XjL(3,:);

N12 = 2.2812e+07*(XjL(1,:)-XiL(1,:))/(L_el);
M12 = 1.5812e+03*(2*c + 6*d*csi);
T12 = 1.5812e+03*(6*d);

figure
%plot FRF axial force EG midspan
nexttile
semilogy(om/(2*pi),abs(N12(:)))
title("FRF magnitude of internal axial force in the midspan of EG")
xlabel('Frequency [Hz]')
ylabel('magnitude(FRF) [N/N]')
grid on
nexttile
plot(om/(2*pi),angle(N12(:)))
title("FRF phase of internal axial force in the midspan of EG")
xlabel('Frequency [Hz]')
ylabel('phase(FRF) [rad]')
grid on
hold off

figure
%plot FRF bending moment EG midspan
nexttile
semilogy(om/(2*pi),abs(M12(:)))
title("FRF magnitude of internal bending moment in the midspan of EG")
xlabel('Frequency [Hz]')
ylabel('magnitude(FRF) [m^-1]')
grid on
nexttile
plot(om/(2*pi),angle(M12(:)))
title("FRF phase of internal bending moment in the midspan of EG")
xlabel('Frequency [Hz]')
ylabel('phase(FRF) [rad]')
grid on
hold off

figure
%plot FRF shear force EG midspan
nexttile
semilogy(om/(2*pi),abs(T12(:)))
title("FRF magnitude of internal shear force in the midspan of EG")
xlabel('Frequency [Hz]')
ylabel('magnitude(FRF) [N/N]')
grid on
nexttile
plot(om/(2*pi),angle(T12(:)))
title("FRF phase of internal shear force in the midspan of EG")
xlabel('Frequency [Hz]')
ylabel('phase(FRF) [rad]')
grid on
hold off
