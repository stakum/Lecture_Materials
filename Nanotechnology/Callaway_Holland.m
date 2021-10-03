clear all; close all; clc;
%% Callaway-Holland model
% Takuma Shiga, The University of Tokyo, 
% shiga@photon.t.u-tokyo.ac.jp

kB=1.380649e-23;    % Boltzmann constant (J/K)
h=6.62607004e-34;   % Plank constant (J s)
hbar=h/(2*pi);      % Dirac constant (J s)

L=0.716e-2;         % Sample size (m)
F=0.8;              % Scaling factor 

vT=5.86e3;          % m/s (omega<omega_1)
vTU=2.0e3;          % m/s (omega>omega_1)
vL=8.48e3;          % m/s 
vb=((2/vT+1/vL)/3)^(-1);

theta1=180;
theta2=210;
theta3=570;

A=1.32e-44/10;      % for impurity scattering (s^3), original parameter was 1.32e44
BT=9.3e-13;         % for transverse three-phonon scattering (K^-3)
BTU=5.5e-18;        % for transverse three-phonon scattering (sec)
BL=2.0e-24;         % for longitudional 3-phonon scattering (s/K^3)

alpha=(kB/hbar)^4*A;
betaT=(kB/hbar)*BT;
betaTU=(kB/hbar)^2*BTU;
betaL=(kB/hbar)^2*BL;
taub_inv=vb/(L*F);

CT=(kB/(2*pi^2*vT))*(kB/hbar)^3;
CTU=(kB/(2*pi^2*vTU))*(kB/hbar)^3;
CL=(kB/(2*pi^2*vL))*(kB/hbar)^3;

%%
T=logspace(0,3.3,50);

kT0=zeros(size(T));
kTU=zeros(size(T));
kL=zeros(size(T));

for ind=1:length(T)
    % Transverse kappa_T0 contribution
    x=linspace(1e-3,theta1/T(ind));
    wBE=weighed_BoseEinstein(x);
    y=CT*wBE./(taub_inv+alpha*x.^4*T(ind)^4+betaT*x*T(ind)^5);
    kT0(ind)=(2/3)*T(ind)^3*trapz(x,y);
    
    % Transverse kappa_TU contribution
    x=linspace(theta1/T(ind),theta2/T(ind));
    wBE=weighed_BoseEinstein(x);
    y=CTU*wBE./(taub_inv+alpha*x.^4*T(ind)^4+betaTU*x.^2*T(ind)^2./sinh(x));
    kTU(ind)=(2/3)*T(ind)^3*trapz(x,y);

    % Longitudinal kappa_L contribution
    x=linspace(1e-3,theta3/T(ind));
    wBE=weighed_BoseEinstein(x);
    y=CL*wBE./(taub_inv+alpha*x.^4*T(ind)^4+betaL*x.^2*T(ind)^5);
    kL(ind)=(1/3)*T(ind)^3*trapz(x,y);
end

k=kT0+kTU+kL;
exp=load('Holland_Experiment_Fig.3.txt');

f=figure;
ax=gca;
loglog(T,k,'k-','LineWidth',2)
hold on
loglog(exp(:,1),exp(:,2)*1e2,'bo')
loglog(T,kT0,'b-','LineWidth',1.5)
loglog(T,kTU,'b--','LineWidth',1.5)
loglog(T,kL,'r-','LineWidth',1.5)
ylim([0.1 70]*1e2)
xlim([1 2000])
xlabel('Temperature (K)')
ylabel('Thermal conductivity (W/m-K)')
set(ax,'FontSize',18);
set(ax,'FontName','Arial')
set(ax,'LineWidth',1.5)


%% function
function wBE=weighed_BoseEinstein(x)
    wBE=x.^4.*exp(x)./(exp(x)-1).^2;
end