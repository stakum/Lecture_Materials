close all; clear all; clc;
%% transmission function for single-atomic junction

omega=linspace(1e-3,2,500);
m1=1; % mass for left-semi-infinite lattice
beta1=1; % spring constant for left semi-infinite lattice
m2=3; % mass for right-semi-infinite lattice
beta2=1; % spring constant for right semi-infinite lattice
beta12=1; % spring constant at the junction

%% Plot transmittance
T1=transmittance_atomic_junction(omega,m1,beta1,m2,beta2,1); %beta12=1
T2=transmittance_atomic_junction(omega,m1,beta1,m2,beta2,0.2); %beta12=0.2

f=figure;ax=gca; box on;
plot(omega/max(omega),T1,'-','LineWidth',2)
hold on
plot(omega/max(omega),T2,'-','LineWidth',2)
xline(sqrt(1/3),'k--','LineWidth',1.5)
grid on
xlabel('$\omega/\omega_0$','Interpreter','latex')
ylabel('Transmittance')
legend('$\beta_{12}/\beta_1=1$','$\beta_{12}/\beta_1=0.2$','Interpreter','latex')
set(ax,'FontSize',18);
set(ax,'FontName','Arial')
set(ax,'LineWidth',2.0)



%% function 
function T=transmittance_atomic_junction(omega,m1,beta1,m2,beta2,beta12)
    omega1=sqrt(4*beta1/m1);
    omega2=sqrt(4*beta2/m2);
    ka1=2*asin(omega/omega1);
    ka2=2*asin(omega/omega2);
    lambda1=exp(1i*ka1);
    lambda2=exp(1i*ka2);
    coefficient=sqrt(4*beta2*m2-omega.^2*m2^2)./sqrt(4*beta1*m1-omega.^2*m1^2);
    t12=(-beta12*beta1*(lambda1-1./lambda1))...
    ./((m1*omega.^2-beta1-beta12+beta1.*lambda1).*(m2*omega.^2-beta2-beta12+beta2.*lambda2)-beta12^2);
    T=coefficient.*abs(t12).^2;
end
