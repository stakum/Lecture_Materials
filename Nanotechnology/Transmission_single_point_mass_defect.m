clear all; close all; clc;
%% Transmission on single-point mass defect one-dimensional lattice 

ka=linspace(0,pi,100); % wavenumber 
beta=1; % Spring constant for one-dimensional lattice
m=1; % mass of a particle in host lattice
M=2; % mass of a single-point defect

omega=sqrt(2*beta/m*(1-cos(ka))); % dispersion relation
omega0=sqrt(4*beta/m); % maximum angular frequency

Tka=transmittance_ka(ka,omega,beta,m,M);
Tom=transmittance_omega(omega,omega0,beta,m,M);

%% Plot transmittance in omega and wavenumber representations
f=figure; 
subplot(1,2,1); ax=gca; box on;
plot(omega/omega0,Tom,'b-','LineWidth',2)
grid on
xlabel('$\omega/\omega_0$','Interpreter','latex')
ylabel('Transmittance')
set(ax,'FontSize',18);
set(ax,'FontName','Arial')
set(ax,'LineWidth',1.5)

subplot(1,2,2); ax=gca; box on;
plot(ka/pi,Tka,'b-','LineWidth',2)
grid on
xlabel('$ka/\pi$','Interpreter','latex')
ylabel('Transmittance')
set(ax,'FontSize',18);
set(ax,'FontName','Arial')
set(ax,'LineWidth',1.5)

%% Plot transmittance with different defect masses
Tom2=transmittance_omega(omega,omega0,beta,m,2);
Tom5=transmittance_omega(omega,omega0,beta,m,5);
Tom10=transmittance_omega(omega,omega0,beta,m,10);

f=figure;ax=gca; box on;
plot(omega/omega0,Tom2,'-','LineWidth',2)
hold on
plot(omega/omega0,Tom5,'-','LineWidth',2)
plot(omega/omega0,Tom10,'-','LineWidth',2)
grid on
xlabel('$\omega/\omega_0$','Interpreter','latex')
ylabel('Transmittance')
legend('$M=2$','$M=5$','$M=10$','Interpreter','latex')
set(ax,'FontSize',18);
set(ax,'FontName','Arial')
set(ax,'LineWidth',1.5)




%% function 
function T=transmittance_ka(ka,omega,beta,m,M) 
    % Transmittance in the wavenumber representation
    T=4*beta^2*sin(ka).^2./((M-m)^2*omega.^4+4*beta^2*sin(ka).^2);
end

function T=transmittance_omega(omega,omega0,beta,m,M)
    % Transmittance in the angular-frequency representation
    T=(omega0^2-omega.^2)./((M/m-1)^2*omega.^2+(omega0^2-omega.^2));
end