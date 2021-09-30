clear all; close all; clc;
%% Transient Fourier's heat conduction: analytical solution
% Takuma Shiga, The University of Tokyo, 
% shiga@photon.t.u-tokyo.ac.jp

Kn=0.1;                   % Knudsen number
tau=[0.1 1 10];         % Normalized time
xi=linspace(0,1,40);    % Normalized position
summax=1000;            % Upper bound of the summation

%% Temperature profile
theta(1:length(tau),1:length(xi))=0; % Normalized temperature

for t=1:length(tau)
    nonsteady(1:length(xi))=0.0;
    for m=1:summax
        nonsteady=nonsteady+...
            sin(m*pi*xi)*exp(-m^2*pi^2*tau(t)*Kn^2/3)/m;
    end
    theta(t,:)=1-xi-(2/pi).*nonsteady;
end

f=figure; ax=gca; box on;
plot(xi,theta(1,:),'bo-','LineWidth',1.5,'MarkerFaceColor','b')
hold on
plot(xi,theta(2,:),'k^-','LineWidth',1.5,'MarkerFaceColor','k')
plot(xi,theta(3,:),'rs-','LineWidth',1.5,'MarkerFaceColor','r')
grid on
xlabel('\xi')
ylabel('\theta')
legend({'\tau=0.1','\tau=1','\tau=10'})
title(['Kn=',num2str(Kn)])
set(ax,'FontSize',18);
set(ax,'FontName','Arial')
set(ax,'LineWidth',1.5)

%% Heat flux profile: this will be normalized by radiative heat flux, beta
Q(1:length(tau),1:length(xi))=0; % Normalized heat flux

for t=1:length(tau)
    nonsteady(1:length(xi))=0.0;
    for m=1:summax
        nonsteady=nonsteady+...
            cos(m*pi*xi)*exp(-m^2*pi^2*tau(t)*Kn^2/3);
    end
    Q(t,:)=1+2*nonsteady;
end

f=figure; ax=gca; box on;
plot(xi,Q(1,:),'bo-','LineWidth',1.5,'MarkerFaceColor','b')
hold on
plot(xi,Q(2,:),'k^-','LineWidth',1.5,'MarkerFaceColor','k')
plot(xi,Q(3,:),'rs-','LineWidth',1.5,'MarkerFaceColor','r')
grid on
xlabel('\xi')
ylabel('$Q$','Interpreter','latex')
legend({'\tau=0.1','\tau=1','\tau=10'})
title(['Kn=',num2str(Kn)])
set(ax,'FontSize',18);
set(ax,'FontName','Arial')
set(ax,'LineWidth',1.5)