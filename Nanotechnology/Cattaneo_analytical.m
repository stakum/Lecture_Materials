clear all; close all; clc;
%% Analytical solution of Cattaneo equation
% Takuma Shiga, The University of Tokyo, 
% shiga@photon.t.u-tokyo.ac.jp

Kn=10;                   % Knudsen number
tau=[0.1 1 10];         % Normalized time
xi=linspace(0,1,40);    % Normalized position
summax=1000;            % Upper bound of the summation

%
% Temperature profile
%
theta(1:length(tau),1:length(xi))=0; % Normalized temperature

for t=1:length(tau)
    for i=1:length(xi)
        for n=1:10000
            deltan=sqrt(1-4*n^2*pi^2*Kn^2/3);
            term1=((1+deltan)/2/deltan)*exp(deltan*tau(t)/2);
            term2=((1-deltan)/2/deltan)*exp(-deltan*tau(t)/2);
            nonsteady=2*(term1-term2)*sin(n*pi*xi(i))/(pi*n)*exp(-tau(t)/2);
            theta(t,i)=theta(t,i)-nonsteady;
        end
        theta(t,i)=theta(t,i)+1-xi(i);
    end
end

f=figure; ax=gca; box on;
plot(xi,theta(1,:),'bo-','LineWidth',1.5,'MarkerFaceColor','b')
hold on
plot(xi,theta(2,:),'k^-','LineWidth',1.5,'MarkerFaceColor','k')
plot(xi,theta(3,:),'rs-','LineWidth',1.5,'MarkerFaceColor','r')
ylim([0 1])
grid on
xlabel('\xi')
ylabel('\theta')
legend({'\tau=0.1','\tau=1','\tau=10'})
title(['Kn=',num2str(Kn)])
set(ax,'FontSize',18);
set(ax,'FontName','Arial')
set(ax,'LineWidth',1.5)