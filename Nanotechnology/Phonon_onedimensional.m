clear all; close all; clc;
%% one-dimensional monoatomic lattice
K=1;
M=1;
N=50;
q=(-N/2+1:N/2)*2*pi/N;
omega=sqrt(4*K/M)*abs(sin(q/2));
qq=(-pi:0.01:pi);
omega_continuum=sqrt(4*K/M)*abs(sin(qq/2));
f=figure; 
ax=gca;
box on 
plot(qq/pi,omega_continuum,'b-','LineWidth',2)
hold on
plot(q/pi,omega,'bo','MarkerFaceColor','b')
plot(q/pi,q,'k--','LineWidth',1.5)
pbaspect([2 1 1])
grid on
set(gca, 'XTick', [-1 0 1]);
set(gca, 'XTickLabel', {'-$\pi/a$','0','$\pi/a$'}, 'TickLabelInterpreter', 'latex');
set(gca, 'YTick', [0 1 2]);
set(gca, 'YTickLabel', {'0','$\omega_0/2$','$\omega_0$'}, 'TickLabelInterpreter', 'latex');
ylim([0 2.2])
set(ax,'FontSize',18);
set(ax,'FontName','Arial')
set(ax,'LineWidth',1.5)

%% one-dimensional diatomic lattice
K=1; % spring constant
m=1;
M=2; % mass
alat=1; % lattice constant
N=50; % number of primitive unit cell
nmin=-N/2+1; nmax=N/2;
n=(nmin:1:nmax); % allowable q-points
q=(2*pi/alat)*n/(N); % q-points
q2=(0:0.01:0.7*pi/alat);

omega_opt=sqrt(K*(1/m+1/M)+K*sqrt((1/m+1/M)^2-4*sin(abs(q*alat/2)).^2/M/m));
omega_ac=sqrt(K*(1/m+1/M)-K*sqrt((1/m+1/M)^2-4*sin(abs(q*alat/2)).^2/M/m));

f=figure; 
ax=gca;
plot(q/(pi/alat),omega_opt,'bo-','MarkerFaceColor','b','LineWidth',2)
hold on
plot(q/(pi/alat),omega_ac,'ro-','MarkerFaceColor','r','LineWidth',2)
plot(q2/(pi/alat),sqrt(K/2/(m+M))*q2*alat,'k--','LineWidth',1.5)

grid on
set(gca, 'XTick', [-1 0 1]);
set(gca, 'XTickLabel', {'-$\pi/a$','0','$\pi/a$'}, 'TickLabelInterpreter', 'latex');
set(gca, 'YTick', [0 1 sqrt(2*K*(1/m+1/M))]);
set(gca, 'YTickLabel', {'0',' ','$\omega_{\max}$'}, 'TickLabelInterpreter', 'latex');
ylim([0 sqrt(2*K*(1/m+1/M)*1.1)])
set(ax,'FontSize',18);
set(ax,'FontName','Arial')
set(ax,'LineWidth',1.5)