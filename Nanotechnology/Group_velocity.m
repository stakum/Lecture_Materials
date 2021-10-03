clear all; close all; clc;
%% group velocity
% Takuma Shiga, The University of Tokyo, 
% shiga@photon.t.u-tokyo.ac.jp

n=1000;
a=1;
x=(1:n)*a*0.25;
k=(0:n-1)*2*pi/a/n;
omega0=2;
T=2*pi/omega0;
dt=T/10;

%% Atomic motion (Beating)
m1=200;
m2=210;
omega1=omega0*abs(sin(k(m1)*a/2));
omega2=omega0*abs(sin(k(m2)*a/2));
u1=sin(k(m1)*x);
u2=sin(k(m2)*x);
u=u1+u2;
[up,lo]=envelope(u);

f=figure;
h1=plot(x,u1,'b:','LineWidth',1.0);
hold on
h2=plot(x,u2,'r--','LineWidth',1.0);
h3=plot(x,u,'k-','LineWidth',2.0);
h4=plot(x,up,'k--','LineWidth',1.5);
h5=plot(x,lo,'k--','LineWidth',1.5);
ax=gca;
xlabel('Coordinate')
ylabel('Displacement')
set(ax,'FontSize',18);
set(ax,'FontName','Arial')
set(ax,'LineWidth',1.5)
ylim([-2 2])

t=0;
for ind=1:round(50*T/dt)
    t=t+dt;
    u1=cos(k(m1)*x-omega1*t);
    u2=cos(k(m2)*x-omega2*t);
    u=u1+u2;
    [up,lo]=envelope(u);
    h1.YData=u1;
    h2.YData=u2;
    h3.YData=u;
    h4.YData=up;
    h5.YData=lo;
    drawnow;
end

