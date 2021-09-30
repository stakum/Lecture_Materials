clear all; close all; clc;  
%% phonon wave-packet dynamics (pwd) for one-dimensional FPU lattice
% Takuma Shiga, The University of Tokyo, 
% shiga@photon.t.u-tokyo.ac.jp

% SI unit is used
amu=1.660539040e-27; % kg
ang=1e-10;    % m
eV=1.602e-19; % J
k0=eV/ang;
omega0=sqrt(4*k0/amu);
dt0=1/omega0;
kb=1.38064852e-23;
temp0=eV/kb;

%% Input parameters
k1=1;              % Harmonic force constant (eV/ang^2)
k12=1;
k2=1;
m=1;              % Atomic mass (kg)
a=1;              % lattice constant (L)
na=4000;          % Number of atoms 
nk=2000;          % Number of k-points
c_q=100;         % Central wavenumber point for WP
c_x=500;          % Central spatial point for WP
sigma=100;         % Gaussian window, coherence length

%% Initialization
filename = 'Animation_PWD.gif';           % File name for movie
x=(0:na-1)*a;                             % Atomic posistion
mass(1:na)=1;                             % Atomic mass
mass(na/2:na)=3;
u(1:na)=0;                                % Atomic displacement  
v(1:na)=0;                                % Atomic velocity 
force(1:na)=0;                            % Atomic force 
q=pi/nk/a*(0:nk-1);                       % Wavenumber
q0=q(c_q);                                % Center of WP in k-space
x0=x(c_x);                                % Center of WP in real-space

%% Dispersion relation and group velocity
for j=1:nk
    omega(j)=sqrt(4*k1/m)*sin(abs(q(j)*a/2));
end
dq=q(2)-q(1);
for j=1:nk
    if(j==1 || j==nk)
        vg(j)=0;
    else
        vg(j)=(omega(j+1)-omega(j-1))/2/dq;
    end
end

f=figure;
subplot(1,2,1)
ax=gca; box on;
plot(q/pi,omega,'b-','LineWidth',2)
xlabel('$ka/\pi$','Interpreter','latex')
ylabel('$\omega/\omega_0$','Interpreter','latex')
grid on
set(ax,'FontSize',18);
set(ax,'FontName','Arial')
set(ax,'LineWidth',1.5)
title('Angular frequency')

subplot(1,2,2)
ax=gca; box on;
plot(q/pi,vg,'r-','LineWidth',2)
xlabel('$ka/\pi$','Interpreter','latex')
ylabel('$v_{g}$','Interpreter','latex')
grid on
set(ax,'FontSize',18);
set(ax,'FontName','Arial')
set(ax,'LineWidth',1.5)
title('Group velocity')

%% Generate an initial condition
upw(1:na)=0; % displacement
vpw(1:na)=0; % velocity
% eigenvector is omitted because longitudinal acoustic wave is considered. 
for j=1:na
    for jj=1:nk
        upw(j)=upw(j)+...
            cos(q(jj)*(x(j)-x0))*exp(-0.5*(q(jj)-q0)^2*sigma^2)/sqrt(mass(j));
        vpw(j)=vpw(j)+...
            omega(jj)*sin(q(jj)*(x(j)-x0))*exp(-0.5*(q(jj)-q0)^2*sigma^2)/sqrt(m);
    end
end
normalize=max(upw);
upw=upw/normalize*0.1;
vpw=vpw/normalize*0.1;
umax=max(abs(upw))*1.05;
vmax=max(abs(vpw))*1.05;
kemax=max(0.5*mass.*vpw.^2)*1.05;

f=figure;
visualize_data(x,upw,vpw,mass,umax,vmax,kemax)

%% Molecular dynamics
estimated_time=1.5*na/2/vg(c_q);
frequency=max(omega);
period=1/(frequency/2/pi);
dt=period/100;
mdlength=round(estimated_time/dt);
nsnap=round(mdlength/100);
u=u+upw;
v=v+vpw;
force=harmonic_force_interfere(na,u,k1,k12,k2);

fig=figure(4);
fig.PaperType       = 'a4';
fig.PaperUnits      = 'centimeters';
fig.PaperPosition   = [2,2,18,24];
fig.Units           = 'centimeters';
fig.Position        = [2,2,18,24];
fig.Color           = 'w';
fig.InvertHardcopy  = 'off';

counter=0;
for md=1:mdlength
    [u,v]=updatexv(na,u,v,force,dt,mass);
    force=harmonic_force_interfere(na,u,k1,k12,k2);
    v=updatev(na,v,force,dt,mass);
    if(mod(md,nsnap)==0)
        visualize_data(x,u,v,mass,umax,vmax,kemax)
        frame=getframe(4);
        im=frame2im(frame);
        [A,map]=rgb2ind(im,32);
        if(counter==0)
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
        else
            imwrite(A,map,filename,'gif','WriteMode','append',...
                'DelayTime',0.1);
        end
        counter=counter+1;
     end
end

kin=0.5*mass.*v.^2;
kin1=sum(kin(1:na/2-1));
kin2=sum(kin(na/2:na));
reflectance=kin1/(kin1+kin2);
transmittance=kin2/(kin1+kin2);
result=['Reflectance=',num2str(reflectance)];disp(result);
result=['Transmittance=',num2str(transmittance)];disp(result);

%% Functions
function [u,v]=updatexv(na,u,v,f,dt,m)
    for j=1:na
        v(j)=v(j)+f(j)*dt/m(j)*0.5;
        u(j)=u(j)+v(j)*dt;
    end
end

function v=updatev(na,v,f,dt,m)
    for j=1:na
        v(j)=v(j)+f(j)*dt/m(j)*0.5;
    end
    p0=mean(m.*v);
    for j=1:na
        v(j)=v(j)-p0/m(j);
    end
end

function f=harmonic_force_interfere(na,u,k1,k12,k2)
    f(1:na)=0;
    for j=1:na
        if(j==1)
            x1=u(1)-u(2);
            x2=u(1)-u(na);
            f(j)=-k1*x1-k12*x2;
        elseif(j>=2 && j<=na/2-2)
            f(j)=-k1*(2*u(j)-u(j-1)-u(j+1));
        elseif(j==na/2-1)
            x1=u(j)-u(j-1);
            x2=u(j)-u(j+1);
            f(j)=-k1*x1-k12*x2;
        elseif(j==na/2)
            x1=u(na/2)-u(na/2-1);
            x2=u(na/2)-u(na/2+1);
            f(j)=-k12*x1-k2*x2;
        elseif(j>=na/2+1 && j<=na-1)
            f(j)=-k2*(2*u(j)-u(j-1)-u(j+1));
        elseif(j==na)
            x1=u(na)-u(na-1);
            x2=u(na)-u(1);
            f(j)=-k2*x1-k1*x2;
        end
    end
end

function visualize_data(x,u,v,mass,umax,vmax,kemax)
    subplot(3,1,1)
    ax=gca; box on;
    plot(x,u)
    grid on
    xlim([0 4000]) ; ylim([-umax umax])
    xlabel('Position') ; ylabel('Displacement')
    set(ax,'FontSize',18);
    set(ax,'FontName','Arial')
    set(ax,'LineWidth',1.5)
    set(gcf, 'Color', 'w'); 
    subplot(3,1,2)
    ax=gca; box on;
    plot(x,v)
    grid on
    xlim([0 4000]) ; ylim([-vmax vmax])
    xlabel('Position') ; ylabel('Velocity')
    set(ax,'FontSize',18);
    set(ax,'FontName','Arial')
    set(ax,'LineWidth',1.5)
    set(gcf, 'Color', 'w'); 
    subplot(3,1,3)
    ke=0.5*mass.*v.^2;
    ax=gca; box on;
    plot(x,ke)
    grid on
    xlim([0 4000]) ; ylim([0 kemax])
    xlabel('Position') ; ylabel('Kinetic energy')
    set(ax,'FontSize',18);
    set(ax,'FontName','Arial')
    set(ax,'LineWidth',1.5)
    set(gcf, 'Color', 'w'); 
    drawnow
end