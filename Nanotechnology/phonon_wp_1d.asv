clear all; close all; clc;  
%% phonon wave-packet dynamics (pwd) for one-dimensional FPU lattice
% E: energy
amu=1.660539040e-27; % kg
ang=1e-10;    % m
eV=1.602e-19; % J
k0=eV/ang;
omega0=sqrt(4*k0/amu);
dt0=1/omega0;
kb=1.38e-23;
temp0=eV/kb;

%% Input parameters
k=1;              % Harmonic force constant (eV/ang^2)
k1=1;
k2=1;
alpha=0;          % Cubic anharmonic force constants (E/L^3)
beta=0;           % Quartic anharmonic force constants (E/L^4)
phi4=0;           % Quartic force constants for onsite potential (E/L^4)
m=1;              % Atomic mass (kg)
a=1;              % lattice constant (L)
na=4002;          % Number of atoms 
nk=2000;          % Number of k-points
c_q=1800;         % Central wavenumber point for WP
c_x=500;          % Central spatial point for WP
sigma=50;         % Gaussian window, coherence length

%% Initialization
filename = 'Animation_PWD.gif';           % File name for movie
x=(0:na-1)*a;                             % Atomic posistion
mass(1:na)=1;                             % Atomic mass
mass(4001)=1;
mass(4002)=3;
u(1:na)=0;                                % Atomic displacement  
v(1:na)=0;                                % Atomic velocity 
force(1:na)=0;                            % Atomic force 
q=pi/nk/a*(0:nk-1);                       % Wavenumber
q0=q(c_q);                                % Center of WP in k-space
x0=x(c_x);                                % Center of WP in real-space

%% calc mode and make wave packet
for j=1:nk
    omega(j)=sqrt(4*k/m)*sin(abs(q(j)*a/2));
end
dq=q(2)-q(1);
for j=1:nk
    if(j==1 || j==nk)
        vg(j)=0;
    else
        vg(j)=(omega(j+1)-omega(j-1))/2/dq;
    end
end
figure(1)
plotyy(q/pi,omega,q/pi,vg)
xlabel('Normalized wave number')
ylabel('Angular frequency')
legend('Dispersion','Group velocity')

upw(1:na)=0;
vpw(1:na)=0;
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
figure(2)
visualize_data(x,upw,vpw,mass,umax,vmax,kemax)

%% Molecular dynamics
estimated_time=na/2/vg(c_q);
frequency=max(omega);
period=1/(frequency/2/pi);
dt=period/100;
mdlength=round(estimated_time/dt);
nsnap=round(mdlength/100);
u=u+upw;
v=v+vpw;
[force]=harmonic_force_interfere(na,u,k,k1,k2);

fig=figure(4);
counter=0;
for md=1:mdlength
    [u,v]=updatexv(na,u,v,force,dt,mass);
    [force]=harmonic_force_interfere(na,u,k,k1,k2);
    [v]=updatev(na,v,force,dt,mass);
    if(mod(md,nsnap)==0)
        fname=['snap-',num2str(md),'.txt'];
        clear out;
        out(:,1)=x;
        out(:,2)=v;
        save(fname,'-ascii','out')
%         visualize_data_for_movie(x,u,v,mass,umax,vmax,kemax,nsta_FPU,nend_FPU,nsta_mass,nend_mass)
        visualize_data(x,u,v,mass,umax,vmax,kemax)
%         F=getframe(fig);
%         writeVideo(writerObj,F);
        frame=getframe(4);
        im=frame2im(frame);
        [A,map]=rgb2ind(im,32);
        if(counter==0)
%             imwrite(A,map,filename);
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
        else
%             imwrite(A,map,filename,'WriteMode','append');
            imwrite(A,map,filename,'gif','WriteMode','append',...
                'DelayTime',0.1);
        end
        counter=counter+1;
     end
end
% close(writerObj);
% figure
% plot(u)
pot1=sum(m/2.0*v(1:na/2).^2)
pot2=sum(m/2.0*v(na/2+1:na).^2)
reflectance=pot1/(pot1+pot2);
transmittance=pot2/(pot1+pot2);
result=['Reflectance=',num2str(reflectance)];disp(result);
result=['Transmittance=',num2str(transmittance)];disp(result);

%% Fourier transform
NFFT=1024;            % Length of velocity data
NH=NFFT/2;            % Half of NFFT
nsnap=50;             % Velocity output interval
f=[0:NH-1]/2/NH...    % Angular frequency for FFT
    /(nsnap*dt)*2*pi;
vsave(na,nsnap)=0;    % Atomic velocity for FFT
usave(na,nsnap)=0;

counter=0;
for md=1:NFFT*nsnap
    [u,v]=updatexv(na,u,v,force,dt,mass);
    [force]=harmonic_force_interfere2(na,u,k,k1,k2);
    [v]=updatev(na,v,force,dt,mass);
    if(mod(md,nsnap)==0)
        counter=counter+1;
        vsave(:,counter)=v(:);
        usave(:,counter)=u(:);
     end
end
Yl(1:NH)=0;
Yr(1:NH)=0;
Ym(1:NH)=0;
for j=1:na-1
    Z=abs(fft(vsave(j,:)))/na/NFFT;
%     Y=Y+Z(1:NH)+Z(NFFT:-1:NFFT-NH+1);
%     Zu=abs(fft(usave(j,:)))/na/NFFT;
    if(j<=1999)
        Yl=Yl+Z(1:NH);
    elseif(j>=2002)
        Yr=Yr+Z(1:NH);
    else
        Ym=Ym+Z(1:NH);
    end
end

% Plot single-sided amplitude spectrum.
figure(5)
plot(f,Yl,'b')
hold on
plot(f,Yr,'r')
plot(f,Ym,'k')
% plot(f,Yu)
title('Amplitude Spectrum of atomic velocity')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
legend('Left lead','Right lead','Middle')
clear out;
out(:,1)=f;
out(:,2)=Yl;
out(:,3)=Ym;
out(:,4)=Yr;
save('pwd_power.txt','-ascii','out')


%% Functions
function [u,v]=updatexv(na,u,v,f,dt,m)
    for j=1:na
        v(j)=v(j)+f(j)*dt/m(j)*0.5;
        u(j)=u(j)+v(j)*dt;
    end
end

function [v]=updatev(na,v,f,dt,m)
    for j=1:na
        v(j)=v(j)+f(j)*dt/m(j)*0.5;
    end
    p0=mean(m.*v);
    for j=1:na
        v(j)=v(j)-p0/m(j);
    end
end