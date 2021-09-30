clear all; close all; clc;
%% Graphene phonon dispersion relation
% Reference: 
% - R. Saito, G. Dresselhaus, M. S. Dresselhaus, Physical Properties of Carbon Nanotube (1998) Chapter 9

% Takuma Shiga, The University of Tokyo, 
% shiga@photon.t.u-tokyo.ac.jp

dyn=1e-5; % Newton (N)
cm=1e-2;  % meter (m)
au=1.6726219e-27;
M=12*au;
phi=[36.50 24.50 9.82;...
     8.80 -3.23 -0.40;...
     3.00 -5.25 0.15;...
    -1.92 2.29 -0.58];
phi=phi*1e4*dyn/cm/M;


acc=1.42;
a=sqrt(3)*acc;
n=100;
omega=zeros(6,3*(n+1));
wavenumber=zeros(1,3*(n+1));
% G-M-K-G
counter=0;
dkx_GM=2*pi/(sqrt(3)*a)/n;
dkx_MK=2*pi/(sqrt(3)*a)/n; dky_MK=2*pi/(3*a)/n;
for ind=0:n % G-M
    counter=counter+1;
    kvec=[ind*dkx_GM 0];    
    [DAA,DAB,DBA,DBB]=dynamat(phi,kvec,a);
    D=[DAA DAB; DBA DBB];
    tmp=sqrt(abs(eig(D)));
    wavenumber(counter)=ind*dkx_GM;
    omega(:,counter)=sort(tmp);
end
for ind=0:n % M-K
    counter=counter+1;
    kvec=[2*pi/(sqrt(3)*a) ind*dky_MK];
    [DAA,DAB,DBA,DBB]=dynamat(phi,kvec,a);
    D=[DAA DAB; DBA DBB];
    tmp=sqrt(abs(eig(D)));
    wavenumber(counter)=wavenumber(n+1)+ind*dky_MK;
    omega(:,counter)=sort(tmp);
end
for ind=0:n % K-G
    counter=counter+1;
    kvec=[2*pi/(sqrt(3)*a)-ind*dkx_MK 2*pi/(3*a)-ind*dky_MK];
    [DAA,DAB,DBA,DBB]=dynamat(phi,kvec,a);
    D=[DAA DAB; DBA DBB];
    tmp=sqrt(abs(eig(D)));
    wavenumber(counter)=wavenumber(2*n+2)+ind*sqrt(dkx_MK^2+dky_MK^2);
    omega(:,counter)=sort(tmp);
end

figure
% plot(wavenumber,omega(1:6,:)/(2*pi*1e12*0.03)) % (cm^{-1})
plot(wavenumber,omega(1:6,:)/(2*pi*1e12),'LineWidth',2.0) % (THz)
hold on
plot([wavenumber(n+1) wavenumber(n+1)],[0 50],'k-','LineWidth',1.5)
plot([wavenumber(2*n+2) wavenumber(2*n+2)],[0 50],'k-','LineWidth',1.5)
xlim([0 max(wavenumber)])
xticks([0 wavenumber(n+1) wavenumber(2*n+2) wavenumber(end)])
xticklabels({'\Gamma','M','K','\Gamma'})
xlabel('Wavenumber, {\it q}')
ylabel('Frequency (THz)')

%%
n=20;
k=2*pi/(sqrt(3)*a);
kx=linspace(-k,k,n); % Uniform grid in the x-direction
ky=linspace(-k,k,n); % Uniform grid in the y-direction

omega=zeros(6,length(kx),length(ky));
for ind=1:length(kx)
    for jnd=1:length(ky)
        kvec=[kx(ind) ky(jnd)];
        [DAA,DAB,DBA,DBB]=dynamat(phi,kvec,a);
        D=[DAA DAB; DBA DBB];
        tmp=sqrt(abs(eig(D)));
        tmp=sort(tmp);
        omega(:,ind,jnd)=tmp/(2*pi*1e12);
    end
end
figure
hold on
for ind=1:6
    surf(kx,ky',reshape(omega(ind,:,:),[length(kx),length(ky)]))
end
grid on
xlabel('Normalized {\it q_{x}}')
ylabel('Normalized {\it q_{y}}')
zlabel('Frequency (THz)')
view(3)
        hold on
        line([k k],[-k/sqrt(3) k/sqrt(3)],[0 0],'Color','k','LineWidth',1.5) 
        line([k 0],[k/sqrt(3) k],[0 0],'Color','k','LineWidth',1.5) 
        line([0 -k],[k k/sqrt(3)],[0 0],'Color','k','LineWidth',1.5) 
        line([-k -k],[k/sqrt(3) -k/sqrt(3)],[0 0],'Color','k','LineWidth',1.5) 
        line([-k 0],[-k/sqrt(3) -k],[0 0],'Color','k','LineWidth',1.5) 
        line([0 k],[-k -k/sqrt(3) ],[0 0],'Color','k','LineWidth',1.5) 
        hold off

figure
for ind=1:6
    subplot(2,3,ind)
    [C,h]=contourf(kx,ky',reshape(omega(ind,:,:),[length(kx),length(ky)]));
    clabel(C,h)
end


%% Functions

function [DAA,DAB,DBA,DBB]=dynamat(phi,kvec,a)
phi1NN=diag(phi(1,1:3));
phi2NN=diag(phi(2,1:3));
phi3NN=diag(phi(3,1:3));
phi4NN=diag(phi(4,1:3));

DAA=zeros(3,3);
DAB=zeros(3,3);
DBA=zeros(3,3);
DBB=zeros(3,3);
% 1NN
for ind=0:2
    theta=ind*(2*pi/3);
    U=rotmat(theta);
    rvec=relative_vector(a/sqrt(3),theta);
    K=inv(U)*phi1NN*U;
    DAB=DAB-K*exp(1i*dot(kvec,rvec));
    DAA=DAA+K;
    
    theta=ind*(2*pi/3)+pi/3;
    U=rotmat(theta);
    rvec=relative_vector(a/sqrt(3),theta);
    K=inv(U)*phi1NN*U;
    DBA=DBA-K*exp(1i*dot(kvec,rvec));
    DBB=DBB+K;
end
% 2NN
for ind=0:5
    theta=ind*(2*pi/6)+pi/6;
    U=rotmat(theta);
    rvec=relative_vector(a,theta);
    K=inv(U)*phi2NN*U;
    DAA=DAA+K*(1-exp(1i*dot(kvec,rvec)));
    DBB=DBB+K*(1-exp(1i*dot(kvec,rvec)));
end
% 3NN
for ind=0:2
    theta=ind*(2*pi/3)+pi/3;
    U=rotmat(theta);
    rvec=relative_vector(2*a/sqrt(3),theta);
    K=inv(U)*phi3NN*U;
    DAB=DAB-K*exp(1i*dot(kvec,rvec));
    DAA=DAA+K;
    theta=ind*(2*pi/3);
    U=rotmat(theta);
    rvec=relative_vector(2*a/sqrt(3),theta);
    K=inv(U)*phi3NN*U;
    DBA=DBA-K*exp(1i*dot(kvec,rvec));
    DBB=DBB+K;
end
% 4NN
theta4NN_AA=[0.3335 1.7609 2.4279];
theta4NN_BB=[2.8081 1.3807 0.7137];
for ind=-1:1
    if ind~=0
        for jnd=1:3
            theta=ind*theta4NN_AA(jnd);
            U=rotmat(theta);
            rvec=relative_vector(sqrt(21)*a/3,theta);
            K=inv(U)*phi4NN*U;
            DAB=DAB-K*exp(1i*dot(kvec,rvec));
            DAA=DAA+K;
        end
        for jnd=1:3
            theta=ind*theta4NN_BB(jnd);
            U=rotmat(theta);
            rvec=relative_vector(sqrt(21)*a/3,theta);
            K=inv(U)*phi4NN*U;
            DBA=DBA-K*exp(1i*dot(kvec,rvec));
            DBB=DBB+K;
        end        
    end
end
end

function U=rotmat(theta)
    U=[cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
end

function Rij=relative_vector(distance,theta)
    Rij=distance*[cos(theta) sin(theta)];
end

