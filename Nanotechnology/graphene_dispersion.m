clear all; close all; clc;
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
n=40;
kx=linspace(-2*pi/(sqrt(3)*a),2*pi/(sqrt(3)*a),n);
ky=linspace(-2*pi/(3*a),2*pi/(3*a),n);
kmax=2*pi/(sqrt(3)*a);
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
    surf(kx/max(kx),ky'/max(ky),reshape(omega(ind,:,:),[length(kx),length(ky)]))
end
grid on
xlabel('Normalized {\it q_{x}}')
ylabel('Normalized {\it q_{y}}')
zlabel('Frequency (THz)')
view(3)

figure
for ind=1:6
    subplot(2,3,ind)
    [C,h]=contourf(kx/kmax,ky'/kmax,reshape(omega(ind,:,:),[length(kx),length(ky)]));
    clabel(C,h)
end
