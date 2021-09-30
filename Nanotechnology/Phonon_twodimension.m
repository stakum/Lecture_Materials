clear all; close all; clc;
%% Phonon dispersion relation of two-dimensional monoatomic lattice
% Takuma Shiga, The University of Tokyo, 
% shiga@photon.t.u-tokyo.ac.jp

alat=1.0; % lattice constant 
K1=1.0; % spring constant for nearest neighbor
K2=0.3; % spring constant for next nearest neighbor 
Nx=20; % number of primitive unit cell along x-direction 
Ny=20; % number of primitive unit cell along y-direction

step=0;
for qqx=-Nx/2+1:Nx/2
    for qqy=-Ny/2+1:Ny/2
        step=step+1
        qx=2*pi*qqx/alat/Nx;
        qy=2*pi*qqy/alat/Ny;
        D(1,1)=2*K1*(1-cos(qx*alat))+2*K2*(1-cos(qx*alat)*cos(qy*alat));
        D(1,2)=2*K2*sin(qx*alat)*sin(qy*alat);
        D(2,1)=D(1,2);
        D(2,2)=2*K1*(1-cos(qy*alat))+2*K2*(1-cos(qx*alat)*cos(qy*alat));
        temp(1,:)=eig(D);
        omega1(qqx+Nx/2,qqy+Ny/2)=sqrt(temp(1,1));
        omega2(qqx+Nx/2,qqy+Ny/2)=sqrt(temp(1,2));
    end
end
qx=2*pi*(-Nx/2+1:Nx/2)/alat/Nx;
qy=2*pi*(-Ny/2+1:Ny/2)/alat/Ny;
[qx,qy]=meshgrid(qx,qy);

f=figure; 
ax=gca;
surf(qx/(pi/alat),qy/(pi/alat),omega1,'LineWidth',1.5)
hold on
surf(qx/(pi/alat),qy/(pi/alat),omega2,'LineWidth',1.5)
grid on
xlim([-1 1]); ylim([-1 1]);zlim([0 3])

set(gca, 'XTick', [-1 -0.5 0 0.5 1]);
set(gca, 'XTickLabel', {'-$\pi/a$','-$\pi/a$/2','0','$\pi/a$/2','$\pi/a$'}, 'TickLabelInterpreter', 'latex');
set(gca, 'YTick', [-1 -0.5 0 0.5 1]);
set(gca, 'YTickLabel', {'-$\pi/a$','-$\pi/a$/2','0','$\pi/a$/2','$\pi/a$'}, 'TickLabelInterpreter', 'latex');
set(ax,'FontSize',18);
set(ax,'FontName','Arial')
set(ax,'LineWidth',1.5)
xlabel('$k_x$','Interpreter','latex')
ylabel('$k_y$','Interpreter','latex')
zlabel('Frequency (arb. unit)')


f=figure; 
subplot(1,2,1)
ax=gca;
contour(qx/(pi/alat),qy/(pi/alat),omega1,'ShowText','on','LineWidth',2)
set(gca, 'XTick', [-1 -0.5 0 0.5 1]);
set(gca, 'XTickLabel', {'-$\pi/a$','-$\pi/a$/2','0','$\pi/a$/2','$\pi/a$'}, 'TickLabelInterpreter', 'latex');
set(gca, 'YTick', [-1 -0.5 0 0.5 1]);
set(gca, 'YTickLabel', {'-$\pi/a$','-$\pi/a$/2','0','$\pi/a$/2','$\pi/a$'}, 'TickLabelInterpreter', 'latex');
set(ax,'FontSize',18);
set(ax,'FontName','Arial')
set(ax,'LineWidth',1.5)
xlabel('$k_x$','Interpreter','latex')
ylabel('$k_y$','Interpreter','latex')

subplot(1,2,2)
ax=gca;
contour(qx/(pi/alat),qy/(pi/alat),omega2,'ShowText','on','LineWidth',2)
set(gca, 'XTick', [-1 -0.5 0 0.5 1]);
set(gca, 'XTickLabel', {'-$\pi/a$','-$\pi/a$/2','0','$\pi/a$/2','$\pi/a$'}, 'TickLabelInterpreter', 'latex');
set(gca, 'YTick', [-1 -0.5 0 0.5 1]);
set(gca, 'YTickLabel', {'-$\pi/a$','-$\pi/a$/2','0','$\pi/a$/2','$\pi/a$'}, 'TickLabelInterpreter', 'latex');
set(ax,'FontSize',18);
set(ax,'FontName','Arial')
set(ax,'LineWidth',1.5)
xlabel('$k_x$','Interpreter','latex')
ylabel('$k_y$','Interpreter','latex')

%% Phonon dispersion relation of two-dimensional monoatomic lattice along hight symmetry line 
clear all
alat=1.0; % lattice constant 
K1=1.0; % spring constant for nearest neighbor
K2=0.3; % spring constant for next nearest neighbor 
Nx=100; % number of primitive unit cell along x-direction 
Ny=100; % number of primitive unit cell along y-direction

hsyml=[0];
% W -> G
step=0;
for ii=Nx/2:-1:0
        step=step+1
        qx=2*pi*ii/alat/Nx;
        qy=2*pi*ii/alat/Ny;
        %qq(step)=sqrt((sqrt(2)*pi*ii/alat/Nx)^2-qx*qx-qy*qy);
        qq(step)=sqrt(2)*pi/alat*step/(Nx/2);
        D(1,1)=2*K1*(1-cos(qx*alat))+2*K2*(1-cos(qx*alat)*cos(qy*alat));
        D(1,2)=2*K2*sin(qx*alat)*sin(qy*alat);
        D(2,1)=D(1,2);
        D(2,2)=2*K1*(1-cos(qy*alat))+2*K2*(1-cos(qx*alat)*cos(qy*alat));
        temp(1,:)=eig(D);
        omega1(step)=sqrt(temp(1,1));
        omega2(step)=sqrt(temp(1,2));
end
hsyml=[hsyml qq(step)];

% G -> X 
stepold=step;
for ii=0:Nx/2
        step=step+1
        qx=2*pi*ii/alat/Nx;
        qy=0;
        qq(step)=sqrt(qx*qx+qy*qy)+qq(stepold);
        D(1,1)=2*K1*(1-cos(qx*alat))+2*K2*(1-cos(qx*alat)*cos(qy*alat));
        D(1,2)=2*K2*sin(qx*alat)*sin(qy*alat);
        D(2,1)=D(1,2);
        D(2,2)=2*K1*(1-cos(qy*alat))+2*K2*(1-cos(qx*alat)*cos(qy*alat));
        temp(1,:)=eig(D);
        omega1(step)=sqrt(temp(1,1));
        omega2(step)=sqrt(temp(1,2));
end
hsyml=[hsyml qq(step)];

% X -> W
stepold=step;
for ii=0:Ny/2
        step=step+1
        qx=pi/alat;
        qy=2*pi*ii/alat/Ny;
        qq(step)=sqrt(qy*qy)+qq(stepold);
        D(1,1)=2*K1*(1-cos(qx*alat))+2*K2*(1-cos(qx*alat)*cos(qy*alat));
        D(1,2)=2*K2*sin(qx*alat)*sin(qy*alat);
        D(2,1)=D(1,2);
        D(2,2)=2*K1*(1-cos(qy*alat))+2*K2*(1-cos(qx*alat)*cos(qy*alat));
        temp(1,:)=eig(D);
        omega1(step)=sqrt(temp(1,1));
        omega2(step)=sqrt(temp(1,2));
end
hsyml=[hsyml qq(step)];

f=figure; 
ax=gca;
plot(qq,omega1,'b-','LineWidth',2)
hold on
plot(qq,omega2,'r-','LineWidth',2)
xline(hsyml)
hold off 
grid on
xlim([0 max(hsyml)])
ylabel('Frequency (arb. unit)')
set(gca, 'XTick', [0    4.5317    7.6733   10.8149]);
set(gca, 'XTickLabel', {'W','\Gamma','X','W'});
set(ax,'FontSize',18);
set(ax,'FontName','Arial')
set(ax,'LineWidth',1.5)





