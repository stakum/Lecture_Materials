clear all; close all; clc;
%% Phonon dispersion relation of van der Waals crystal
% Takuma Shiga, The University of Tokyo, 
% shiga@photon.t.u-tokyo.ac.jp

an  = 6.022e23;     % Avogadro number
wm  = 40.0e-3/an;   % MassÅiatomic weight/an)
sigma = 3.40e-10;   % Sigma (m)
eps = 2.67e-21;     % Epsilon (J)
alat=5.315e-10;     % Lattice constant (m)

%% Primitive lattice vectors and their reciprocal vectors
a1=alat/2*[1 1 0];
a2=alat/2*[0 1 1];
a3=alat/2*[1 0 1];
Omega=abs(dot(a1,cross(a2,a3)));
b1=2*pi/Omega*cross(a2,a3);
b2=2*pi/Omega*cross(a3,a1);
b3=2*pi/Omega*cross(a1,a2);
bb=2*pi/alat;

%% Identifying interaction pairs with Verlet's neighbering list
ncell=8;            % isotropic cells
natom=(2*ncell+1)^3;
r=zeros(3,natom);
iatom=0;
for ind=-ncell:ncell
    for jnd=-ncell:ncell
        for knd=-ncell:ncell
            iatom=iatom+1;
            r(:,iatom)=ind*a1+jnd*a2+knd*a3;
            if(ind==0 && jnd==0 && knd==0)
                center=iatom;
            end
        end
    end
end
npair=0;
vvlist=0;
for ind=1:natom
    if(ind~=center)
        distance=abs(norm(r(:,center)-r(:,ind)));
        if(distance<=2.5*sigma)
            npair=npair+1;
            vvlist(npair)=ind;
        end
    end
end

%% Phonon dispersion relations along G-X-K-G-L 
n=100;
omega=zeros(3,3*(n+1));
wavenumber=zeros(1,3*(n+1));
counter=0;
dGX=2*pi/alat/n;
dGL=pi/alat/n;
for ind=0:n % G-X
    counter=counter+1;
    kvec=[ind*dGX 0 0]; 
    Dyn=dynmat(eps,sigma,r,npair,center,vvlist,wm,kvec);
    tmp=sqrt(abs(eig(Dyn)));
    omega(:,counter)=sort(tmp)/(2*pi*1e12);
    wavenumber(counter)=ind*dGX;
end
for ind=0:n %X-K-G
    counter=counter+1;
    kvec=[2*pi/alat-ind*dGX 2*pi/alat-ind*dGX 0]; 
    Dyn=dynmat(eps,sigma,r,npair,center,vvlist,wm,kvec);
    tmp=sqrt(abs(eig(Dyn)));
    omega(:,counter)=sort(tmp)/(2*pi*1e12);
    wavenumber(counter)=wavenumber(n+1)+ind*sqrt(2)*dGX;
end
for ind=0:n %G-L
    counter=counter+1;
    kvec=(ind*dGL)*[1 1 1]; 
    Dyn=dynmat(eps,sigma,r,npair,center,vvlist,wm,kvec);
    tmp=sqrt(abs(eig(Dyn)));
    omega(:,counter)=sort(tmp)/(2*pi*1e12);
    wavenumber(counter)=wavenumber(2*n+2)+ind*sqrt(3)*dGL;
end


f=figure; ax=gca;
plot(wavenumber,omega(1:3,:),'LineWidth',2.0) % (THz)
hold on
tmp=wavenumber(n+1)+0.25*(wavenumber(2*n+2)-wavenumber(n+1));
xline(wavenumber(n+1),'k-','LineWidth',1.5)
xline(tmp,'k--','LineWidth',1.5)
xline(wavenumber(2*n+2),'k-','LineWidth',1.5)
xlim([0 max(wavenumber)])
xticks([0 wavenumber(n+1) tmp wavenumber(2*n+2) wavenumber(end)])
xticklabels({'\Gamma','X','K','\Gamma','L'})
ylabel('Frequency (THz)')
set(ax,'FontSize',18);
set(ax,'FontName','Arial')
set(ax,'LineWidth',1.5)

%% Contour plot
frequency=1.1;
n=10;
kx=linspace(-2*pi/alat,2*pi/alat,n);
ky=linspace(-2*pi/alat,2*pi/alat,n);
kz=linspace(-2*pi/alat,2*pi/alat,n);
omega1=zeros(n,n,n);
omega2=zeros(n,n,n);
omega3=zeros(n,n,n);

for branch=1:3    
    for ind=1:n
        for jnd=1:n
            for knd=1:n
                if(abs(kx(ind))+abs(ky(jnd))+abs(kz(knd))<= 1.5*(2*pi/alat))
                    kvec=[kx(ind) ky(jnd) kz(knd)];
                    Dyn=dynmat(eps,sigma,r,npair,center,vvlist,wm,kvec);
                    tmp=sort(sqrt(abs(eig(Dyn))),'ascend')/(2*pi*1e12);
                    if(branch==1)
                        omega1(ind,jnd,knd)=tmp(branch);
                    elseif(branch==2)
                        omega2(ind,jnd,knd)=tmp(branch);
                    elseif(branch==3)
                        omega3(ind,jnd,knd)=tmp(branch);
                    end
                else
                    if(branch==1)
                        omega1(ind,jnd,knd)=nan;
                    elseif(branch==2)
                        omega2(ind,jnd,knd)=nan;
                    elseif(branch==3)
                        omega3(ind,jnd,knd)=nan;
                    end
                end
            end
        end
    end
end

fig=figure;
fig.PaperType       = 'a4';
fig.PaperUnits      = 'centimeters';
fig.PaperPosition   = [5,30,50,14];
fig.Units           = 'centimeters';
fig.Position        = [5,30,50,14];
fig.Color           = 'w';
fig.InvertHardcopy  = 'off';

subplot(1,3,1)
ax=gca; box on; grid on;
p1=patch(isosurface(kx/(2*pi/alat),ky/(2*pi/alat),kz/(2*pi/alat),omega1,frequency));
isonormals(kx/(2*pi/alat),ky/(2*pi/alat),kz/(2*pi/alat),omega1,p1);
p1.FaceColor='yellow';
p1.EdgeColor='none';
draw_FBZ
set(ax,'FontSize',18);
set(ax,'FontName','Arial')
set(ax,'LineWidth',1.5)
xlabel('$k_x$','Interpreter','latex')
ylabel('$k_y$','Interpreter','latex')
zlabel('$k_z$','Interpreter','latex')
view(3);
camlight
lighting gouraud

subplot(1,3,2)
ax=gca; box on; grid on;
p2=patch(isosurface(kx/(2*pi/alat),ky/(2*pi/alat),kz/(2*pi/alat),omega2,frequency));
isonormals(kx/(2*pi/alat),ky/(2*pi/alat),kz/(2*pi/alat),omega2,p2);
p2.FaceColor='green';
p2.EdgeColor='none';
draw_FBZ
set(ax,'FontSize',18);
set(ax,'FontName','Arial')
set(ax,'LineWidth',1.5)
xlabel('$k_x$','Interpreter','latex')
ylabel('$k_y$','Interpreter','latex')
zlabel('$k_z$','Interpreter','latex')
view(3);
camlight
lighting gouraud

subplot(1,3,3)
ax=gca; box on; grid on;
p3=patch(isosurface(kx/(2*pi/alat),ky/(2*pi/alat),kz/(2*pi/alat),omega3,frequency));
isonormals(kx/(2*pi/alat),ky/(2*pi/alat),kz/(2*pi/alat),omega3,p3);
p3.FaceColor='magenta';
p3.EdgeColor='none';
draw_FBZ
set(ax,'FontSize',18);
set(ax,'FontName','Arial')
set(ax,'LineWidth',1.5)
xlabel('$k_x$','Interpreter','latex')
ylabel('$k_y$','Interpreter','latex')
zlabel('$k_z$','Interpreter','latex')
view(3);
camlight
lighting gouraud


%% Generate dynamical matrix
function Dyn=dynmat(eps,sigma,r,npair,center,vvlist,wm,kvec)
    Dyn=zeros(3,3);
    for xyz1=1:3
        for xyz2=1:3
            for pair=1:npair
                rij_v=r(:,center)-r(:,vvlist(pair));
                rij=abs(norm(rij_v));
                prefactor=(96*eps/rij^2)*(7*(sigma/rij)^12-2*(sigma/rij)^6);
                chi=prefactor*rij_v(xyz1)*rij_v(xyz2)/rij^2;
                if(xyz1==xyz2)
                    chi=chi-(24*eps/rij^2)*(2*(sigma/rij)^12-(sigma/rij)^6);
                end
                Dyn(xyz1,xyz2)=Dyn(xyz1,xyz2)...
                    +(chi/wm)*(1-exp(-1i*dot(kvec,rij_v)));
            end
        end
    end
end

%% Draw the first Brillouin zone
function draw_FBZ_first_quadrant
    line([0 0],[0 0],[0 1],'Color','k','LineWidth',2);
    line([0 0],[0 1],[0 0],'Color','k','LineWidth',2);
    line([0 1],[0 0],[0 0],'Color','k','LineWidth',2);
    line([0 0.5],[0 0],[1 1],'Color','k','LineWidth',2);
    line([0 0],[0 0.5],[1 1],'Color','k','LineWidth',2);
    line([0 0.5],[1 1],[0 0],'Color','k','LineWidth',2);
    line([0 0],[1 1],[0 0.5],'Color','k','LineWidth',2);
    line([1 1],[0 0.5],[0 0],'Color','k','LineWidth',2);
    line([1 1],[0 0],[0 0.5],'Color','k','LineWidth',2);
    line([0.5 1],[1 0.5],[0 0],'Color','k','LineWidth',2);
    line([0 0],[0.5 1],[1 0.5],'Color','k','LineWidth',2);
    line([0.5 1],[0 0],[1 0.5],'Color','k','LineWidth',2);
    line([0.5 0],[0 0.5],[1 1],'Color','k','LineWidth',2);
    line([1 1],[0.5 0],[0 0.5],'Color','k','LineWidth',2);
    line([0.5 0],[1 1],[0 0.5],'Color','k','LineWidth',2);
end
function draw_FBZ

    line([0.5 1],[1 0.5],[0 0],'Color','k','LineWidth',2);
    line([0 0],[0.5 1],[1 0.5],'Color','k','LineWidth',2);
    line([0.5 1],[0 0],[1 0.5],'Color','k','LineWidth',2);
    line([0.5 0],[0 0.5],[1 1],'Color','k','LineWidth',2);
    line([1 1],[0.5 0],[0 0.5],'Color','k','LineWidth',2);
    line([0.5 0],[1 1],[0 0.5],'Color','k','LineWidth',2);
    
    
    line([0.5 1],[1 0.5],[0 0],'Color','k','LineWidth',2);
    line([0 0],[0.5 1],[-1 -0.5],'Color','k','LineWidth',2);
    line([0.5 1],[0 0],[-1 -0.5],'Color','k','LineWidth',2);
    line([0.5 0],[0 0.5],[-1 -1],'Color','k','LineWidth',2);
    line([1 1],[0.5 0],[0 -0.5],'Color','k','LineWidth',2);
    line([0.5 0],[1 1],[0 -0.5],'Color','k','LineWidth',2);
    
    
    line([0.5 1],[1 0.5],[0 0],'Color','k','LineWidth',2);
    line([0 0],[0.5 1],[1 0.5],'Color','k','LineWidth',2);
    line([0.5 1],[0 0],[1 0.5],'Color','k','LineWidth',2);
    line([0.5 0],[0 0.5],[1 1],'Color','k','LineWidth',2);
    line([1 1],[0.5 0],[0 0.5],'Color','k','LineWidth',2);
    line([0.5 0],[1 1],[0 0.5],'Color','k','LineWidth',2);
    
    
    line([0.5 1],[-1 -0.5],[0 0],'Color','k','LineWidth',2);
    line([0 0],[-0.5 -1],[1 0.5],'Color','k','LineWidth',2);
    line([0.5 1],[0 0],[1 0.5],'Color','k','LineWidth',2);
    line([0.5 0],[0 -0.5],[1 1],'Color','k','LineWidth',2);
    line([1 1],[-0.5 0],[0 0.5],'Color','k','LineWidth',2);
    line([0.5 0],[-1 -1],[0 0.5],'Color','k','LineWidth',2);
    
    
    line([0.5 1],[-1 -0.5],[0 0],'Color','k','LineWidth',2);
    line([0 0],[-0.5 -1],[-1 -0.5],'Color','k','LineWidth',2);
    line([0.5 1],[0 0],[-1 -0.5],'Color','k','LineWidth',2);
    line([0.5 0],[0 -0.5],[-1 -1],'Color','k','LineWidth',2);
    line([1 1],[-0.5 0],[0 -0.5],'Color','k','LineWidth',2);
    line([0.5 0],[-1 -1],[0 -0.5],'Color','k','LineWidth',2);
    
    
    line([-0.5 -1],[-1 -0.5],[0 0],'Color','k','LineWidth',2);
    line([0 0],[-0.5 -1],[1 0.5],'Color','k','LineWidth',2);
    line([-0.5 -1],[0 0],[1 0.5],'Color','k','LineWidth',2);
    line([-0.5 0],[0 -0.5],[1 1],'Color','k','LineWidth',2);
    line([-1 -1],[-0.5 0],[0 0.5],'Color','k','LineWidth',2);
    line([-0.5 0],[-1 -1],[0 0.5],'Color','k','LineWidth',2);
    
    line([-0.5 -1],[-1 -0.5],[0 0],'Color','k','LineWidth',2);
    line([0 0],[-0.5 -1],[-1 -0.5],'Color','k','LineWidth',2);
    line([-0.5 -1],[0 0],[-1 -0.5],'Color','k','LineWidth',2);
    line([-0.5 0],[0 -0.5],[-1 -1],'Color','k','LineWidth',2);
    line([-1 -1],[-0.5 0],[0 -0.5],'Color','k','LineWidth',2);
    line([-0.5 0],[-1 -1],[0 -0.5],'Color','k','LineWidth',2);
    
    
    
    line([-0.5 -1],[1 0.5],[0 0],'Color','k','LineWidth',2);
    line([0 0],[0.5 1],[1 0.5],'Color','k','LineWidth',2);
    line([-0.5 -1],[0 0],[1 0.5],'Color','k','LineWidth',2);
    line([-0.5 0],[0 0.5],[1 1],'Color','k','LineWidth',2);
    line([-1 -1],[0.5 0],[0 0.5],'Color','k','LineWidth',2);
    line([-0.5 0],[1 1],[0 0.5],'Color','k','LineWidth',2);
    
    
    
    line([-0.5 -1],[1 0.5],[0 0],'Color','k','LineWidth',2);
    line([0 0],[0.5 1],[-1 -0.5],'Color','k','LineWidth',2);
    line([-0.5 -1],[0 0],[-1 -0.5],'Color','k','LineWidth',2);
    line([-0.5 0],[0 0.5],[-1 -1],'Color','k','LineWidth',2);
    line([-1 -1],[0.5 0],[0 -0.5],'Color','k','LineWidth',2);
    line([-0.5 0],[1 1],[0 -0.5],'Color','k','LineWidth',2);
    
end
