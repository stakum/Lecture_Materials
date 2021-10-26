x1=-16:1:6; % 変数x1の範囲 
x2=-12:1:6; % 変数x2の範囲
tmax=10; % 時刻の上限
tspan=linspace(0,tmax,50); % 時刻の範囲
x0=[-5 -2]; % 平衡点
x_ini=[-5 -1; -5 0; -5 1; -5 2;]; % t=0における複数の位置ベクトル

ODE=@(t,x)[x(1)-2*x(2)+1;
           x(1)-x(2)+3];

sol=zeros(length(tspan),2,length(x0));
for ind=1:length(x_ini) 
    [t,tmp]=ode45(ODE,tspan,[x_ini(ind,1);x_ini(ind,2)]); 
    if(max(t)~=tmax)
        tspan_tmp=linspace(0,max(t),50);
        [t,tmp]=ode45(ODE,tspan_tmp,[x_ini(ind,1);x_ini(ind,2)]); 
    end
    sol(:,:,ind)=tmp;
end

[x1,x2]=meshgrid(x1,x2); 
dx1=x1-2*x2+1; 
dx2=x1-x2+3;
dx=sqrt(dx1.^2+dx2.^2);
dx1=dx1./dx; dx2=dx2./dx; % 長さ1に正規化

fig=figure;
ax=gca; 
box on; 
quiver(x1,x2,dx1,dx2) % ベクトル場のプロット
hold on
for ind=1:length(x_ini)
    plot(sol(:,1,ind),sol(:,2,ind),'r-','LineWidth',2) % ode45で解いた軌跡のプロット
end
plot(x0(:,1),x0(:,2),'k*','MarkerFaceColor','k','MarkerSize',12) % 平衡点のプロット
grid on
xlabel('$$x_1$$','InterPreter','latex')
ylabel('$$x_2$$','InterPreter','latex')
set(ax,'FontSize',18);
set(ax,'FontName','Arial')
set(ax,'LineWidth',1.5)
xlim([-14 4])
ylim([-10 4])
