clear all; close all; clc;

x=linspace(-1,1,100);
figure
plot(x,sqrt(1-x.^2))
hold on
plot(x,-sqrt(1-x.^2))
plot(x,sin(x))
pbaspect([1 1 1])

eps=1e-6;
err=1;
x0=[1;1];
iter=1;
x=x0;
x1(iter)=x(1);
x2(iter)=x(2);
while err>eps
    J=Jfun(x);
    F=Ffun(x);
    delta = - J\F;
    x=x+delta;
    err=norm(delta);
    iter=iter+1;
    x1(iter)=x(1);
    x2(iter)=x(2);
end

for ind=1:iter-1
    u(ind)=x1(ind+1)-x1(ind);
    v(ind)=x2(ind+1)-x2(ind);
end
u(iter)=0;
v(iter)=0;
plot(x1,x2,'bo-')
quiver(x1,x2,u,v)




function F=Ffun(x)
    F(1,1)=x(1)^2+x(2)^2-1;
    F(2,1)=x(2)-sin(x(1));
    return
end

function J=Jfun(x)
    pi2=0.5*pi;
    J(1,1)=2*x(1);
    J(1,2)=2*x(2);
    J(2,1)=-cos(x(1));
    J(2,2)=1;
    return
end