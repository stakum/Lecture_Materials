clear all; close all; clc;
x=linspace(0,3,100);
figure
plot(x,f(x),'b-','LineWidth',2)
hold on
plot([min(x) max(x)],[0 0],'k--','LineWidth',2)

%%
clear all; 

eps=1e-6;
% First solution 
% a=0; b=1;
% [x,iter]=bisection(a,b,eps);
% x0=(a+b)/2;
% [x,iter]=newton(x0,eps);
% Second solution
a=1; b=3;
[xx1,iter]=bisection_verbose(a,b,eps);
x0=(a+b)/2;
[xx2,iter]=newton_verbose(x0,eps);

figure
plot(xx1,'bo-')
hold on
plot(xx2,'rs-')

%%
function y=f(x)
y=(x-1).*log(x)-1;
% y=tanh(x)-x/2;
end

function y=g(x)
y=log(x)+(x-1)./x;
% y=1./(cosh(x).^2)-1/2;
end

function [x,iter]=bisection(a,b,eps)
iter=0;
x=(a+b)/2;
fx=f(x); fa=f(a); fb=f(b);
while (b-a)>eps && iter<100
    if(fx*fa<0)
        a=a; b=x;
    else
        a=x; b=b;
    end
    x=(a+b)/2;
    fx=f(x); fa=f(a); fb=f(b);
    iter=iter+1;
end
sol=['Bisection: x = ',num2str(round(x,4,'significant')),' # of iteration = ',num2str(iter)];
disp(sol);
end

function [xx,iter]=bisection_verbose(a,b,eps)
iter=0;
x=(a+b)/2;
fx=f(x); fa=f(a); fb=f(b);
while (b-a)>eps && iter<100
    if(fx*fa<0)
        a=a; b=x;
    else
        a=x; b=b;
    end
    x=(a+b)/2;
    fx=f(x); fa=f(a); fb=f(b);
    iter=iter+1;
    xx(iter)=x;
end
sol=['Bisection: x = ',num2str(round(x,4,'significant')),' # of iteration = ',num2str(iter)];
disp(sol);
end

function [x,iter]=newton(x0,eps)
diff=1;
iter=0;
x=x0;
while abs(diff)>eps && iter<100
    fx=f(x);
    gx=g(x);
    diff=fx/gx;
    x=x-diff;
    iter=iter+1;
end
sol=['Newton: x = ',num2str(round(x,4,'significant')),' # of iteration = ',num2str(iter)];
disp(sol);
end

function [xx,iter]=newton_verbose(x0,eps)
diff=1;
iter=0;
x=x0;
while abs(diff)>eps && iter<100
    fx=f(x);
    gx=g(x);
    diff=fx/gx;
    x=x-diff;
    iter=iter+1;
    xx(iter)=x;
end
sol=['Newton: x = ',num2str(round(x,4,'significant')),' # of iteration = ',num2str(iter)];
disp(sol);
end