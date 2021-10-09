clear all; close all; clc;
%% Bisection method
a0=-0.25;
b0=1.25;

[alpha,xx,aa,bb]=bisection_verbose(a0,b0);
tmp=linspace(0,3,50);

subplot(1,2,1)
plot(xx,'bo-')
xlabel('Iteration step')
ylabel('x')
subplot(1,2,2)
plot(tmp,f(tmp)) 
hold on
plot(xx,f(xx),'bo')
xlabel('x')
ylabel('f(x)')

x=linspace(-1,1.5,50);

figure
plot(x,f(x))
hold on
plot(aa(1),f(aa(1)),'bo')
plot(bb(1),f(bb(1)),'bo')
plot(xx(1:4),f(xx(1:4)),'rs')


function [alpha,xx,aa,bb]=bisection_verbose(a,b)
    eps=1e-6;
    delta=1e-10;
    iter=1;
    x=(a+b)/2;
    aa(iter)=a;
    bb(iter)=b;
    xx(iter)=x;
    fx=f(x);
    fa=f(a);
    fb=f(b);
    while (b-a)>eps || abs(fx)>delta
        if(fx*fa<0)
            a=a;
            b=x;
        elseif(fx*fb<0)
            a=x;
            b=b;
        end
        x=(a+b)/2;
        fx=f(x);
        fa=f(a);
        fb=f(b);
        iter=iter+1;
        xx(iter)=x;
        aa(iter)=a;
        bb(iter)=b;
    end
    alpha=x;
end

function [alpha,xx]=bisection(a,b)
    eps=1e-6;
    delta=1e-10;
    iter=1;
    x=(a+b)/2;
    xx(iter)=x;
    fx=f(x);
    fa=f(a);
    fb=f(b);
    while (b-a)>eps || abs(fx)>delta
        if(fx*fa<0)
            a=a;
            b=x;
        elseif(fx*fb<0)
            a=x;
            b=b;
        end
        x=(a+b)/2;
        fx=f(x);
        fa=f(a);
        fb=f(b);
        iter=iter+1;
        xx(iter)=x;
    end
    alpha=x;
end

function y=f(x)
    %y=tanh(x)-x/2;
    y=x.^2-1;
end