a0=-0.25; b0=1.25;
eps=1e-6; delta=1e-10;
iter=1;
x=(a+b)/2;
fx=f(x); fa=f(a); fb=f(b);
while (b-a)>eps || abs(fx)>delta
    if(fx*fa<0)
        a=a; b=x;
    elseif(fx*fb<0)
        a=x; b=b;
    end
    x=(a+b)/2;
    fx=f(x);
    fa=f(a);
    fb=f(b);
    iter=iter+1;
end
alpha=x;

function y=f(x)
y=x.^2-1;
end