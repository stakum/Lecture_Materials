a=1; b=3;
x0=(a+b)/2;
eps=1e-6;
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

function y=f(x)
 y=tanh(x)-x/2;
end
function y=g(x)
 y=1./(cosh(x).^2)-1/2;
end
