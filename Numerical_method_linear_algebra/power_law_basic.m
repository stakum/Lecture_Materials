n=10;
A=pascal(n); 
x0=ones(n,1);
lambda0=dot(x0,A*x0);
iter=0;

while iter<50
    x=A*x0;
    lambda=dot(x,x)/dot(x,x0);
    if abs(lambda-lambda0)<1e-6
        return
    end
    lambda0=lambda;
    x0=x;
    iter=iter+1;
end
e=x/norm(x);