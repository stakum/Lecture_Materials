[n m]=size(A);
D=zeros(n,n);
Dinv=zeros(n,n);
for ind=1:n
    D(ind,ind)=A(ind,ind);
    Dinv(ind,ind)=1/D(ind,ind);
end
N=A-D;
H=-Dinv*N;
c=Dinv*b;
iter_max=50;
eps=1e-6;
iter=0;
while iter<iter_max
    iter=iter+1;
    x=H*x0+c;
    if(abs(norm(x-x0))<eps)
        return
    end
    x0=x;
end