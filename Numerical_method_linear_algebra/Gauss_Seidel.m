[n m]=size(A);
D=zeros(n,n);
L=zeros(n,n);
U=zeros(n,n);
for ind=1:n
    D(ind,ind)=A(ind,ind);
    for jnd=ind+1:n
        U(ind,jnd)=A(ind,jnd);
    end
end
L=A-D-U;
H=-inv(D+L)*U;
c=inv(D+L)*b;
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