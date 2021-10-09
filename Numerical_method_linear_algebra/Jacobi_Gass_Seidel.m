clear all; close all; clc; 
%%
A=[4 1 1;
   2 6 3;
   -1 4 7];
b=[15;12;21];

x0=[0;0;0];
[x_J,residual_J,iter_J]=Jacobi(A,b,x0);
[x_GS,residual_GS,iter_GS]=Gauss_Seidel(A,b,x0);

f=figure; ax=gca; box on; grid on
semilogy(residual_J,'b--','LineWidth',2.0)
hold on
semilogy(residual_GS,'r-','LineWidth',2.0)
yline(1e-6,'k-.','LineWidth',1.5)
legend('Jacobi','Gauss-Seidel')
xlabel('Number of iterations')
ylim([1e-8 10])
ylabel('Residual')
set(ax,'FontSize',18);
set(ax,'FontName','Arial')
set(ax,'LineWidth',1.5)

function [x,residual,iter]=Jacobi(A,b,x0)
    [n m]=size(A);
    if n~=m
        disp('input matrix should be square one')
        return
    end

    D=zeros(n,n);
    Dinv=zeros(n,n);
    for ind=1:n
        D(ind,ind)=A(ind,ind);
        Dinv(ind,ind)=1/D(ind,ind);
    end
    N=A-D;
    H=-Dinv*N;
    res=['Spectral radius of matrix H =',num2str(max(abs(eig(H))))];
    disp(res);
    c=Dinv*b;
    iter_max=50;
    eps=1e-6;
    iter=0;
    while iter<iter_max
        iter=iter+1;
        x=H*x0+c;
        residual(iter)=abs(norm(x-x0));
        if(abs(norm(x-x0))<eps)
            return
        end
        x0=x;
    end

end

function [x,residual,iter]=Gauss_Seidel(A,b,x0)
    [n m]=size(A);
    if n~=m
        disp('input matrix should be square one')
        return
    end
    
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
    res=['Spectral radius of matrix H =',num2str(max(abs(eig(H))))];
    disp(res);
    c=inv(D+L)*b;
    iter_max=50;
    eps=1e-6;
    iter=0;
    while iter<iter_max
        iter=iter+1;
        x=H*x0+c;
        residual(iter)=abs(norm(x-x0));
        if(abs(norm(x-x0))<eps)
            return
        end
        x0=x;
    end

end


