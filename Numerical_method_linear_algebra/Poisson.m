clear all; close all; clc;
%% 

m=30;
F=zeros(m,m);
F(2*m/5:3*m/5,2*m/5:3*m/5)=1;
K=1e6;
tol=1e-8;
VV=fastpoisson(F);
[k_J,V_J]=jdp(F,K,tol);
w=2/(1+sin(pi/(m+1)));
[k_SOR,V_SOR]=sordp(F,K,w,tol);

k_J
k_SOR

N=[10 30 50 100];
N_J=[342 2684 7226 28207];
N_SOR=[31 86 142 279];
f=figure; ax=gca; box on; grid on
semilogy(N,N_J,'b--','LineWidth',2.0)
hold on
semilogy(N,N_SOR,'r-','LineWidth',2.0)
legend('Jacobi','SOR')
ylabel('Number of iterations')
xlabel('Number of bins')
set(ax,'FontSize',18);
set(ax,'FontName','Arial')

f=figure; ax=gca; box on; grid on
surf(V_SOR)
xlim([1 m])
ylim([1 m])
ylabel('{\it Y}')
xlabel('{\it X}')
set(ax,'FontSize',18);
set(ax,'FontName','Arial')

%% Fast Poisson
function V=fastpoisson(F)
    %function V=fastpoisson(F)
    m=length(F); h=1/(m+1); hv=pi*h*(1:m)';
    sigma=sin(hv/2).^2;
    S=sin(hv*(1:m));
    G=S*F*S;
    X=h^4*G./(sigma*ones(1,m)+ ones(m,1)*sigma');
    V=zeros(m+2,m+2);
    V(2:m+1,2:m+1)=S*X*S;
end

%% Jacobi method
function [k,V]=jdp(F,K,tol)
    % k=jdp(F,K,tol)
    m=length(F); U=fastpoisson(F);
    V=zeros(m+2,m+2); E=F/(m+1)^2;
    for k=1:K
        V(2:m+1,2:m+1)=(V(1:m,2:m+1)+V(3:m+2,2:m+1)...
            +V(2:m+1,1:m)+V(2:m+1,3:m+2)+E)/4;
        if max(max(abs(V-U)))<tol, return
        end
    end
    k=K+1;
end

%% SOR method
function [k,V]=sordp(F,K,w,tol)
    % k=sordp(F,K,w,tol)
    m=length(F); U=fastpoisson(F); V=zeros(m+2,m+2); E=F/(m+1)^2;
    for k=1:K
        for j=2:m+1
            for i=2:m+1
                V(i,j)=w*(V(i-1,j)+V(i+1,j)+V(i,j-1)...
                    +V(i,j+1)+E(i-1,j-1))/4+(1-w)*V(i,j);
            end
        end
        if max(max(abs(V-U)))<tol, return
        end
    end
k=K+1;
end