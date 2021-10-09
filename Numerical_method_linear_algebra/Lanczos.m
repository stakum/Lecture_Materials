clear all; close all; clc;
%% Lanczos method
n=10;
A=pascal(n);
u=zeros(n,n);
u(1,1)=-1;
alpha(1)=dot(u(:,1),A*u(:,1));
v=A*u(:,1)-alpha(1)*u(:,1);
beta(1)=norm(v);
u(:,2)=v/beta(1);

for k=2:n-1
    alpha(k)=dot(u(:,k),A*u(:,k));
    v=A*u(:,k)-beta(k-1)*u(:,k-1)-alpha(k)*u(:,k);
    beta(k)=norm(v);
    u(:,k+1)=v/beta(k);
end

alpha(n)=dot(u(:,n),A*u(:,n));

B=diag(alpha)+diag(beta(1:n-1),-1)+diag(beta(1:n-1),1);