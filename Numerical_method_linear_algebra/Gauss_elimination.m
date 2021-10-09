A=[4 1 1; 1 3 1; 2 1 4];
b=[9; 10; 16];
%% Gauss elimination process
[n,m]=size(A);
for k = 1:n-1
    for i = k+1:n
        mik=A(i,k)/A(k,k);
        for j = 1:n % k+1:n
            A(i,j)=A(i,j)-mik*A(k,j);
        end
        b(i)=b(i)-mik*b(k);
    end
end

%% Backward substition
x=b;
x(n)=b(n)/A(n,n);
for k = n-1:-1:1
%     x(k)=(b(k)-A(k,k+1:n)*x(k+1:n))/A(k,k);
    x(k)=b(k)/A(k,k);
    for l= k+1:n
        x(k)=x(k)-A(k,l)*x(l)/A(k,k);
    end        
end
