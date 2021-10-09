A=[4 1 1; 1 3 1; 2 1 4];
%% LU decomposition
L=eye(size(A));
U=A;
[n,m]=size(A);
for k = 1:n-1
    for i = k+1:n
        L(i,k)=U(i,k)/U(k,k);
        for j = 1:n % k+1:n
            U(i,j)=U(i,j)-L(i,k)*U(k,j);
        end
    end
end