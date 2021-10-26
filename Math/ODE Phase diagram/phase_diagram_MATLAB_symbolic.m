Z0=[(-5:5); (-5:5)];
x0=[-5 -2];

syms x(t) y(t)
A=[1 -2; 1 -1];B=[1; 3];
Z = [x; y];sol = diff(Z) == A*Z + B;
[xSol(t), ySol(t)] = dsolve(sol);
xSol(t) = simplify(xSol(t))
ySol(t) = simplify(ySol(t))

f=figure;
ax=gca; 
box on; 
hold on
for ind=1:length(Z0)
    C = Z(0) == Z0(:,ind);
    [xSol(t),ySol(t)]=dsolve(sol,C);
    fplot(xSol,ySol,'LineWidth',2)
end
plot(x0(:,1),x0(:,2),'k*','MarkerFaceColor','k','MarkerSize',12) % 平衡点のプロット
grid on
xlabel('$$x_1$$','InterPreter','latex')
ylabel('$$x_2$$','InterPreter','latex')
set(ax,'FontSize',18);
set(ax,'FontName','Arial')
set(ax,'LineWidth',1.5)
xlim([-16 6]);
ylim([-10 6])
