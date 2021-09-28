clear all; close all; clc;
%%

x=linspace(1e-3,2,100);
Bose=1./(exp(x)-1);
Maxwell=exp(-x);

f=figure; 
ax=gca;
semilogy(x,Bose,'LineWidth',2)
hold on
semilogy(x,Maxwell,'LineWidth',2)
grid on
set(ax,'FontSize',18);
set(ax,'FontName','Arial')
set(ax,'LineWidth',1.5)
legend('Bose-Einstein','Maxwell')
xlabel('$\hbar\omega_{k}/k_\mathrm{B}T$','Interpreter','latex')
ylabel('$n_k^0$','Interpreter','latex')