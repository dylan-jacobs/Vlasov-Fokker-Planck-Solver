

figure(2);clf;loglog(Lambdavals,E1,'k-','linewidth',1.5);
hold on;loglog(Lambdavals,E2,'b-','linewidth',1.5);
hold on;loglog(Lambdavals,E3,'m-','linewidth',1.5);
hold on;loglog(Lambdavals(ceil(0.35*end):ceil(0.75*end)),.7*Lambdavals(ceil(0.35*end):ceil(0.75*end)).^1,'k-.','linewidth',1.5);
hold on;loglog(Lambdavals(ceil(0.35*end):ceil(0.75*end)),0.01*Lambdavals(ceil(0.35*end):ceil(0.75*end)).^2,'b-.','linewidth',1.5);
hold on;loglog(Lambdavals(ceil(0.35*end):ceil(0.75*end)),0.0012*Lambdavals(ceil(0.35*end):ceil(0.75*end)).^3,'m-.','linewidth',1.5);
xlabel('\lambda');ylabel('L^1 error');
%title(['Rigid body rotation solved to T_f = ',num2str(Tf)]);
legend('IMEX111','IMEX222','IMEX443','Order 1','Order 2','Order 3','location','eastoutside');
axis([Lambdavals(1),Lambdavals(end),1.0e-6,1.0e1]);
fontsize(18,"points");
set(gcf,'Units','pixels','Position',[100 100 800 500]);

% exportgraphics(gcf,'test9_error.pdf','ContentType','vector')

