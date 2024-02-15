% thetaB v h1min, h2min, cmax and tbut

close all
clear
clc

% load files
a = load("LC_3times5CB_thetaDepEvap_extrema_v_theta_dim.mat");
[Max,Ind] = max(a.h2min_dim) % find thetaB val of break up

% AL min at final time
% figure();
% plot(a.thetaBval,a.h1min_dim,'linewidth',2)
% xlabel('$\theta_B$','Interpreter','LaTex')
% ylabel('Min AL thickness $(\mu m)$','Interpreter','LaTex')
% set(gca,'FontSize',20)
% set(gca,'XTick',[0,pi/6,pi/3,pi/2]) 
% set(gca,'XTickLabel',{'0','pi/6','pi/3','pi/2'})

% AL min at final time
h=figure();
hold on
ax=gca;
set(ax,'XMinorTick','on','YMinorTick','on','FontSize',32,'Units','pixels','Position', [80 100 375 375])
% set(h,'defaultAxesColorOrder',[[0.4940 0.1840 0.5560]; [0 0.4470 0.7410]]);
yyaxis right
plot(a.thetaBval,a.h1min_dim,'linewidth',3.5,'Color',[0 0.4470 0.7410])
xlabel('$\theta_B$','FontSize',32,'Interpreter','LaTeX')
ylabel('$h_{1min}$ ($\mu$m)','FontSize',32,'Interpreter','LaTeX')
ylim([0.45,2.1])
yyaxis left
plot(a.TBUTthetaB,a.TBUT_dim,'linewidth',3.5,'Color',[0.8500 0.3250 0.0980])
ylabel('TBUT (s)','FontSize',32,'Interpreter','LaTeX')
ylim([0,60])
set(ax,'XTick',[0,pi/6,pi/3,pi/2]) 
set(ax,'XTickLabel',{'$0$','$\pi/6$','$\pi/3$','$\pi/2$'})
set(ax,'TickLabelInterpreter', 'latex')
set(h, 'Position', [0 0 550 500])
box on
ax.YAxis(1).Color = [0.8500 0.3250 0.0980];
ax.YAxis(2).Color = [0 0.4470 0.7410];
xline(a.thetaBval(Ind),'--k','linewidth',3.5)


% LL min at final time

h=figure();
hold on
ax=gca;
set(ax,'XMinorTick','on','YMinorTick','on','FontSize',32,'Units','pixels','Position', [140 100 375 375])
plot(a.thetaBval,a.h2min_dim,'linewidth',3.5)
xlabel('$\theta_B$','FontSize',32,'Interpreter','LaTeX')
ylabel('$h_{2min}$ (nm)','FontSize',32,'Interpreter','LaTeX')
ylim([6 18])
set(h, 'Position', [0 0 550 500])
box on
set(gca,'XTick',[0,pi/6,pi/3,pi/2]) 
set(ax,'XTickLabel',{'$0$','$\pi/6$','$\pi/3$','$\pi/2$'})
set(ax,'TickLabelInterpreter', 'latex')
xline(a.thetaBval(Ind),'--k','linewidth',3.5)



% figure();
% yyaxis right
% plot(a.thetaBval,a.h2min_dim,'linewidth',2)
% xlabel('\theta_B')
% ylabel('Min LL thickness (nm)')
% ylim([6,22])
% yyaxis left
% plot(a.TBUTthetaB,a.TBUT_dim,'linewidth',2)
% ylabel('TBUT (s)')
% ylim([0,60])
% set(gca,'FontSize',20)
% set(gca,'XTick',[0,pi/6,pi/3,pi/2]) 
% set(gca,'XTickLabel',{'0','\pi/6','\pi/3','\pi/2'})

% max osmolarity at final time
h=figure();
hold on
ax=gca;
set(ax,'XMinorTick','on','YMinorTick','on','FontSize',32,'Units','pixels','Position', [140 100 375 375])
plot(a.thetaBval,a.cmax_dim,'linewidth',3.5)
xlabel('$\theta_B$','FontSize',32,'Interpreter','LaTeX')
ylabel('$c_{max}$ (mOsM)','FontSize',32,'Interpreter','LaTeX')
ylim([400 1800])
set(h, 'Position', [0 0 550 500])
box on
set(gca,'XTick',[0,pi/6,pi/3,pi/2]) 
% set(ax,'XTickLabel',{'$0$','$\frac{\pi}{6}$','$\frac{\pi}{3}$','$\frac{\pi}{2}$'})
set(ax,'XTickLabel',{'$0$','$\pi/6$','$\pi/3$','$\pi/2$'})
set(ax,'TickLabelInterpreter', 'latex')
xline(a.thetaBval(Ind),'--k','linewidth',3.5)
text(25,1500,'TBUT > 60 s')

% figure();
% yyaxis right
% plot(a.thetaBval,a.cmax_dim,'linewidth',2)
% xlabel('\theta_B')
% ylabel('Max Osmolarity (mOsM)')
% ylim([500,1700])
% yyaxis left
% plot(a.TBUTthetaB,a.TBUT_dim,'linewidth',2)
% ylabel('TBUT (s)')
% ylim([0,60])
% set(gca,'FontSize',20)
% set(gca,'XTick',[0,pi/6,pi/3,pi/2]) 
% set(gca,'XTickLabel',{'0','\pi/6','\pi/3','\pi/2'})
