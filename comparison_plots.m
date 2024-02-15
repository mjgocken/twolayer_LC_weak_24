%% Comparison of minimum values through time for various values of thetaB
% Using theta-dependent evaporation function
% for dimensional data, the following units are used:
%   t_dim (s), h1_dim (micrometers), h2_dim (nanometers), p1_dim (Pa)
%   u2_dim (micrometers/min), g_dim (mol/m^2), c_dim (mOsM), y_dim (mm)

close all
clear
clc

% other files
% theta0 = load("LC_53timesCB_thetaDepEvap_thetaB0.mat");
% thetapi12 = load("LC_53timesCB_thetaDepEvap_thetaBpi12.mat");
% thetapi4 = load("LC_53timesCB_thetaDepEvap_thetaBpi4.mat");
% thetapi2 = load("LC_53timesCB_thetaDepEvap_thetaBpi2.mat");
% a = load("LC_53timesCB_thetaDepEvap_thetaB0_dim.mat");
% b = load("LC_53timesCB_thetaDepEvap_thetaBpi12_dim.mat");
% c = load("LC_53timesCB_thetaDepEvap_thetaBpi4_dim.mat");
% d = load("LC_53timesCB_thetaDepEvap_thetaBpi2_dim.mat");

% for osmolarity table
a = load("LC_3times5CB_thetaDepEvap_thetaB0_dim_finer.mat");
b = load("LC_3times5CB_thetaDepEvap_thetaBpi18_dim_finer.mat");
c = load("LC_3times5CB_thetaDepEvap_thetaBpi12_dim_finer.mat");
d = load("LC_3times5CB_thetaDepEvap_thetaBpi4_dim_finer.mat");
e = load("LC_3times5CB_thetaDepEvap_thetaBpi2_dim_finer.mat");
C = {a,b,c,d,e};

% for everything else
% a = load("LC_3times5CB_thetaDepEvap_thetaB0_dim_finer.mat");
% b = load("LC_3times5CB_thetaDepEvap_thetaBpi12_dim_finer.mat");
% c = load("LC_3times5CB_thetaDepEvap_thetaBpi4_dim_finer.mat");
% d = load("LC_3times5CB_thetaDepEvap_thetaBpi2_dim_finer.mat");
% e = load("newtonian_thetaDepEvap_thetaBpi2_dim.mat");
% C = {a,b,c,d};
half=floor(length(a.y_dim)/2); % for plotting half on yyplot

 % Minimum values through time
mk = {'--','-','-.',':','--','-','-.'};
fpath = '/Users/mjgocken/Desktop/figures';

% %min values of h1 thru time
% h=figure();
% hold on
% ax=gca;
% set(ax,'XMinorTick','on','YMinorTick','on','FontSize',32,'Units','pixels','Position', [150 100 383 385])
% for k = 1:numel(C)
%     v = C{k};
%     for j = 1:length(v.t_dim)
%         h1_min(j) = min(v.h1_dim(j,:));
%     end
%     fprintf("at 30s, h1min is %f, at 60s it is %f\n",h1_min(round(length(v.t_dim)/2)),h1_min(end))
%     plot(v.t_dim,h1_min,'LineWidth',3.5,'LineStyle',mk{k})
% end
% % lgd = legend('$0$','$\pi/12$','$\pi/4$','$\pi/2$','Interpreter','LaTeX','FontSize',29);
% % title(lgd,'$\theta_B$','Interpreter','LaTeX','FontSize',29);
% %     legend('$\theta_B=0$','$\theta_B=\pi/12$','$\theta_B=\pi/4$','$\theta_B=\pi/2$','Interpreter','LaTex')
% xlabel('$t$ (s)','FontSize',32,'Interpreter','LaTeX')
% ylabel('$h_{1min}$ ($\mu$m)','FontSize',32,'Interpreter','LaTeX')
% box on
% set(h, 'Position', [0 0 550 500])
% fname = 'h1min_v_time';
% saveas(gcf,fullfile(fpath,fname),'epsc');
% % % saveas(gcf,fullfile(fpath,fname),'fig');
% 
% 
% 
% %min values of h2 thru time
% h=figure();
% hold on
% ax=gca;
% set(ax,'XMinorTick','on','YMinorTick','on','FontSize',32,'Units','pixels','Position', [150 100 383 385])
% for k = 1:numel(C)
%     v = C{k};
%     for j = 1:length(v.t_dim)
%         h2_min(j) = min(v.h2_dim(j,:));
%     end
%     fprintf("at 30s, h2min is %f, at 60s it is %f\n",h2_min(round(length(v.t_dim)/2)),h2_min(end))
%     plot(v.t_dim,h2_min,'LineWidth',3.5,'LineStyle',mk{k})
% end
% lgd = legend('$0$','$\pi/18$','$\pi/12$','$\pi/4$','$\pi/2$','Interpreter','LaTeX','FontSize',29,'location','northwest');
% title(lgd,'$\theta_B$','Interpreter','LaTeX','FontSize',29);
% %     legend('$\theta_B=0$','$\theta_B=\pi/12$','$\theta_B=\pi/4$','$\theta_B=\pi/2$','Interpreter','LaTex')
% xlabel('$t$ (s)','FontSize',32,'Interpreter','LaTeX')
% ylabel('$h_{2min}$ (nm)','FontSize',32,'Interpreter','LaTeX')
% box on
% set(h, 'Position', [0 0 550 500])
% fname = 'h2min_v_time';
% saveas(gcf,fullfile(fpath,fname),'epsc');
% % % saveas(gcf,fullfile(fpath,fname),'fig');
% 
% % %peak values of u2 thru time
% % h=figure();
% % hold on
% % ax=gca;
% % set(ax,'XMinorTick','on','YMinorTick','on','FontSize',32,'Units','pixels','Position', [150 100 375 375])
% % for k = 1:numel(C)
% %     v = C{k};
% %     for j = 1:length(v.t_dim)
% %         u2_max(j) = max(v.u2_dim(j,:));
% %     end
% %     u2_max(1) = 62.1036;
% %     plot(v.t_dim,u2_max,'LineWidth',2.5,'LineStyle',mk{k})
% % end
% % ylim([0 65])
% % xlabel('$t$ (s)','FontSize',32,'Interpreter','LaTeX')
% % ylabel('$u_{2max}$ ($\mu$m/min)','FontSize',32,'Interpreter','LaTeX')
% % box on
% % set(h, 'Position', [0 0 550 500])
% % fname = 'u2max_v_time';
% % % saveas(gcf,fullfile(fpath,fname),'epsc');
% % % % saveas(gcf,fullfile(fpath,fname),'fig');
% % axes('Position',[.24 .755 .05 .05]); % to show break in y axis
% % px=[1 5];
% % py1=[1 2];
% % height=1;
% % py2=py1+height;
% % plot(px,py1,'k','LineWidth',2);hold all;
% % plot(px,py2,'k','LineWidth',2);hold all;
% % fill([px flip(px)],[py1 flip(py2)],'w','EdgeColor','none');
% % box off;
% % axis off;
% 
% %peak values of u2 thru time
% h=figure();
% hold on
% ax=gca;
% set(ax,'XMinorTick','on','YMinorTick','on','FontSize',32,'Units','pixels','Position', [150 100 383 385])
% for k = 1:numel(C)
%     v = C{k};
%     for j = 1:length(v.t_dim)
%         u2_max(j) = max(v.u2_dim(j,:));
%     end
%     plot(v.t_dim(2:end),u2_max(2:end),'LineWidth',3.5,'LineStyle',mk{k})
% end
% % ylim([0 65])
% xlabel('$t$ (s)','FontSize',32,'Interpreter','LaTeX')
% ylabel('$u_{2max}$ ($\mu$m/min)','FontSize',32,'Interpreter','LaTeX')
% box on
% set(h, 'Position', [0 0 550 500])
% fname = 'u2max_v_time';
% saveas(gcf,fullfile(fpath,fname),'epsc');
% % % saveas(gcf,fullfile(fpath,fname),'fig');
% 
% %min values of p1 thru time
% h=figure();
% hold on
% ax=gca;
% set(ax,'XMinorTick','on','YMinorTick','on','FontSize',32,'Units','pixels','Position', [150 100 383 385])
% for k = 1:numel(C)
%     v = C{k};
%     for j = 1:length(v.t_dim)
%         p1_min(j) = min(v.p1_dim(j,:));
%     end
%     plot(v.t_dim,p1_min,'LineWidth',3.5,'LineStyle',mk{k})
% end
% xlabel('$t$ (s)','FontSize',32,'Interpreter','LaTeX')
% ylabel('$p_{1min}$ (Pa)','FontSize',32,'Interpreter','LaTeX')
% box on
% set(h, 'Position', [0 0 550 500])
% fname = 'p1min_v_time';
% saveas(h,fullfile(fpath,fname),'epsc');
% % % saveas(gcf,fullfile(fpath,fname),'fig');
% 
% %peak values of c thru time
% h=figure();
% hold on
% ax=gca;
% set(ax,'XMinorTick','on','YMinorTick','on','FontSize',32,'Units','pixels','Position', [150 100 383 385])
% ax.FontSize = 32
% for k = 1:numel(C)
%     v = C{k};
%     for j = 1:length(v.t_dim)
%         c_max(j) = max(v.c_dim(j,:));
%     end
%     fprintf("at 30s, cmax is %f, at 60s it is %f\n",c_max(round(length(v.t_dim)/2)),c_max(end))
%     plot(v.t_dim,c_max,'LineWidth',3.5,'LineStyle',mk{k})
% end
% ylim([0 1700])
% % yticks([0 250 500 750 1000 1250 1500])
% % yticklabels()
% % % lgd = legend('$0$','$\pi/12$','$\pi/4$','$\pi/2$','Interpreter','LaTeX','FontSize',29,'location','northeast');
% % % title(lgd,'$\theta_B$','Interpreter','LaTeX','FontSize',29);
% % %     legend('$\theta_B=0$','$\theta_B=\pi/12$','$\theta_B=\pi/4$','$\theta_B=\pi/2$','Interpreter','LaTex')
% xlabel('$t$ (s)','FontSize',32,'Interpreter','LaTeX')
% ylabel('$c_{max}$ (mOsM)','FontSize',32,'Interpreter','LaTeX')
% box on
% set(h, 'Position', [0 0 550 500])
% fname = 'cmax_v_time.eps';
% % exportgraphics(h,fname,'Resolution',1200)
% saveas(gcf,fullfile(fpath,fname),'epsc');
% % % saveas(gcf,fullfile(fpath,fname),'fig');
% 
% %peak values of gamma thru time
% h=figure();
% hold on
% ax=gca;
% set(ax,'XMinorTick','on','YMinorTick','on','FontSize',32,'Units','pixels','Position', [150 100 383 385])
% for k = 1:numel(C)
%     v = C{k};
%     for j = 1:length(v.t_dim)
%         g_max(j) = max(v.g_dim(j,:))/4e-7; % nondimensionalize to give relative value
%     end
%     plot(v.t_dim,g_max,'LineWidth',3.5,'LineStyle',mk{k})
% end
% % lgd = legend('$0$','$\pi/12$','$\pi/4$','$\pi/2$','Interpreter','LaTeX','FontSize',29,'location','northeast');
% % title(lgd,'$\theta_B$','Interpreter','LaTeX','FontSize',29);
% %     legend('$\theta_B=0$','$\theta_B=\pi/12$','$\theta_B=\pi/4$','$\theta_B=\pi/2$','Interpreter','LaTex')
% xlabel('$t$ (s)','FontSize',32,'Interpreter','LaTeX')
% ylabel('$\Gamma_{max}$','FontSize',32,'Interpreter','LaTeX')
% box on
% set(h, 'Position', [0 0 550 500])
% fname = 'gmax_v_time.eps';
% % exportgraphics(h,fname,'Resolution',1200)
% saveas(gcf,fullfile(fpath,fname),'epsc');
% % saveas(gcf,fullfile(fpath,fname),'fig');
% 
% % FigWidth = 750;
% % FigHeight = 750;
% % BufferWidth = 200;
% % BufferHeight = 120;
% % F=figure('PaperPositionMode','manual');
% % AX = gca;
% % F.Resize = 'off';
% % F.Units = 'points';
% % F.PaperUnits = 'points';
% % AX.Units = 'points';
% % AX.Clipping = 'on';
% % AX.PositionConstraint = 'innerposition';
% % AX.InnerPosition = [BufferWidth BufferHeight FigWidth-BufferWidth-50 FigHeight-BufferHeight-120]*72/96; % converting from pixels to points
% % F.OuterPosition = [0 0 FigWidth FigHeight]*72/96; % converting from pixels to points
% % F.PaperPosition = [0 0 FigWidth FigHeight]*72/96;% converting from pixels to points
% % hold on
% % set(AX,'XMinorTick','on','YMinorTick','on','FontSize',32)
% % for k = 1:numel(C)
% %     v = C{k};
% %     for j = 1:length(v.t_dim)
% %         g_max(j) = max(v.g_dim(j,:))/4e-7; % nondimensionalize to give relative value
% %     end
% %     plot(v.t_dim,g_max,'LineWidth',2.5,'LineStyle',mk{k})
% % end
% % xlabel('$t$ (s)','FontSize',32,'Interpreter','LaTeX')
% % ylabel('$\Gamma_{max}$','FontSize',32,'Interpreter','LaTeX')
% % box on
% % fname = 'gmax_v_time';
% % % exportgraphics(F,fname,'Resolution',1200)
% % saveas(F,fullfile(fpath,fname),'epsc');

% final time plots
%h1
h=figure();
hold on
ax=gca;
set(ax,'XMinorTick','on','YMinorTick','on','FontSize',32,'Units','pixels','Position', [150 100 383 385])
mk = {'--','-','-.',':','--','-','-.'};
% plot(a.y_dim(1,:),a.h1_dim(1,:),'-k','LineWidth',2)
for k = 1:numel(C)
    v = C{k};
    plot(v.y_dim(end,:),v.h1_dim(end,:),'LineWidth',3.5,'LineStyle',mk{k})
end
ylim([0 4]);
xticks([-1 -0.5 0 0.5 1])
xlabel('$x$ (mm)','FontSize',32,'Interpreter','LaTeX')
ylabel('$h_1$ ($\mu$m)','FontSize',32,'Interpreter','LaTeX')
set(h, 'Position', [0 0 550 500])
box on
fname = 'h1_finaltime.eps';
saveas(gcf,fullfile(fpath,fname),'epsc');

%h2
h=figure();
hold on
ax=gca;
set(ax,'XMinorTick','on','YMinorTick','on','FontSize',32,'Units','pixels','Position', [150 100 383 385])
mk = {'--','-','-.',':','--','-','-.'};
% plot(a.y_dim(1,:),a.h1_dim(1,:),'-k','LineWidth',2)
for k = 1:numel(C)
    v = C{k};
    plot(v.y_dim(end,:),v.h2_dim(end,:),'LineWidth',3.5,'LineStyle',mk{k})
end
lgd = legend('$0$','$\pi/18$','$\pi/12$','$\pi/4$','$\pi/2$','Interpreter','LaTeX','FontSize',29,'location','southeast');
title(lgd,'$\theta_B$','Interpreter','LaTeX','FontSize',29);
legend('Interpreter','LaTeX','FontSize',29);
xticks([-1 -0.5 0 0.5 1])
% ylim([0 50]);
xlabel('$x$ (mm)','FontSize',32,'Interpreter','LaTeX')
ylabel('$h_2$ (nm)','FontSize',32,'Interpreter','LaTeX')
set(h, 'Position', [0 0 550 500])
box on
fname = 'h2_finaltime.eps';
saveas(gcf,fullfile(fpath,fname),'epsc');

% %h1 and h2
% figure();
% ax=gca;
% ax.FontSize = 32;
% mk = {'--','-','-.',':','--','-','-.'};
% yyaxis left
% hold on
% for k = 1:numel(C)
%     v = C{k};
%     plot(v.y_dim(end,1:half),v.h1_dim(end,1:half),'LineWidth',2.5,'LineStyle',mk{k})
% end
% xlabel('$x$ (mm)','FontSize',32,'Interpreter','LaTeX')
% ylabel('$h_1$ ($\mu$m)','FontSize',32,'Interpreter','LaTeX')
% colororder('default');
% xticks([-1 -0.5 0 0.5 1])
% ylim([0 3])
% yyaxis right
% hold on
% for k = 1:numel(C)
%     v = C{k};
%     plot(v.y_dim(end,half:end),v.h2_dim(end,half:end),'LineWidth',2.5,'LineStyle',mk{k})
% end
% lgd = legend('$0$','$\pi/12$','$\pi/4$','$\pi/2$','Interpreter','LaTeX','FontSize',29,'location','southeast');
% title(lgd,'$\theta_B$','Interpreter','LaTeX','FontSize',29);
% legend('Interpreter','LaTeX','FontSize',29);
% % ylim([0 50]);
% xlabel('$x$ (mm)','FontSize',32,'Interpreter','LaTeX')
% ylabel('$h_2$ (nm)','FontSize',32,'Interpreter','LaTeX')
% hold off
% colororder('default');
% ylim([0 60])
% ax.YAxis(1).Color = 'k';
% ax.YAxis(2).Color = 'k';

%p1
h=figure();
hold on
ax=gca;
set(ax,'XMinorTick','on','YMinorTick','on','FontSize',32,'Units','pixels','Position', [150 100 383 385])
mk = {'--','-','-.',':','--','-','-.'};
for k = 1:numel(C)
    v = C{k};
    plot(v.y_dim(end,:),v.p1_dim(end,:),'LineWidth',3.5,'LineStyle',mk{k})
end
ylim([-7 2]);
xticks([-1 -0.5 0 0.5 1])
xlabel('$x$ (mm)','FontSize',32,'Interpreter','LaTeX')
ylabel('$p_1$ (Pa)','FontSize',32,'Interpreter','LaTeX')
set(h, 'Position', [0 0 550 500])
box on
fname = 'p1_finaltime.eps';
saveas(gcf,fullfile(fpath,fname),'epsc');

%u2
h=figure();
hold on
ax=gca;
set(ax,'XMinorTick','on','YMinorTick','on','FontSize',32,'Units','pixels','Position', [150 100 383 385])
mk = {'--','-','-.',':','--','-','-.'};
for k = 1:numel(C)
    v = C{k};
    plot(v.y_dim(end,:),v.u2_dim(end,:),'LineWidth',3.5,'LineStyle',mk{k})
end
% ylim([0 4]);
xticks([-1 -0.5 0 0.5 1])
xlabel('$x$ (mm)','FontSize',32,'Interpreter','LaTeX')
ylabel('$u_2$ ($\mu$m/min)','FontSize',32,'Interpreter','LaTeX')
set(h, 'Position', [0 0 550 500])
box on
fname = 'u2_finaltime.eps';
saveas(gcf,fullfile(fpath,fname),'epsc');

%c
h=figure();
hold on
ax=gca;
set(ax,'XMinorTick','on','YMinorTick','on','FontSize',32,'Units','pixels','Position', [150 100 383 385])
mk = {'--','-','-.',':','--','-','-.'};
for k = 1:numel(C)
    v = C{k};
    plot(v.y_dim(end,:),v.c_dim(end,:),'LineWidth',3.5,'LineStyle',mk{k})
    fprintf('max osmolarity is %f, min is %f, and diff is %f\n',max(v.c_dim(end,:)),min(v.c_dim(end,:)),max(v.c_dim(end,:))-min(v.c_dim(end,:)))
end
ylim([0 1800]);
xticks([-1 -0.5 0 0.5 1])
xlabel('$x$ (mm)','FontSize',32,'Interpreter','LaTeX')
ylabel('$c$ (mOsM)','FontSize',32,'Interpreter','LaTeX')
set(h, 'Position', [0 0 550 500])
box on
fname = 'c_finaltime.eps';
saveas(gcf,fullfile(fpath,fname),'epsc');


%g
% h=figure();
% hold on
% ax=gca;
% set(ax,'XMinorTick','on','YMinorTick','on','FontSize',32,'Units','pixels','Position', [150 100 375 375])
% mk = {'--','-','-.',':','--','-','-.'};
% for k = 1:numel(C)
%     v = C{k};
%     plot(v.y_dim(end,:),v.g_dim(end,:),'LineWidth',2.5,'LineStyle',mk{k})
% end
% lgd = legend('$0$','$\pi/12$','$\pi/4$','$\pi/2$','Interpreter','LaTeX','FontSize',29,'location','northeast');
% title(lgd,'$\theta_B$','Interpreter','LaTeX','FontSize',29);
% legend('Interpreter','LaTeX','FontSize',29);
% % ylim([0 4]);
% xlabel('$x$ (mm)','FontSize',32,'Interpreter','LaTeX')
% ylabel('$\Gamma$ (mol/m$^2$)','FontSize',32,'Interpreter','LaTeX')
% set(h, 'Position', [0 0 550 500])

%g
h=figure();
hold on
ax=gca;
set(ax,'XMinorTick','on','YMinorTick','on','FontSize',32,'Units','pixels','Position', [150 100 383 385])
mk = {'--','-','-.',':','--','-','-.'};
for k = 1:numel(C)
    v = C{k};
    plot(v.y_dim(end,:),v.g_dim(end,:)/4e-7,'LineWidth',3.5,'LineStyle',mk{k})%nondimensional
end
ylim([0.999 1.002]);
xticks([-1 -0.5 0 0.5 1])
xlabel('$x$ (mm)','FontSize',32,'Interpreter','LaTeX')
ylabel('$\Gamma$','FontSize',32,'Interpreter','LaTeX')
set(h, 'Position', [0 0 550 500])
box on
fname = 'g_finaltime.eps';
saveas(gcf,fullfile(fpath,fname),'epsc');


%% max osm table

for k=1:numel(C)
    v=C{k};
    if length(v.t_dim)==160
        fprintf("max osmolarity at: t = %d is %.1f, t = %d is %.1f\n",v.t_dim(151),max(v.c_dim(151,:)),v.t_dim(end),max(v.c_dim(end,:)));
    else
        fprintf("max osmolarity at: t = %d is %.1f, t = %d is %.1f, t = %d is %.1f\n",v.t_dim(151),max(v.c_dim(151,:)),v.t_dim(301),max(v.c_dim(301,:)),v.t_dim(end),max(v.c_dim(end,:)));
    end
end


