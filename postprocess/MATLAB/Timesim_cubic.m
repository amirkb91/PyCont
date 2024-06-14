clear all
close all
%% Figure Properties
fsizex  = 23;
ratio = 1;
fsizey  = fsizex*ratio;
afont = 28;
aline = 2;
plotlinew = 2.5;

f1 = figure('Units','centimeters','PaperUnits','centimeters',...
    'Position',[2 2 fsizex fsizey],'PaperPosition',[0 0 fsizex fsizey],...
    'PaperSize',[fsizex fsizey],'PaperPositionMode','manual');
a1 = axes(f1,'FontSize',afont,'LineWidth',aline,...
    'TickLabelInterpreter','latex','Box','on');
a1.XLabel.Interpreter = "latex";
a1.YLabel.Interpreter = "latex";
% legend(a1,'Interpreter','latex','Location','northeast');
grid(a1, 'on'); hold(a1, 'on');
% a1.YLimMode = 'auto';
a1.XLimitMethod = 'tight';
xlabel(a1, 'Time (s)');
ylabel(a1, 'Position (m)');
a1.YTick = [-50 -40 -30 -20 -10 0 10 20 30 40 50];

f2 = figure('Units','centimeters','PaperUnits','centimeters',...
    'Position',[2 2 fsizex fsizey],'PaperPosition',[0 0 fsizex fsizey],...
    'PaperSize',[fsizex fsizey],'PaperPositionMode','manual');
a2 = axes(f2,'FontSize',afont,'LineWidth',aline,...
    'TickLabelInterpreter','latex','Box','on');
a2.XLabel.Interpreter = "latex";
a2.YLabel.Interpreter = "latex";
% legend(a1,'Interpreter','latex','Location','northeast');
grid(a2, 'on'); hold(a2, 'on');
% a2.YLimMode = 'auto';
a2.XLimitMethod = 'tight';
xlabel(a2, 'Position $x_1$ (m)');
ylabel(a2, 'Position $x_2$ (m)');
a2.YTick = [-50 -40 -30 -20 -10 0 10 20 30 40 50];

f3 = figure('Units','centimeters','PaperUnits','centimeters',...
    'Position',[2 2 fsizex fsizey],'PaperPosition',[0 0 fsizex fsizey],...
    'PaperSize',[fsizex fsizey],'PaperPositionMode','manual');
a3 = axes(f3,'FontSize',afont,'LineWidth',aline,...
    'TickLabelInterpreter','latex','Box','on');
a3.XLabel.Interpreter = "latex";
a3.YLabel.Interpreter = "latex";
% legend(a1,'Interpreter','latex','Location','northeast');
grid(a3, 'on'); hold(a3, 'on');
% a3.YLimMode = 'auto';
a3.XLimitMethod = 'tight';
xlabel(a3, 'Time (s)');
ylabel(a3, 'Velocity (m/s)');
% a3.YTick = [-50 -40 -30 -20 -10 0 10 20 30 40 50];


%%
t0 = readmatrix('partition_0_t.txt');
t1 = readmatrix('partition_1_t.txt');
t2 = readmatrix('partition_2_t.txt');
p0 = readmatrix('partition_0_pose_time.txt');
p1 = readmatrix('partition_1_pose_time.txt');
p2 = readmatrix('partition_2_pose_time.txt');
v0 = readmatrix('partition_0_vel_time.txt');
v1 = readmatrix('partition_1_vel_time.txt');
v2 = readmatrix('partition_2_vel_time.txt');

%%

plot(a1,t0,p0(:,1),'-','LineWidth',plotlinew,'Color',"#0072BD");
plot(a1,t1,p1(:,1),'-','LineWidth',plotlinew,'Color',"#D95319");
plot(a1,t2,p2(:,1),'-','LineWidth',plotlinew,'Color',"#77AC30");
plot(a1,t0,p0(:,2),'--','LineWidth',plotlinew,'Color',"#0072BD");
plot(a1,t1,p1(:,2),'--','LineWidth',plotlinew,'Color',"#D95319");
plot(a1,t2,p2(:,2),'--','LineWidth',plotlinew,'Color',"#77AC30");

plot(a2,p0(:,1),p0(:,2),'-','LineWidth',plotlinew,'Color',"#0072BD");
plot(a2,p1(:,1),p1(:,2),'-','LineWidth',plotlinew,'Color',"#D95319");
plot(a2,p2(:,1),p2(:,2),'-','LineWidth',plotlinew,'Color',"#77AC30");

plot(a3,t0,v0(:,1),'-','LineWidth',plotlinew,'Color',"#0072BD");
plot(a3,t1,v1(:,1),'-','LineWidth',plotlinew,'Color',"#D95319");
plot(a3,t2,v2(:,1),'-','LineWidth',plotlinew,'Color',"#77AC30");
plot(a3,t0,v0(:,2),'--','LineWidth',plotlinew,'Color',"#0072BD");
plot(a3,t1,v1(:,2),'--','LineWidth',plotlinew,'Color',"#D95319");
plot(a3,t2,v2(:,2),'--','LineWidth',plotlinew,'Color',"#77AC30");
