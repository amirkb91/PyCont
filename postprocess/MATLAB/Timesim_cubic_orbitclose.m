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
a1.YLimMode = 'auto';
a1.XLimitMethod = 'tight';
xlabel(a1, '$x_2$ (m)');
ylabel(a1, '$\dot{x}_1$ (m/s)');

f2 = figure('Units','centimeters','PaperUnits','centimeters',...
    'Position',[2 2 fsizex fsizey],'PaperPosition',[0 0 fsizex fsizey],...
    'PaperSize',[fsizex fsizey],'PaperPositionMode','manual');
a2 = axes(f2,'FontSize',afont,'LineWidth',aline,...
    'TickLabelInterpreter','latex','Box','on');
a2.XLabel.Interpreter = "latex";
a2.YLabel.Interpreter = "latex";
% legend(a1,'Interpreter','latex','Location','northeast');
grid(a2, 'on'); hold(a2, 'on');
a2.YLimMode = 'auto';
a2.XLimitMethod = 'tight';
xlabel(a2, '$x_2$ (m)');
ylabel(a2, '$\dot{x}_1$ (m/s)');
% 

%%
t0a = readmatrix('corr0_partition_0_t.txt');
t1a = readmatrix('corr0_partition_1_t.txt');
t2a = readmatrix('corr0_partition_2_t.txt');
p0a = readmatrix('corr0_partition_0_pose_time.txt');
p1a = readmatrix('corr0_partition_1_pose_time.txt');
p2a = readmatrix('corr0_partition_2_pose_time.txt');
v0a = readmatrix('corr0_partition_0_vel_time.txt');
v1a = readmatrix('corr0_partition_1_vel_time.txt');
v2a = readmatrix('corr0_partition_2_vel_time.txt');

t0b = readmatrix('corr2_partition_0_t.txt');
t1b = readmatrix('corr2_partition_1_t.txt');
t2b = readmatrix('corr2_partition_2_t.txt');
p0b = readmatrix('corr2_partition_0_pose_time.txt');
p1b = readmatrix('corr2_partition_1_pose_time.txt');
p2b = readmatrix('corr2_partition_2_pose_time.txt');
v0b = readmatrix('corr2_partition_0_vel_time.txt');
v1b = readmatrix('corr2_partition_1_vel_time.txt');
v2b = readmatrix('corr2_partition_2_vel_time.txt');

%%

plot(a1,p0a(:,2),v0a(:,1),'-','LineWidth',plotlinew,'Color',"#0072BD");
plot(a1,p1a(:,2),v1a(:,1),'-','LineWidth',plotlinew,'Color',"#D95319");
plot(a1,p2a(:,2),v2a(:,1),'-','LineWidth',plotlinew,'Color',"#77AC30");
a1.XLim=[-24 -20];

plot(a2,p0b(:,2),v0b(:,1),'-','LineWidth',plotlinew,'Color',"#0072BD");
plot(a2,p1b(:,2),v1b(:,1),'-','LineWidth',plotlinew,'Color',"#D95319");
plot(a2,p2b(:,2),v2b(:,1),'-','LineWidth',plotlinew,'Color',"#77AC30");
a2.XLim=[-24 -20];
