clear all; 
close all;

%% PLOT SWEEP RESULTS

%% Figure Properties
fsizex  = 25;
ratio = 3/4;
fsizey  = fsizex*ratio;
afont = 23;
aline = 2;
plotlinew = 2.5;

f1 = figure('Units','centimeters','PaperUnits','centimeters',...
    'Position',[2 2 fsizex fsizey],'PaperPosition',[0 0 fsizex fsizey],...
    'PaperSize',[fsizex fsizey],'PaperPositionMode','manual');
a1 = axes(f1,'FontSize',afont,'LineWidth',aline,...
    'TickLabelInterpreter','latex','Box','on');
a1.XLabel.Interpreter = "latex";
a1.YLabel.Interpreter = "latex";
grid(a1, 'on'); hold(a1, 'on');
a1.YLimMode = 'auto';
a1.XLimitMethod = 'padded';
axis(a1,'square')
xlabel(a1, '$\Omega$');
ylabel(a1, '$x_p / L$');

f2 = figure('Units','centimeters','PaperUnits','centimeters',...
    'Position',[2 2 fsizex fsizey],'PaperPosition',[0 0 fsizex fsizey],...
    'PaperSize',[fsizex fsizey],'PaperPositionMode','manual');
a2 = axes(f2,'FontSize',afont,'LineWidth',aline,...
    'TickLabelInterpreter','latex','Box','on');
a2.XLabel.Interpreter = "latex";
a2.YLabel.Interpreter = "latex";
grid(a2, 'on'); hold(a2, 'on');
a2.YLimMode = 'auto';
a2.XLimitMethod = 'padded';
axis(a2,'square')
xlabel(a2, '$\Omega$');
ylabel(a2, '$y_p / L$');

f3 = figure('Units','centimeters','PaperUnits','centimeters',...
    'Position',[2 2 fsizex fsizey],'PaperPosition',[0 0 fsizex fsizey],...
    'PaperSize',[fsizex fsizey],'PaperPositionMode','manual');
a3 = axes(f3,'FontSize',afont,'LineWidth',aline,...
    'TickLabelInterpreter','latex','Box','on');
a3.XLabel.Interpreter = "latex";
a3.YLabel.Interpreter = "latex";
grid(a3, 'on'); hold(a3, 'on');
a3.YLimMode = 'auto';
a3.XLimitMethod = 'padded';
axis(a3,'square')
xlabel(a3, '$\Omega$');
ylabel(a3, '$\theta$ ($\pi$ rad)');

%% Data Load
folder = 'C:\Users\akb110\OneDrive - Imperial College London\PhD Files\Simulation Results\LieC++-PyCont_Thesis2024\Vertical_Cantilever\Sweep_Final\';
file = 'beam_sim_300.h5';
pose = h5read([folder file],'/dynamic_analysis/FEModel/POSE/MOTION').';
% vel = h5read(file,'/dynamic_analysis/FEModel/VELOCITY/MOTION').';
% time = h5read(file,'/dynamic_analysis/FEModel/time');

node_number = 23; % first node is 0
xconfig = node_number*4+3;
yconfig = node_number*4+4;
thetaconfig_c = node_number*4+1;
thetaconfig_s = node_number*4+2;
xdof = node_number*3+2;
ydof = node_number*3+1;
beam_length = 72e-3;
nat_freq = 6.33914;

%% Sort out single outlier
index = find(pose(yconfig,:)==0);
pose(xconfig,index) = 0.5*(pose(xconfig,index-1)+pose(xconfig,index+1));
pose(yconfig,index) = 0.5*(pose(yconfig,index-1)+pose(yconfig,index+1));

%% Find Theta
cos_theta_half = pose(thetaconfig_c, :);
sin_theta_half = pose(thetaconfig_s, :);
theta = 2 * atan2(sin_theta_half, cos_theta_half);

%% Sweep Function
f0 = 10;
f1 = 20;
sweep_rate = 1.0;
sampling_freq = 4000;
tend = abs(f1-f0)*(60/sweep_rate);
k = (f1-f0)/tend;
time = 0:1/sampling_freq:tend ;
finst = k*time + f0 ;
dt = time(2);

%% Plots
plot(a1,finst, pose(xconfig,:)/beam_length,'-','LineWidth',plotlinew);
plot(a2,finst, (pose(yconfig,:)-pose(8,:))/beam_length,'-','LineWidth',plotlinew);
plot(a3,finst, theta/pi,'-','LineWidth',plotlinew);

% set(a1, 'Children', flipud(get(gca, 'Children')));

%%
% lines = findall(gcf, 'Type', 'Line');
% numLines = length(lines);
% colors = turbo(numLines);
% for i = 1:numLines
%     lines(i).Color = colors(i, :);
% end