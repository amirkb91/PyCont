clear all; close all;

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
xlabel(a1, 'Frequency (Hz)');
ylabel(a1, 'X position (m)');

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
xlabel(a2, 'Frequency (Hz)');
ylabel(a2, 'Y position (m)');

%% Data Load
file = 'C:\Users\akb110\OneDrive - Imperial College London\PhD Files\Simulation Results\LieC++-PyCont_Thesis2024\Vertical_Cantilever\Sweep_New\beam_sim_sweep_amp10.h5';
pose = h5read(file,'/dynamic_analysis/FEModel/POSE/MOTION').';
vel = h5read(file,'/dynamic_analysis/FEModel/VELOCITY/MOTION').';
time = h5read(file,'/dynamic_analysis/FEModel/time');

node_number = 32; % first node is 0
xconfig = node_number*4+3;
yconfig = node_number*4+4;
xdof = node_number*3+1;
ydof = node_number*3+2;

%% Sweep Function
f0 = 11.5;
f1 = 15;
sweep_rate = 0.5;
sampling_freq = 1400;
    
tend = abs(f1-f0)*(60/sweep_rate);
k = (f1-f0)/tend ;
time2 = 0:1/sampling_freq:tend ;
finst = k*time2 + f0 ;
    
%% Plots
plot(a1,finst, pose(xconfig,:),'-','LineWidth',plotlinew);
plot(a2,finst, pose(yconfig,:),'-','LineWidth',plotlinew);

% set(a1, 'Children', flipud(get(gca, 'Children')));