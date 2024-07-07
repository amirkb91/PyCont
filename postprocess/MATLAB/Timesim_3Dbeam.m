clear all; close all;

%% GIVE BEAM TIME SIMULATION FOR SINGLE POINT ON FEP
% Need to have a h5 file with the time series of the single solution we
% want to plot. Need to have a beam_sim from cpp for this to run. Can copy
% the file into the python folder and read from there

% file names should have _timesim

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
xlabel(a1, 'X position (m)');
ylabel(a1, 'X velocity (m/s)');

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
xlabel(a2, 'Y position (m)');
ylabel(a2, 'Y velocity (m/s)');

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
xlabel(a3, 'Z position (m)');
ylabel(a3, 'Z velocity (m/s)');

%% Data Load
file = '//wsl$/ubuntu/home/akb110/Codes/PyCont/examples/beam_cpp/right_sing_NNM1d_timesim.h5';
pose = h5read(file,'/dynamic_analysis/FEModel/POSE/MOTION').';
vel = h5read(file,'/dynamic_analysis/FEModel/VELOCITY/MOTION').';
vel_spatial = h5read(file,'/dynamic_analysis/FEModel/VELOCITY_SPATIAL/MOTION').';
time = h5read(file,'/dynamic_analysis/FEModel/time');

node_number = 42; % first node is 0
xconfig = node_number*7+5;
yconfig = node_number*7+6;
zconfig = node_number*7+7;
xdof = node_number*6+1;
ydof = node_number*6+2;
zdof = node_number*6+3;
thetaxdof = node_number*6+4;
thetaydof = node_number*6+5;
thetazdof = node_number*6+6;

%% Plots
plot(a1,pose(xconfig,:) ,vel(xdof,:),'-','LineWidth',plotlinew);
plot(a2,pose(yconfig,:) ,vel(ydof,:),'-','LineWidth',plotlinew);
plot(a3,pose(zconfig,:) ,vel(zdof,:),'-','LineWidth',plotlinew);
