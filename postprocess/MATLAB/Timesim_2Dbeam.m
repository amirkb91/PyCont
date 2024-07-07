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

%% Data Load
file = '//wsl$/ubuntu/home/akb110/Codes/PyCont/examples/beam_cpp/FRF1_2000_timesim.h5';
pose = h5read(file,'/dynamic_analysis/FEModel/POSE/MOTION').';
vel = h5read(file,'/dynamic_analysis/FEModel/VELOCITY/MOTION').';
% vel_spatial = h5read(file,'/dynamic_analysis/FEModel/VELOCITY_SPATIAL/MOTION').';
time = h5read(file,'/dynamic_analysis/FEModel/time');

node_number = 21; % first node is 0
xconfig = node_number*4+3;
yconfig = node_number*4+4;
xdof = node_number*3+1;
ydof = node_number*3+2;

%% Plots
plot(a1,pose(xconfig,:) ,vel(xdof,:),'-','LineWidth',plotlinew);
plot(a2,pose(yconfig,:) ,vel(ydof,:),'-','LineWidth',plotlinew);
