clear all; close all;

%% GIVE BEAM SHAPE PLOT
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
plotlinew = 2;

f1 = figure('Units','centimeters','PaperUnits','centimeters',...
    'Position',[2 2 fsizex fsizey],'PaperPosition',[0 0 fsizex fsizey],...
    'PaperSize',[fsizex fsizey],'PaperPositionMode','manual');
a1 = axes(f1,'FontSize',afont,'LineWidth',aline,...
    'TickLabelInterpreter','latex','Box','on');
a1.XLabel.Interpreter = "latex";
a1.YLabel.Interpreter = "latex";
a1.ZLabel.Interpreter = "latex";
grid(a1, 'on'); hold(a1, 'on');
a1.YLimMode = 'auto';
a1.XLimitMethod = 'padded';
a1.YLimitMethod = 'padded';
a1.ZLimitMethod = 'padded';
a1.PlotBoxAspectRatio=[1,1,1];
% a1.XLim = [-0.5 1.5];
% a1.YLim = [-1 1];
xlabel(a1, '$X$');
ylabel(a1, '$Y$');
zlabel(a1, '$Z$');
a1.View=[40 13];


%% Data Load
file = '//wsl$/ubuntu/home/akb110/Codes/PyCont/examples/beam_cpp/right_sing_NNM1d_timesim.h5';
pose = h5read(file,'/dynamic_analysis/FEModel/POSE/MOTION').';

xconfig = (0:42)*7+5;
yconfig = (0:42)*7+6;
zconfig = (0:42)*7+7;
color = "#D95319";

%% Plot Pose
n_snaps_grey = 50;
n_snaps = 15;
n_time = (length(pose)-1)/2;

index = ceil(linspace(0,n_time,n_snaps_grey));
for i=1:n_snaps_grey
    plot3(a1,pose(xconfig,index(i)+1),pose(yconfig,index(i)+1),pose(zconfig,index(i)+1),'o','LineWidth',plotlinew,'MarkerSize',0.5,'Color','#bcbcbc');
end

index = ceil(linspace(0,n_time,n_snaps));
for i=1:n_snaps
    plot3(a1,pose(xconfig,index(i)+1),pose(yconfig,index(i)+1),pose(zconfig,index(i)+1),'o-','LineWidth',plotlinew,'MarkerSize',3,'Color',color);
end

axis equal
