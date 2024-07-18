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
% a1.ZLabel.Interpreter = "latex";
grid(a1, 'on'); hold(a1, 'on');
a1.YLimMode = 'auto';
a1.XLimitMethod = 'padded';
a1.YLimitMethod = 'padded';
% a1.ZLimitMethod = 'padded';
% a1.PlotBoxAspectRatio=[1,1,1];
% a1.XLim = [-0.5 1.5];
% a1.YLim = [-1 1];
xlabel(a1, '$X$');
ylabel(a1, '$Y$');
% zlabel(a1, '$Z$');
% a1.View=[40 13];


%% Data Load
folder = 'C:\Users\akb110\OneDrive - Imperial College London\PhD Files\Simulation Results\LieC++-PyCont_Thesis2024\Vertical_Cantilever\Sweep_New\wide\';
file = 'beam_sim_sweep_amp180.h5';
pose = h5read([folder file],'/dynamic_analysis/FEModel/POSE/MOTION').';

% file = '//wsl$/ubuntu/home/akb110/Codes/PyCont/examples/beam_cpp/FRF2_12000c_NSlargest_steadystate_timesim.h5';
% pose = h5read(file,'/dynamic_analysis/FEModel/POSE/MOTION').';
pose(1:4,:)=[];
pose(5:8,:)=[];

xconfig = (0:21)*4+3;
yconfig = (0:21)*4+4;
% zconfig = (0:21)*7+7;
x_beamtip = 21*4+4;
color = "#D95319";

%% Plot Pose
pose = pose(:,900229:900657);
n_snaps_grey = 0;
n_snaps = 30;
n_time = (length(pose));
beam_length = 72e-3;
index = 1:floor(n_time/n_snaps_grey):n_time;
for i=1:n_snaps_grey
    plot(a1,pose(xconfig,index(i))/beam_length,(pose(yconfig,index(i))-pose(4,index(i)))/beam_length,'o','LineWidth',plotlinew,'MarkerSize',0.5,'Color','#bcbcbc');
end

index = 1:floor(n_time/n_snaps):n_time;
abs_pose = abs(pose(x_beamtip, index)-1);
abs_pose_normalized = (abs_pose - min(abs_pose)) / (max(abs_pose) - min(abs_pose));
cmap = custom_viridis(n_snaps);
colors = interp1(linspace(0, 1, n_snaps), cmap, abs_pose_normalized);
for i=1:n_snaps
    plot(a1,pose(xconfig,index(i))/beam_length,(pose(yconfig,index(i))-pose(4,index(i)))/beam_length,'o-','LineWidth',plotlinew,'MarkerSize',3,'Color',colors(i,:));
end

axis equal

%%
function c = custom_viridis(m)

if nargin < 1, m = size(get(gcf,'colormap'),1); end

% Define the RGB values for the custom Viridis colormap
cm_data = [
    68, 1, 84
    71, 44, 122
    59, 81, 139
    44, 113, 142
    33, 144, 140
    39, 173, 129
    92, 200, 99
    170, 220, 50
    253, 231, 37
    ] / 255;

% Interpolate the colormap to the desired number of colors
x = linspace(1, length(cm_data), m);
c = interp1(1:length(cm_data), cm_data, x);

end