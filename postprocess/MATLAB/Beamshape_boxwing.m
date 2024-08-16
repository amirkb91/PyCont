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
afont = 15;
aline = 2;
plotlinew = 2;

f1 = figure('Units','centimeters','PaperUnits','centimeters',...
    'Position',[2 2 fsizex fsizey],'PaperPosition',[0 0 fsizex fsizey],...
    'PaperSize',[fsizex fsizey],'PaperPositionMode','manual');
a1 = axes(f1,'FontSize',afont,'LineWidth',aline,...
    'TickLabelInterpreter','latex','Box','off');
a1.XLabel.Interpreter = "latex";
a1.YLabel.Interpreter = "latex";
a1.ZLabel.Interpreter = "latex";
grid(a1, 'off'); hold(a1, 'on');
% a1.YLimMode = 'auto';
a1.XLimitMethod = 'tight';
a1.YLimitMethod = 'tight';
a1.ZLimitMethod = 'tight';
a1.PlotBoxAspectRatio=[1,1,1];
% xlabel(a1, '$X$');
ylabel(a1, '$Y$');
zlabel(a1, '$Z$');
a1.View=[-73 10];
a1.XTick = [];
% a1.YTick = [];
% a1.ZTick = [];
% a1.ZLim = [-1 1];
% a1.YLim = [-1 1];
% a1.XLim = [0 1];


%% Data Load
folder = 'C:\Users\akb110\OneDrive - Imperial College London\PhD Files\Simulation Results\LieC++-PyCont_Thesis2024\Boxwing\Rigid_Body\';
% file_static = {};
% file_static{end+1} = 'beam_ss_01.h5';
% file_static{end+1} = 'beam_ss_04.h5';
% file_static{end+1} = 'beam_ss_09.h5';

file_dynamic = 'FRC_lat_amp5b_peak_a_timesim.h5';

model = 'boxwing_XZsupport.h5';
nodes = {};
nodes{1} = h5read([folder model],'/Components/Component_0/Nodes/NodalFrame<3, (Fields::type)0, FrameParametersT_ER>')';
nodes{2} = h5read([folder model],'/Components/Component_1/Nodes/NodalFrame<3, (Fields::type)0, FrameParametersT_ER>')';
nodes{3} = h5read([folder model],'/Components/Component_2/Nodes/NodalFrame<3, (Fields::type)0, FrameParametersT_ER>')';
nodes{4} = h5read([folder model],'/Components/Component_3/Nodes/NodalFrame<3, (Fields::type)0, FrameParametersT_ER>')';
nodes{5} = h5read([folder model],'/Components/Component_4/Nodes/NodalFrame<3, (Fields::type)0, FrameParametersT_ER>')';
nodes{6} = h5read([folder model],'/Components/Component_5/Nodes/NodalFrame<3, (Fields::type)0, FrameParametersT_ER>')';
nodes{7} = h5read([folder model],'/Components/Component_6/Nodes/NodalFrame<3, (Fields::type)0, FrameParametersT_ER>')';
nodes{8} = h5read([folder model],'/Components/Component_7/Nodes/NodalFrame<3, (Fields::type)0, FrameParametersT_ER>')';
nodes{9} = h5read([folder model],'/Components/Component_8/Nodes/NodalFrame<3, (Fields::type)0, FrameParametersT_ER>')';
nodes{10} = h5read([folder model],'/Components/Component_9/Nodes/NodalFrame<3, (Fields::type)0, FrameParametersT_ER>')';
nodes{11} = h5read([folder model],'/Components/Component_10/Nodes/NodalFrame<3, (Fields::type)0, FrameParametersT_ER>')';
nodes{12} = h5read([folder model],'/Components/Component_11/Nodes/NodalFrame<3, (Fields::type)0, FrameParametersT_ER>')';
nodes_config = h5read([folder model],'/FEModel/Nodes_config/NodalFrame<3, (Fields::type)0, FrameParametersT_ER>')';

%% Plot Pose Static
% for j=1:length(file_static)
%     for i=1:12
%         pose = h5read([folder file_static{j}],'/dynamic_analysis/FEModel/POSE/MOTION').';
%         xconfig = nodes_config(6,nodes{i}+1)+1;
%         yconfig = nodes_config(7,nodes{i}+1)+1;
%         zconfig = nodes_config(8,nodes{i}+1)+1;
%         plot3(a1,pose(xconfig,end),pose(yconfig,end),pose(zconfig,end),'o-','LineWidth',plotlinew,'MarkerSize',3,'Color','k');
%     end
% end
% axis equal

%% Plot Pose Dynamic Snaps
n_snaps = 6;

for i=1:12
    pose = h5read([folder file_dynamic],'/dynamic_analysis/FEModel/POSE/MOTION').';
    n_time = size(pose,2);
    xconfig = nodes_config(6,nodes{i}+1)+1;
    yconfig = nodes_config(7,nodes{i}+1)+1;
    zconfig = nodes_config(8,nodes{i}+1)+1;
    index = 1:floor(n_time/n_snaps):n_time;
    cmap = custom_viridis(n_snaps);
%     if i == 1
%         offset = pose(zconfig(1),index(1))*0;
%     end
    for j=1:n_snaps
        plot3(a1,pose(xconfig,index(j))+j-1,pose(yconfig,index(j)),pose(zconfig,index(j)),'-','LineWidth',plotlinew,'MarkerSize',3,'Color',cmap(j,:));
    end
end
axis equal
a1.View=[-86 8];
a1.DataAspectRatio=[1 1 1];
% a1.XAxis.Visible = 'off';
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