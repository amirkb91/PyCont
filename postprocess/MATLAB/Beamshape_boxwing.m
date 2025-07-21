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
% a1.XLimitMethod = 'tight';
% a1.YLimitMethod = 'tight';
% a1.ZLimitMethod = 'tight';
a1.PlotBoxAspectRatio=[1,1,1];
xlabel(a1, '$X$');
ylabel(a1, '$Y$');
zlabel(a1, '$Z$');
% a1.View=[-73 10];
% a1.XTick = [];
% a1.YTick = [];
% a1.ZTick = [];
% a1.ZLim = [-1 1];
% a1.YLim = [-1 1];
% a1.XLim = [0 1];


%% Data Load
folder = '/home/akb110/Codes/PyCont/examples/beam_cpp/';
% file_static = {};
% file_static{end+1} = 'beam_ss_01.h5';
% file_static{end+1} = 'beam_ss_04.h5';
% file_static{end+1} = 'beam_ss_09.h5';

file_dynamic = 'beam_sim.h5';

model = 'boxwing.h5';
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
% % n_snaps = 6;
% % 
% % for i=1:12
% %     pose = h5read([folder file_dynamic],'/dynamic_analysis/FEModel/POSE/MOTION').';
% %     n_time = size(pose,2);
% %     xconfig = nodes_config(6,nodes{i}+1)+1;
% %     yconfig = nodes_config(7,nodes{i}+1)+1;
% %     zconfig = nodes_config(8,nodes{i}+1)+1;
% %     index = 1:floor(n_time/n_snaps):n_time;
% %     cmap = custom_viridis(n_snaps);
% % %     if i == 1
% % %         offset = pose(zconfig(1),index(1))*0;
% % %     end
% %     for j=1:n_snaps
% %         plot3(a1,pose(xconfig,index(j))+j-1,pose(yconfig,index(j)),pose(zconfig,index(j)),'-','LineWidth',plotlinew,'MarkerSize',3,'Color',cmap(j,:));
% %     end
% % end
% % axis equal
% % a1.View=[-86 8];
% % a1.DataAspectRatio=[1 1 1];
% % % a1.XAxis.Visible = 'off';


%% Animate Pose Dynamic and Create GIF
% Read pose data once for the full time series
pose = h5read([folder file_dynamic],'/dynamic_analysis/FEModel/POSE/MOTION').';
n_time = size(pose,2);

% Pre-compute node indices for each beam component
xconfigs = cell(12,1);
yconfigs = cell(12,1);
zconfigs = cell(12,1);
for i = 1:12
    xconfigs{i} = nodes_config(6, nodes{i}+1)+1;
    yconfigs{i} = nodes_config(7, nodes{i}+1)+1;
    zconfigs{i} = nodes_config(8, nodes{i}+1)+1;
end

% Compute overall min and max for X, Y, and Z across all components and time steps
allXmin = inf; allXmax = -inf;
allYmin = inf; allYmax = -inf;
allZmin = inf; allZmax = -inf;
for i = 1:12
    curX = pose(xconfigs{i}, :);
    curY = pose(yconfigs{i}, :);
    curZ = pose(zconfigs{i}, :);
    allXmin = min(allXmin, min(curX(:)));
    allXmax = max(allXmax, max(curX(:)));
    allYmin = min(allYmin, min(curY(:)));
    allYmax = max(allYmax, max(curY(:)));
    allZmin = min(allZmin, min(curZ(:)));
    allZmax = max(allZmax, max(curZ(:)));
end

% Define padding fraction and compute padded limits
padding_fraction = 0.05;
X_range = allXmax - allXmin;
Y_range = allYmax - allYmin;
Z_range = allZmax - allZmin;
paddedXmin = allXmin - padding_fraction * X_range;
paddedXmax = allXmax + padding_fraction * X_range;
paddedYmin = allYmin - padding_fraction * Y_range;
paddedYmax = allYmax + padding_fraction * Y_range;
paddedZmin = allZmin - padding_fraction * Z_range;
paddedZmax = allZmax + padding_fraction * Z_range;

% Clear the axes and set fixed axes properties with computed limits
cla(a1);
axis equal;
a1.XLim = [paddedXmin, paddedXmax];
a1.YLim = [paddedYmin, paddedYmax];
a1.ZLim = [paddedZmin, paddedZmax];
a1.DataAspectRatio = [1 1 1];
a1.View = [-40, 33];

% Initialize plot objects for each component (one per beam) in black
h_anim = gobjects(12,1);
for i = 1:12
    h_anim(i) = plot3(a1, pose(xconfigs{i}, 1), pose(yconfigs{i}, 1), pose(zconfigs{i}, 1), ...
        'LineWidth', plotlinew, 'Color', 'k');
end
drawnow;

% Option for frame skipping for quicker animation
m_skip = 5;  % Recommended skip value (adjust as needed for smoothness)
frame_indices = 1:m_skip:n_time;

% Set up GIF file parameters
gif_filename = 'beam_animation.gif';
delay_time = 0.05; % Delay time in seconds per frame

% Loop over each selected time step to update the pose and capture frames
for t = frame_indices
    for i = 1:12
        % Update the plotted beam for component i
        set(h_anim(i), 'XData', pose(xconfigs{i}, t), ...
                       'YData', pose(yconfigs{i}, t), ...
                       'ZData', pose(zconfigs{i}, t));
    end
    drawnow;
    
    % Capture the frame
    frame = getframe(f1);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    % Write to the GIF file
    if t == frame_indices(1)
        imwrite(imind, cm, gif_filename, 'gif', 'LoopCount', Inf, 'DelayTime', delay_time);
    else
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
    end
end



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