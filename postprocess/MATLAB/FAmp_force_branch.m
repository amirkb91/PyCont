clear all; close all;

%% GIVE AMPLITUDE PLOT FRF.
% Need to have a h5 file with the time series of all the solutions on the
% NNM branch. run python timesim_branch to create such file

% file names should have _withtime
% pick if you want X Y or theta by commenting out

%% Figure Properties
fsizex  = 25;
ratio = 3/4;
fsizey  = fsizex*ratio;
afont = 23;
aline = 1;
plotlinew = 2.0;

f = figure('Units','centimeters','PaperUnits','centimeters',...
    'Position',[2 2 fsizex fsizey],'PaperPosition',[0 0 fsizex fsizey],...
    'PaperSize',[fsizex fsizey],'PaperPositionMode','manual');
a = axes(f,'FontSize',afont,'LineWidth',aline,...
    'TickLabelInterpreter','latex','Box','on','BoxStyle','full');
a.XLabel.Interpreter = "latex";
a.YLabel.Interpreter = "latex";
a.ZLabel.Interpreter = "latex";
grid(a, 'on'); hold(a, 'on');
% a.BoxStyle="full";

a.XLimitMethod = 'padded';
a.YLimitMethod = 'padded';
a.ZLimitMethod = 'padded';

a.View=[30,40];
axis(a,'square')
zlabel(a, '$y_p$');
xlabel(a, '$\Omega/\omega_1$');
ylabel(a, '$z_p$');
% a.YTick = [];

%% Geometries Undeformed Positions
normalise_amp = 0.25;
normalise_frq = 5.789e+0;

%% DoF Index (+1 for MATLAB starting index)
node_number = 21;

configpernode = 7;
xconfig = configpernode*node_number + 5;  % x direction
yconfig = configpernode*node_number + 6;  % y direction
zconfig = configpernode*node_number + 7;  % y direction

%% Data Load
files = {};
folder = 'C:\Users\akb110\OneDrive - Imperial College London\PhD Files\Simulation Results\LieC++-PyCont_Thesis2024\3D_Cantilever/';

files{end+1} = [folder 'FRF1_0_005_withtime.h5'];
files{end+1} = [folder 'FRF1_0_005_branch_withtime.h5'];


%% Plot
for i=1:length(files)
    pose = h5read(files{i},'/Config_Time/POSE');
    pose = permute(pose,[3,2,1]);
    T = h5read(files{i},'/T');
    f = 1./T/normalise_frq;
    
    if i == 2
        pose = pose(1:end,1:end,1:196);
        T = T(1:196);
        f = f(1:196);
    end
    stability = h5read(files{i},'/Bifurcation/Stability');
    stable_index = find(diff(stability))+1;
    stable_index = [1, stable_index', length(T)];
    
    [maxpose_Y, ~] = max(abs(pose(yconfig,:,:))/normalise_amp);
    maxpose_Y = reshape(maxpose_Y, size(T));
    
    [maxpose_Z, ~] = max(abs(pose(zconfig,:,:))/normalise_amp);
    maxpose_Z = reshape(maxpose_Z, size(T));
    
    if i == 1
        color = "#0072BD";
    else
        color = "#D95319";
    end
    
    for j=1:length(stable_index)-1
        is_stable = stability(stable_index(j+1)-1);
        if is_stable
            linestyle = '-';
        else
            linestyle = ':';
        end
        plot3(a,f(stable_index(j):stable_index(j+1)),maxpose_Z(stable_index(j):stable_index(j+1)),maxpose_Y(stable_index(j):stable_index(j+1)),...
            'LineWidth',plotlinew,'LineStyle',linestyle,'Color',color);       
    end
    
end
