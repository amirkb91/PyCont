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

a.View=[-47,21];
axis(a,'square')
zlabel(a, '$x_p/L$');
xlabel(a, '$\Omega/\omega_1$');
ylabel(a, 'F (N)');
a.YTick = [];
a.ZDir = 'reverse';

%% Geometries Undeformed Positions
normalise_amp = 1.0;
normalise_frq = 41.7182203819309;

%% DoF Index (+1 for MATLAB starting index)
node_number = 21;

configpernode = 4;
xconfig = configpernode*node_number + 3;  % x direction
yconfig = configpernode*node_number + 4;  % y direction

%% Data Load
files = {};
folder = '//wsl$/ubuntu/home/akb110/Codes/PyCont/examples/beam_cpp/';

files{end+1} = [folder 'FRF1_200_withtime.h5'];
files{end+1} = [folder 'FRF1_400_withtime.h5'];
files{end+1} = [folder 'FRF1_600_withtime.h5'];
files{end+1} = [folder 'FRF1_2000_withtime.h5'];
files{end+1} = [folder 'FRF1_4000_withtime.h5'];
files{end+1} = [folder 'FRF1_8000_withtime.h5'];
% files{end+1} = [folder 'FRF1_12000_withtime.h5'];

force_amplitudes = [1 2 3 4 5 6];
% force_amplitudes = [200 400 600 2000 4000 8000];

%% Plot
for i=1:length(files)
    pose = h5read(files{i},'/Config_Time/POSE');
    pose = permute(pose,[3,2,1]);
    T = h5read(files{i},'/T');
    f = 1./T/normalise_frq;
    
    % find normalised max pose for each cont solution
    % -1 for x becuase that's the X position of tip undeformed
    [maxpose_X, ~] = max(abs(pose(xconfig,:,:)-1)/normalise_amp);
    maxpose_X = 1-reshape(maxpose_X, size(T));
    
    [maxpose_Y, ~] = max(abs(pose(yconfig,:,:))/normalise_amp);
    maxpose_Y = reshape(maxpose_Y, size(T));
    
    color = "#0072BD";
    
    plot3(a,f,force_amplitudes(i)*ones(length(f),1),maxpose_X,...
        '-','LineWidth',plotlinew,'Color',color);
 
end
