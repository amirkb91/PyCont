clear all; close all;

%% Figure Properties
fsizex  = 25;
ratio = 3/4;
fsizey  = fsizex*ratio;
afont = 23;
aline = 2;
plotlinew = 1.5;

f = figure('Units','centimeters','PaperUnits','centimeters',...
    'Position',[2 2 fsizex fsizey],'PaperPosition',[0 0 fsizex fsizey],...
    'PaperSize',[fsizex fsizey],'PaperPositionMode','manual');
a = axes(f,'FontSize',afont,'LineWidth',aline,...
    'TickLabelInterpreter','latex','Box','on');
a.XLabel.Interpreter = "latex";
a.YLabel.Interpreter = "latex";
a.ZLabel.Interpreter = "latex";
hold(a, 'on');
a.View=[-60,20];
axis(a,'square')

a.XLim = [0.97 1.1];
a.YScale= 'log';
% a.YLimMode = 'auto';
% a.XLimitMethod = 'tight';
xlabel(a, 'F/Omega');
zlabel(a, 'Normalised Max Y Amp');

%% Geometries Undeformed Positions
normalise_amp = 1.0;
normalise_frq = 41.82280070074808;

%% DoF Index (+1 for MATLAB starting index)
node = 21;
dof = 4*node + 4;

%% Data Load
files = {};
folder = '//wsl$/ubuntu/home/akb110/Codes/PyCont/examples/beam_2D/';
files{end+1} = [folder 'FRF1_amp10000.h5'];
files{end+1} = [folder 'FRF1_amp8000.h5'];
files{end+1} = [folder 'FRF1_amp6000.h5'];
files{end+1} = [folder 'FRF1_amp4000.h5'];
files{end+1} = [folder 'FRF1_amp3000.h5'];
files{end+1} = [folder 'FRF1_amp2200.h5'];
files{end+1} = [folder 'FRF1_amp1800.h5'];
files{end+1} = [folder 'FRF1_amp1200.h5'];
files{end+1} = [folder 'FRF1_amp1000.h5'];
files{end+1} = [folder 'FRF1_amp880b.h5'];
files{end+1} = [folder 'FRF1_amp880.h5'];
files{end+1} = [folder 'FRF1_amp660.h5'];
files{end+1} = [folder 'FRF1_amp440.h5'];
files{end+1} = [folder 'FRF1_amp220.h5'];

%% Plot
for i=1:length(files)
    pose = h5read(files{i},'/Config_Time/POSE');
    pose = permute(pose,[3,2,1]);
    E = h5read(files{i},'/Energy');
    T = h5read(files{i},'/T');
    f = 1./T;
    
    force = jsondecode(h5read(files{i},'/Parameters'));
    force = force.forcing.amplitude;
    force = repmat(force,[length(T),1]);
    
    % find normalised max pose for each cont solution
    [amp, ~] = max(abs(pose(dof,:,:))/normalise_amp);
    amp = reshape(amp, size(T));
    
%     plot3(a,f/normalise_frq,E,amp,'-','LineWidth',plotlinew,'Color','k'); ylabel(a, 'Energy (J)');
    plot3(a,f/normalise_frq,force,amp,'-','LineWidth',plotlinew,'Color','k'); ylabel(a, 'Force Amp (N)');
end

