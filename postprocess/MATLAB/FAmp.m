clear all; close all;

%% GIVE AMPLITUDE PLOT.
% Need to have a h5 file with the time series of all the solutions on the
% NNM branch. run python timesim_branch to create such file

% file names should have _withtime 

%% Figure Properties
fsizex  = 25;
ratio = 3/4;
fsizey  = fsizex*ratio;
afont = 23;
aline = 2;
plotlinew = 2.5;

f = figure('Units','centimeters','PaperUnits','centimeters',...
    'Position',[2 2 fsizex fsizey],'PaperPosition',[0 0 fsizex fsizey],...
    'PaperSize',[fsizex fsizey],'PaperPositionMode','manual');
a = axes(f,'FontSize',afont,'LineWidth',aline,...
    'TickLabelInterpreter','latex','Box','on');
a.XLabel.Interpreter = "latex";
a.YLabel.Interpreter = "latex";
% legend(a,'Interpreter','latex','Location','northwest');
grid(a, 'on'); hold(a, 'on');

a.XLim = [0 inf];
a.YLimMode = 'auto';
a.XLimitMethod = 'tight';
ylabel(a, 'Frequency (Hz)');

xlabel(a, 'Z Position (m)');

%% DoF Index (+1 for MATLAB starting index)
node_number = 42; % first node is 0
xconfig = node_number*7+5;
yconfig = node_number*7+6;
zconfig = node_number*7+7;

normalise = 1.0;

%% Data Load
files1 = {};
folder = '//wsl$/ubuntu/home/akb110/Codes/PyCont/examples/beam_cpp/';
files1{end+1} = [folder 'right_sing_NNM1_withtime.h5'];
files1{end+1} = [folder 'right_sing_NNM1b_withtime.h5'];
files1{end+1} = [folder 'right_sing_NNM1c_withtime.h5'];
files1{end+1} = [folder 'right_sing_NNM1d_withtime.h5'];

%% Plot
for i=1:length(files1)
    pose = h5read(files1{i},'/Config_Time/POSE');
    pose = permute(pose,[3,2,1]);
    T = h5read(files1{i},'/T');
    f = 1./T;
    % find normalised max pose for each cont solution
    [maxpose, r] = max(abs(pose(zconfig,:,:))/normalise);
    maxpose = (reshape(maxpose, size(T)));
    r = reshape(r, size(T));
    
    color = [0.8500 0.3250 0.0980];
    if i==1
        plot(a,maxpose,f,'-','LineWidth',plotlinew,'Color',color);
    else
        plot(a,maxpose,f,'-','LineWidth',plotlinew,'Color',color,...
            'HandleVisibility','off');
    end
end

