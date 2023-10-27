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
% legend(a,'Interpreter','latex','Location','northwest');
grid(a, 'on'); hold(a, 'on');

a.XLim = [0.97 1.1];
a.YLimMode = 'auto';
a.XLimitMethod = 'tight';
ylabel(a, 'Normalised Max Y');
xlabel(a, 'F/Omega');
% a.Title.String = 'x/L = 1.0';

%% Geometries Undeformed Positions
normalise_amp = 1.0;
normalise_frq = 41.82280070074808;

%% DoF Index (+1 for MATLAB starting index)
node = 21;
dof = 4*node + 3+1;

%% Data Load
files1 = {}; files2 = {};
folder = '//wsl$/ubuntu/home/akb110/Codes/PyCont/examples/beam_2D/';
name1 = 'NNM';
files1{end+1} = [folder 'NNM1.h5'];
files1{end+1} = [folder 'NNM1b.h5'];
files1{end+1} = [folder 'NNM1c.h5'];
files1{end+1} = [folder 'NNM1d.h5'];
files1{end+1} = [folder 'NNM1e.h5'];

name2 = 'FRF';
% files2{end+1} = [folder 'FRF1_amp10000.h5'];
files2{end+1} = [folder 'FRF1_amp8000.h5'];
files2{end+1} = [folder 'FRF1_amp6000.h5'];
% files2{end+1} = [folder 'FRF1_amp4000.h5'];
% files2{end+1} = [folder 'FRF1_amp3000.h5'];
% files2{end+1} = [folder 'FRF1_amp2200.h5'];
% files2{end+1} = [folder 'FRF1_amp1800.h5'];
% files2{end+1} = [folder 'FRF1_amp1200.h5'];
% files2{end+1} = [folder 'FRF1_amp1000.h5'];
% files2{end+1} = [folder 'FRF1_amp880b.h5'];
% files2{end+1} = [folder 'FRF1_amp880.h5'];
% files2{end+1} = [folder 'FRF1_amp660.h5'];
% files2{end+1} = [folder 'FRF1_amp440.h5'];
% files2{end+1} = [folder 'FRF1_amp220.h5'];

%% Plot
for i=1:length(files1)
    pose = h5read(files1{i},'/Config_Time/POSE');
    pose = permute(pose,[3,2,1]);
    T = h5read(files1{i},'/T');
    f = 1./T;
    
    % find normalised max pose for each cont solution
    % following finds the max over T in one go, no need for for loop
    [amp, ~] = max(abs(pose(dof,:,:))/normalise_amp);
    amp = reshape(amp, size(T));
    
    color = 'k';
    if i==1
        plot(a,f/normalise_frq,amp,'-','LineWidth',plotlinew,'Color',color,'DisplayName',name1);
    else
        plot(a,f/normalise_frq,amp,'-','LineWidth',plotlinew,'Color',color,'DisplayName',name1,...
            'HandleVisibility','off');
    end
end

for i=1:length(files2)
    pose = h5read(files2{i},'/Config_Time/POSE');
    pose = permute(pose,[3,2,1]);
    T = h5read(files2{i},'/T');
    f = 1./T/normalise_frq;
    stability = h5read(files2{i},'/Bifurcation/Stability');
    stable_index = find(diff(stability))+1;
    stable_index = [1, stable_index', length(T)];
    
    % find normalised max pose for each cont solution
    [amp, ~] = max(abs(pose(dof,:,:))/normalise_amp);
    amp = reshape(amp, size(T));
    
    color = [0.8500 0.3250 0.0980];
    
    for j=1:length(stable_index)-1
        is_stable = stability(stable_index(j+1)-1);
        if is_stable
            linestyle = '-';
        else
            linestyle = '--';
        end
        plot(a,f(stable_index(j):stable_index(j+1)),amp(stable_index(j):stable_index(j+1)),...
            'LineStyle',linestyle,'LineWidth',plotlinew,'Color',color,'DisplayName',name2,'HandleVisibility','off');
    end
    
end
