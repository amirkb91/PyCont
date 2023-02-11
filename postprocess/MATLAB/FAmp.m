clear all; close all;

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
legend(a,'Interpreter','latex','Location','northwest');
grid(a, 'on'); hold(a, 'on');

a.XLim = [0 inf];
a.YLimMode = 'auto';
a.XLimitMethod = 'tight';
ylabel(a, 'Frequency (Hz)');

xlabel(a, 'Normalised Max Y Displacement');
a.Title.String = 'x/L = 1.0';

%% Geometries Undeformed Positions
beamthickness = 0.01;
% beam = linspace(0,1,31)';
% beam = [0;1;beam(2:end-1)];
% beam = [(0:30)',beam,zeros(31,1)];
% 
% R = 12.5;
% beta_max = 2. * asin(0.5 / R);
% arch = [R*cos(0.5*pi+beta_max*(0.5-linspace(0,1,31)))',...
%         R*sin(0.5*pi+beta_max*(0.5-linspace(0,1,31)))'];
% arch = [arch(1,:);arch(end,:);arch(2:end-1,:)];
% arch = [(0:30)',arch];

%% DoF Index (+1 for MATLAB starting index) 
% node10 = 0.3 / node12 = 0.36
node = 1;
% X
dof_VK   = 3*node + 1;  dof_SE23 = 4*node + 3;
% Y
% dof_VK   = 3*node + 2;  dof_SE23 = 4*node + 4;

%% Data Load
load('pose0.mat');
files1 = {}; files2 = {};
folder = '//wsl$/ubuntu/home/akb110/Codes/PyCont/examples/beam_2D/Results/cantilever/';
name1 = 'vonK';
files1{end+1} = [folder 'NNM2/cant_NNM2_VK_400.h5'];
name2 = 'SE(2)';
files2{end+1} = [folder 'NNM2/cant_NNM2_SE_400.h5'];
pose0_VK = pose0.VK_cclamped;
pose0_SE23 = pose0.SE23_cclamped;

%% Plot
for i=1:length(files1)
    pose = h5read(files1{i},'/Config/POSE_time');
    pose = permute(pose,[3,2,1]);
    T = h5read(files1{i},'/T');
    f = 1./T;
    % find normalised max pose for each cont solution
    [maxpose, r] = max(abs(pose(dof_VK,:,:)-pose0_VK(dof_VK))/beamthickness);
    maxpose = (reshape(maxpose, size(T)));
    r = reshape(r, size(T));
    
    color = 'k';
    if i==1
        plot(a,maxpose,f,'-','LineWidth',plotlinew,'Color',color,'DisplayName',name1);
    else
        plot(a,maxpose,f,'-','LineWidth',plotlinew,'Color',color,'DisplayName',name1,...
            'HandleVisibility','off');
    end
end

for i=1:length(files2)
    pose = h5read(files2{i},'/Config/POSE_time');
    pose = permute(pose,[3,2,1]);
    T = h5read(files2{i},'/T');
    f = 1./T;
    % find normalised max pose for each cont solution
    [maxpose, r] = max(abs(pose(dof_SE23,:,:)-pose0_SE23(dof_SE23))/beamthickness);
    maxpose = (reshape(maxpose, size(T)));
    r = reshape(r, size(T));
    
%     maxpose = (abs(pose(dof_SE23,1,:)-pose0_SE23(dof_SE23))/beamthickness);
%     maxpose = (reshape(maxpose, size(T)));
    
    color = [0 0.4470 0.7410];
    if i==1
        plot(a,maxpose,f,'-.','LineWidth',plotlinew,'Color',color,'DisplayName',name2);
    else
        plot(a,maxpose,f,'-.','LineWidth',plotlinew,'Color',color,'DisplayName',name2,...
            'HandleVisibility','off');
    end
end
