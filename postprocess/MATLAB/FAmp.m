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
a.Title.String = 'x/L = 0.5';

%% Geometries Undeformed Positions
beamthickness = 0.01;
load('pose0.mat');
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
node = 16;
% X
dof_VK   = 3*node+0 + 1;  dof_SE23 = 4*node+2 + 1;
% Y
% dof_VK   = 3*node+1 + 1;  dof_SE23 = 4*node+3 + 1;

%% Data Load
files1 = {}; files2 = {};

% folder = 'cclamped/NNM1/';
% name1 = 'vonK';
% files1{end+1} = [folder 'VK/NNM1_VK.h5'];
% name2 = 'SE(2)';
% files2{end+1} = [folder 'SE23/NNM1_SE23_1.h5'];
% files2{end+1} = [folder 'SE23/NNM1_SE23_2.h5'];
% files2{end+1} = [folder 'SE23/NNM1_SE23_3.h5'];
% pose0_VK = pose0.VK_cclamped;
% pose0_SE23 = pose0.SE23_cclamped;
%===============================================
% folder = 'cclamped/NNM2/';
% name1 = 'vonK';
% files1{end+1} = [folder 'VK/NNM2_VK_1.h5'];
% files1{end+1} = [folder 'VK/NNM2_VK_2.h5'];
% files1{end+1} = [folder 'VK/NNM2_VK_3.h5'];
% files1{end+1} = [folder 'VK/NNM2_VK_4_trunc.h5'];
% name2 = 'SE(2)';
% files2{end+1} = [folder 'SE23/NNM2_SE23_1.h5'];
% files2{end+1} = [folder 'SE23/NNM2_SE23_2.h5'];
% files2{end+1} = [folder 'SE23/NNM2_SE23_3.h5'];
% files2{end+1} = [folder 'SE23/NNM2_SE23_4.h5'];
% files2{end+1} = [folder 'SE23/NNM2_SE23_5.h5'];
% pose0_VK = pose0.VK_cclamped;
% pose0_SE23 = pose0.SE23_cclamped;
%===============================================
% folder = 'arch/NNM1/';
% name1 = 'vonK';
% files1{end+1} = [folder 'VK/NNM1_VK_1.h5'];
% files1{end+1} = [folder 'VK/NNM1_VK_2.h5'];
% files1{end+1} = [folder 'VK/NNM1_VK_3.h5'];
% files1{end+1} = [folder 'VK/NNM1_VK_4_trunc.h5'];
% name2 = 'SE(2)';
% files2{end+1} = [folder 'SE23/NNM1_SE23_1.h5'];
% files2{end+1} = [folder 'SE23/NNM1_SE23_2.h5'];
% files2{end+1} = [folder 'SE23/NNM1_SE23_3.h5'];
% files2{end+1} = [folder 'SE23/NNM1_SE23_4.h5'];
% pose0_VK = pose0.VK_arch;
% pose0_SE23 = pose0.SE23_arch;
%===============================================
folder = 'arch/NNM2/';
name1 = 'vonK';
files1{end+1} = [folder 'VK/NNM2_VK_1_trunc.h5'];
files1{end+1} = [folder 'VK/NNM2_VK_2.h5'];
files1{end+1} = [folder 'VK/NNM2_VK_3.h5'];
files1{end+1} = [folder 'VK/NNM2_VK_4.h5'];
name2 = 'SE(2)';
files2{end+1} = [folder 'SE23/NNM2_SE23_1.h5'];
files2{end+1} = [folder 'SE23/NNM2_SE23_2.h5'];
files2{end+1} = [folder 'SE23/NNM2_SE23_3.h5'];
files2{end+1} = [folder 'SE23/NNM2_SE23_4.h5'];
files2{end+1} = [folder 'SE23/NNM2_SE23_5.h5'];
files2{end+1} = [folder 'SE23/NNM2_SE23_6.h5'];
files2{end+1} = [folder 'SE23/NNM2_SE23_7.h5'];
files2{end+1} = [folder 'SE23/NNM2_SE23_8.h5'];
files2{end+1} = [folder 'SE23/NNM2_SE23_9.h5'];
files2{end+1} = [folder 'SE23/NNM2_SE23_10.h5'];
pose0_VK = pose0.VK_arch;
pose0_SE23 = pose0.SE23_arch;

%% Plot
for i=1:length(files1)
    pose = h5read(files1{i},'/POSE');
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
    pose = h5read(files2{i},'/POSE');
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
