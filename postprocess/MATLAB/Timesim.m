clear all; close all;

%% Figure Properties
fsizex  = 25;
ratio = 3/4;
fsizey  = fsizex*ratio;
afont = 23;
aline = 2;
plotlinew = 2.5;

f1 = figure('Units','centimeters','PaperUnits','centimeters',...
    'Position',[2 2 fsizex fsizey],'PaperPosition',[0 0 fsizex fsizey],...
    'PaperSize',[fsizex fsizey],'PaperPositionMode','manual');
a1 = axes(f1,'FontSize',afont,'LineWidth',aline,...
    'TickLabelInterpreter','latex','Box','on');
a1.XLabel.Interpreter = "latex";
a1.YLabel.Interpreter = "latex";
legend(a1,'Interpreter','latex','Location','northeast');
grid(a1, 'on'); hold(a1, 'on');
a1.YLimMode = 'auto';
a1.XLimitMethod = 'tight';
xlabel(a1, 'Time (s)');
ylabel(a1, 'Normalised Displacement');

f2 = figure('Units','centimeters','PaperUnits','centimeters',...
    'Position',[2 2 fsizex fsizey],'PaperPosition',[0 0 fsizex fsizey],...
    'PaperSize',[fsizex fsizey],'PaperPositionMode','manual');
a2 = axes(f2,'FontSize',afont,'LineWidth',aline,...
    'TickLabelInterpreter','latex','Box','on');
a2.XLabel.Interpreter = "latex";
a2.YLabel.Interpreter = "latex";
legend(a2,'Interpreter','latex','Location','northeast');
grid(a2, 'on'); hold(a2, 'on');
xlabel(a2, 'Frequency Ratio');
ylabel(a2, 'Normalised Displacement');
a2.XLim = [0 12];
a2.XTick = (0:1:12);

f3 = figure('Units','centimeters','PaperUnits','centimeters',...
    'Position',[2 2 fsizex fsizey],'PaperPosition',[0 0 fsizex fsizey],...
    'PaperSize',[fsizex fsizey],'PaperPositionMode','manual');
a3 = axes(f3,'FontSize',afont,'LineWidth',aline,...
    'TickLabelInterpreter','latex','Box','on');
a3.XLabel.Interpreter = "latex";
a3.YLabel.Interpreter = "latex";
legend(a3,'Interpreter','latex','Location','northeast');
grid(a3, 'on'); hold(a3, 'on');
a3.YLimMode = 'auto';
a3.XLimitMethod = 'padded';
xlabel(a3, 'Normalised Displacement');
ylabel(a3, 'Velocity (m/s)');

%% Geometries Undeformed Positions
beamthickness = 0.01;
load('pose0.mat');

%% DoF Index (+1 for MATLAB starting index) 
node = 12;
% X
% dof_VK   = 3*node+0 + 1;  dof_SE23 = 4*node+2 + 1;
% Y
dof_VK   = 3*node+1 + 1;  dof_SE23 = 4*node+3 + 1;

%% Fundamental Frequency FFT
% fundamental = 220; %cclamped NNM2
% fundamental = 74.8; %arch NNM1
% fundamental = 170; %arch NNM2 a
% fundamental = 172; %arch NNM2 b
fundamental = 169; %arch NNM2 c
% fundamental = 176.4; %arch NNM2 d
% fundamental = 189; %arch NNM2 e

%% Data Load
%===============================================
% folder = 'cclamped/NNM2/';
% name1 = 'vonK';
% file1 = [folder 'VK/TIME_NNM2_VK_a.h5'];
% name2 = 'SE(2)';
% file2 = [folder 'SE23/TIME_NNM2_SE23_a.h5'];
% pose0_VK = pose0.VK_cclamped;
% pose0_SE23 = pose0.SE23_cclamped;
%===============================================
% folder = 'arch/NNM1/';
% name1 = 'vonK';
% % file1 = [folder 'VK/TIME_NNM1_VK_tip.h5'];
% file1 = [folder 'VK/TIME_NNM1_VK_mid.h5'];
% name2 = 'SE(2)';
% % file2 = [folder 'SE23/TIME_NNM1_SE23_tip.h5'];
% file2 = [folder 'SE23/TIME_NNM1_SE23_mid.h5'];
% pose0_VK = pose0.VK_arch;
% pose0_SE23 = pose0.SE23_arch;
%===============================================
folder = 'arch/NNM2/';
name1 = 'vonK';
file1 = [folder 'VK/TIME_NNM2_VK_c.h5'];
name2 = 'SE(2)';
file2 = [folder 'SE23/TIME_NNM2_SE23_c.h5'];
pose0_VK = pose0.VK_arch;
pose0_SE23 = pose0.SE23_arch;

%% Plot Pose
pose = h5read(file1,'/Config/POSE').';
p1 = (pose(dof_VK,:)-pose0_VK(dof_VK))/beamthickness;
t1 = h5read(file1,'/time');
color = 'k';
plot(a1,t1,p1,'-','LineWidth',plotlinew,'Color',color,'DisplayName',name1);

pose = h5read(file2,'/Config/POSE').';
p2 = (pose(dof_SE23,:)-pose0_SE23(dof_SE23))/beamthickness;
t2 = h5read(file2,'/time');
color = 'r';
plot(a1,t2,p2,'-.','LineWidth',plotlinew,'Color',color,'DisplayName',name2);

%% Plot Config Space
pose = h5read(file1,'/Config/POSE').';
vel = h5read(file1,'/Config/VELOCITY').';
p1 = (pose(dof_VK,:)-pose0_VK(dof_VK))/beamthickness;
v1 = vel(dof_VK,:);
color = 'k';
plot(a3,p1,v1,'-','LineWidth',plotlinew,'Color',color,'DisplayName',name1);

pose = h5read(file2,'/Config/POSE').';
vel = h5read(file2,'/Config/VELOCITY').';
p2 = (pose(dof_SE23,:)-pose0_SE23(dof_SE23))/beamthickness;
v2 = vel(dof_VK,:);  % dof VK is correct as velocity has same size as VK
color = 'r';
plot(a3,p2,v2,'-.','LineWidth',plotlinew,'Color',color,'DisplayName',name2);

%% FFT
Fs = 1/(t1(2)-t1(1));
L = length(t1)-1;
L2 = floor(L/2);
ZF = fft(p1);
ZM = abs(ZF)/L;
ZP = angle(ZF);
P = ZM(1:L2+1);
P(2:end-1) = 2*P(2:end-1);
P2 = ZP(1:L2+1);
f = Fs*(0:(L2))/L;
color = 'k';
plot(a2,f/fundamental,P,'-','LineWidth',plotlinew,'Color',color,'DisplayName',name1);

Fs = 1/(t2(2)-t2(1));
L = length(t2)-1;
L2 = floor(L/2);
ZF = fft(p2);
ZM = abs(ZF)/L;
ZP = angle(ZF);
P = ZM(1:L2+1);
P(2:end-1) = 2*P(2:end-1);
P2 = ZP(1:L2+1);
f = Fs*(0:(L2))/L;
color = [0.9290 0.6940 0.1250];
plot(a2,f/fundamental,P,'-.','LineWidth',plotlinew,'Color',color,'DisplayName',name2);
