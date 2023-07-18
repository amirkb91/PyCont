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
% legend(a1,'Interpreter','latex','Location','northeast');
grid(a1, 'on'); hold(a1, 'on');
a1.YLimMode = 'auto';
a1.XLimitMethod = 'tight';
xlabel(a1, 'Time (s)');
ylabel(a1, 'Normalised Amplitude');

f2 = figure('Units','centimeters','PaperUnits','centimeters',...
    'Position',[2 2 fsizex fsizey],'PaperPosition',[0 0 fsizex fsizey],...
    'PaperSize',[fsizex fsizey],'PaperPositionMode','manual');
a2 = axes(f2,'FontSize',afont,'LineWidth',aline,...
    'TickLabelInterpreter','latex','Box','on');
a2.XLabel.Interpreter = "latex";
a2.YLabel.Interpreter = "latex";
% legend(a2,'Interpreter','latex','Location','northeast');
grid(a2, 'on'); hold(a2, 'on');
xlabel(a2, 'Frequency Ratio');
ylabel(a2, 'Normalised Amplitude');
a2.XLim = [0 12];
a2.XTick = (0:1:12);

f3 = figure('Units','centimeters','PaperUnits','centimeters',...
    'Position',[2 2 fsizex fsizey],'PaperPosition',[0 0 fsizex fsizey],...
    'PaperSize',[fsizex fsizey],'PaperPositionMode','manual');
a3 = axes(f3,'FontSize',afont,'LineWidth',aline,...
    'TickLabelInterpreter','latex','Box','on');
a3.XLabel.Interpreter = "latex";
a3.YLabel.Interpreter = "latex";
% legend(a3,'Interpreter','latex','Location','northeast');
grid(a3, 'on'); hold(a3, 'on');
a3.YLimMode = 'auto';
a3.XLimitMethod = 'padded';
xlabel(a3, 'Normalised Amplitude');
ylabel(a3, 'Velocity (m/s)');

f4 = figure('Units','centimeters','PaperUnits','centimeters',...
    'Position',[2 2 fsizex fsizey],'PaperPosition',[0 0 fsizex fsizey],...
    'PaperSize',[fsizex fsizey],'PaperPositionMode','manual');
a4 = axes(f4,'FontSize',afont,'LineWidth',aline,...
    'TickLabelInterpreter','latex','Box','on');
a4.XLabel.Interpreter = "latex";
a4.YLabel.Interpreter = "latex";
% legend(a4,'Interpreter','latex','Location','northeast');
grid(a4, 'on'); hold(a4, 'on');
xlabel(a4, 'Frequency Ratio');
ylabel(a4, 'Normalised Amplitude');
a4.XLim = [-inf 25];
a4.XTick = (0:1:25);

%% Geometries Undeformed Positions
beamthickness = 0.01;
load('pose0_new-node-order.mat');

%% DoF Index 
% node10 = 0.3 / node12 = 0.36 / node16 = 0.5 / node30 = 1.0
node = 30;
dof_VK   = 3*node+1;  dof_SE = 4*node+3;  % X
% dof_VK   = 3*node+2;  dof_SE = 4*node+4;  % Y

%% Data Load
folder = '//wsl$/ubuntu/home/akb110/Codes/PyCont/examples/beam_2D/Results/cantilever/NNM2/';
name1 = 'vonK';
file1 = [folder 'VK/Time_BB.h5'];
fileeig1 = [folder 'VK/beam_eig.h5'];
pose0_VK = pose0.VK_cclamped;
name2 = 'SE(2)';
file2 = [folder 'SE/Time_BB.h5'];
fileeig2 = [folder 'SE/beam_eig.h5'];
pose0_SE = pose0.SE23_cclamped;

pose1 = h5read(file1,'/dynamic_analysis/FEModel/POSE/MOTION').';
vel1 = h5read(file1,'/dynamic_analysis/FEModel/VELOCITY/MOTION').';
eig1 = h5read(fileeig1,'/eigen_analysis/Eigenvectors/MOTION').';
p1 = (pose1(dof_VK,:)-pose0_VK(dof_VK))/beamthickness;
t1 = h5read(file1,'/dynamic_analysis/FEModel/time');
color1 = 'k';

pose2 = h5read(file2,'/dynamic_analysis/FEModel/POSE/MOTION').';
vel2 = h5read(file2,'/dynamic_analysis/FEModel/VELOCITY/MOTION').';
eig2 = h5read(fileeig2,'/eigen_analysis/Eigenvectors/MOTION').';
p2 = (pose2(dof_SE,:)-pose0_SE(dof_SE))/beamthickness;
t2 = h5read(file2,'/dynamic_analysis/FEModel/time');
color2 = "#D95319";

%% Plot Pose
plot(a1,t1,p1,'-','LineWidth',plotlinew,'Color',color1,'DisplayName',name1);
plot(a1,t2,p2,'-.','LineWidth',plotlinew,'Color',color2,'DisplayName',name2);

%% Plot Config Space
v1 = vel1(dof_VK,:);
plot(a3,p1,v1,'-','LineWidth',plotlinew,'Color',color1,'DisplayName',name1);

% Rotate SE velocities
for j=1:length(t2)
    for i=0:30
        q1 = pose2(4*i+1,j);  %cos
        q2 = pose2(4*i+2,j);  %sin
        % either theta calculation works fine
        theta = 2*atan(q2/q1);
%         theta = sign(q2)*2*acos(q1);
        vx = vel2(3*i+1,j);
        vy = vel2(3*i+2,j);
        V = [cos(theta),-sin(theta);sin(theta),cos(theta)]*[vx;vy];
        vel2(3*i+1,j) = V(1);
        vel2(3*i+2,j) = V(2);
    end
end
v2 = vel2(dof_VK,:);  % dof VK is correct as velocity has same size as VK
plot(a3,p2,v2,'-.','LineWidth',plotlinew,'Color',color2,'DisplayName',name2);

%% FFT
Fs = 1/(t1(2)-t1(1));
L = length(t1)-1;
L2 = floor(L/2);
ZF = fft(p1(1:end-1));
ZM = abs(ZF)/L;
ZP = angle(ZF);
P = ZM(1:L2+1);
P(2:end-1) = 2*P(2:end-1);
P2 = ZP(1:L2+1);
f = Fs*(0:(L2))/L;
fund_VK = 1/t1(end);  % fundamental frequencies are the end times
plot(a2,f/fund_VK,P,'-','LineWidth',plotlinew,'Color',color1,'DisplayName',name1);
% fitobj1 = fit(t1,p1','sin5');
% fitobj1.b1/(2*pi*fund_VK)
Q1 = P; F1 = f/fund_VK;

Fs = 1/(t2(2)-t2(1));
L = length(t2)-1;
L2 = floor(L/2);
ZF = fft(p2(1:end-1));
ZM = abs(ZF)/L;
ZP = angle(ZF);
P = ZM(1:L2+1);
P(2:end-1) = 2*P(2:end-1);
P2 = ZP(1:L2+1);
f = Fs*(0:(L2))/L;
fund_SE = 1/t2(end);  % fundamental frequencies are the end times
plot(a2,f/fund_SE,P,'-.','LineWidth',plotlinew,'Color',color2,'DisplayName',name2);
% fitobj2 = fit(t2,p2','sin5');
% fitobj2.b1/(2*pi*fund_SE)
% fitobj2.b2/(2*pi*fund_SE)
% fitobj2.b3/(2*pi*fund_SE)
% fitobj2.b4/(2*pi*fund_SE)
% fitobj2.b5/(2*pi*fund_SE)
Q2 = P; F2 = f/fund_SE;

Q = [Q1; Q2];
h = bar(a4,F2,Q);

%% Project onto Modes
def1 = pose1(:,1) - pose0_VK;
proj1=abs(eig1\def1);
figure;bar(proj1,'FaceColor',color1);

% Using Relative Inc but make 3d first
pose0_SE_3d = reshape(pose0_SE,[4,31]);
pose0_SE_3d = [pose0_SE_3d(1,:);zeros(2,31);pose0_SE_3d(2,:);pose0_SE_3d(3:end,:);zeros(1,31)];
pose0_SE_3d = pose0_SE_3d(:);
pose2_3d = reshape(pose2(:,1),[4,31]);
pose2_3d = [pose2_3d(1,:);zeros(2,31);pose2_3d(2,:);pose2_3d(3:end,:);zeros(1,31)];
pose2_3d = pose2_3d(:);
def2 = relative_inc(pose0_SE_3d,pose2_3d);
def2 = reshape(def2,[6,31]);
def2 = [def2(1:2,:);def2(6,:)];
def2 = def2(:);
proj2=abs(eig2\def2);
figure;bar(proj2,'FaceColor',color2);

% Alternatively, taking absolute difference between poses. convert to 3dof first
pose0_SE_row = reshape(pose0_SE,[4,31]);
theta = 2*atan(pose0_SE_row(2,:)./pose0_SE_row(1,:));
pose0_SE_3dof = [pose0_SE_row(3:4,:);theta];
pose0_SE_3dof = pose0_SE_3dof(:);
pose2_row = reshape(pose2(:,1),[4,31]);
theta = 2*atan(pose2_row(2,:)./pose2_row(1,:));
pose2_3dof = [pose2_row(3:4,:);theta];
pose2_3dof = pose2_3dof(:);
def2 = pose2_3dof - pose0_SE_3dof;
proj2=abs(eig2\def2);
figure;bar(proj2,'FaceColor','m');