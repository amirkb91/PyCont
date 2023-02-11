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
% a2.XLim = [0 12];
% a2.XTick = (0:1:12);

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

%% DoF Index (+1 for MATLAB starting index) 
% node10 = 0.3 / node12 = 0.36
node = 1;
% X
% dof_VK   = 3*node + 1;  dof_SE23 = 4*node + 3;
% Y
dof_VK   = 3*node + 2;  dof_SE23 = 4*node + 4;

%% Fundamental Frequency for FFT
% read value from FEP corresponding to point for which time plot is taken
fund_VK = 45.7;
fund_SE = 48.83;

%% Data Load
load('pose0.mat');
%===============================================
folder = '//wsl$/ubuntu/home/akb110/Codes/PyCont/examples/beam_2D/Results/cantilever/';
name1 = 'vonK';
file1 = [folder 'NNM2/Time_NNM2_VK.h5'];
file1_eig = [folder 'NNM2/beam_eig_VK.h5'];
name2 = 'SE(2)';
file2 = [folder 'NNM2/Time_NNM2_SE.h5'];
pose0_VK = pose0.VK_cclamped;
pose0_SE23 = pose0.SE23_cclamped;

%% Read data
pose1 = h5read(file1,'/Config/POSE').';
vel1 = h5read(file1,'/Config/VELOCITY').';
p1 = (pose1(dof_VK,:)-pose0_VK(dof_VK))/beamthickness;
v1 = vel1(dof_VK,:);
t1 = h5read(file1,'/time');
eig1 = h5read(file1_eig,'/eigen_analysis/Eigenvectors').';

%%%%
% mode2plot = 8;
% figure;plot(eig1(3*(0:30)+1,mode2plot),'o');
def1 = pose1(:,1) - pose0_VK;
proj=abs(eig1\def1);
figure;bar(proj);

pose2 = h5read(file2,'/Config/POSE').';
vel2 = h5read(file2,'/Config/VELOCITY').';
p2 = (pose2(dof_SE23,:)-pose0_SE23(dof_SE23))/beamthickness;
t2 = h5read(file2,'/time');
vel2_rot = vel2;
% Rotate velocities
for j=1:length(t2)
    for i=0:30
        q1 = pose2(4*i+1,j);  %cos
        q2 = pose2(4*i+2,j);  %sin
        % either theta calculation works fine
%         theta = 2*atan(q2/q1);
        theta = sign(q2)*2*acos(q1);
        vx = vel2(3*i+1,j);
        vy = vel2(3*i+2,j);
        V = [cos(theta),-sin(theta);sin(theta),cos(theta)]*[vx;vy];
        vel2_rot(3*i+1,j) = V(1);
        vel2_rot(3*i+2,j) = V(2);
    end
end
v2 = vel2_rot(dof_VK,:);

%% Plot Pose
plot(a1,t1,p1,'-','LineWidth',plotlinew,'Color','k','DisplayName',name1);
plot(a1,t2,p2,'-.','LineWidth',plotlinew,'Color','r','DisplayName',name2);

%% Plot Config Space
plot(a3,p1,v1,'-','LineWidth',plotlinew,'Color','k','DisplayName',name1);
plot(a3,p2,v2,'-.','LineWidth',plotlinew,'Color','r','DisplayName',name2);
% plot(a3,p2,vel2(dof_VK,:),'-.','LineWidth',plotlinew,'Color','g','DisplayName',name2);
%% Plot FFT
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
plot(a2,f/fund_VK,P,'-','LineWidth',plotlinew,'Color',color,'DisplayName',name1);

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
plot(a2,f/fund_SE,P,'-.','LineWidth',plotlinew,'Color',color,'DisplayName',name2);
