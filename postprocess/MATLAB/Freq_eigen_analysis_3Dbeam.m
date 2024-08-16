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
grid(a1, 'on'); hold(a1, 'on');
a1.YLimMode = 'auto';
a1.XLimitMethod = 'tight';
legend(a1,'Interpreter','latex','Location','northwest');
axis(a1,'square')
xlabel(a1, 'Continuation Step');
ylabel(a1, 'Normalised Modal Amplitude');

f2 = figure('Units','centimeters','PaperUnits','centimeters',...
    'Position',[2 2 fsizex fsizey],'PaperPosition',[0 0 fsizex fsizey],...
    'PaperSize',[fsizex fsizey],'PaperPositionMode','manual');
a2 = axes(f2,'FontSize',afont,'LineWidth',aline,...
    'TickLabelInterpreter','latex','Box','on');
a2.XLabel.Interpreter = "latex";
a2.YLabel.Interpreter = "latex";
grid(a2, 'on'); hold(a2, 'on');
a2.YLimitMethod = 'tight';
a2.XLimitMethod = 'tight';
axis(a2,'square')
xlabel(a2, '$\Omega/\omega_4$');
ylabel(a2, 'Mode Number');
colormap(a2, 'turbo');
colorbar(a2);
caxis(a2, [0 1]);
a2.YDir = 'reverse';

f3 = figure('Units','centimeters','PaperUnits','centimeters',...
    'Position',[2 2 fsizex fsizey],'PaperPosition',[0 0 fsizex fsizey],...
    'PaperSize',[fsizex fsizey],'PaperPositionMode','manual');
a3 = axes(f3,'FontSize',afont,'LineWidth',aline,...
    'TickLabelInterpreter','latex','Box','on');
a3.XLabel.Interpreter = "latex";
a3.YLabel.Interpreter = "latex";
a3.ZLabel.Interpreter = "latex";
grid(a3, 'on'); hold(a3, 'on');
xlabel(a3, 'Harmonic');
ylabel(a3, 'Step');
zlabel(a3, 'Amplitude');
a3.View=[18 25];

f4 = figure('Units','centimeters','PaperUnits','centimeters',...
    'Position',[2 2 fsizex fsizey],'PaperPosition',[0 0 fsizex fsizey],...
    'PaperSize',[fsizex fsizey],'PaperPositionMode','manual');
a4 = axes(f4,'FontSize',afont,'LineWidth',aline,...
    'TickLabelInterpreter','latex','Box','on');
a4.XLabel.Interpreter = "latex";
a4.YLabel.Interpreter = "latex";
grid(a4, 'on'); hold(a4, 'on');
a4.YLimitMethod = 'tight';
axis(a4,'square')
xlabel(a4, 'Harmonic');
ylabel(a4, 'Step');
colormap(a4, 'turbo');
colorbar(a4);
a4.YDir = 'normal';

%% Data Load
files1 = {};
folder = 'C:\Users\akb110\OneDrive - Imperial College London\PhD Files\Simulation Results\LieC++-PyCont_Thesis2024\Boxwing/';
files1{end+1} = [folder 'FRC4_amp08_withtime.h5'];
% files1{end+1} = [folder 'right_sing_NNM1b_withtime.h5'];
% files1{end+1} = [folder 'right_sing_NNM1c_withtime.h5'];
% files1{end+1} = [folder 'right_sing_NNM1d_trunc_withtime.h5'];

file_eig = 'C:\Users\akb110\OneDrive - Imperial College London\PhD Files\Simulation Results\LieC++-PyCont_Thesis2024\Boxwing\beam_eig.h5';
eig = h5read(file_eig,'/eigen_analysis/Eigenvectors/MOTION').';

%% Selections
num_modes = 15;


% node_number = 36; % right wing tip front
% node_number = 48; % right wing tip rear 
% node_number = 63; % vtail top
node_number = 79; % fuselage rear

configpernode = 7;
xconfig = configpernode*node_number + 5;  % x direction
yconfig = configpernode*node_number + 6;  % y direction
zconfig = configpernode*node_number + 7;  % y direction


% node_number = 42; % first node is 0
% xconfig = node_number*7+5;
% yconfig = node_number*7+6;
% zconfig = node_number*7+7;

normalise_frq = 7.616919386305014;

%% Modal amplitude and MAC plot
modal_amp = zeros(size(eig,2),0);
MAC = zeros(size(eig,2),0);
INC = {};
POSE = {};
T = {};
for i=1:length(files1)
    inc = h5read(files1{i},'/Config_Time/INC');
    pose = h5read(files1{i},'/Config_Time/POSE');
    inc = permute(inc,[3,2,1]);
    pose = permute(pose,[3,2,1]);
    time = h5read(files1{i},'/T');
    INC{i} = inc;
    POSE{i} = pose;
    T{i} = time;
    f = 1./T{i}/normalise_frq;
    for j=1:size(inc,3)
        proj = abs(eig\inc(:,1,j));
        proj = proj / norm(proj);
        modal_amp = [modal_amp,proj];
    end
    for i2=1:size(eig,2)
        for j2=1:size(inc,3)
        mac(i2,j2) = (eig(:,i2).'*inc(:,1,j2))^2 / ((eig(:,i2).'*eig(:,i2))*(inc(:,1,j2).'*inc(:,1,j2)));
        end
    end
    MAC = [MAC, mac];
    mac = [];
end
INC = cat(3, INC{:});
POSE = cat(3, POSE{:});
T = cat(1,T{:});

for i=1:num_modes
    plot(a1,f, modal_amp(i,:),'-','LineWidth',plotlinew,'DisplayName',sprintf("Mode %i",i));
end
imagesc(a2, f,(1:num_modes),MAC(1:num_modes,:))
colormap(a2,flipud(bone(20)));
caxis(a2, [0.2 1])
a2.YTick = (1:num_modes);
%% Fourier Analysis
for i=1:size(INC,3)
    t1 = linspace(0,T(i),size(INC,2));
    p1 = POSE(zconfig,:,i);
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
    fund = 1/t1(end);  % fundamental frequencies are the end times
%     f = f/fund;
    p_all(i,:) = P;
    f_all(i,:) = f;
end

%% WATERFALL PLOT & HEATMAP
num_signals = size(p_all, 1); % Number of signals
num_frequencies = size(p_all, 2); % Number of frequency points per signal

% Plot each spectrum with an offset in the y-axis
for i = 1:num_signals
    plot3(a3, f_all(i, :), i * ones(1, num_frequencies), p_all(i, :), 'LineWidth', plotlinew);
end
xlim(a3, [0 num_modes*2]);

imagesc(a4, f_all(1, :), 1:size(p_all, 1), p_all);
xlim(a4, [0.5 num_modes+0.5]);