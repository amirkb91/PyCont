clear all; close all;

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
xlabel(a1, '$\Omega$ (Hz)');
ylabel(a1, 'Instantaneous Frequency (Hz)');

%% Data Load
folder = 'C:\Users\akb110\OneDrive - Imperial College London\PhD Files\Simulation Results\LieC++-PyCont_Thesis2024\Vertical_Cantilever\Sweep_Final\';
file = 'beam_sim_180.h5';
pose = h5read([folder file],'/dynamic_analysis/FEModel/ACCELERATION/MOTION').';

node_number = 23; % first node is 0
xconfig = node_number*4+3;
yconfig = node_number*4+4;
beam_length = 72e-3;
nat_freq = 6.374;
xdof = node_number*3+2;

%% Sort out single outlier and select xconfig pose
% index = find(pose(yconfig,:)==0);
% pose(xconfig,index) = 0.5*(pose(xconfig,index-1)+pose(xconfig,index+1));
% pose(yconfig,index) = 0.5*(pose(yconfig,index-1)+pose(yconfig,index+1));
X = pose(xdof,:);

%% Sweep parameters which created the sine sweep signal
f0 = 10;
f1 = 20;
sweep_rate = 1;
sampling_freq = 4000;
tend = abs(f1-f0)*(60/sweep_rate);
k = (f1-f0)/tend;
time_signal = 0:1/sampling_freq:tend;
finst_signal = k*time_signal + f0;

%% Wavelet

dec = 2;
fac2plot = 20;
nsamples  = length(X);
X = X(1:dec:end);

% Wavelet parameters
fs = 1000;
time = 0:1/fs:(nsamples-1)/fs ;
fi = 0;
ff = 60;
nf = 500; % number of frequency lines
f_morlet = 5; % Morlet frequency

% perform wavelet analysis
[~, wtinst, timeWT, y, freqWT] = myWT(X', fs/dec, fi, ff, nf, f_morlet ,1);

[~,index_time] = min(abs(timeWT - time(end)));
if timeWT(index_time) > time(end), index_time = index_time-1; end
finst = interp1(time, finst_signal, timeWT(1:index_time));

% plot results
figure;
plot(finst_signal(1:dec:end), X); xlim([finst_signal(1) finst_signal(end)])
xlabel('Sweep Frequency (Hz)');
ylabel('Normalised Position');

% figure;
% imagesc(finst(1:fac2plot:end), freqWT, (abs(y(1:fac2plot:index_time,:))'));
% axis xy;
% xlabel('Sweep Frequency (Hz)');
% ylabel('Instantaneous Frequency (Hz)');
% c = colorbar;
% colormap jet;

max_abs_y = max(max(abs(y(1:fac2plot:index_time,:))));
levels = 20;
% levels = linspace(20,100,10);
% contourf(a1, finst(1:fac2plot:end), freqWT, abs(y(1:fac2plot:index_time,:))', levels, 'LineColor', 'none');
contourf(a1, finst(1:fac2plot:end), freqWT, 20*log10(abs(y(1:fac2plot:index_time,:))'), levels, 'LineColor', 'none');
colorbar(a1);
colormap(a1, jet(20)); % Use fewer colors for stronger contrast

set(a1, 'FontSize', afont, 'LineWidth', aline, 'TickLabelInterpreter', 'latex', 'Box', 'on');
a1.XLabel.Interpreter = "latex";
a1.YLabel.Interpreter = "latex";
xlabel(a1, '$\Omega$ (Hz)');
ylabel(a1, 'Instantaneous Frequency (Hz)');
% xlim([11 15])
% ylim([0 30])

% caxis([0.01, 0.05]); % Set color limits
% hold on;
% contour(finst(1:fac2plot:end), freqWT, abs(y(1:fac2plot:index_time,:))', levels, 'LineColor', 'k', 'LineWidth', 1);
% ylim([20 40])

