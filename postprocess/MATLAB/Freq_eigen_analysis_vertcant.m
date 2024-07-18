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
legend(a1,'Interpreter','latex','Location','northwest');
xlabel(a1, '$\Omega$ (Hz)');
ylabel(a1, 'Normalised Modal Amplitude');


%% Data Load
folder = 'C:\Users\akb110\OneDrive - Imperial College London\PhD Files\Simulation Results\LieC++-PyCont_Thesis2024\Vertical_Cantilever\Sweep_New\wide\';
file = 'beam_sim_sweep_amp180_INC.h5';
inc = h5read([folder file],'/dynamic_analysis/FEModel/INC/MOTION')';

file_eig = [folder 'beam_eig.h5'];
eig = h5read(file_eig,'/eigen_analysis/Eigenvectors/MOTION').';

% remove auxiliary beams from inc and eig
inc(1:3,:) = [];
inc(4:6,:) = [];
eig(1:3,:) = [];
eig(4:6,:) = [];

%% Selections
num_modes = 10;
node_number = 21; % first node is 0, side beams removed
xconfig = node_number*4+3;
yconfig = node_number*4+4;
xdof = node_number*3+1;
ydof = node_number*3+2;

%% Modal amplitude and MAC plot
dec = 20;
num_iterations = ceil(size(inc, 2) / dec);
modal_amp = zeros(size(eig,2),num_iterations);

k = 1;
for j=1:dec:size(inc,2)
    proj = abs(eig\inc(:,j));
    modal_amp(:,k) = proj;
    k = k+1;
end

% Normalise such that max peak of mode 1 is one
factor = max(modal_amp(1,:));
modal_amp = modal_amp / factor;


%%
f0 = 9;
f1 = 20;
sweep_rate = 1;
sampling_freq = 3000;

% f0 = 10;
% f1 = 18;
% sweep_rate = 1;
% sampling_freq = 4000;

tend = abs(f1-f0)*(60/sweep_rate);
k = (f1-f0)/tend;
time = 0:1/sampling_freq:tend ;
finst = k*time + f0 ;
finst = finst(1:dec:end);

% for i=1:num_modes
%     plot(a1,finst, modal_amp(i,:),'-','LineWidth',plotlinew,'DisplayName',sprintf("Mode %i",i-1));
% end

[max_values, indices] = max(modal_amp(1:num_modes,:), [], 2);
[~, sorted_indices] = sort(max_values, 'descend');
% Plot the lines in the sorted order
customColors = [
    0.0 0.4470 0.7410; % blue
    0.8500 0.3250 0.0980; % orange
    0.9290 0.6940 0.1250; % yellow
    0.4940 0.1840 0.5560; % purple
    0.4660 0.6740 0.1880; % green
    0.3010 0.7450 0.9330; % light blue
    0.6350 0.0780 0.1840; % dark red
    0.75 0.75 0; % olive
    0.0 0.5 0.0; % dark green
    0.5 0 0.5; % dark purple
    0.25 0.25 0.25 % dark grey
    ];
for i = 1:length(sorted_indices)
    idx = sorted_indices(i);
    plot(a1, finst, modal_amp(idx, :), '-', 'LineWidth', plotlinew, 'Color', customColors(i, :), 'DisplayName', sprintf("Mode %i", idx));
end

xlim([10 18]);
ylim([0 1.05]);
