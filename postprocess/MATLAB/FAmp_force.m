clear all; close all;

%% GIVE AMPLITUDE PLOT FRF.
% Need to have a h5 file with the time series of all the solutions on the
% NNM branch. run python timesim_branch to create such file

% file names should have _withtime
% pick if you want X Y or theta by commenting out

%% Figure Properties
fsizex  = 25;
ratio = 3/4;
fsizey  = fsizex*ratio;
afont = 23;
aline = 2;
plotlinew = 2.0;

f = figure('Units','centimeters','PaperUnits','centimeters',...
    'Position',[2 2 fsizex fsizey],'PaperPosition',[0 0 fsizex fsizey],...
    'PaperSize',[fsizex fsizey],'PaperPositionMode','manual');
a = axes(f,'FontSize',afont,'LineWidth',aline,...
    'TickLabelInterpreter','latex','Box','on');
a.XLabel.Interpreter = "latex";
a.YLabel.Interpreter = "latex";
grid(a, 'off'); hold(a, 'on');

a.XLim = [0.98 1.08];
a.YLimMode = 'auto';
a.XLimitMethod = 'tight';
ylabel(a, '$y_P$ (m)');
xlabel(a, '$\Omega/\omega_1$');

%% Geometries Undeformed Positions
normalise_amp = 1.0;
normalise_frq = 41.7182203819309;

%% DoF Index (+1 for MATLAB starting index)
node = 21;

configpernode = 4;
config2plot_X = configpernode*node + 3;  % x direction
config2plot_Y = configpernode*node + 4;  % y direction

%% Data Load
files1 = {}; files2 = {};
folder = '//wsl$/ubuntu/home/akb110/Codes/PyCont/examples/beam_cpp/';
name1 = 'NNM';
files1{end+1} = [folder 'NNM1a_withtime.h5'];
files1{end+1} = [folder 'NNM1b_withtime.h5'];
files1{end+1} = [folder 'NNM1c_withtime.h5'];

name2 = 'FRF';
files2{end+1} = [folder 'FRF1_200_withtime.h5'];
files2{end+1} = [folder 'FRF1_400_withtime.h5'];
files2{end+1} = [folder 'FRF1_600_withtime.h5'];
files2{end+1} = [folder 'FRF1_2000_withtime.h5'];

%% Plot
for i=1:length(files1)
    pose = h5read(files1{i},'/Config_Time/POSE');
    pose = permute(pose,[3,2,1]);
    T = h5read(files1{i},'/T');
    f = 1./T;
    
    nnodes = size(pose,1)/4;
    theta = zeros(size(pose,1)/configpernode,size(pose,2),size(pose,3));
    
    for sol_idx = 1:size(pose, 3)
        for time_idx = 1:size(pose, 2)
            for node_idx = 1:size(pose,1)/configpernode
                cos_theta_half_idx = (node_idx - 1) * 4 + 1;
                sin_theta_half_idx = cos_theta_half_idx + 1;
                cos_theta_half = pose(cos_theta_half_idx, time_idx, sol_idx);
                sin_theta_half = pose(sin_theta_half_idx, time_idx, sol_idx);
                theta(node_idx, time_idx, sol_idx) = ...
                    2 * atan2(sin_theta_half, cos_theta_half);
            end
        end
    end
    
    % find normalised max pose for each cont solution
    % following finds the max over T in one go, no need for for loop
    % -1 for x becuase that's the X position of tip undeformed
    [maxpose_X, ~] = max(abs(pose(config2plot_X,:,:)-1)/normalise_amp);
    maxpose_X = 1-reshape(maxpose_X, size(T));
    
    [maxpose_Y, ~] = max(abs(pose(config2plot_Y,:,:))/normalise_amp);
    maxpose_Y = reshape(maxpose_Y, size(T));
    
    [maxtheta, ~] = max(abs(theta(node+1,:,:)));
    maxtheta = reshape(maxtheta, size(T))/pi;
    
    color = 'k';
%     if i==1
%         plot(a,f/normalise_frq,maxpose_X,'-','LineWidth',plotlinew,'Color',color,'DisplayName',name1);
%     else
%         plot(a,f/normalise_frq,maxpose_X,'-','LineWidth',plotlinew,'Color',color,'DisplayName',name1,...
%             'HandleVisibility','off');
%     end
    
    if i==1
        plot(a,f/normalise_frq,maxpose_Y,'-','LineWidth',plotlinew,'Color',color,'DisplayName',name1);
    else
        plot(a,f/normalise_frq,maxpose_Y,'-','LineWidth',plotlinew,'Color',color,'DisplayName',name1,...
            'HandleVisibility','off');
    end
    
%     if i==1
%         plot(a,f/normalise_frq,maxtheta,'-','LineWidth',plotlinew,'Color',color,'DisplayName',name1);
%     else
%         plot(a,f/normalise_frq,maxtheta,'-','LineWidth',plotlinew,'Color',color,'DisplayName',name1,...
%             'HandleVisibility','off');
%     end    
end

for i=1:length(files2)
    pose = h5read(files2{i},'/Config_Time/POSE');
    pose = permute(pose,[3,2,1]);
    T = h5read(files2{i},'/T');
    f = 1./T/normalise_frq;
    
    stability = h5read(files2{i},'/Bifurcation/Stability');
    stable_index = find(diff(stability))+1;
    stable_index = [1, stable_index', length(T)];
    
    nnodes = size(pose,1)/4;
    theta = zeros(size(pose,1)/configpernode,size(pose,2),size(pose,3));
    
    for sol_idx = 1:size(pose, 3)
        for time_idx = 1:size(pose, 2)
            for node_idx = 1:size(pose,1)/configpernode
                cos_theta_half_idx = (node_idx - 1) * 4 + 1;
                sin_theta_half_idx = cos_theta_half_idx + 1;
                cos_theta_half = pose(cos_theta_half_idx, time_idx, sol_idx);
                sin_theta_half = pose(sin_theta_half_idx, time_idx, sol_idx);
                theta(node_idx, time_idx, sol_idx) = ...
                    2 * atan2(sin_theta_half, cos_theta_half);
            end
        end
    end
    
    % find normalised max pose for each cont solution
    % -1 for x becuase that's the X position of tip undeformed
    [maxpose_X, ~] = max(abs(pose(config2plot_X,:,:)-1)/normalise_amp);
    maxpose_X = 1-reshape(maxpose_X, size(T));
    
    [maxpose_Y, ~] = max(abs(pose(config2plot_Y,:,:))/normalise_amp);
    maxpose_Y = reshape(maxpose_Y, size(T));
    
    [maxtheta, ~] = max(abs(theta(node+1,:,:)));
    maxtheta = reshape(maxtheta, size(T))/pi;
    
    color = "#0072BD";
    
%     for j=1:length(stable_index)-1
%         is_stable = stability(stable_index(j+1)-1);
%         if is_stable
%             linestyle = '-';
%         else
%             linestyle = '--';
%         end
%         plot(a,f(stable_index(j):stable_index(j+1)),maxpose_X(stable_index(j):stable_index(j+1)),...
%             'LineStyle',linestyle,'LineWidth',plotlinew,'Color',color,'DisplayName',name2,'HandleVisibility','off');
%     end
    
    for j=1:length(stable_index)-1
        is_stable = stability(stable_index(j+1)-1);
        if is_stable
            linestyle = '-';
        else
            linestyle = '--';
        end
        plot(a,f(stable_index(j):stable_index(j+1)),maxpose_Y(stable_index(j):stable_index(j+1)),...
            'LineStyle',linestyle,'LineWidth',plotlinew,'Color',color,'DisplayName',name2,'HandleVisibility','off');
    end
    
%     for j=1:length(stable_index)-1
%         is_stable = stability(stable_index(j+1)-1);
%         if is_stable
%             linestyle = '-';
%         else
%             linestyle = '--';
%         end
%         plot(a,f(stable_index(j):stable_index(j+1)),maxtheta(stable_index(j):stable_index(j+1)),...
%             'LineStyle',linestyle,'LineWidth',plotlinew,'Color',color,'DisplayName',name2,'HandleVisibility','off');
%     end    
end
axis square
