clear all; close all;

%% Figure Properties
fsizex  = 25;
ratio = 3/4;
fsizey  = fsizex*ratio;
afont = 23;
aline = 2;
plotlinew = 2;

f1 = figure('Units','centimeters','PaperUnits','centimeters',...
    'Position',[2 2 fsizex fsizey],'PaperPosition',[0 0 fsizex fsizey],...
    'PaperSize',[fsizex fsizey],'PaperPositionMode','manual');
a1 = axes(f1,'FontSize',afont,'LineWidth',aline,...
    'TickLabelInterpreter','latex','Box','on');
a1.XLabel.Interpreter = "latex";
a1.YLabel.Interpreter = "latex";
grid(a1, 'on'); hold(a1, 'on');
a1.YLimMode = 'auto';
a1.XLimitMethod = 'padded';
a1.YLimitMethod = 'padded';
a1.PlotBoxAspectRatio=[1,1,1];
a1.XLim = [-0.5 1.5];
a1.YLim = [-1 1];
xlabel(a1, 'x/L');
ylabel(a1, 'y/L');

%% Data Load
% file = '//wsl$/ubuntu/home/akb110/Codes/PyCont/examples/beam_2D/Results/cantilever/NNM1/VK_reducedintegration/Time_VK.h5';
file = '//wsl$/ubuntu/home/akb110/Codes/PyCont/examples/beam_2D/Results/cantilever/NNM1/SE/Time_SE_BB.h5';
pose = h5read(file,'/dynamic_analysis/FEModel/POSE/MOTION').';
model = "SE";
if model == "VK"
    xdof = (0:30)*3+1;
    ydof = (0:30)*3+2;
    color = 'k';
    n_snaps_grey = 50;
elseif model == "SE"
    xdof = (0:30)*4+3;
    ydof = (0:30)*4+4;
    color = "#D95319";
    n_snaps_grey = 200;
end
%% Plot Pose
n_snaps = 13;
n_time = (length(pose)-1)/2;

index = ceil(linspace(0,n_time,n_snaps_grey));
for i=1:n_snaps_grey
    plot(a1,pose(xdof,index(i)+1),pose(ydof,index(i)+1),'o','LineWidth',plotlinew,'MarkerSize',0.5,'Color','#bcbcbc');
end

index = ceil(linspace(0,n_time,n_snaps));
for i=1:n_snaps
    plot(a1,pose(xdof,index(i)+1),pose(ydof,index(i)+1),'o-','LineWidth',plotlinew,'MarkerSize',3,'Color',color);
end

