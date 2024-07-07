clear all; close all;

%% GIVE FLOQUET and AMPLITUDE PLOT FRF.
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
plotlinew = 1.5;

f2 = figure('Units','centimeters','PaperUnits','centimeters',...
    'Position',[2 2 fsizex fsizey],'PaperPosition',[0 0 fsizex fsizey],...
    'PaperSize',[fsizex fsizey],'PaperPositionMode','manual');
a2 = axes(f2,'FontSize',afont,'LineWidth',aline,...
    'TickLabelInterpreter','latex','Box','on');
a2.XLabel.Interpreter = "latex";
a2.YLabel.Interpreter = "latex";
grid(a2, 'on'); hold(a2, 'on');

a2.XLimitMethod = 'padded';
a2.YLimitMethod = 'padded';
axis(a2, 'square');
ylabel(a2, '${Im}$');
xlabel(a2, '${Re}$');

%% Data Load
files1 = {}; 
folder = '//wsl$/ubuntu/home/akb110/Codes/PyCont/examples/beam_cpp/';
files1{end+1} = [folder 'NNM1a_withtime.h5'];
files1{end+1} = [folder 'NNM1b_withtime.h5'];
files1{end+1} = [folder 'NNM1c_withtime.h5'];

%%
for i=1:length(files1)
    floq = h5read(files1{i},'/Bifurcation/Floquet');
    floq2plot_r = floq.r(:,:)';
    floq2plot_i = floq.i(:,:)';
    plot(a2,floq2plot_r,floq2plot_i,'o','MarkerFaceColor','none','MarkerEdgeColor',"#4DBEEE",'LineWidth',plotlinew);
    
    circle = linspace(0,2*pi,1000);
    plot(a2,cos(circle),sin(circle),'LineWidth',2.0,'LineStyle','--','Color',"#D95319");
end