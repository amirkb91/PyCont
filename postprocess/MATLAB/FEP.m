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

a.XScale = 'log';
a.YLimMode = 'auto';
a.XLimitMethod = 'padded';
a.XTick = [1e-3 1e-2 1e-1 1 1e1 1e2 1e3];
xlabel(a, 'Energy (J)');
ylabel(a, 'Frequency (Hz)');

%% Data Load
files1 = {}; files2 = {};
folder = '//wsl$/ubuntu/home/akb110/Codes/PyCont/examples/beam_rightangle/Results/';
name1 = 'VK';
files1{end+1} = [folder 'nnm1_VK/nnm1.h5'];
files1{end+1} = [folder 'nnm1_VK/nnm1b.h5'];
files1{end+1} = [folder 'nnm1_VK/nnm1c.h5'];
files1{end+1} = [folder 'nnm1_VK/nnm1d.h5'];
files1{end+1} = [folder 'nnm1_VK/nnm1e.h5'];
name2 = 'SE';
files2{end+1} = [folder 'nnm1_SE/nnm1.h5'];
files2{end+1} = [folder 'nnm1_SE/nnm1b.h5'];
files2{end+1} = [folder 'nnm1_SE/nnm1c.h5'];
files2{end+1} = [folder 'nnm1_SE/nnm1d.h5'];
files2{end+1} = [folder 'nnm1_SE/nnm1e.h5'];

%% Plot
for i=1:length(files1)
    E = h5read(files1{i},'/Energy');
    T = h5read(files1{i},'/T');
    f = 1./T;
    color = 'k';
    if i==1
        plot(a,E,f,'-','LineWidth',plotlinew,'Color',color,'DisplayName',name1);
    else
        plot(a,E,f,'-','LineWidth',plotlinew,'Color',color,'DisplayName',name1,...
            'HandleVisibility','off');
    end
end

for i=1:length(files2)
    E = h5read(files2{i},'/Energy');
    T = h5read(files2{i},'/T');
    f = 1./T;
    color = [0.8500 0.3250 0.0980];
    if i==1
        plot(a,E,f,'-.','LineWidth',plotlinew,'Color',color,'DisplayName',name2);
    else
        plot(a,E,f,'-.','LineWidth',plotlinew,'Color',color,'DisplayName',name2,...
            'HandleVisibility','off');
    end
end
