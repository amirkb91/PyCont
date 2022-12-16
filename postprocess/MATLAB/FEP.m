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

% folder = 'cclamped/NNM1/';
% name1 = 'vonK';
% files1{end+1} = [folder 'VK/NNM1_VK.h5'];
% name2 = 'SE(2)';
% files2{end+1} = [folder 'SE23/NNM1_SE23_1.h5'];
% files2{end+1} = [folder 'SE23/NNM1_SE23_2.h5'];
% files2{end+1} = [folder 'SE23/NNM1_SE23_3.h5'];
%===============================================
% folder = 'cclamped/NNM2/';
% name1 = 'vonK';
% files1{end+1} = [folder 'VK/NNM2_VK_1.h5'];
% files1{end+1} = [folder 'VK/NNM2_VK_2.h5'];
% files1{end+1} = [folder 'VK/NNM2_VK_3.h5'];
% files1{end+1} = [folder 'VK/NNM2_VK_4_trunc.h5'];
% name2 = 'SE(2)';
% files2{end+1} = [folder 'SE23/NNM2_SE23_1.h5'];
% files2{end+1} = [folder 'SE23/NNM2_SE23_2.h5'];
% files2{end+1} = [folder 'SE23/NNM2_SE23_3.h5'];
% files2{end+1} = [folder 'SE23/NNM2_SE23_4.h5'];
% files2{end+1} = [folder 'SE23/NNM2_SE23_5.h5'];
%===============================================
% folder = 'arch/NNM1/';
% name1 = 'vonK';
% files1{end+1} = [folder 'VK/NNM1_VK_1.h5'];
% files1{end+1} = [folder 'VK/NNM1_VK_2.h5'];
% files1{end+1} = [folder 'VK/NNM1_VK_3.h5'];
% files1{end+1} = [folder 'VK/NNM1_VK_4_trunc.h5'];
% name2 = 'SE(2)';
% files2{end+1} = [folder 'SE23/NNM1_SE23_1.h5'];
% files2{end+1} = [folder 'SE23/NNM1_SE23_2.h5'];
% files2{end+1} = [folder 'SE23/NNM1_SE23_3.h5'];
% files2{end+1} = [folder 'SE23/NNM1_SE23_4.h5'];
%===============================================
folder = 'arch/NNM2/';
name1 = 'vonK';
files1{end+1} = [folder 'VK/NNM2_VK_1_trunc.h5'];
files1{end+1} = [folder 'VK/NNM2_VK_2.h5'];
files1{end+1} = [folder 'VK/NNM2_VK_3.h5'];
files1{end+1} = [folder 'VK/NNM2_VK_4.h5'];
name2 = 'SE(2)';
files2{end+1} = [folder 'SE23/NNM2_SE23_1.h5'];
files2{end+1} = [folder 'SE23/NNM2_SE23_2.h5'];
files2{end+1} = [folder 'SE23/NNM2_SE23_3.h5'];
files2{end+1} = [folder 'SE23/NNM2_SE23_4.h5'];
files2{end+1} = [folder 'SE23/NNM2_SE23_5.h5'];
files2{end+1} = [folder 'SE23/NNM2_SE23_6.h5'];
files2{end+1} = [folder 'SE23/NNM2_SE23_7.h5'];
files2{end+1} = [folder 'SE23/NNM2_SE23_8.h5'];
files2{end+1} = [folder 'SE23/NNM2_SE23_9.h5'];
files2{end+1} = [folder 'SE23/NNM2_SE23_10.h5'];

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
