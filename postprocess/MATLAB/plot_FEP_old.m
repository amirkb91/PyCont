clear all; close all;
files = {};

files{end+1} = 'NNM2_VK_1.h5';
files{end+1} = 'NNM2_VK_2.h5';
files{end+1} = 'NNM2_VK_3.h5';
files{end+1} = 'NNM2_VK_4.h5';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ref = 'paper_NNM2.fig';
open(ref);
a = gca;
grid(a, 'on');
hold(a, 'on');
for i=1:length(files)
    E = h5read(files{i},'/Energy');
    T = h5read(files{i},'/T');
    color = 'r';
    plot(a,E,2*pi./T,'.-','LineWidth',1.5,'Color',color);
end
xlim auto; ylim auto;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ref = 'paper_NNM2.h5';
% f = figure;
% a = axes(f);
% a.XScale = 'log';
% grid(a, 'on');
% hold(a, 'on');
% xlabel(a, 'Energy (J)');
% ylabel(a, 'Frequency (Hz)');
% 
% for i=1:length(files)
%     E = h5read(files{i},'/Energy');
%     T = h5read(files{i},'/T');
% %     color = rand(1,3);
%     color = 'r';
%     plot(a,E,1./T,'.-','LineWidth',1.5,'Color',color);
%     plot(a,E(1),1./T(1),'*b','LineWidth',1);
% end
% 
% if false
%     E = h5read(ref,'/Energy');
%     T = h5read(ref,'/T');    
%     plot(a,E,1./T,'-k','LineWidth',1.5);
% end
