%% Generate initial condition vector for my MATLAB vonK code using the NNM solution
% want to check the vector of nonlinear forces from the initial condition
clear all; close all;

%% Data Load
folder = '//wsl$/ubuntu/home/akb110/Codes/PyCont/examples/beam_rightangle/Results/';
file_VK = [folder 'nnm9_VK/nnm9.h5'];
file_SE = [folder 'nnm9_SE/nnm9.h5'];

analysis = "VK";
sol_number = 16;  % NNM solution to do FFT on
node_number = 21;  % elbow (starting from 0) elbow 21, tip 43

if analysis == "VK"
    file = file_VK;
    vec_len = 6;
    ix = 1; iy = 2; iz = 3;
elseif analysis == "SE"
    file = file_SE;
    vec_len = 7;
    ix = 5; iy = 6; iz = 7;
end

T = h5read(file,'/T',sol_number,1);
pose = permute(h5read(file,'/Config/POSE_time'),[3,2,1]);
pose = squeeze(pose(:,:,sol_number));
pose_node = pose(vec_len*node_number + iz, :)';

%% FFT of single DoF for a particular solution
npts = length(pose_node);
time = linspace(0,T,npts);
Fs = 1/(time(2)-time(1));
L = npts-1;
L2 = floor(L/2);
ZF = fft(pose_node);
ZM = abs(ZF)/L;
ZP = angle(ZF);
P = ZM(1:L2+1);
P(2:end-1) = 2*P(2:end-1);
P2 = ZP(1:L2+1);
f = Fs*(0:(L2))/L;
subplot(1,2,1);
plot(time,pose_node,'o-m');grid;
subplot(1,2,2);
bar(f*T,P);
xlim([0 10])

%% Generate IC vector.
load pose0_beamrightangle.mat
if analysis == "VK"
    % pose - pose0 is already equal to inc
    pose0 = pose0.VK;
    ic = pose(:,1) - pose0;
   
elseif analysis == "SE"
    pose0 = pose0.SE;
    ic = se2vk(pose0, pose(:,1));    
end

% remove repeated elbow node and re-order (base/elbow/tip/rest) remove fix dof
% % elbow_node = 21;
% % elbow_ind = 6*elbow_node+1:6*elbow_node+6;
% % ic(elbow_ind) = [];
% % ic = [ic(elbow_ind);ic(end-6+1:end);ic(6+1:elbow_ind(1)-1);ic(elbow_ind(end)+1:end-6)];

% store solution for my MATLAB code
q0 = ic;
timeint.T     = T;
timeint.dt    = time(2)-time(1);
timeint.nstep = length(time);
timeint.time  = time;
save("ic_NNM.mat","q0","timeint");

% h5 file
delete 'ic.h5';
h5create('ic.h5','/dynamic_analysis/FEModel/POSE/MOTION',[1 length(ic)]);
h5write('ic.h5','/dynamic_analysis/FEModel/POSE/MOTION',ic');

%%
% load('del.mat');
% close all
% I=1;nexttile;plot(q0((0:43)*6+I),'r');hold on;plot(q0_VK((0:43)*6+I),'b'); title('1');
% I=2;nexttile;plot(q0((0:43)*6+I),'r');hold on;plot(q0_VK((0:43)*6+I),'b'); title('2');
% I=3;nexttile;plot(q0((0:43)*6+I),'r');hold on;plot(q0_VK((0:43)*6+I),'b'); title('3');
% I=4;nexttile;plot(q0((0:43)*6+I),'r');hold on;plot(q0_VK((0:43)*6+I),'b'); title('4');
% I=5;nexttile;plot(q0((0:43)*6+I),'r');hold on;plot(q0_VK((0:43)*6+I),'b'); title('5');
% I=6;nexttile;plot(q0((0:43)*6+I),'r');hold on;plot(q0_VK((0:43)*6+I),'b'); title('6');
