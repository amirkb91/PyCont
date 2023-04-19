clear all; close all;

% quick script to test some pose to inc and compare

% pose_A = h5read("beam_eig.h5",'/eigen_analysis/POSE/MOTION')';
% INC = h5read("beam_eig.h5",'/eigen_analysis/Eigenvectors/MOTION')';
% INC = INC(:,1);
% pose_B = h5read("beam_sim.h5",'/dynamic_analysis/FEModel/POSE/MOTION')';
% pose_B = pose_B(:,1);

load poseinc.mat

%%
nnodes = length(pose_A)/7;
pose_A = reshape(pose_A,7,nnodes);
pose_B = reshape(pose_B,7,nnodes);
INC = reshape(INC,6,nnodes)';

Frame_A = zeros(4,4,nnodes);
Frame_B  = zeros(4,4,nnodes);
HAm1_HB = zeros(4,4,nnodes);
inc = zeros(nnodes,6);

for i=1:nnodes
    Frame_A(:,:,i) = [SO3_SE3.quat2R(pose_A(1:4,i)) pose_A(5:7,i); 0 0 0 1];
    Frame_B(:,:,i) = [SO3_SE3.quat2R(pose_B(1:4,i)) pose_B(5:7,i); 0 0 0 1];
    HAm1_HB(:,:,i) = SO3_SE3.InvSE3(Frame_A(:,:,i))*Frame_B(:,:,i);
    % in cpp code, logse3 isn't used, parametrisation used instead
%     inc(i,:) = SO3_SE3.LogSE3(HAm1_HB(:,:,i));
    inc(i,:) = SO3_SE3.get_parameters_from_frame(HAm1_HB(:,:,i));
end

%%
imagesc(abs(inc-INC));
colorbar;
xticklabels({'X', 'Y', 'Z', 'Theta_X', 'Theta_Y', 'Theta_Z'});
ylabel('Node')
title("Absolute error of LOG vs INC")

figure;
imagesc(abs(inc-INC)./abs(INC));
colorbar;
xticklabels({'X', 'Y', 'Z', 'Theta_X', 'Theta_Y', 'Theta_Z'});
ylabel('Node')
title("Relative error of LOG vs INC")