function [pose] = vk2se(pose0, pose_VK)
% input: pose0 is SE reference undeformed pose, pose_VK is deformed VK pose
% finds the initial conditions for an SE sim, which corresponds to the
% deformation of VK

nnodes = length(pose0)/7;
pose0 = reshape(pose0,7,nnodes);
pose_VK = reshape(pose_VK,6,nnodes);

pose = zeros(7,nnodes);
pose(5:7,:) = pose_VK(1:3,:);

for i=1:nnodes
    q_A = pose0(1:4,i);
    R = SO3_SE3.quat2R(q_A);
    
    % rotate thetas back first
    theta = pose_VK(4:6,i);
    theta = R' * theta;
    
    q_rel = SO3_SE3.euler2quat(theta);
    q_B = SO3_SE3.quatmult(q_A, q_rel);
    pose(1:4,i) = q_B;
end

pose = pose(:);

end