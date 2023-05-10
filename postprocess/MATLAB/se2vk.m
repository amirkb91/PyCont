function [pose_VK] = se2vk(pose0, pose)
% input: pose0 is SE reference undeformed pose, pose is SE deformed pose
% finds the initial conditions for a VK sim, which corresponds to the
% deformation of pose

nnodes = length(pose0)/7;
pose0 = reshape(pose0,7,nnodes);
pose = reshape(pose,7,nnodes);

% for VK MATLAB code, XYZ are in terms of displacements, not positions
% xyz = pose_B(5:7,:) - pose_A(5:7,:);
% for VK cpp code, XYZ are positions.
xyz = pose(5:7,:);

theta = zeros(size(xyz));

for i=1:nnodes
    q_A = pose0(1:4,i);
    q_B = pose(1:4,i);
    
    q_A_m1 = [q_A(1);-q_A(2);-q_A(3);-q_A(4)];
    q_rel = SO3_SE3.quatmult(q_A_m1, q_B);
    
    psi = 2*atan2(sqrt(q_rel(2)^2+q_rel(3)^2+q_rel(4)^2),q_rel(1));
    theta(:,i) = [q_rel(2);q_rel(3);q_rel(4)] * psi/sin(psi/2);
%     theta(:,i) = SO3_SE3.quat2euler(q_rel);  % does same job as above
    % rotate to local frame
    R = SO3_SE3.quat2R(q_A);
    theta(:,i) = R * theta(:,i);
end
theta(isnan(theta))=0;

pose_VK = [xyz; theta];
pose_VK = pose_VK(:);

end