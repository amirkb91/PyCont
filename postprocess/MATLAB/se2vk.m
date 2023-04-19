function [ic_vk] = se2vk(pose_A,pose_B)

nnodes = length(pose_A)/7;
pose_A = reshape(pose_A,7,nnodes);
pose_B = reshape(pose_B,7,nnodes);

% for VK MATLAB code, XYZ are in terms of displacements, not positions
xyz_disp = pose_B(5:7,:) - pose_A(5:7,:);
theta = zeros(size(xyz_disp));

for i=1:nnodes
    q_A = pose_A(1:4,i);
    q_B = pose_B(1:4,i);
    
    q_A_m1 = [q_A(1);-q_A(2);-q_A(3);-q_A(4)];
    q_rel = SO3_SE3.quatmult(q_A_m1, q_B);
    
    psi = 2*atan2(sqrt(q_rel(2)^2+q_rel(3)^2+q_rel(4)^2),q_rel(1));
    theta(:,i) = [q_rel(2);q_rel(3);q_rel(4)] * psi/sin(psi/2);
    
    % rotate to local frame
    R = SO3_SE3.quat2R(q_A);
    theta(:,i) = R * theta(:,i);
end
theta(isnan(theta))=0;

ic_vk = [xyz_disp; theta];
ic_vk = ic_vk(:);

end