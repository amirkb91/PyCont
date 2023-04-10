function [inc] = relative_inc(pose_A,pose_B)

nnodes = length(pose_A)/7;
pose_A = reshape(pose_A,7,nnodes);
pose_B = reshape(pose_B,7,nnodes);

Frame_A = zeros(4,4,nnodes);
Frame_B  = zeros(4,4,nnodes);
HAm1_HB = zeros(4,4,nnodes);
inc = zeros(nnodes*6,1);

for i=1:nnodes
    Frame_A(:,:,i) = [SO3_SE3.quat2R(pose_A(1:4,i)) pose_A(5:7,i); 0 0 0 1];
    Frame_B(:,:,i) = [SO3_SE3.quat2R(pose_B(1:4,i)) pose_B(5:7,i); 0 0 0 1];
    HAm1_HB(:,:,i) = SO3_SE3.InvSE3(Frame_A(:,:,i))*Frame_B(:,:,i);
    inc((i-1)*6+1:(i-1)*6+6) = SO3_SE3.LogSE3(HAm1_HB(:,:,i));
%     inc((i-1)*6+1:(i-1)*6+6) = SO3_SE3.get_parameters_from_frame(HAm1_HB(:,:,i));
    
    R = SO3_SE3.quat2R(pose_A(1:4,i));
    inc((i-1)*6+1:(i-1)*6+3) = R * inc((i-1)*6+1:(i-1)*6+3);
    inc((i-1)*6+4:(i-1)*6+6) = R * inc((i-1)*6+4:(i-1)*6+6);
end

end