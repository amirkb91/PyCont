function [phi,theta,psi] = quat2euler(q)
% Converts a quaternion to Euler angles (ZYX convention)
% Input:
% - q: unit quaternion [q1;q2;q3;q4]
% Output:
% - phi: roll angle in radians
% - theta: pitch angle in radians
% - psi: yaw angle in radians

% Extract quaternion elements
q1 = q(1);
q2 = q(2);
q3 = q(3);
q4 = q(4);

% Compute Euler angles
phi = atan2(2*(q3*q4 + q1*q2), q1^2 - q2^2 - q3^2 + q4^2);
theta = asin(2*(q1*q3 - q2*q4));
psi = atan2(2*(q2*q3 + q1*q4), q1^2 + q2^2 - q3^2 - q4^2);

end