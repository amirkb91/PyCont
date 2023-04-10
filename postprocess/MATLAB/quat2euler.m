function [angles] = quat2euler(q)

qw = q(1); qx = q(2); qy = q(3); qz = q(4);

sinr_cosp = 2 * (qw * qx + qy * qz);
cosr_cosp = 1 - 2 * (qx * qx + qy * qy);
roll = atan2(sinr_cosp, cosr_cosp);

% alternative pitch formula (asin) is more stable
% sinp = sqrt(1 + 2 * (qw * qy - qx * qz));
% cosp = sqrt(1 - 2 * (qw * qy - qx * qz));
% pitch = 2 * atan2(sinp, cosp) - pi / 2;
pitch = asin(2 * (qw * qy - qx * qz));

siny_cosp = 2 * (qw * qz + qx * qy);
cosy_cosp = 1 - 2 * (qy * qy + qz * qz);
yaw = atan2(siny_cosp, cosy_cosp);

angles = [roll;pitch;yaw];
end

