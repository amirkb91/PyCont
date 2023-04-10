function [q] = euler2quat(angles)

roll = angles(1); pitch = angles(2); yaw = angles(3);

cr = cos(roll * 0.5);
sr = sin(roll * 0.5);
cp = cos(pitch * 0.5);
sp = sin(pitch * 0.5);
cy = cos(yaw * 0.5);
sy = sin(yaw * 0.5);

qw = cr * cp * cy + sr * sp * sy;
qx = sr * cp * cy - cr * sp * sy;
qy = cr * sp * cy + sr * cp * sy;
qz = cr * cp * sy - sr * sp * cy;

q = [qw;qx;qy;qz];
end

