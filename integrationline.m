function [angle] = integrationline(angular_vel)

sum = 0;
for i = 1:length(angular_vel)-1
    sum = sum + ((angular_vel(i)+angular_vel(i+1))/2);
    angle(i) = -sum/100;
end    