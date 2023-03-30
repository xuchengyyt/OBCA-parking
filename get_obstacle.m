function [ret] = get_obstacle(Veh,state)
length = Veh.LB+Veh.LF;
rear_dist = length/2 - Veh.LB;
vertices = get_vehicle_vertices(length,Veh.W,rear_dist,state);
left = min(vertices(:,1));right = max(vertices(:,1));
down = min(vertices(:,2));up = max(vertices(:,2));

% A1 = [-1,0;0,-1;1,0;0,1];
% b = [-left;-down;right;up]; OBCA changed in 2022/3/20
A1 = [1,0;0,1; -1,0;0,-1];
b = [right;up;-left;-down];

ret = {A1;b};
end
function vertices = get_vehicle_vertices(length,width,rear_dist,state)
x = state(1);y = state(2);theta = state(3);
xc = x + (rear_dist) * cos(theta); % 后轴中心到几何中心
yc = y + rear_dist * sin(theta);
vertices = [xc + length / 2 * cos(theta) - width / 2 * sin(theta),yc + length / 2 * sin(theta) + width / 2 * cos(theta);
            xc + length / 2 * cos(theta) + width / 2 * sin(theta),yc + length / 2 * sin(theta) - width / 2 * cos(theta);
            xc - length / 2 * cos(theta) + width / 2 * sin(theta),yc - length / 2 * sin(theta) - width / 2 * cos(theta);
            xc - length / 2 * cos(theta) - width / 2 * sin(theta),yc - length / 2 * sin(theta) + width / 2 * cos(theta)];
end