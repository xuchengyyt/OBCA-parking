%经过Hybrid A Star 搜索的轨迹最终通过H-OBCA进行优化以及避障；
%% 侧方位泊车
clc
close all
clear
addpath('D:\Program Files\Polyspace\R2020a\bin\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*
opti = casadi.Opti();
%% get the reference
load('hybrid_tra2.mat')
anchorPoint = getAnchorPoint(Deriction);
point1 = [Xref(anchorPoint(1));Yref(anchorPoint(1));Thref(anchorPoint(1));0];
point2 = [Xref(anchorPoint(2));Yref(anchorPoint(2));Thref(anchorPoint(2));0];

    
%% vehicle
Veh = struct;
Veh.LF = 3.8;
Veh.LB = 1.0;
Veh.L = Veh.LF + Veh.LB;
Veh.W = 1.8;
Veh.wheel_base = 2.9;
%% parameter
Ps = 1000;
t = 0.1; % sampling time 10ms or 20ms
N = length(Xref) - 1; %预测时域；第一个位置就是pos，所以要减一
Ncbf = N;

v_max = 3;
v_min = -3;
acc_max = 2;
acc_min = -2;
steer_max = 0.6;
steer_min = -0.6;
pos = [-3.9;5;0;0]; % the initial state;
target = [1.2;1.2;0;0];
control_stash = [];
pos_stash = [];
dist_stash = [];
open_loop_stash = cell(100,1);
d_min = 0.1;% the min dist to obs
%% vehicle and obs
% rear-axis-rep
ego_A = [1,0;0,1;-1,0;0,-1]; %4*2  % ego_a is constant but ego_b will change 
ego_b = [Veh.LF;Veh.W/2;Veh.LB ;Veh.W/2]; %[-left;-down;right;up]; %[-left;-down;right;up]
offset = (ego_b(1) + ego_b(3)) / 2 - ego_b(3);
% obstacle
obs(1,:) = [7.5,1.2];obs(2,:) = [-4.5,1.2,]; % 平行泊车
num = 2;
method  = 1; % 1 平行泊车 2 垂直泊车
mat_A = [];mat_b = [];
for i = 1 : num
    ret = get_obstacle(Veh,[obs(i,1);obs(i,2);0]);
    mat_A = [mat_A;ret{1}];
    mat_b = [mat_b;ret{2}];
end
constant_v = 1;




%% state function
x = SX.sym('x');
y = SX.sym('y');
yaw = SX.sym('yaw');
vel = SX.sym('vel');

states = [x;y;yaw;vel];
acc = SX.sym('acc');
steer = SX.sym('steer');
controls = [steer;acc];
rhs = [vel * cos(yaw);vel*sin(yaw);vel*tan(steer)/Veh.wheel_base;acc];
f = Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u) % 后轴中心
n_states = length(states);
n_controls = length(controls);



Q = zeros(4,4); Q(1,1) = 1000;Q(2,2) = 1000;Q(3,3) = 1000;Q(4,4) =  10; % weighing matrices (states)
R = zeros(2,2); R(1,1) = 10; R(2,2) = 10; % weighing matrices (controls)
R1 = zeros(2,2);R1(1,1) = 100;R1(2,2) = 100; % smooth matrices

u0 = zeros(2,N);
x0 = zeros(n_states,N+1);

%% start simulation
mpciter = 1;

pos_stash = [pos_stash;pos'];%存储第一次的位置

P = zeros(n_states,N+1);
P(:,1) = pos;
    
for i = 2:(N+1)
    P(1,i) = Xref(i);
    P(2,i) = Yref(i);
    P(3,i) = Thref(i);
    if (i == anchorPoint(1) || i == anchorPoint(2) )
         P(4,i) = 0;
    elseif (Deriction(i)>0)
          P(4,i) =constant_v;
    else
          P(4,i) =-constant_v;
    end
end
X = opti.variable(4,N+1);opti.set_initial(X,P)
U = opti.variable(2,N);opti.set_initial(U,u0)
T = opti.variable(1,N);opti.set_initial(T,t)
% initial constrints
opti.subject_to(X(:,1)-pos == 0);
% terminal constraints
opti.subject_to(0.1 >= X(1,end) - target(1) >= -0.1);
opti.subject_to(0.1 >= X(2,end) - target(2) >= -0.1);
opti.subject_to(0.01 >= X(3,end) - target(3) >= -0.01);
opti.subject_to(X(4,end) - target(4) == 0);
opti.subject_to(U(:,end) == 0);
% anchor point constraints
% opti.subject_to(X(:,anchorPoint(1)) == point1);
% opti.subject_to(X(:,anchorPoint(2)) == point2);
% kinematic constraints
for i = 1:N
    opti.subject_to(X(:,i+1) - T(i) * f(X(:,i),U(:,i)) - X(:,i) == 0);
    opti.subject_to(t *0.7 <= T(i) <= t*1.3);
end
% variable constraints
for i =1:N
    opti.subject_to(steer_max>=U(1,i)>=steer_min)
    opti.subject_to(acc_max>=U(2,i)>=acc_min)
end
% state constarints
for i = 1:N
%         opti.subject_to(x_max>=X(1,i+1)>=x_min);
%         opti.subject_to(y_max>=X(2,i+1)>=y_min);
%         opti.subject_to(yaw_max>=X(3,i+1)>=yaw_min);
    opti.subject_to(v_max>=X(4,i+1)>=v_min);
end
    % get the state of vehicle
matrix_vehicle = get_vehicle_matrix(pos,ego_A,ego_b);


% s = opti.variable(Ncbf,1);opti.set_initial(s,zeros(Ncbf,1));
 lamb_ = opti.variable(length(ego_A) * 2,Ncbf);
 mu_ = opti.variable(length(ego_A) * 2,Ncbf);
 [lamb_ini,mu_ini,obj_val] = DualMultWS(ego_A,ego_b,obs,Xref',Yref',Thref',N,num);
 opti.set_initial(lamb_,lamb_ini);
 opti.set_initial(mu_,mu_ini);
 obj_val
% geo-rep
g = [Veh.L/2;Veh.W/2;Veh.L/2;Veh.W/2];
G = [1,0;0,1;-1,0;0,-1]; 
for i = 1:Ncbf
    for j = 1:num
        %% dual obs constraints
        obs_A = mat_A((4*(j-1) + 1):4*j,:);
        obs_B = mat_b((4*(j-1) + 1):4*j,:);
        lamb = lamb_((4*(j-1) + 1):4*j,:);
        mu = mu_((4*(j-1) + 1):4*j,:);
        opti.subject_to(lamb(:,i)>=0);
        opti.subject_to(mu(:,i)>=0);
        % rear axle center to geo center
        % geo_center(1) = X(1,i+1) + offset * cos(X(3,i+1));
        % geo_center(2) = X(2,i+1) + offset * sin(X(3,i+1));
        % -g'*mu + (A*t - b)*lambda > 0 
        opti.subject_to(-g'* mu(:,i) + (obs_A * [X(1,i+1) + offset * cos(X(3,i+1));X(2,i+1) + offset * sin(X(3,i+1))] - obs_B)' * lamb(:,i)> d_min ); % without slack variable
        % G'*mu + R'*A*lambda = 0
        opti.subject_to(G'* mu(:,i) + rotate(X(:,i+1))' * obs_A' * lamb(:,i) == 0)      % 2 * 4  *   4 * 1 + 2 * 2 * 2 * 4 * 4 * 1 = 0
        % norm(A'*lambda) <= 1
        temp = obs_A'*lamb(:,i); %  （2*4） * （4*1）
        opti.subject_to(temp' * temp <=1)
    end
end
    
    % cost function
obj = 0; % Objective function
for k = 1:N
obj = obj+(X(:,k+1)-P(:,k+1))'*Q*(X(:,k+1)-P(:,k+1)) + U(:,k)'*R*U(:,k); % calculate obj
end
% smooth function
for k = 1:N-1
obj = obj +(U(:,k+1) - U(:,k))' * R1 * (U(:,k+1) - U(:,k));
end
% time optimal
for k = 1: N
    obj = obj + T(i) * 100 + T(i) * T(i) * 100;
end
% slack penalty
% for k = 1:Ncbf
%     obj = obj + (s(i) -1)^2* Ps ;
% end
    
    % solve optimization
%     n = 20;
%     T2 = zeros(1,n);

tic
opti.minimize(obj)
opts = struct;
opts.ipopt.max_iter = 1000;
opts.ipopt.print_level =1;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-3;
opts.ipopt.acceptable_obj_change_tol = 1e-2;
opti.solver('ipopt',opts);
opt_sol = opti.solve();
u0 = opt_sol.value(U);
x0 = opt_sol.value(X);
% s = opt_sol.value(s);
open_loop = x0;
open_loop_stash{mpciter} = open_loop;
pos = open_loop(:,2);%位置更新
pos_stash = [pos_stash;pos'];
toc
WS2 = {x0,u0};
%save('WS2.mat','WS2');

% save('T2.mat','T2');
mpciter
% get the time
t1 = zeros(1,N);
for i = 2:N
    t1(i) = t1(i-1) + 0.1;
end
%% result show
if 1 
figure(1)
% title('基于TDR-OBCA的泊车轨迹优化')
set(gca,'FontSize',16)
hold on 
plot(x0(1,:),x0(2,:),'r-','LineWidth',1);
plot(Xref,Yref,'b--','LineWidth',1)


xlim([-8,14])
ylim([0,7])
line([-8,14],[7,7],'Color','b','Linewidth',4)
line([-8,14],[0,0],'Color','b','Linewidth',4)
line([-6,12],[2.5,2.5],'Color','k','Linewidth',2,'linestyle','--')
line([-6,-6],[0,2.5],'Color','k','Linewidth',2,'linestyle','--')
line([0,0],[0,2.5],'Color','k','Linewidth',2,'linestyle','--')
line([6,6],[0,2.5],'Color','k','Linewidth',2,'linestyle','--')

line([12,12],[0,2.5],'Color','k','Linewidth',2,'linestyle','--')
for i = 1: length(x0)
   [vehx,vehy] = vehicleMotion(x0(1,i),x0(2,i),x0(3,i),Veh); % 根据后轴中心的位姿计算车辆边框的位姿
   plot(vehx,vehy,'Color','g','Linewidth',1); % 车辆边框
end
[vehx,vehy] = vehicleMotion(x0(1,1),x0(2,1),x0(3,1),Veh); % 根据后轴中心的位姿计算车辆边框的位姿
plot(vehx,vehy,'k'); % 车辆边框
[vehx,vehy] = vehicleMotion(x0(1,end),x0(2,end),x0(3,end),Veh); % 根据后轴中心的位姿计算车辆边框的位姿
plot(vehx,vehy,'k'); % 车辆边框
% obs plot

[vehx,vehy] = vehicleMotion(obs(1,1),obs(1,2),0,Veh); % 根据后轴中心的位姿计算车辆边框的位姿
plot(vehx,vehy,'k'); % 车辆边框
[vehx,vehy] = vehicleMotion(obs(2,1),obs(2,2),0,Veh); % 根据后轴中心的位姿计算车辆边框的位姿
plot(vehx,vehy,'k'); % 车辆边框
set(gcf,'unit','centimeters','position',[1,2,22*1.3,7*1.3])
legend('优化轨迹','参考路径')

figure(2)
plot(t1,u0(1,:),'Linewidth',2);
xlabel('时间（s）','FontSize',16)
ylabel('前轮转角（rad）','FontSize',16)
set(gca,'FontSize',20)
set(gcf,'unit','centimeters','position',[1,2,20,7])

figure(3)
plot(t1,u0(2,:),'Linewidth',2);
xlabel('时间（s）','FontSize',16)
ylabel('加速度（m/s^2）','FontSize',16)
set(gca,'FontSize',20)
set(gcf,'unit','centimeters','position',[1,2,20,7])

figure(4)
plot([t1';0.1 * (N+1)]',x0(4,:),'Linewidth',2)
xlabel('时间（s）','FontSize',16)
ylabel('速度（m/s）','FontSize',16)
set(gca,'FontSize',20)
set(gcf,'unit','centimeters','position',[1,2,20,7])
end


function distance = dis_pos(pos,target_pos)
distance = sqrt((pos(1) - target_pos(1))^2 + (pos(2) - target_pos(2))^2);
end

function matrix_a_b = get_vehicle_matrix(state,mat_A,vec_b)
R = rotate(state);T = translation(state);
matrix_A = mat_A * R';
matrix_B = mat_A*R'* T + vec_b;
matrix_a_b = {matrix_A,matrix_B};
end
function heading = rotate(state)
    heading = [cos(state(3)),-sin(state(3));
               sin(state(3)),cos(state(3))];
end
function trans = translation(state)
    trans = [state(1);state(2)];
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
function dist_dual_variable =  get_dist_region_to_region(mat_A, vec_a,mat_B,vec_b)
    % A 是障碍物，b是自车
    % Return distance between a point and a convex region
    opti = casadi.Opti(); 
    % variables and cost
    point_in_region1 = opti.variable(length(mat_A(1,:)), 1);
    point_in_region2 = opti.variable(length(mat_B(1,:)), 1);
    cost = 0;
    % constraints
    con1 = mat_A * point_in_region1 <= vec_a;
    con2 = mat_B * point_in_region2 <= vec_b;
    opti.subject_to(con1)
    opti.subject_to(con2)
    dist_vec = point_in_region1 - point_in_region2;
    cost  = cost +  dist_vec' * dist_vec;
    % solve optimization
    opti.minimize(cost)
    p_opts = struct('expand',true);
    s_opts = struct('max_iter',300);
    opti.solver('ipopt',p_opts,s_opts);
    opt_sol = opti.solve();
    % minimum distance & dual variables
    dist = opt_sol.value(cost);
    dist = sqrt(dist);
    
    if dist > 0
        lamb = opt_sol.value(opti.dual(con1)) / (2 * dist);
        mu = opt_sol.value(opti.dual(con2)) / (2 * dist);
    else
        lamb = zeros(length(mat_A),1);
        mu = zeros(length(mat_B),1);
    end
    dist_dual_variable = {dist,lamb,mu};
    
end

