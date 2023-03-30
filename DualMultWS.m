function [lambda,mul,obj_val] = DualMultWS(ego_A,ego_b,obs,xref,yref,thref,N,nobs)
addpath('D:\Program Files\Polyspace\R2020a\bin\casadi-windows-matlabR2016a-v3.5.5')
Veh = struct;
Veh.LF = 3.8;
Veh.LB = 1.0;
Veh.L = Veh.LF + Veh.LB;
Veh.W = 1.8;
Veh.wheel_base = 2.9;
opti = casadi.Opti();
X = zeros(3,N+1);
X(1,:) = xref;
X(2,:) = yref;
X(3,:) = thref;
lamb_ = opti.variable(8,N);
mu_ = opti.variable(8,N);
mat_A = [];mat_b = [];
num = nobs;
for i = 1 : num
    ret = get_obstacle(Veh,[obs(i,1);obs(i,2);0]);
    mat_A = [mat_A;ret{1}];
    mat_b = [mat_b;ret{2}];
end
dmin = 0;
offset = (ego_b(1) + ego_b(3)) / 2 - ego_b(3);
obj = 0; % Objective function


for i = 1:N
    for j = 1:num
        %% dual obs constraints
        obs_A = mat_A((4*(j-1) + 1):4*j,:);
        obs_B = mat_b((4*(j-1) + 1):4*j,:);
        lamb = lamb_((4*(j-1) + 1):4*j,:);
        mu = mu_((4*(j-1) + 1):4*j,:);
        opti.subject_to(lamb(:,i)>=0);
        opti.subject_to(mu(:,i)>=0);
        % -g'*mu + (A*t - b)*lambda > 0 
        % rear axle center to geo center
        geo_center = X(1:2,i+1); 
        geo_center(1) = X(1,i+1) + offset * cos(X(3,i+1));
        geo_center(2) = X(2,i+1) + offset * sin(X(3,i+1));
        opti.subject_to(-ego_b'* mu(:,i) + (obs_A * geo_center - obs_B)' * lamb(:,i)> dmin); % without slack variable
        % G'*mu + R'*A*lambda = 0
        opti.subject_to(ego_A'* mu(:,i) + rotate(X(:,i+1))' * obs_A' * lamb(:,i) == 0)      % 2 * 4  *   4 * 1 + 2 * 2 * 2 * 4 * 4 * 1 = 0
        % norm(A'*lambda) <= 1
        temp = obs_A'*lamb(:,i); %  （2*4） * （4*1）
        opti.subject_to(temp' * temp <=1)
        
        obj = obj + -ego_b'* mu(:,i) + (obs_A * geo_center - obs_B)' * lamb(:,i);
    end
end
obj = -obj;
opti.minimize(obj)
opts = struct;
opts.ipopt.max_iter = 100000;
opts.ipopt.print_level =1;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-5;
opts.ipopt.acceptable_obj_change_tol = 1e-5;
opti.solver('ipopt',opts);
opt_sol = opti.solve();
lambda = opt_sol.value(lamb_);
mul = opt_sol.value(mu_);
obj_val = opt_sol.value(obj);
end
function heading = rotate(state)
    heading = [cos(state(3)),-sin(state(3));
               sin(state(3)),cos(state(3))];
end

