function state_dot = UAV_dynamic_function(state, control_inputs)
% UAV_dynamic_function - 修复版UAV动力学模型
% 确保正确的坐标系和力的方向

% 参数定义
g = 9.81; % 重力加速度
m = 1.5;  % 质量
Ixx = 0.03; Iyy = 0.03; Izz = 0.06; % 转动惯量

% 状态分解
pos = state(1:3);      % [x, y, z]
angles = state(4:6);   % [phi, theta, psi]
vel_body = state(7:9); % [u, v, w] - 机体坐标系速度
rates = state(10:12);  % [p, q, r] - 角速度

% 控制输入
U1 = control_inputs(1); % 总推力 (向上为正)
U2 = control_inputs(2); % 滚转力矩
U3 = control_inputs(3); % 俯仰力矩  
U4 = control_inputs(4); % 偏航力矩

% 姿态角
phi = angles(1); theta = angles(2); psi = angles(3);
p = rates(1); q = rates(2); r = rates(3);

% === 关键修复：确保正确的坐标系转换 ===
% 机体到惯性系的旋转矩阵
R_body_to_inertial = [cos(psi)*cos(theta), cos(psi)*sin(theta)*sin(phi)-sin(psi)*cos(phi), cos(psi)*sin(theta)*cos(phi)+sin(psi)*sin(phi);
                      sin(psi)*cos(theta), sin(psi)*sin(theta)*sin(phi)+cos(psi)*cos(phi), sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi);
                      -sin(theta),         cos(theta)*sin(phi),                            cos(theta)*cos(phi)];

% 位置导数 = 惯性系速度
vel_inertial = R_body_to_inertial * vel_body;
pos_dot = vel_inertial;

% 姿态导数
angles_dot = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);
              0, cos(phi),            -sin(phi);
              0, sin(phi)/cos(theta), cos(phi)/cos(theta)] * rates;

% === 关键修复：力的方向和重力补偿 ===
% 机体坐标系下的力 (推力沿机体z轴向上)
F_body = [0; 0; U1]; % U1应该是正值向上推力

% 转换到惯性坐标系
F_inertial = R_body_to_inertial * F_body;

% 重力在惯性系下 (向下为负)
F_gravity = [0; 0; -m*g];

% 总的惯性力
F_total_inertial = F_inertial + F_gravity;

% 转换回机体坐标系进行动力学计算
F_total_body = R_body_to_inertial' * F_total_inertial;

% 机体坐标系下的速度导数 (考虑科氏力)
vel_body_dot = F_total_body/m - cross(rates, vel_body);

% 角加速度 (转动动力学)
rates_dot = [U2/Ixx; U3/Iyy; U4/Izz] - [((Iyy-Izz)/Ixx)*q*r; ((Izz-Ixx)/Iyy)*p*r; ((Ixx-Iyy)/Izz)*p*q];

% 组装状态导数
state_dot = [pos_dot; angles_dot; vel_body_dot; rates_dot];

end