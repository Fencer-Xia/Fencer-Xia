function [U2, U3, U4, l_hat_att_dot, att_state_out] = UAV_attitude_controller_AFBC(state, desired_angles, l_hat_att, att_state_in, dt, params)
% UAV_attitude_controller_AFBC - 无人机姿态环自适应分数阶反步控制器
% 针对小惯量模型引入了防奇异正则化、低通微分滤波与sigma-修正

% --- 1. 状态分解 ---
phi = state(4); theta = state(5); psi = state(6);
p = state(10); q = state(11); r = state(12);
Omega = [p; q; r];          % 实际机体角速度
Theta = [phi; theta; psi];  % 实际姿态角

% 解析参数
Ixx = params.Ixx; Iyy = params.Iyy; Izz = params.Izz;
J = diag([Ixx, Iyy, Izz]);
K3 = params.K3_att;         % 运动学环反馈增益
K4 = params.K4_att;         % 动力学环反馈增益
K_frac = params.K_frac_att; % 分数阶增益
beta = params.beta_att;
Gamma = params.Gamma_att;
sigma = params.sigma_att;   % 鲁棒泄漏项系数（防爆关键）

% 读取历史状态
prev_Theta_d = att_state_in.Theta_d;

% --- 2. 期望姿态导数计算 (低通滤波器防高频噪声爆炸) ---
tau_f = 0.05; % 滤波时间常数
raw_Theta_d_dot = (desired_angles - prev_Theta_d) / dt;
% 处理偏航角跨越 pi/-pi 时的数值跳变
raw_Theta_d_dot(3) = atan2(sin(raw_Theta_d_dot(3)*dt), cos(raw_Theta_d_dot(3)*dt)) / dt;
% 一阶低通滤波
Theta_d_dot = (dt / (tau_f + dt)) * raw_Theta_d_dot + (tau_f / (tau_f + dt)) * att_state_in.Theta_d_dot;

% --- 3. 角度误差 (z3) ---
z3 = Theta - desired_angles;
z3(3) = atan2(sin(z3(3)), cos(z3(3))); % 偏航角归一化

% 姿态运动学逆矩阵 (欧拉角速率转机体角速度)
W_inv = [1, 0, -sin(theta);
         0, cos(phi), sin(phi)*cos(theta);
         0, -sin(phi), cos(phi)*cos(theta)];

% --- 4. 虚拟控制量 (期望机体角速度 Omega_d) ---
Omega_d = W_inv * (-K3 * z3 + Theta_d_dot);

% --- 5. 期望角速度导数计算 (滤波) ---
raw_Omega_d_dot = (Omega_d - att_state_in.Omega_d) / dt;
Omega_d_dot = (dt / (tau_f + dt)) * raw_Omega_d_dot + (tau_f / (tau_f + dt)) * att_state_in.Omega_d_dot;

% --- 6. 角速度误差 (z4) ---
z4 = Omega - Omega_d;

% --- 7. AFBC 控制力矩计算 ---
% (A) 陀螺力矩交叉项 (模型非线性补偿)
gyro_torque = cross(Omega, J * Omega);

% (B) 分数阶项 (加入 epsilon=1e-3 防奇异)
epsilon = 1e-3;
frac_term = K_frac * diag((abs(z4) + epsilon).^(2*beta - 1)) * sign(z4);

% (C) 自适应鲁棒项 (tanh代替sign消除高频抖振)
robust_term = diag(l_hat_att) * tanh(z4 / 0.05);

% (D) 物理力矩合成
tau = gyro_torque + J * Omega_d_dot - K4 * z4 - robust_term - frac_term;

U2 = tau(1);
U3 = tau(2);
U4 = tau(3);

% --- 8. 内环自适应律更新 ---
% 使用带 sigma-修正 的更新律，防止在存在持续外扰时自适应参数发散到无穷大
l_hat_att_dot = Gamma * abs(z4) - sigma .* l_hat_att;

% --- 9. 更新状态字典供下一步使用 ---
att_state_out.Theta_d = desired_angles;
att_state_out.Theta_d_dot = Theta_d_dot;
att_state_out.Omega_d = Omega_d;
att_state_out.Omega_d_dot = Omega_d_dot;
end