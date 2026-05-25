function [U2, U3, U4, l_hat_att_dot, att_state_out] = UAV_attitude_controller_AFBC(state, desired_angles, l_hat_att, att_state_in, dt, params)
% UAV_attitude_controller_AFBC - 无模型平滑正则化自适应姿态环 (紧凑标量版)

% --- 1. 状态分解 ---
phi = state(4); theta = state(5); psi = state(6);
p = state(10); q = state(11); r = state(12);
Omega = [p; q; r];          
Theta = [phi; theta; psi];  

% ==========================================================
% 【核心替换】：剥离真实转动惯量矩阵，使用紧凑虚拟惯量标量
% ==========================================================
J_uav_virtual = 0.065; % 略大于最大的偏航轴真实惯量 0.06

K3 = params.K3_att; K4 = params.K4_att;
K_frac = params.K_frac_att; beta = params.beta_att;
Gamma = params.Gamma_att; sigma = params.sigma_att;

% --- 2. 第一级 跟踪微分器 ---
R1 = 25; 
v1_theta = att_state_in.prev_Theta_d;
v2_theta = att_state_in.Theta_d_dot;

err_theta = v1_theta - desired_angles;
err_theta(3) = wrapToPi(err_theta(3));

v1_theta_dot = v2_theta;
v2_theta_dot = -R1^2 * err_theta - 2 * R1 * v2_theta;

Theta_d_smooth = v1_theta + v1_theta_dot * dt;
Theta_d_smooth(3) = wrapToPi(Theta_d_smooth(3));
Theta_d_dot = v2_theta + v2_theta_dot * dt;

% --- 3. 运动学环 ---
z3 = Theta - Theta_d_smooth;
z3(3) = wrapToPi(z3(3));

W_inv = [1, 0, -sin(theta);
         0, cos(phi), sin(phi)*cos(theta);
         0, -sin(phi), cos(phi)*cos(theta)];

Omega_d_target = W_inv * (-K3 * z3 + Theta_d_dot);

% --- 4. 第二级 跟踪微分器 ---
R2 = 30;
v1_omega = att_state_in.Omega_d;
v2_omega = att_state_in.Omega_d_dot;

err_omega = v1_omega - Omega_d_target;
v1_omega_dot = v2_omega;
v2_omega_dot = -R2^2 * err_omega - 2 * R2 * v2_omega;

Omega_d_smooth = v1_omega + v1_omega_dot * dt;
Omega_d_dot = v2_omega + v2_omega_dot * dt;

% --- 5. 动力学环 ---
z4 = Omega - Omega_d_smooth;

norm_z4 = norm(z4);
l_hat_att_dot = diag(Gamma) * norm_z4 - sigma .* l_hat_att;
l_hat_att_dot = max(0, l_hat_att_dot); 

% 正则化处理
epsilon = 0.02; 
reg_term_z4 = reg_sig(z4, beta, epsilon);

% ==========================================================
% 纯无模型架构计算：利用标量乘法替代原有的矩阵前馈
% ==========================================================
tau_robust = - l_hat_att .* tanh(z4 / 0.05); 
tau = J_uav_virtual * (Omega_d_dot - K4 * z4 - K_frac * reg_term_z4) + tau_robust;

U2 = max(-5, min(tau(1), 5));
U3 = max(-5, min(tau(2), 5));
U4 = max(-2, min(tau(3), 2));

% --- 6. 状态缓存更新 ---
att_state_out.prev_Theta_d = Theta_d_smooth;
att_state_out.Theta_d_dot = Theta_d_dot;
att_state_out.Omega_d = Omega_d_smooth;
att_state_out.Omega_d_dot = Omega_d_dot;
end

function y = reg_sig(x, alpha, epsilon)
    % 修复后的严格连续的正则化函数
    y = zeros(size(x));
    for i = 1:length(x)
        if abs(x(i)) >= epsilon
            y(i) = (abs(x(i))^alpha) * sign(x(i)); % 外部使用标准的符号函数，保证理论收敛
        else
            y(i) = x(i) * (epsilon^(alpha - 1));   % 内部通过线性斜率完美连接边界，消除奇点断层
        end
    end
end