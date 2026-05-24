function [U, ctrl_state] = UAV_controller_BSC_DSC(state, desired_state, ctrl_state, params, dt)
% 传统反步法 + 动态面控制 (UAV) - 复数爆炸修复版

    pos = state(1:3); angles = state(4:6); vel = state(7:9); omega = state(10:12);
    m = params.m; g = params.g; J = diag([params.Ixx, params.Iyy, params.Izz]);
    phi = angles(1); theta = angles(2); psi = angles(3);
    R_b2i = [cos(psi)*cos(theta), cos(psi)*sin(theta)*sin(phi)-sin(psi)*cos(phi), cos(psi)*sin(theta)*cos(phi)+sin(psi)*sin(phi);
             sin(psi)*cos(theta), sin(psi)*sin(theta)*sin(phi)+cos(psi)*cos(phi), sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi);
             -sin(theta),         cos(theta)*sin(phi),                            cos(theta)*cos(phi)];
    vel_i = R_b2i * vel;

    % --- 反步与滤波器参数 ---
    % --- 反步与滤波器参数 (修复版：温和基准参数) ---
    % 1. 调低位置环增益，防止起步索要过大的姿态倾角
    K1_p = diag([1.0, 1.0, 1.5]);   K2_p = diag([2.0, 2.0, 3.0]);
    
    % 2. 调低姿态环增益，降低机动刚度
    K1_a = diag([4.0, 4.0, 4.0]);   K2_a = diag([3.0, 3.0, 3.0]);
    
    % 3. 【核心修复】增大 DSC 滤波时间常数
    % dt=0.01，tau设为0.15意味着有15个采样点来进行平滑过渡，彻底抹平微分尖峰
    tau_p = 0.15; 
    tau_a = 0.15; 
    tau_Theta = 0.15;

    % ================= 位置环 =================
    z1_p = pos - desired_state.pos;
    
    % 【修复2】对送入线性控制器的位置误差进行安全限幅 (防止起步巨大冲激)
    max_err = 3.0;
    z1_p_sat = max(min(z1_p, max_err), -max_err);
    
    alpha_p = desired_state.vel - K1_p * z1_p_sat;
    
    alpha_p_f = ctrl_state.alpha_p_f;
    dot_alpha_p_f = (alpha_p - alpha_p_f) / tau_p;
    ctrl_state.alpha_p_f = alpha_p_f + dot_alpha_p_f * dt;
    
    z2_p = vel_i - alpha_p_f;
    F_d = m * dot_alpha_p_f - z1_p - K2_p * z2_p;
    F_total = F_d + [0; 0; m*g];
    
    % 【核心修复 1】严格保障力学推力限幅与单位几何向量的代数一致性
    if norm(F_total) > params.max_thrust
        F_total = (F_total / norm(F_total)) * params.max_thrust;
    end
    U1 = norm(F_total);
    if U1 < 1e-6
        z_b = [0; 0; 1];
    else
        z_b = F_total / U1; % 现在绝对保证是单位向量
    end
    
    psi_d = atan2(desired_state.vel(2), desired_state.vel(1));
    phi_d = asin(max(-0.95, min(0.95, z_b(1)*sin(psi_d) - z_b(2)*cos(psi_d))));
    theta_d = atan2(z_b(1)*cos(psi_d) + z_b(2)*sin(psi_d), z_b(3));
    desired_angles = [phi_d; max(-0.95, min(0.95, theta_d)); psi_d];

    % ================= 姿态环 =================
    Theta_d_f = ctrl_state.prev_Theta_d;
    dot_Theta_d_f = (desired_angles - Theta_d_f) / tau_Theta;
    dot_Theta_d_f(3) = atan2(sin(dot_Theta_d_f(3)*tau_Theta), cos(dot_Theta_d_f(3)*tau_Theta)) / tau_Theta;
    ctrl_state.prev_Theta_d = Theta_d_f + dot_Theta_d_f * dt;
    
    z1_a = angles - Theta_d_f;
    z1_a(3) = wrapToPi(z1_a(3));
    
    W_inv = [1, 0, -sin(theta); 0, cos(phi), sin(phi)*cos(theta); 0, -sin(phi), cos(phi)*cos(theta)];
    alpha_a = W_inv * (dot_Theta_d_f - K1_a * z1_a);
    
    alpha_a_f = ctrl_state.alpha_a_f;
    dot_alpha_a_f = (alpha_a - alpha_a_f) / tau_a;
    ctrl_state.alpha_a_f = alpha_a_f + dot_alpha_a_f * dt;
    
    z2_a = omega - alpha_a_f;
    tau_cmd = cross(omega, J * omega) + J * dot_alpha_a_f - W_inv' * z1_a - K2_a * z2_a;
    
    U = [U1; tau_cmd(1); tau_cmd(2); tau_cmd(3)];
end