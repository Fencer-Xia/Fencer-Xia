function [U, ctrl_state] = UAV_controller_PID(state, desired_state, ctrl_state, params, dt)
% 带抗积分饱和的 PID 控制器 (UAV) - 最终临界阻尼版

    pos = state(1:3); angles = state(4:6); vel = state(7:9); omega = state(10:12);
    m = params.m; g = params.g; % 注意：这里不再需要 J 矩阵去缩放力矩了
    phi = angles(1); theta = angles(2); psi = angles(3);
    
    R_b2i = [cos(psi)*cos(theta), cos(psi)*sin(theta)*sin(phi)-sin(psi)*cos(phi), cos(psi)*sin(theta)*cos(phi)+sin(psi)*sin(phi);
             sin(psi)*cos(theta), sin(psi)*sin(theta)*sin(phi)+cos(psi)*cos(phi), sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi);
             -sin(theta),         cos(theta)*sin(phi),                            cos(theta)*cos(phi)];
    vel_i = R_b2i * vel;

    % ==========================================================
    % 【核心修复】：PID 增益重构
    % 1. 位置环：提供稳定、不超调的加速度引导
    Kp_p = diag([1.2, 1.2, 2.0]); 
    Ki_p = diag([0.1, 0.1, 0.2]); 
    Kd_p = diag([2.0, 2.0, 3.0]);
    
    % 2. 姿态环：使用直接扭矩增益，匹配临界阻尼 (Zeta=1.0)
    Kp_a = diag([8.0, 8.0, 5.0]);   % 每偏离1弧度，输出8N.m的纠偏力矩
    Ki_a = diag([0.5, 0.5, 0.5]); 
    Kd_a = diag([1.2, 1.2, 1.0]);   % 每1rad/s的角速度误差，提供1.2N.m的刹车阻尼
    % ==========================================================

    max_int_p = [2.0; 2.0; 2.0]; 
    max_int_a = [0.5; 0.5; 0.5];

    % ================= 位置 PID =================
    e_p = desired_state.pos - pos;
    max_e_p = 3.0;
    e_p_sat = max(min(e_p, max_e_p), -max_e_p); % 误差截断
    
    e_v = desired_state.vel - vel_i;
    
    ctrl_state.int_p = ctrl_state.int_p + e_p_sat * dt;
    ctrl_state.int_p = max(min(ctrl_state.int_p, max_int_p), -max_int_p); 
    
    a_cmd = Kp_p * e_p_sat + Ki_p * ctrl_state.int_p + Kd_p * e_v;
    a_cmd(1) = max(min(a_cmd(1), 6.0), -6.0); % 加速度限幅，防止大倾角
    a_cmd(2) = max(min(a_cmd(2), 6.0), -6.0);
    
    F_total = m * (a_cmd + [0; 0; g]);
    
    if norm(F_total) > params.max_thrust
        F_total = (F_total / norm(F_total)) * params.max_thrust;
    end
    U1 = norm(F_total);
    if U1 < 1e-6; z_b = [0;0;1]; else; z_b = F_total / U1; end
    
    psi_d = atan2(desired_state.vel(2), desired_state.vel(1));
    phi_d = asin(max(-0.95, min(0.95, z_b(1)*sin(psi_d) - z_b(2)*cos(psi_d))));
    theta_d = atan2(z_b(1)*cos(psi_d) + z_b(2)*sin(psi_d), z_b(3));
    desired_angles = [phi_d; max(-0.95, min(0.95, theta_d)); psi_d];

    % ================= 姿态 PID =================
    tau_Theta = 0.08; % 稍微增加平滑度
    Theta_d_f = ctrl_state.prev_Theta_d;
    dot_Theta_d_f = (desired_angles - Theta_d_f) / tau_Theta;
    dot_Theta_d_f(3) = atan2(sin(dot_Theta_d_f(3)*tau_Theta), cos(dot_Theta_d_f(3)*tau_Theta)) / tau_Theta;
    
    next_Theta_d_f = Theta_d_f + dot_Theta_d_f * dt;
    next_Theta_d_f(3) = wrapToPi(next_Theta_d_f(3)); 
    ctrl_state.prev_Theta_d = next_Theta_d_f;

    e_a = Theta_d_f - angles;
    e_a(3) = wrapToPi(e_a(3));
    
    ctrl_state.int_a = ctrl_state.int_a + e_a * dt;
    ctrl_state.int_a = max(min(ctrl_state.int_a, max_int_a), -max_int_a);
    
    W_inv = [1, 0, -sin(theta);
             0, cos(phi), sin(phi)*cos(theta);
             0, -sin(phi), cos(phi)*cos(theta)];
             
    e_a_body = W_inv * e_a;
    int_a_body = W_inv * ctrl_state.int_a;
    e_omega_body = W_inv * dot_Theta_d_f - omega; 
    
    % 【绝杀修复】：移除 J 矩阵，使用计算好的纯扭矩增益输出
    tau_cmd = Kp_a * e_a_body + Ki_a * int_a_body + Kd_a * e_omega_body;
    
    U = [U1; tau_cmd(1); tau_cmd(2); tau_cmd(3)];
end