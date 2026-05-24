function [U, ctrl_state] = UAV_controller_BSC(state, desired_state, ctrl_state, params, dt)
% 传统反步法 (UAV) - 移除DSC动态面，使用数值差分计算导数

    pos = state(1:3); angles = state(4:6); vel = state(7:9); omega = state(10:12);
    m = params.m; g = params.g; J = diag([params.Ixx, params.Iyy, params.Izz]);
    phi = angles(1); theta = angles(2); psi = angles(3);
    R_b2i = [cos(psi)*cos(theta), cos(psi)*sin(theta)*sin(phi)-sin(psi)*cos(phi), cos(psi)*sin(theta)*cos(phi)+sin(psi)*sin(phi);
             sin(psi)*cos(theta), sin(psi)*sin(theta)*sin(phi)+cos(psi)*cos(phi), sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi);
             -sin(theta),         cos(theta)*sin(phi),                            cos(theta)*cos(phi)];
    vel_i = R_b2i * vel;

    % --- 反步参数 ---
    K1_p = diag([1.0, 1.0, 1.5]);   K2_p = diag([2.0, 2.0, 3.0]);
    K1_a = diag([4.0, 4.0, 4.0]);   K2_a = diag([3.0, 3.0, 3.0]);

    % ================= 位置环 =================
    z1_p = pos - desired_state.pos;
    
    % 位置误差安全限幅
    max_err = 3.0;
    z1_p_sat = max(min(z1_p, max_err), -max_err);
    
    % 虚拟控制律 alpha_p
    alpha_p = desired_state.vel - K1_p * z1_p_sat;
    
    % 【经典反步核心】: 直接求导取代 DSC 滤波
    if isfield(ctrl_state, 'prev_alpha_p')
        dot_alpha_p = (alpha_p - ctrl_state.prev_alpha_p) / dt;
    else
        dot_alpha_p = zeros(3,1); % 初始状态处理
    end
    ctrl_state.prev_alpha_p = alpha_p; % 记录当前状态供下一步求导
    
    z2_p = vel_i - alpha_p;
    F_d = m * dot_alpha_p - z1_p - K2_p * z2_p;
    F_total = F_d + [0; 0; m*g];
    
    % 力学推力限幅
    if norm(F_total) > params.max_thrust
        F_total = (F_total / norm(F_total)) * params.max_thrust;
    end
    U1 = norm(F_total);
    if U1 < 1e-6
        z_b = [0; 0; 1];
    else
        z_b = F_total / U1; 
    end
    
    psi_d = atan2(desired_state.vel(2), desired_state.vel(1));
    phi_d = asin(max(-0.95, min(0.95, z_b(1)*sin(psi_d) - z_b(2)*cos(psi_d))));
    theta_d = atan2(z_b(1)*cos(psi_d) + z_b(2)*sin(psi_d), z_b(3));
    desired_angles = [phi_d; max(-0.95, min(0.95, theta_d)); psi_d];

    % ================= 姿态环 =================
    % 期望姿态求导取代 DSC
    if isfield(ctrl_state, 'prev_desired_angles')
        dot_Theta_d = (desired_angles - ctrl_state.prev_desired_angles) / dt;
        % 处理偏航角跨越 pi/-pi 时的数值跳变
        dot_Theta_d(3) = atan2(sin(dot_Theta_d(3)*dt), cos(dot_Theta_d(3)*dt)) / dt;
    else
        dot_Theta_d = zeros(3,1);
    end
    ctrl_state.prev_desired_angles = desired_angles;
    
    z1_a = angles - desired_angles;
    z1_a(3) = wrapToPi(z1_a(3));
    
    W_inv = [1, 0, -sin(theta); 0, cos(phi), sin(phi)*cos(theta); 0, -sin(phi), cos(phi)*cos(theta)];
    alpha_a = W_inv * (dot_Theta_d - K1_a * z1_a);
    
    % 姿态虚拟控制律求导取代 DSC
    if isfield(ctrl_state, 'prev_alpha_a')
        dot_alpha_a = (alpha_a - ctrl_state.prev_alpha_a) / dt;
    else
        dot_alpha_a = zeros(3,1);
    end
    ctrl_state.prev_alpha_a = alpha_a;
    
    z2_a = omega - alpha_a;
    tau_cmd = cross(omega, J * omega) + J * dot_alpha_a - W_inv' * z1_a - K2_a * z2_a;
    
    U = [U1; tau_cmd(1); tau_cmd(2); tau_cmd(3)];
end