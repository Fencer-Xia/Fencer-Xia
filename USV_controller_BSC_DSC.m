function [tau, ctrl_state] = USV_controller_BSC_DSC(state, z1, desired_speed, ctrl_state, dt)
% 传统反步法 + 动态面控制 (USV)

    u = state(1); r = state(5);
    xe = z1(1); ye = z1(2); psie = z1(3);
    
    % --- 增益参数 (提升反馈刚度以应对风浪) ---
    k_xe = 2.0; ky = 1.0; k_psi = 2.0;
    k_u = 5.0;  k_r = 5.0;
    tau_f = 0.1; % DSC 滤波时间常数
    m_matrix = diag([198, 181]);
    
    % 1. 运动学虚拟控制律
    alpha_u = desired_speed * cos(psie) + k_xe * xe;
    alpha_r = ky * ye + k_psi * psie;
    alpha = [alpha_u; alpha_r];
    
    % 2. DSC 滤波
    alpha_f = ctrl_state.alpha_f;
    dot_alpha_f = (alpha - alpha_f) / tau_f;
    ctrl_state.alpha_f = alpha_f + dot_alpha_f * dt;
    
    % 3. 动力学控制律
    z2 = [u; r] - alpha_f;
    tau_cmd = m_matrix * (dot_alpha_f - z1([1,3]) - diag([k_u, k_r]) * z2);
    
    tau_max = [500; 300];
    tau_cmd = max(min(tau_cmd, tau_max), -tau_max);
    
    % 维度对齐补偿 (推力; 横漂设为0; 偏航)
    tau = [tau_cmd(1); 0; tau_cmd(2)];
end