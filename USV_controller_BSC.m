function [tau, ctrl_state] = USV_controller_BSC(state, z1, desired_speed, ctrl_state, dt)
% 传统反步法 (USV) - 移除DSC动态面，使用数值差分

    u = state(1); r = state(5);
    xe = z1(1); ye = z1(2); psie = z1(3);
    
    % --- 增益参数 ---
    k_xe = 2.0; ky = 1.0; k_psi = 2.0;
    k_u = 5.0;  k_r = 5.0;
    m_matrix = diag([198, 181]);
    
    % 1. 运动学虚拟控制律
    alpha_u = desired_speed * cos(psie) + k_xe * xe;
    alpha_r = ky * ye + k_psi * psie;
    alpha = [alpha_u; alpha_r];
    
    % 2. 【经典反步核心】直接求导取代 DSC 滤波
    if isfield(ctrl_state, 'prev_alpha')
        dot_alpha = (alpha - ctrl_state.prev_alpha) / dt;
    else
        dot_alpha = zeros(2,1); % 初始状态处理
    end
    ctrl_state.prev_alpha = alpha;
    
    % 3. 动力学控制律
    z2 = [u; r] - alpha;
    tau_cmd = m_matrix * (dot_alpha - z1([1,3]) - diag([k_u, k_r]) * z2);
    
    tau_max = [500; 300];
    tau_cmd = max(min(tau_cmd, tau_max), -tau_max);
    
    % 维度对齐补偿 (推力; 横漂设为0; 偏航)
    tau = [tau_cmd(1); 0; tau_cmd(2)];
end