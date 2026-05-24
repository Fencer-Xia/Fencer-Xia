function [tau, ctrl_state] = USV_controller_PID(state, z1, desired_speed, ctrl_state, dt)
% 带抗积分饱和的 PID 控制器 (USV) - 解除封印的工业级暴力版

    u = state(1); r = state(5);
    xe = z1(1); ye = z1(2); psie = z1(3);
    
    % ==========================================================
    % 【核心修复 1】：提供极高的系统刚度，彻底碾压 200N 的海风扰动
    kp_u = 300;  ki_u = 50;   kd_u = 50;
    kp_r = 300;  ki_r = 50;   kd_r = 50;
    
    % 【核心修复 2】：放开积分限幅！允许积分项提供高达 2500N 的力
    % 这样 PID 就能轻松抵消海风，让 P 项纯粹用来追踪速度
    max_int = [50; 50]; 
    % ==========================================================
    
    % --- 外环：运动学视线制导 (LOS) ---
    k_xe = 1.5;
    % 【核心修复 3】：解除 1.5m/s 的龟速限制！允许最大 4.0m/s 的额外追赶速度
    u_catchup = max(min(k_xe * xe, 4.0), -4.0); 
    
    % 保持 cos(psie) 的几何投影，防止横向偏差剧烈拉大
    target_u = desired_speed * cos(psie) + u_catchup;
    e_u = target_u - u;
    
    % 保持完美的 atan 视线角映射，防止方向盘反折
    k_ye = 0.8; 
    e_steer = wrapToPi(psie + atan(k_ye * ye)); 
    
    % --- 内环：PID 执行 ---
    ctrl_state.int_u = ctrl_state.int_u + e_u * dt;
    ctrl_state.int_r = ctrl_state.int_r + e_steer * dt; 
    ctrl_state.int_u = max(min(ctrl_state.int_u, max_int(1)), -max_int(1));
    ctrl_state.int_r = max(min(ctrl_state.int_r, max_int(2)), -max_int(2));
    
    de_u = (e_u - ctrl_state.prev_e_u) / dt;
    de_steer = wrapToPi(e_steer - ctrl_state.prev_e_psi) / dt;
    
    ctrl_state.prev_e_u = e_u; 
    ctrl_state.prev_e_psi = e_steer;
    
    tau_u = kp_u * e_u + ki_u * ctrl_state.int_u + kd_u * de_u;
    tau_r = kp_r * e_steer + ki_r * ctrl_state.int_r + kd_r * de_steer;
    
    % 物理推力硬约束 (最终防线)
    tau_max = [500; 300];
    tau_cmd = [tau_u; tau_r];
    tau_cmd = max(min(tau_cmd, tau_max), -tau_max);
    
    tau = [tau_cmd(1); 0; tau_cmd(2)];
end