function [U1, desired_angles, l_hat_dot] = UAV_position_controller_AFBC(current_state, desired_state, l_hat, params)
% UAV_position_controller_AFBC - 无模型纯净反步法位置环 (紧凑标量版)

    g = params.g; 
    
    % ==========================================================
    % 【核心替换】：剥离真实物理质量，使用紧凑虚拟质量标量
    % ==========================================================
    m_uav_virtual = 1.55; 

    beta = params.beta; Gamma = params.Gamma;

    K1_kin = diag([1.5, 1.5, 2.0]); 
    K2_kin = diag([0.5, 0.5, 0.8]); 
    K1_dyn = diag([8.0, 8.0, 12.0]);
    K2_dyn = diag([4.0, 4.0, 6.0]); 

    pos = current_state(1:3); angles = current_state(4:6); vel_body = current_state(7:9);
    phi = angles(1); theta = angles(2); psi = angles(3);
    R_body_to_inertial = [cos(psi)*cos(theta), cos(psi)*sin(theta)*sin(phi)-sin(psi)*cos(phi), cos(psi)*sin(theta)*cos(phi)+sin(psi)*sin(phi);
                          sin(psi)*cos(theta), sin(psi)*sin(theta)*sin(phi)+cos(psi)*cos(phi), sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi);
                          -sin(theta),         cos(theta)*sin(phi),                            cos(theta)*cos(phi)];
    vel_inertial = R_body_to_inertial * vel_body;

    z1 = pos - desired_state.pos;
    z1_dot = vel_inertial - desired_state.vel;
    max_err = 3.0;
    z1_sat = max(min(z1, max_err), -max_err);

    % 剥离出反馈纠偏速度
    feedback_v = K1_kin * tanh(z1_sat) + K2_kin * (abs(z1_sat).^(2*beta-1)) .* tanh(z1_sat / 0.1);
    feedback_v = max(min(feedback_v, 2.0), -2.0); 
    
    alpha_p = desired_state.vel - feedback_v;

    z2 = vel_inertial - alpha_p;

    sech_z1_sq = 1 - tanh(z1_sat).^2;
    epsilon_deriv = 1e-3;
    frac_deriv = (2*beta-1) * (abs(z1_sat) + epsilon_deriv).^(2*beta-2);
    raw_alpha_p_dot = desired_state.acc - K1_kin * diag(sech_z1_sq) * z1_dot - K2_kin * diag(frac_deriv) * z1_dot;
    alpha_p_dot = max(min(raw_alpha_p_dot, 5.0), -5.0); 

    norm_vel_inertial = norm(vel_inertial);
    sigma = l_hat(1) + l_hat(2) * norm_vel_inertial + l_hat(3) * norm_vel_inertial^2;

    % ==========================================================
    % 纯无模型架构计算：全部使用 m_uav_virtual 标量计算前馈与补偿
    % ==========================================================
    F_d = m_uav_virtual * alpha_p_dot - z1 - K1_dyn * tanh(z2) - sigma * tanh(z2/0.1) - K2_dyn * diag(abs(z2).^(2*beta-1)) * tanh(z2 / 0.1);
    Force_total_desired = F_d + [0; 0; m_uav_virtual * g];

    max_thrust_local = m_uav_virtual * g * 2.5; 
    if norm(Force_total_desired) > max_thrust_local
        Force_total_desired = (Force_total_desired / norm(Force_total_desired)) * max_thrust_local;
    end

    U1 = norm(Force_total_desired);
    if U1 < 1e-6; z_b_desired = [0; 0; 1]; else; z_b_desired = Force_total_desired / U1; end

    desired_vel_xy = desired_state.vel(1:2);
    if norm(desired_vel_xy) > 0.1; psi_d = atan2(desired_vel_xy(2), desired_vel_xy(1)); else; psi_d = angles(6); end

    asin_arg = z_b_desired(1)*sin(psi_d) - z_b_desired(2)*cos(psi_d);
    asin_arg = max(-0.95, min(0.95, asin_arg)); 
    phi_d   = asin(asin_arg);
    theta_d = atan2(z_b_desired(1)*cos(psi_d) + z_b_desired(2)*sin(psi_d), z_b_desired(3));
    theta_d = max(-0.95, min(0.95, theta_d));

    desired_angles = [phi_d; theta_d; psi_d];

    k1 = min(diag(K1_dyn)); k2 = min(diag(K2_dyn)); 
    l_hat_non_negative = max(0, l_hat);

    l0_hat_dot = Gamma(1) * norm(z2) - Gamma(1) * k1 * sign(l_hat(1)) - 2*beta*(Gamma(1)/k2)^(1-beta)*k2*l_hat_non_negative(1)^(2*beta-1);
    l1_hat_dot = Gamma(2) * norm(z2) * norm_vel_inertial - Gamma(2) * k1 * sign(l_hat(2)) - 2*beta*(Gamma(2)/k2)^(1-beta)*k2*l_hat_non_negative(2)^(2*beta-1);
    l2_hat_dot = Gamma(3) * norm(z2) * norm_vel_inertial^2 - Gamma(3) * k1 * sign(l_hat(3)) - 2*beta*(Gamma(3)/k2)^(1-beta)*k2*l_hat_non_negative(3)^(2*beta-1);

    l_hat_dot = [l0_hat_dot; l1_hat_dot; l2_hat_dot];
end