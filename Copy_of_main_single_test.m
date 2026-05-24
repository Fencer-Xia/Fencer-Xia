% =========================================================================
% 单一算法测试台 (Testbench) - 用于标定不同控制器的参数
% 选项: 'AFBC' (本文核心), 'BSC' (经典反步), 'PID' (工业基准)
% =========================================================================
clc; clear; close all;

TEST_ALGO = 'BSC';  % 填写 'AFBC', 'BSC', 或 'PID'
fprintf('当前运行的算法架构为: %s\n', TEST_ALGO);

rng(0); 

dt = 0.01; T_sim = 300; time_vec = 0:dt:T_sim; N = length(time_vec);
num_uavs = 3; uav_spacing_default = 25; num_usvs = 8; usv_spacing_default = 10;
maneuver.v_c = 2.5; maneuver.R_c = 65; maneuver.usv_v_max = 3.0;

s_usv_new = max(3.0, (maneuver.R_c / (num_usvs/2 - 0.5)) * (maneuver.usv_v_max / maneuver.v_c - 1));
s_uav_new = s_usv_new * (uav_spacing_default / usv_spacing_default);

uav_params.m = 1.5; uav_params.g = 9.81;
uav_params.Ixx = 0.03; uav_params.Iyy = 0.03; uav_params.Izz = 0.06;
uav_params.max_thrust = 2 * uav_params.m * uav_params.g;
uav_params.K1 = diag([20, 20, 30]); uav_params.K2 = diag([16, 16, 12]);
uav_params.beta = 1.2; uav_params.Gamma = [0.01; 0.01; 0.01];
uav_params.K3_att = diag([10, 10, 15]); uav_params.K4_att = diag([1.5, 1.5, 2.0]);
uav_params.K_frac_att = diag([0.5, 0.5, 0.5]); uav_params.beta_att = 0.8;
uav_params.Gamma_att = diag([1e-3, 1e-3, 1e-3]); uav_params.sigma_att = [0.01; 0.01; 0.01];

usv_params.gamma_afbc = [0.5; 0.05; 0.02]; usv_params.l_hat_max = [10.0; 1.0; 0.5];

uav_state_history = zeros(12, N, num_uavs);
uav_l_hat = repmat([0.1; 0.1; 0.1], 1, num_uavs);
uav_l_hat_att = repmat([0.1; 0.01; 0.01], 1, num_uavs);

uav_ctrl_state = cell(1, num_uavs);
for k = 1:num_uavs
    uav_state_history(1:3, 1, k) = [(-2+4*rand()); (k - 2) * uav_spacing_default; (10+5*rand())];
    
    uav_ctrl_state{k}.Theta_d = [0;0;0]; 
    uav_ctrl_state{k}.Theta_d_dot = [0;0;0]; 
    uav_ctrl_state{k}.Omega_d = [0;0;0]; 
    uav_ctrl_state{k}.Omega_d_dot = [0;0;0];
    
    uav_ctrl_state{k}.prev_Theta_d = [0;0;0]; 
    
    % BSC 专用状态
    uav_ctrl_state{k}.prev_alpha_p = [0;0;0]; 
    uav_ctrl_state{k}.prev_alpha_a = [0;0;0];
    uav_ctrl_state{k}.prev_desired_angles = [0;0;0];
    
    uav_ctrl_state{k}.int_p = [0;0;0]; 
    uav_ctrl_state{k}.int_a = [0;0;0];
end

usv_state.x = (rand(1, num_usvs) - 0.5) * 20;
usv_state.psi = (rand(1, num_usvs) - 0.5) * (pi/2);
for j = 1:num_usvs; usv_state.y(j) = (j - 4.5) * usv_spacing_default + (rand - 0.5) * 10; end
usv_state.u = ones(1, num_usvs) * maneuver.v_c;
usv_state.r = zeros(1, num_usvs);
usv_state.l_hat = repmat([0.1; 0.01; 0.01], 1, num_usvs);

usv_state_history = cell(1, num_usvs);
usv_ctrl_state = cell(1, num_usvs);
for j = 1:num_usvs
    usv_state_history{j} = zeros(2, N); usv_state_history{j}(:, 1) = [usv_state.x(j); usv_state.y(j)];
    % BSC 专用
    usv_ctrl_state{j}.prev_alpha = [maneuver.v_c; 0];
    
    usv_ctrl_state{j}.int_u = 0; usv_ctrl_state{j}.int_r = 0;
    usv_ctrl_state{j}.prev_e_u = 0; usv_ctrl_state{j}.prev_e_psi = 0;
    usv_ctrl_state{j}.integral_err_kin = [0;0]; 
end

start_pos = [10; 0; 20]; start_psi = 0;
plan = build_trajectory_plan(T_sim, maneuver.v_c, maneuver.R_c, start_pos, start_psi);

comm_params.rho_max = 0.2; comm_params.rho_min = 0.1; comm_params.rho_psi = 0.05;
comm_params.T_relax = 8.0; comm_params.k_accel = 0.1; comm_params.wakeup_time = 0.2;
smooth_params.Kp = 25.0; smooth_params.Kd = 10.0;

usv_comm_state = cell(1, num_usvs);
for j = 1:num_usvs
    if j<=3; l_idx=1; elseif j<=5; l_idx=2; else; l_idx=3; end
    p_init = uav_state_history(1:2, 1, l_idx); v_init = [maneuver.v_c; 0];
    usv_comm_state{j}.p_uav_last = p_init; usv_comm_state{j}.v_uav_last = v_init;
    usv_comm_state{j}.psi_uav_last = 0; usv_comm_state{j}.r_uav_last = 0;
    usv_comm_state{j}.t_last = 0; usv_comm_state{j}.rho = comm_params.rho_max;
    usv_comm_state{j}.p_smooth = p_init; usv_comm_state{j}.v_smooth = v_init;
    usv_comm_state{j}.psi_smooth = 0; usv_comm_state{j}.r_smooth = 0;
    usv_comm_state{j}.comm_events = []; usv_comm_state{j}.error_history = zeros(1, N);
end

uav_error_history = zeros(1, N, num_uavs);
usv_error_history = zeros(1, N, num_usvs);
uav_desired_accel = zeros(3, num_uavs);
traj_leader_actual = zeros(3, N); 
traj_leader_actual(:, 1) = start_pos;
uav_U_history = zeros(4, N, num_uavs);

fprintf('启动仿真...\n');
for i = 1:N-1
    t = time_vec(i);
    d_uav = [120+70*sin(0.4*t)+30*cos(1.5*t); 40+50*cos(0.5*t)+5*sin(1.8*t); 0] * 0.005;
    d_usv = [120+70*sin(0.4*t)+30*cos(1.5*t); 40+50*cos(0.5*t)+5*sin(1.8*t); 0; 0; 50+50*sin(0.3*t)+15*cos(1.2*t)];

    [pos_d_leader, vel_d_leader, p_ratio] = get_trajectory_state(t, plan, maneuver.v_c, maneuver.R_c);
    traj_leader_actual(:, i+1) = pos_d_leader;
    target_spacing_uav = uav_spacing_default + (s_uav_new - uav_spacing_default) * p_ratio;
    target_spacing_usv = usv_spacing_default + (s_usv_new - usv_spacing_default) * p_ratio;
    
    psi_form = atan2(vel_d_leader(2), vel_d_leader(1));
    R_z_form = [cos(psi_form), -sin(psi_form); sin(psi_form), cos(psi_form)];
    
    des_state_l.pos = pos_d_leader; des_state_l.vel = vel_d_leader; des_state_l.acc = [0;0;0];

    % --------- UAV ---------
    for k = 1:num_uavs
        curr_uav = uav_state_history(:, i, k);
        offset_i = R_z_form * [0; (k - 2) * target_spacing_uav];
        des_uav = des_state_l; des_uav.pos(1:2) = des_state_l.pos(1:2) + offset_i;
        
        switch TEST_ALGO
            case 'AFBC'
                [U1, des_a, l_dot] = UAV_position_controller_AFBC(curr_uav, des_uav, uav_l_hat(:,k), uav_params);
                [U2, U3, U4, latt_dot, uav_ctrl_state{k}] = UAV_attitude_controller_AFBC(curr_uav, des_a, uav_l_hat_att(:,k), uav_ctrl_state{k}, dt, uav_params);
                U = [U1; U2; U3; U4];
                uav_l_hat(:,k) = max(0, uav_l_hat(:,k) + l_dot * dt);
                uav_l_hat_att(:,k) = max(0, uav_l_hat_att(:,k) + latt_dot * dt);
            case 'BSC'
                [U, uav_ctrl_state{k}] = UAV_controller_BSC(curr_uav, des_uav, uav_ctrl_state{k}, uav_params, dt);
            case 'PID'
                [U, uav_ctrl_state{k}] = UAV_controller_PID(curr_uav, des_uav, uav_ctrl_state{k}, uav_params, dt);
        end
        
        U(1) = max(0, min(U(1), uav_params.max_thrust));
        U(2) = max(-5, min(U(2), 5));
        U(3) = max(-5, min(U(3), 5));
        U(4) = max(-2, min(U(4), 2));
        uav_U_history(:, i+1, k) = U;
        s_dot = UAV_dynamic_function(curr_uav, U);
        s_dot(7:9) = s_dot(7:9) + d_uav / uav_params.m;
        uav_state_history(:, i+1, k) = curr_uav + s_dot * dt;
        uav_error_history(1, i+1, k) = norm(curr_uav(1:3) - des_uav.pos);
        uav_desired_accel(:,k) = (des_uav.vel - curr_uav(7:9)); 
    end

    % --------- USV ---------
    is_wakeup = any(abs(t - ([plan.t_start] - comm_params.wakeup_time)) < dt/2);
    for j = 1:num_usvs
        if j<=3; l_idx=1; elseif j<=5; l_idx=2; else; l_idx=3; end
        uav_full = uav_state_history(:, i, l_idx);
        p_real = uav_full(1:2); psi_real = uav_full(6); r_real = uav_full(12);
        v_real = [uav_full(7)*cos(psi_real)-uav_full(8)*sin(psi_real); uav_full(7)*sin(psi_real)+uav_full(8)*cos(psi_real)];
        
        p_pred = usv_comm_state{j}.p_uav_last + usv_comm_state{j}.v_uav_last * (t - usv_comm_state{j}.t_last);
        psi_pred = wrapToPi(usv_comm_state{j}.psi_uav_last + usv_comm_state{j}.r_uav_last * (t - usv_comm_state{j}.t_last));
        
        err_p = norm(p_real - p_pred);
        usv_comm_state{j}.rho = max(comm_params.rho_min, min(comm_params.rho_max, usv_comm_state{j}.rho + ((comm_params.rho_max - usv_comm_state{j}.rho)/comm_params.T_relax - comm_params.k_accel*norm(uav_desired_accel(1:2,l_idx))) * dt));
        
        if is_wakeup || (err_p >= usv_comm_state{j}.rho)
            usv_comm_state{j}.p_uav_last = p_real; usv_comm_state{j}.v_uav_last = v_real;
            usv_comm_state{j}.psi_uav_last = psi_real; usv_comm_state{j}.r_uav_last = r_real;
            usv_comm_state{j}.t_last = t; usv_comm_state{j}.comm_events = [usv_comm_state{j}.comm_events, t];
            p_pred = p_real; psi_pred = psi_real;
        end
        
        p_smooth = usv_comm_state{j}.p_smooth; v_smooth = usv_comm_state{j}.v_smooth;
        a_sm = smooth_params.Kp*(p_pred - p_smooth) + smooth_params.Kd*(usv_comm_state{j}.v_uav_last - v_smooth);
        usv_comm_state{j}.v_smooth = v_smooth + a_sm*dt; usv_comm_state{j}.p_smooth = p_smooth + usv_comm_state{j}.v_smooth*dt;
        
        psi_sm = usv_comm_state{j}.psi_smooth; r_sm = usv_comm_state{j}.r_smooth;
        a_psi_sm = smooth_params.Kp*wrapToPi(psi_pred - psi_sm) + smooth_params.Kd*(usv_comm_state{j}.r_uav_last - r_sm);
        usv_comm_state{j}.r_smooth = r_sm + a_psi_sm*dt; usv_comm_state{j}.psi_smooth = wrapToPi(psi_sm + usv_comm_state{j}.r_smooth*dt);
        
        [~, ~, next_p] = get_trajectory_state(t+dt, plan, maneuver.v_c, maneuver.R_c);
        next_s_usv = usv_spacing_default + (s_usv_new - usv_spacing_default) * next_p;
        if j<=3; s_idx=j-2; elseif j<=5; s_idx=j-4.5; else; s_idx=j-7; end
        
        c_Rz = [cos(psi_sm), -sin(psi_sm); sin(psi_sm), cos(psi_sm)];
        n_Rz = [cos(usv_comm_state{j}.psi_smooth), -sin(usv_comm_state{j}.psi_smooth); sin(usv_comm_state{j}.psi_smooth), cos(usv_comm_state{j}.psi_smooth)];
        
        pd_curr = usv_comm_state{j}.p_smooth + c_Rz * [0; s_idx * target_spacing_usv];
        pd_next = (usv_comm_state{j}.p_smooth + usv_comm_state{j}.v_smooth*dt) + n_Rz * [0; s_idx * next_s_usv];
        
        vd_usv = (pd_next - pd_curr)/dt;
        psi_d = atan2(vd_usv(2), vd_usv(1)); speed_d = norm(vd_usv);
        
        x_b = usv_state.x(j); y_b = usv_state.y(j); psi_b = usv_state.psi(j); u_b = usv_state.u(j); r_b = usv_state.r(j);
        err_i = pd_curr - [x_b; y_b];
        ye = -err_i(1)*sin(psi_b) + err_i(2)*cos(psi_b);
        z1 = [err_i(1)*cos(psi_b) + err_i(2)*sin(psi_b); ye; wrapToPi(psi_d - psi_b)];
        
        dyn_s = [u_b; 0; 0; 0; r_b];
        
        switch TEST_ALGO
            case 'AFBC'
                if abs(ye) < 2.0; usv_ctrl_state{j}.integral_err_kin(1) = usv_ctrl_state{j}.integral_err_kin(1) + ye*dt; end
                [tau, z2] = control(dyn_s, z1, usv_state.l_hat(:,j), usv_ctrl_state{j}.integral_err_kin, speed_d);
                nv = norm([u_b; r_b]);
                l_dot = [usv_params.gamma_afbc(1)*norm(z2); usv_params.gamma_afbc(2)*norm(z2)*nv; usv_params.gamma_afbc(3)*norm(z2)*(nv^2)];
                usv_state.l_hat(:,j) = min(usv_params.l_hat_max, max(1e-6, usv_state.l_hat(:,j) + l_dot*dt));
            case 'BSC'
                [tau, usv_ctrl_state{j}] = USV_controller_BSC(dyn_s, z1, speed_d, usv_ctrl_state{j}, dt);
            case 'PID'
                [tau, usv_ctrl_state{j}] = USV_controller_PID(dyn_s, z1, speed_d, usv_ctrl_state{j}, dt);
        end
        
        acc = dynamic_function(tau, dyn_s, d_usv, [x_b; y_b; 0; 0]);
        usv_state.u(j) = u_b + acc(1)*dt; usv_state.r(j) = r_b + acc(5)*dt;
        usv_state.psi(j) = wrapToPi(psi_b + usv_state.r(j)*dt);
        usv_state.x(j) = x_b + usv_state.u(j)*cos(usv_state.psi(j))*dt; usv_state.y(j) = y_b + usv_state.u(j)*sin(usv_state.psi(j))*dt;
        
        usv_state_history{j}(:, i+1) = [usv_state.x(j); usv_state.y(j)];
        usv_error_history(1, i+1, j) = norm(pd_curr - [usv_state.x(j); usv_state.y(j)]);
    end
end
fprintf('仿真完成 [%s]\n', TEST_ALGO);

figure('Name', [TEST_ALGO ' 平均误差验证'], 'Color', 'w');
subplot(2,1,1); plot(time_vec, mean(uav_error_history, 3), 'b', 'LineWidth', 1.5); title('UAV 集群平均位置误差 (m)'); grid on;
subplot(2,1,2); plot(time_vec, mean(usv_error_history, 3), 'r', 'LineWidth', 1.5); title('USV 集群平均位置误差 (m)'); grid on; xlabel('Time (s)');

target_uav = 1; 

figure('Name', [TEST_ALGO ' 诊断 1: 3D 轨迹'], 'Color', 'w');
plot3(traj_leader_actual(1,:), traj_leader_actual(2,:), traj_leader_actual(3,:), 'k--', 'LineWidth', 2); hold on; grid on;
plot3(uav_state_history(1,:,target_uav), uav_state_history(2,:,target_uav), uav_state_history(3,:,target_uav), 'b-', 'LineWidth', 2);
title(['UAV ' num2str(target_uav) ' 3D 空间轨迹']);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
legend('期望领航者轨迹', 'UAV 实际轨迹');
view(3);

figure('Name', [TEST_ALGO ' 诊断 2: 姿态角监控'], 'Color', 'w', 'Position', [100, 100, 800, 600]);
subplot(3,1,1);
plot(time_vec, rad2deg(uav_state_history(4,:,target_uav)), 'r', 'LineWidth', 1.5);
title('滚转角 Roll (\phi)'); ylabel('Deg (°)'); grid on;
subplot(3,1,2);
plot(time_vec, rad2deg(uav_state_history(5,:,target_uav)), 'g', 'LineWidth', 1.5);
title('俯仰角 Pitch (\theta)'); ylabel('Deg (°)'); grid on;
subplot(3,1,3);
plot(time_vec, rad2deg(uav_state_history(6,:,target_uav)), 'b', 'LineWidth', 1.5);
title('偏航角 Yaw (\psi)'); xlabel('Time (s)'); ylabel('Deg (°)'); grid on;

figure('Name', [TEST_ALGO ' 诊断 3: 物理控制输入饱和监控'], 'Color', 'w', 'Position', [900, 100, 800, 800]);
subplot(4,1,1);
plot(time_vec, uav_U_history(1,:,target_uav), 'k', 'LineWidth', 1.5);
yline(uav_params.max_thrust, 'r--'); title('总推力 U1 (N)'); grid on;
subplot(4,1,2);
plot(time_vec, uav_U_history(2,:,target_uav), 'r', 'LineWidth', 1.5);
yline(5, 'r--'); yline(-5, 'r--'); title('滚转力矩 U2 (N·m)'); grid on;
subplot(4,1,3);
plot(time_vec, uav_U_history(3,:,target_uav), 'g', 'LineWidth', 1.5);
yline(5, 'r--'); yline(-5, 'r--'); title('俯仰力矩 U3 (N·m)'); grid on;
subplot(4,1,4);
plot(time_vec, uav_U_history(4,:,target_uav), 'b', 'LineWidth', 1.5);
yline(2, 'r--'); yline(-2, 'r--'); title('偏航力矩 U4 (N·m)'); xlabel('Time (s)'); grid on;

function plan = build_trajectory_plan(T_max, v_c, R_c, start_pos, start_psi)
    T_turn = (R_c * pi / 2) / v_c;  
    block = [
        T_turn,  1, 1, 1; 5,       0, 1, 1; T_turn,  1, 1, 1; 5,       0, 1, 1; 
        27,      0, 1, 0; 100,     0, 0, 0; 27,      0, 0, 1; T_turn, -1, 1, 1; 
        5,       0, 1, 1; T_turn, -1, 1, 1; 5,       0, 1, 1; 27,      0, 1, 0; 
        100,     0, 0, 0; 27,      0, 0, 1  
    ];
    seq = [132, 0, 0, 0; 27,  0, 0, 1];
    while sum(seq(:,1)) < T_max
        seq = [seq; block]; %#ok<AGROW> 
    end
    num_segments = size(seq, 1);
    plan(num_segments).t_start = 0; 
    current_t = 0; current_pos = start_pos; current_psi = start_psi;
    for idx = 1:num_segments
        dur = seq(idx, 1); type = seq(idx, 2);
        plan(idx).t_start = current_t; plan(idx).t_end = current_t + dur;
        plan(idx).duration = dur;      plan(idx).type = type;
        plan(idx).p_start = seq(idx, 3); plan(idx).p_end = seq(idx, 4);
        plan(idx).x0 = current_pos(1); plan(idx).y0 = current_pos(2);
        plan(idx).z0 = current_pos(3); plan(idx).psi0 = current_psi;
        if type == 0
            current_pos(1) = current_pos(1) + v_c * dur * cos(current_psi);
            current_pos(2) = current_pos(2) + v_c * dur * sin(current_psi);
        elseif type == 1
            cx = current_pos(1) - R_c * sin(current_psi); cy = current_pos(2) + R_c * cos(current_psi);
            current_psi = current_psi + pi/2;
            current_pos(1) = cx + R_c * sin(current_psi); current_pos(2) = cy - R_c * cos(current_psi);
        elseif type == -1
            cx = current_pos(1) + R_c * sin(current_psi); cy = current_pos(2) - R_c * cos(current_psi);
            current_psi = current_psi - pi/2;
            current_pos(1) = cx - R_c * sin(current_psi); current_pos(2) = cy + R_c * cos(current_psi);
        end
        current_t = current_t + dur;
    end
end

function [pos_d, vel_d, p] = get_trajectory_state(t, plan, v_c, R_c)
    idx = find(t >= [plan.t_start] & t < [plan.t_end], 1, 'last');
    if isempty(idx)
        idx = length(plan); t = plan(idx).t_end; 
    end
    seg = plan(idx); dt = t - seg.t_start;
    if seg.p_start == seg.p_end
        p = seg.p_start;
    else
        tau = dt / seg.duration;
        p_trans = 6*tau^5 - 15*tau^4 + 10*tau^3; 
        p = seg.p_start + (seg.p_end - seg.p_start) * p_trans;
    end
    if seg.type == 0
        pos_d = [seg.x0 + v_c * dt * cos(seg.psi0); seg.y0 + v_c * dt * sin(seg.psi0); seg.z0];
        vel_d = [v_c * cos(seg.psi0); v_c * sin(seg.psi0); 0];
    elseif seg.type == 1
        cx = seg.x0 - R_c * sin(seg.psi0); cy = seg.y0 + R_c * cos(seg.psi0);
        psi_t = seg.psi0 + (v_c / R_c) * dt;
        pos_d = [cx + R_c * sin(psi_t); cy - R_c * cos(psi_t); seg.z0];
        vel_d = [v_c * cos(psi_t); v_c * sin(psi_t); 0];
    elseif seg.type == -1
        cx = seg.x0 + R_c * sin(seg.psi0); cy = seg.y0 - R_c * cos(seg.psi0);
        psi_t = seg.psi0 - (v_c / R_c) * dt;
        pos_d = [cx - R_c * sin(psi_t); cy + R_c * cos(psi_t); seg.z0];
        vel_d = [v_c * cos(psi_t); v_c * sin(psi_t); 0];
    end
end