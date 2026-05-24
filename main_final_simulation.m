% FILE: main_final_simulation.m (S形平行扫视搜索版本)
% =========================================================================
% 【最终融合方案: 带全机动预判的主动平滑事件触发通信机制】
%   第1层: 主动唤醒包
%   第2层: 自适应事件触发
%   第3层: 目标平滑器
%   附  加: 标准S形海上平行扫视搜索航迹 (严格几何对齐)
% =========================================================================
clc; clear; close all;
fprintf('开始执行【UAV姿态环AFBC + USV协同编队仿真 (S形搜索)】...\n');

% 仿真与编队参数 (延长仿真时间以展示完整S形)
dt = 0.01; T_sim = 400; time_vec = 0:dt:T_sim; N = length(time_vec);
num_uavs = 3; uav_spacing_default = 25; num_usvs = 8; usv_spacing_default = 10;
maneuver.v_c = 2.5; maneuver.R_c = 65; maneuver.usv_v_max = 3.0;

% 计算安全间距
fprintf('正在根据USV性能上限计算转弯所需的安全间距...\n');
usv_outermost_factor = (num_usvs/2 - 0.5);
s_usv_new_safe = (maneuver.R_c / usv_outermost_factor) * (maneuver.usv_v_max / maneuver.v_c - 1);
s_usv_new = max(3.0, s_usv_new_safe);
spacing_ratio = uav_spacing_default / usv_spacing_default;
s_uav_new = s_usv_new * spacing_ratio;
fprintf('  - 计算得出USV安全间距为: %.2f m\n', s_usv_new);
fprintf('  - 按比例UAV应收缩至间距: %.2f m\n', s_uav_new);

% ========================= 智能体模型参数 =========================
% UAV 物理参数
uav_params.m = 1.5;
uav_params.g = 9.81;
uav_params.Ixx = 0.03;
uav_params.Iyy = 0.03;
uav_params.Izz = 0.06;
uav_params.max_thrust = 2 * uav_params.m * uav_params.g;

% UAV 外环 AFBC 参数 (位置环)
uav_params.K1 = diag([20, 20, 30]);
uav_params.K2 = diag([16, 16, 12]);
uav_params.beta = 1.2;
uav_params.Gamma = [0.01; 0.01; 0.01];

% UAV 内环 AFBC 参数 (姿态环)
uav_params.K3_att = diag([10.0, 10.0, 15.0]);      % 运动学环增益
uav_params.K4_att = diag([1.5, 1.5, 2.0]);        % 动力学环增益
uav_params.K_frac_att = diag([0.5, 0.5, 0.5]);    % 分数阶非线性增益
uav_params.beta_att = 0.8;                        % 分数阶指数
uav_params.Gamma_att = diag([1e-3, 1e-3, 1e-3]);  % 内环自适应学习率
uav_params.sigma_att = [0.01; 0.01; 0.01];        % 鲁棒泄漏项系数

% USV 控制参数
usv_params.gamma_afbc = [0.5; 0.05; 0.02];
usv_params.l_hat_max = [10.0; 1.0; 0.5];
usv_params.sigma_usv = [0.02; 0.02; 0.02];

% ========================= 状态初始化 =========================
uav_state_history = zeros(12, N, num_uavs);
uav_l_hat = repmat([0.1; 0.1; 0.1], 1, num_uavs);
uav_center_idx = (num_uavs + 1) / 2;

rng(0);
x_range = [-2, 2]; z_range = [15, 20];
for k = 1:num_uavs
    y_base = (k - uav_center_idx) * uav_spacing_default;
    x_init = x_range(1) + (x_range(2) - x_range(1)) * rand();
    z_init = z_range(1) + (z_range(2) - z_range(1)) * rand();
    uav_state_history(1:3, 1, k) = [x_init; y_base; z_init];
    uav_state_history(4:12, 1, k) = 0;
end

uav_l_hat_att = repmat([0.1; 0.01; 0.01], 1, num_uavs);
uav_att_state = cell(1, num_uavs);
for k = 1:num_uavs
    uav_att_state{k}.Theta_d = uav_state_history(4:6, 1, k);
    uav_att_state{k}.prev_Theta_d = uav_state_history(4:6, 1, k);
    uav_att_state{k}.Theta_d_dot = [0; 0; 0];
    uav_att_state{k}.Omega_d = [0; 0; 0];
    uav_att_state{k}.Omega_d_dot = [0; 0; 0];
end

% USV 状态
usv_state.x = (rand(1, num_usvs) - 0.5) * 20;
usv_state.psi = (rand(1, num_usvs) - 0.5) * (pi/2);
usv_center_idx = (num_usvs + 1) / 2;
for j = 1:num_usvs
    usv_state.y(j) = (j - usv_center_idx) * usv_spacing_default + (rand - 0.5) * 10;
end
usv_state.u = ones(1, num_usvs) * maneuver.v_c;
usv_state.r = zeros(1, num_usvs);
usv_state.l_hat = repmat([0.1; 0.01; 0.01], 1, num_usvs);
usv_state_history = cell(1, num_usvs);
for j = 1:num_usvs
    usv_state_history{j} = zeros(2, N);
    usv_state_history{j}(:, 1) = [usv_state.x(j); usv_state.y(j)];
end

% ========================= 航线生成 (状态机式S形搜索) =========================
fprintf('正在生成S形平行扫视全局航线 (确保X坐标严格对齐)...\n');
start_pos = [10; 0; 20];
start_psi = 0; % 初始向东
plan = build_trajectory_plan(T_sim, maneuver.v_c, maneuver.R_c, start_pos, start_psi);
traj_leader_actual = zeros(3, N);

uav_error_history = zeros(4, N, num_uavs);
usv_error_history = zeros(2, N, num_usvs);
usv_integral_error = zeros(2, num_usvs);
uav_heading_history = zeros(N, num_uavs);
usv_heading_history = zeros(N, num_usvs);

for k = 1:num_uavs; uav_heading_history(1, k) = uav_state_history(6, 1, k); end
for j = 1:num_usvs; usv_heading_history(1, j) = usv_state.psi(j); end

% ========================= 通信与平滑器参数 =========================
comm_params.rho_max = 0.2;
comm_params.rho_min = 0.1;
comm_params.T_relax = 8.0;
comm_params.k_accel = 0.1;
comm_params.wakeup_time = 0.2;


smooth_params.Kp = 25.0;
smooth_params.Kd = 10.0;

usv_comm_state = cell(1, num_usvs);
for j = 1:num_usvs
    if j <= 3; leader_uav_idx = 1; elseif j <= 5; leader_uav_idx = 2; else; leader_uav_idx = 3; end
    p_uav_init = uav_state_history(1:2, 1, leader_uav_idx);
    v_uav_init = [maneuver.v_c; 0];
    usv_comm_state{j}.p_uav_last = p_uav_init;
    usv_comm_state{j}.v_uav_last = v_uav_init;
    usv_comm_state{j}.psi_uav_last = 0; 
    usv_comm_state{j}.r_uav_last = 0;
    usv_comm_state{j}.t_last = 0;
    usv_comm_state{j}.rho = comm_params.rho_max;
    usv_comm_state{j}.p_smooth = p_uav_init;
    usv_comm_state{j}.v_smooth = v_uav_init;
    usv_comm_state{j}.psi_smooth = 0;   
    usv_comm_state{j}.r_smooth = 0;
    usv_comm_state{j}.error_history = zeros(1, N);
    usv_comm_state{j}.rho_history = zeros(1, N);
    usv_comm_state{j}.comm_events = [];
    usv_comm_state{j}.p_predicted_hist = zeros(2, N);
    usv_comm_state{j}.p_smooth_hist = zeros(2, N);
end

uav_desired_accel = zeros(3, num_uavs);
usv_tau_history = zeros(2, N, num_usvs);

% ========================= 主仿真循环 =========================
fprintf('开始主仿真循环...\n');
for i = 1:N-1
    t = time_vec(i);
    wind_force_x = 120 + 70 * sin(0.4 * t) + 30 * cos(1.5 * t);
    wind_force_y = 40 + 50 * cos(0.5 * t) + 5 * sin(1.8 * t);
    wind_moment_z = 50 + 50 * sin(0.3 * t) + 15 * cos(1.2 * t);

    % 获取当前目标状态与队形缩放系数 (摒弃复杂if-else)
    
    [pos_d_leader, vel_d_leader, p_ratio] = get_trajectory_state(t, plan, maneuver.v_c, maneuver.R_c);
    
    % 重新组装为控制器所需的结构体格式
    desired_state_leader.pos = pos_d_leader;
    desired_state_leader.vel = vel_d_leader;
    desired_state_leader.acc = [0; 0; 0]; % 前馈加速度设为0即可
    
    traj_leader_actual(:, i+1) = desired_state_leader.pos;
    
    target_spacing_uav = uav_spacing_default + (s_uav_new - uav_spacing_default) * p_ratio;
    target_spacing_usv = usv_spacing_default + (s_usv_new - usv_spacing_default) * p_ratio;
    
    psi_form = atan2(desired_state_leader.vel(2), desired_state_leader.vel(1));
    R_z_form = [cos(psi_form), -sin(psi_form); sin(psi_form), cos(psi_form)];

    % ==================== UAV 状态更新 ====================
    for k = 1:num_uavs
        current_uav_state = uav_state_history(:, i, k);
        offset_b = [0; (k - uav_center_idx) * target_spacing_uav];
        offset_i = R_z_form * offset_b;
        desired_state_uav = desired_state_leader;
        desired_state_uav.pos(1:2) = desired_state_leader.pos(1:2) + offset_i;
        
        [U1, des_a, l_hat_dot] = UAV_position_controller_AFBC(current_uav_state, desired_state_uav, uav_l_hat(:, k), uav_params);
        uav_desired_accel(:, k) = des_a;

        [U2, U3, U4, l_hat_att_dot, att_state_out] = UAV_attitude_controller_AFBC(...
            current_uav_state, des_a, uav_l_hat_att(:, k), uav_att_state{k}, dt, uav_params);
        
        uav_l_hat_att(:, k) = max(0, uav_l_hat_att(:, k) + l_hat_att_dot * dt);
        uav_att_state{k} = att_state_out;
        
        control_input = [U1; U2; U3; U4];
        control_input(1) = max(0, min(control_input(1), uav_params.max_thrust));
        control_input(2) = max(-5, min(control_input(2), 5));
        control_input(3) = max(-5, min(control_input(3), 5));
        control_input(4) = max(-2, min(control_input(4), 2));
        
        state_dot = UAV_dynamic_function(current_uav_state, control_input);
        uav_wind_scale = 0.005;
        dis_uav = [wind_force_x; wind_force_y; 0] * uav_wind_scale;
        state_dot(7:9) = state_dot(7:9) + dis_uav / uav_params.m;
        uav_state_history(:, i+1, k) = current_uav_state + state_dot * dt;
        
        uav_l_hat(:, k) = max(0, uav_l_hat(:, k) + l_hat_dot * dt);
        
        err_vec = uav_state_history(1:3, i+1, k) - desired_state_uav.pos;
        uav_error_history(1:3, i+1, k) = err_vec;
        uav_error_history(4, i+1, k) = norm(err_vec);
        uav_heading_history(i+1, k) = uav_state_history(6, i+1, k);
    end

    % ==================== USV 状态更新 ====================
    % 动态判断主动唤醒时刻 (预判所有航段切换点)
    is_wakeup_time = any(abs(t - ([plan.t_start] - comm_params.wakeup_time)) < dt/2);

    for j = 1:num_usvs
        if j <= 3; leader_uav_idx = 1; elseif j <= 5; leader_uav_idx = 2; else; leader_uav_idx = 3; end
        
        p_uav_real = uav_state_history(1:2, i, leader_uav_idx);
        current_uav_full_state = uav_state_history(:, i, leader_uav_idx);
        psi_uav_real = current_uav_full_state(6); % 【修改】提取真实航向角
        r_uav_real = current_uav_full_state(12);
        u_uav = current_uav_full_state(7);
        v_uav = current_uav_full_state(8);
        v_uav_real = [u_uav * cos(psi_uav_real) - v_uav * sin(psi_uav_real);
                      u_uav * sin(psi_uav_real) + v_uav * cos(psi_uav_real)];
        a_uav_intent = uav_desired_accel(1:2, leader_uav_idx);
        
        p_last = usv_comm_state{j}.p_uav_last;
        v_last = usv_comm_state{j}.v_uav_last;
        t_last = usv_comm_state{j}.t_last;
        rho_current = usv_comm_state{j}.rho;
        
        p_uav_predicted = p_last + v_last * (t - t_last);
        prediction_error = p_uav_real - p_uav_predicted;
        error_norm = norm(prediction_error);
        
        psi_last = usv_comm_state{j}.psi_uav_last;
        r_last = usv_comm_state{j}.r_uav_last;
        psi_uav_predicted = wrapToPi(psi_last + r_last * (t - t_last));

        rho_dot = (comm_params.rho_max - rho_current) / comm_params.T_relax - comm_params.k_accel * norm(a_uav_intent);
        rho_new = max(comm_params.rho_min, min(comm_params.rho_max, rho_current + rho_dot * dt));
        usv_comm_state{j}.rho = rho_new;
        
        is_adaptive_triggered = (error_norm >= rho_new);
        
        if is_wakeup_time || is_adaptive_triggered
            usv_comm_state{j}.p_uav_last = p_uav_real;
            usv_comm_state{j}.v_uav_last = v_uav_real;
            usv_comm_state{j}.psi_uav_last = psi_uav_real; % 【修改】通信下发真实航向角
            usv_comm_state{j}.r_uav_last = r_uav_real;
            usv_comm_state{j}.t_last = t;
            usv_comm_state{j}.comm_events = [usv_comm_state{j}.comm_events, t];
            p_uav_predicted = p_uav_real;
        end
        v_uav_predicted = usv_comm_state{j}.v_uav_last;
        
        p_smooth_current = usv_comm_state{j}.p_smooth;
        v_smooth_current = usv_comm_state{j}.v_smooth;        
        a_smooth = smooth_params.Kp * (p_uav_predicted - p_smooth_current) + ...
                   smooth_params.Kd * (v_uav_predicted - v_smooth_current);
        v_smooth_new = v_smooth_current + a_smooth * dt;
        p_smooth_new = p_smooth_current + v_smooth_current * dt;
        usv_comm_state{j}.p_smooth = p_smooth_new;
        usv_comm_state{j}.v_smooth = v_smooth_new;
        p_leader_for_usv = p_smooth_new;
        p_leader_for_usv_next = p_smooth_new + v_smooth_new * dt;
        
        % 【修改】航向角/角速度的独立平滑（消除原先由 XY 速度反推造成的极点奇点和严重滞后）
        psi_smooth_current = usv_comm_state{j}.psi_smooth;
        r_smooth_current = usv_comm_state{j}.r_smooth;
        % 使用 wrapToPi 处理角度跳变
        a_smooth_psi = smooth_params.Kp * wrapToPi(psi_uav_predicted - psi_smooth_current) + ...
                       smooth_params.Kd * (usv_comm_state{j}.r_uav_last - r_smooth_current);
        r_smooth_new = r_smooth_current + a_smooth_psi * dt;
        psi_smooth_new = wrapToPi(psi_smooth_current + r_smooth_current * dt);
        usv_comm_state{j}.psi_smooth = psi_smooth_new;
        usv_comm_state{j}.r_smooth = r_smooth_new;

        usv_comm_state{j}.error_history(i+1) = error_norm;
        usv_comm_state{j}.rho_history(i+1) = rho_new;
        usv_comm_state{j}.p_predicted_hist(:, i+1) = p_uav_predicted;
        usv_comm_state{j}.p_smooth_hist(:, i+1) = p_smooth_new;
        
        % ================= 替换开始 (修复参数不足报错) =================
        % 1. 获取当前与下一时刻的理想轨迹状态 (用于提取前馈角速度 r_form)
        [~, vel_d_leader_curr, ~] = get_trajectory_state(t, plan, maneuver.v_c, maneuver.R_c);
        [~, vel_d_leader_next, next_p] = get_trajectory_state(t+dt, plan, maneuver.v_c, maneuver.R_c);
        
        % 2. 计算理想路径航向与前馈角速度 r_form (捍卫稳定性证明的核心项)
        psi_form_curr = atan2(vel_d_leader_curr(2), vel_d_leader_curr(1));
        psi_form_next = atan2(vel_d_leader_next(2), vel_d_leader_next(1));
        r_form = wrapToPi(psi_form_next - psi_form_curr) / dt;

        % 3. 计算编队间距与索引
        next_spacing_usv = usv_spacing_default + (s_usv_new - usv_spacing_default) * next_p;
        if j <= 3; sub_idx = j - 2; elseif j <= 5; sub_idx = j - 4.5; else; sub_idx = j - 7; end
        
        % 4. 计算目标基准点 (使用理想航向旋转，彻底切断 UAV 姿态抖动对 USV 的污染)
        ideal_Rz = [cos(psi_form_curr), -sin(psi_form_curr); sin(psi_form_curr), cos(psi_form_curr)];
        pos_d_usv_current = p_leader_for_usv + ideal_Rz * [0; sub_idx * target_spacing_usv];
        
        % 5. 期望航向锁定几何切线，期望速度进行转弯补偿
        psi_d = psi_form_curr;
        speed_d_usv = maneuver.v_c - r_form * (sub_idx * target_spacing_usv);
        
        % 6. 状态解析
        x_b = usv_state.x(j); y_b = usv_state.y(j); psi_b = usv_state.psi(j);
        u_b = usv_state.u(j); r_b = usv_state.r(j);
        
        err_ix = pos_d_usv_current(1) - x_b; err_iy = pos_d_usv_current(2) - y_b;
        ye = -err_ix * sin(psi_b) + err_iy * cos(psi_b);
        psie = wrapToPi(psi_d - psi_b);
        z1 = [err_ix * cos(psi_b) + err_iy * sin(psi_b); ye; psie];
        
        if abs(ye) < 2.0; usv_integral_error(1, j) = usv_integral_error(1, j) + ye * dt; end
        dyn_s = [u_b; 0; 0; 0; r_b];
        l_hat_k = usv_state.l_hat(:, j);
        
        % 【核心修正】：传入第 6 个参数 r_form 以匹配修改后的 control.m 接口
        [tau, z2] = control(dyn_s, z1, l_hat_k, usv_integral_error(:, j), speed_d_usv, r_form);
        
        usv_tau_history(:, i+1, j) = tau([1, 3]);
        % ================= 替换结束 =================
        
        dis_usv = [120 + 70 * sin(0.4 * t) + 30 * cos(1.5 * t);
                   40 + 50 * cos(0.5 * t) + 5 * sin(1.8 * t);
                   0; 0; 
                   50 + 50 * sin(0.3 * t) + 15 * cos(1.2 * t)];
        
        acc = dynamic_function(tau, dyn_s, dis_usv, [x_b; y_b; 0; 0]);
        usv_state.u(j) = u_b + acc(1) * dt;
        usv_state.r(j) = r_b + acc(5) * dt;
        usv_state.psi(j) = wrapToPi(psi_b + usv_state.r(j) * dt);
        usv_state.x(j) = x_b + usv_state.u(j) * cos(usv_state.psi(j)) * dt;
        usv_state.y(j) = y_b + usv_state.u(j) * sin(usv_state.psi(j)) * dt;
        usv_state_history{j}(:, i+1) = [usv_state.x(j); usv_state.y(j)];
        usv_error_history(1, i+1, j) = norm(pos_d_usv_current - [x_b; y_b]);
        usv_error_history(2, i+1, j) = z1(3);
        
        nv = norm([u_b; r_b]);
        l_hat_dot = [usv_params.gamma_afbc(1) * norm(z2);
                     usv_params.gamma_afbc(2) * norm(z2) * nv;
                     usv_params.gamma_afbc(3) * norm(z2) * (nv^2)];
                     
        % 替换原本的一行 min-max 代码为 final1 的循环限幅
        usv_state.l_hat(:, j) = usv_state.l_hat(:, j) + l_hat_dot * dt;
        for l_idx = 1:3
            if usv_state.l_hat(l_idx, j) > usv_params.l_hat_max(l_idx)
                usv_state.l_hat(l_idx, j) = usv_params.l_hat_max(l_idx);
            end
        end
        usv_state.l_hat(:, j) = max(usv_state.l_hat(:, j), 1e-6);
        usv_heading_history(i+1, j) = usv_state.psi(j);
        
    end
end
fprintf('仿真完成！正在生成结果...\n');

% ========================= 绘图 =========================


% 【优化】：在原色相基础上大幅压暗，提升白底对比度
usv_colors = {'#D96B6B', '#5C96B8', '#71B5B5', '#9686BA', '#7AA3D1', '#4EA6D9', '#D494BA', '#DFA070'};
uav_colors = {'#E6A861', '#D98C3A', '#5FA35F'};

figure('Name', 'UAV-USV S形平行搜索(3D视图)', 'Color', 'w');
hold on; grid on; axis equal; view(30, 40);
plot3(traj_leader_actual(1,:), traj_leader_actual(2,:), traj_leader_actual(3,:), 'k--', 'LineWidth', 2, 'DisplayName', '虚拟领航者');
for k = 1:num_uavs
    plot3(uav_state_history(1,:,k), uav_state_history(2,:,k), uav_state_history(3,:,k), '-', 'Color', uav_colors{k}, 'LineWidth', 2);
end
for j = 1:num_usvs
    plot3(usv_state_history{j}(1,:), usv_state_history{j}(2,:), zeros(1, N), '-', 'Color', usv_colors{j}, 'LineWidth', 1.5);
end
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)'); title('UAV-USV S形平行搜索 3D轨迹');

figure('Name', 'UAV-USV S形平行搜索(2D顶视图)', 'Color', 'w');
hold on; grid on; axis equal;
plot(traj_leader_actual(1,:), traj_leader_actual(2,:), 'k--', 'LineWidth', 2, 'DisplayName', '虚拟领航者');
for k = 1:num_uavs
    plot(uav_state_history(1,:,k), uav_state_history(2,:,k), '-', 'Color', uav_colors{k}, 'LineWidth', 2);
end
for j = 1:num_usvs
    plot(usv_state_history{j}(1,:), usv_state_history{j}(2,:), '-.', 'Color', usv_colors{j}, 'LineWidth', 1.5);
end
xlabel('X (m)'); ylabel('Y (m)'); title('严格对齐的海上S形平行扫视搜寻轨迹');

% ========================= 误差与性能分析绘图 =========================
figure('Name', '全智能体跟踪误差监控', 'Color', 'w', 'Position', [100, 50, 1200, 900]);
for k = 1:num_uavs
    subplot(4, 3, k); hold on; grid on;
    plot(time_vec, uav_error_history(4, :, k), '-', 'Color', uav_colors{k}, 'LineWidth', 1.5);
    title(sprintf('UAV %d Error Norm', k)); ylabel('Error (m)');
end
for j = 1:num_usvs
    subplot(4, 3, 3 + j); hold on; grid on;
    plot(time_vec, usv_error_history(1, :, j), '-', 'Color', usv_colors{j}, 'LineWidth', 1.5);
    title(sprintf('USV %d Pos Error Norm', j)); ylabel('Error (m)');
    if j > 5; xlabel('Time (s)'); end
end

% 事件触发性能图 (示例USV 1,4,8)
figure('Name', '事件触发机制性能分析 (USV 1, 4, 8)', 'Color', 'w');
usv_to_plot = [1, 4, 8];
for plot_idx = 1:length(usv_to_plot)
    j = usv_to_plot(plot_idx);
    subplot(length(usv_to_plot), 1, plot_idx);
    hold on; grid on;
    plot(time_vec, usv_comm_state{j}.rho_history, 'k--', 'LineWidth', 1.5, 'DisplayName', '\rho(t) (性能边界)');
    plot(time_vec, usv_comm_state{j}.error_history, '-', 'Color', usv_colors{j}, 'LineWidth', 1.5, 'DisplayName', '|e(t)| (预测误差)');
    comm_times = usv_comm_state{j}.comm_events;
    if ~isempty(comm_times)
        h_stem = stem(comm_times, ones(size(comm_times)) * comm_params.rho_max * 0.95);
        set(h_stem, 'Color', 'k', 'Marker', '^', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'LineStyle', 'none', 'DisplayName', '通信时刻');
    end
    title(sprintf('USV %d 通信性能', j)); xlabel('Time (s)'); ylabel('距离 (m)'); legend('show'); ylim([0, comm_params.rho_max * 1.1]);
end

% 平滑器性能分析 
figure('Name', '目标平滑器核心价值：消除阶跃跳变', 'Color', 'w', 'Position', [300, 200, 900, 400]);
j = 2; leader_idx = 1; 
hold on; grid on; box on;
plot_idx = find(time_vec > 120 & time_vec < 160);
t_plot = time_vec(plot_idx);

raw_diff = usv_comm_state{j}.p_predicted_hist(1:2, plot_idx) - uav_state_history(1:2, plot_idx, leader_idx);
raw_deviation = vecnorm(raw_diff);
smooth_diff = usv_comm_state{j}.p_smooth_hist(1:2, plot_idx) - uav_state_history(1:2, plot_idx, leader_idx);
smooth_deviation = vecnorm(smooth_diff);

plot(t_plot, raw_deviation, 'k-.', 'LineWidth', 1.5, 'DisplayName', '未平滑目标跳变 (零阶保持器锯齿)');
plot(t_plot, smooth_deviation, '-', 'Color', usv_colors{j}, 'LineWidth', 2.5, 'DisplayName', '平滑后目标过渡 (PD平滑器响应)');

comm_times = usv_comm_state{j}.comm_events;
comm_in_range = comm_times(comm_times > 120 & comm_times < 160);
if ~isempty(comm_in_range)
    plot(comm_in_range, zeros(size(comm_in_range)), 'k^', 'MarkerFaceColor', 'k', 'MarkerSize', 6, 'DisplayName', '通信触发时刻');
end
title(sprintf('USV %d 接收虚拟领航者位置的平滑过渡分析 (2D 欧氏距离偏差)', j));
xlabel('Time (s)'); ylabel('相对于真实UAV的位置误差 (m)');
legend('Location', 'best');

% 控制输入分析
figure('Name', '关键USV控制输入分析', 'Color', 'w', 'Position', [500, 100, 800, 800]);
usv_to_plot_tau = [1, 4, 8];
tau_r_max = 300;
for plot_idx = 1:length(usv_to_plot_tau)
    j = usv_to_plot_tau(plot_idx);
    subplot(length(usv_to_plot_tau), 1, plot_idx);
    hold on; grid on; box on;
    plot(time_vec, squeeze(usv_tau_history(2, :, j)), '-', 'Color', usv_colors{j}, 'LineWidth', 1.5);
    plot(time_vec, ones(1, N) * tau_r_max, 'r-.', 'LineWidth', 1);
    plot(time_vec, -ones(1, N) * tau_r_max, 'r-.', 'LineWidth', 1);
    title(sprintf('USV %d 控制输入 (偏航力矩 \\tau_r)', j));
    xlabel('Time (s)'); ylabel('Torque (N·m)');
    xlim([0 T_sim]);
end

% 【关键修复】：安全调用的 stem 绘制事件触发间隔
figure('Name', '事件触发间隔分析 (USV 1, 4, 8)', 'Color', 'w', 'Position', [600, 150, 900, 500]);
hold on; grid on; box on;

target_usvs = [1, 4, 8];
markers = {'o', 's', '^'};

% --- 主图：全周期散点分布 (采用安全的句柄赋值法) ---
for idx = 1:length(target_usvs)
    j = target_usvs(idx);
    comm_times = usv_comm_state{j}.comm_events;
    if length(comm_times) > 1
        inter_event_times = diff(comm_times);
        h_stem = stem(comm_times(2:end), inter_event_times);
        % 单独设置属性，完全规避键值对报错
        set(h_stem, 'Color', usv_colors{j}, 'Marker', markers{idx}, ...
            'MarkerFaceColor', usv_colors{j}, 'MarkerEdgeColor', usv_colors{j}, ...
            'MarkerSize', 5, 'LineWidth', 1.2, 'DisplayName', sprintf('USV %d', j));
    end
end
title('通信拓扑节点 (USV 1, 4, 8) 事件触发通信间隔时间分布');
xlabel('Time (s)'); ylabel('Inter-event Time (s)');
legend('Location', 'northeast');

if exist('inter_event_times', 'var')
    ylim([0, max(inter_event_times)*1.2]); 
end

% --- 局部放大图 (Inset)：180s - 240s ---
axes('Position', [0.35, 0.55, 0.35, 0.3]); 
hold on; grid on; box on;

for idx = 1:length(target_usvs)
    j = target_usvs(idx);
    comm_times = usv_comm_state{j}.comm_events;
    if length(comm_times) > 1
        inter_event_times = diff(comm_times);
        zoom_idx = find(comm_times(2:end) > 180 & comm_times(2:end) < 240);
        if ~isempty(zoom_idx)
            h_zoom = stem(comm_times(zoom_idx+1), inter_event_times(zoom_idx));
            set(h_zoom, 'Color', usv_colors{j}, 'Marker', markers{idx}, ...
                'MarkerFaceColor', usv_colors{j}, 'MarkerEdgeColor', usv_colors{j}, ...
                'MarkerSize', 4, 'LineWidth', 1.2);
        end
    end
end
xlim([180, 240]);
title('180s-240s 局部放大', 'FontSize', 9);

% ========================= 航线生成核心引擎 =========================
function plan = build_trajectory_plan(T_max, v_c, R_c, start_pos, start_psi)
    % 基于状态机的完美S形对齐算法
    T_turn = (R_c * pi / 2) / v_c;  % 90度圆弧转弯时间
    
    % 标准动作块: [持续时间, 动作类型, 起始间距比, 结束间距比]
    % 动作类型: 0=直行, 1=左转90度, -1=右转90度
    block = [
        T_turn,  1, 1, 1; % 左转90 (横行)
        5,       0, 1, 1; % 保持收缩直行 5s
        T_turn,  1, 1, 1; % 左转90 (完成U-turn, 反向)
        5,       0, 1, 1; % 保持收缩直行 5s
        27,      0, 1, 0; % 展开队形 27s
        100,     0, 0, 0; % 正常间距长直航线 100s
        27,      0, 0, 1; % 收缩队形 27s (准备下次转弯)
        T_turn, -1, 1, 1; % 右转90 (横行)
        5,       0, 1, 1; % 保持收缩直行 5s
        T_turn, -1, 1, 1; % 右转90 (完成U-turn, 正向)
        5,       0, 1, 1; % 保持收缩直行 5s
        27,      0, 1, 0; % 展开队形 27s
        100,     0, 0, 0; % 正常间距长直航线 100s
        27,      0, 0, 1  % 收缩队形 27s (准备下次转弯)
    ];

    % 第一条长边的初段 (总长要确保与后面的长直段对齐 159s = 5+27+100+27)
    seq = [
        132, 0, 0, 0; % 起步正常间距直行 (替代了5+27+100)
        27,  0, 0, 1  % 准备第一次左转的收缩
    ];

    % 循环拼接直到塞满仿真时间
    while sum(seq(:,1)) < T_max
        seq = [seq; block]; %#ok<AGROW> 
    end

    num_segments = size(seq, 1);
    plan(num_segments).t_start = 0; % 预分配
    current_t = 0; current_pos = start_pos; current_psi = start_psi;

    for idx = 1:num_segments
        dur = seq(idx, 1); type = seq(idx, 2);
        
        plan(idx).t_start = current_t; plan(idx).t_end = current_t + dur;
        plan(idx).duration = dur;      plan(idx).type = type;
        plan(idx).p_start = seq(idx, 3); plan(idx).p_end = seq(idx, 4);
        plan(idx).x0 = current_pos(1); plan(idx).y0 = current_pos(2);
        plan(idx).z0 = current_pos(3); plan(idx).psi0 = current_psi;

        % 步进推算终点几何位置
        if type == 0
            current_pos(1) = current_pos(1) + v_c * dur * cos(current_psi);
            current_pos(2) = current_pos(2) + v_c * dur * sin(current_psi);
        elseif type == 1
            cx = current_pos(1) - R_c * sin(current_psi);
            cy = current_pos(2) + R_c * cos(current_psi);
            current_psi = current_psi + pi/2;
            current_pos(1) = cx + R_c * sin(current_psi);
            current_pos(2) = cy - R_c * cos(current_psi);
        elseif type == -1
            cx = current_pos(1) + R_c * sin(current_psi);
            cy = current_pos(2) - R_c * cos(current_psi);
            current_psi = current_psi - pi/2;
            current_pos(1) = cx - R_c * sin(current_psi);
            current_pos(2) = cy + R_c * cos(current_psi);
        end
        current_t = current_t + dur;
    end
end

function [pos_d, vel_d, p] = get_trajectory_state(t, plan, v_c, R_c)
    % 解析任意时刻 t 的绝对物理位置与队形缩放比 p
    idx = find(t >= [plan.t_start] & t < [plan.t_end], 1, 'last');
    if isempty(idx)
        idx = length(plan); t = plan(idx).t_end; 
    end
    seg = plan(idx); dt = t - seg.t_start;

    if seg.p_start == seg.p_end
        p = seg.p_start;
    else
        tau = dt / seg.duration;
        p_trans = 6*tau^5 - 15*tau^4 + 10*tau^3; % 5阶平滑多项式
        p = seg.p_start + (seg.p_end - seg.p_start) * p_trans;
    end

    if seg.type == 0
        pos_d = [seg.x0 + v_c * dt * cos(seg.psi0);
                 seg.y0 + v_c * dt * sin(seg.psi0); seg.z0];
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