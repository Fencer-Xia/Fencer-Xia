% =========================================================================
% 消融实验主程序：对比 平滑前(Ori-AFBC) 与 平滑后(Sm-AFBC) 的控制性能
% 核心评估指标：控制力矩抖振方差 (Chattering Variance) 与 姿态角响应平滑度
% 注意：未对任何仿真环境（初始误差、扰动、S形轨迹、DETC通信）进行简化
% =========================================================================
clc; clear; close all;
fprintf('开始执行【平滑算法消融实验：未平滑 vs 平滑后】...\n');

% 算法配置：
% Algo 1 (未平滑): UAV_position_sm, UAV_attitude_s, control (USV无后缀)
% Algo 2 (平滑后): UAV_position (无后缀), UAV_attitude (无后缀), control_s (USV带后缀)
algo_list = {'Ori-AFBC (未平滑)', 'Sm-AFBC (本文平滑)'};
colors = {'#D95319', '#0072BD'}; 

% 仿真参数 (完全一致)
dt = 0.01; T_sim = 400; time_vec = 0:dt:T_sim; N = length(time_vec);
num_uavs = 3; uav_spacing_default = 25; 
num_usvs = 8; usv_spacing_default = 10;
maneuver.v_c = 2.5; maneuver.R_c = 65; maneuver.usv_v_max = 3.0;

usv_outermost_factor = (num_usvs/2 - 0.5);
s_usv_new_safe = (maneuver.R_c / usv_outermost_factor) * (maneuver.usv_v_max / maneuver.v_c - 1);
s_usv_new = max(3.0, s_usv_new_safe);
spacing_ratio = uav_spacing_default / usv_spacing_default;
s_uav_new = s_usv_new * spacing_ratio;

Results = struct();

for algo_idx = 1:2
    fprintf('正在运行消融实验 [%d/2]: %s ...\n', algo_idx, algo_list{algo_idx});
    
    % 强制重置随机种子，保证两次仿真的初始位置与扰动序列100%一致
    rng(0); 
    uav_state_history = zeros(12, N, num_uavs);
    uav_l_hat = repmat([0.1; 0.1; 0.1], 1, num_uavs);
    uav_l_hat_att = repmat([0.1; 0.01; 0.01], 1, num_uavs); % 必须重置权重
    comm_params.rho_max = 0.2; 
    comm_params.rho_min = 0.1; 
    comm_params.T_relax = 8.0; 
    comm_params.k_accel = 0.1; 
    comm_params.wakeup_time = 0.2;
    uav_att_state = cell(1, num_uavs);
    for k = 1:num_uavs
        % 重置所有状态机组件
        uav_att_state{k}.Theta_d = [0; 0; 0]; % 或初始姿态
        uav_att_state{k}.prev_Theta_d = [0; 0; 0];
        uav_att_state{k}.Theta_d_dot = [0; 0; 0];
        uav_att_state{k}.Omega_d = [0; 0; 0];
        uav_att_state{k}.Omega_d_dot = [0; 0; 0];
        
        % 恢复初始位置 (如果之前是在循环外随机生成的，请在此处重新生成)
        % y_base = (k - uav_center_idx) * uav_spacing_default;
        % uav_state_history(1:3, 1, k) = ...
    end
    usv_comm_state = cell(1, num_usvs);
        for j = 1:num_usvs
            leader_uav_idx = 1; if j>3; leader_uav_idx=2; end; if j>5; leader_uav_idx=3; end
            p_uav_init = uav_state_history(1:2, 1, leader_uav_idx);
            usv_comm_state{j}.p_uav_last = p_uav_init; 
            usv_comm_state{j}.v_uav_last = [maneuver.v_c; 0];
            usv_comm_state{j}.psi_uav_last = 0; 
            usv_comm_state{j}.r_uav_last = 0; 
            usv_comm_state{j}.t_last = 0;
            usv_comm_state{j}.rho = comm_params.rho_max;
            usv_comm_state{j}.p_smooth = p_uav_init; 
            usv_comm_state{j}.v_smooth = [maneuver.v_c; 0];
            usv_comm_state{j}.psi_smooth = 0; 
            usv_comm_state{j}.r_smooth = 0;
            usv_comm_state{j}.comm_events = []; % 必须重置触发记录
        end
    % --- 物理与控制参数 (保持不变) ---
    uav_params.m = 1.5; uav_params.g = 9.81;
    uav_params.Ixx = 0.03; uav_params.Iyy = 0.03; uav_params.Izz = 0.06;
    uav_params.max_thrust = 2 * uav_params.m * uav_params.g;
    uav_params.K1 = diag([20, 20, 30]); uav_params.K2 = diag([16, 16, 12]);
    uav_params.beta = 1.2; uav_params.Gamma = [0.01; 0.01; 0.01];
    uav_params.K3_att = diag([10.0, 10.0, 15.0]); uav_params.K4_att = diag([1.5, 1.5, 2.0]);
    uav_params.K_frac_att = diag([0.5, 0.5, 0.5]); uav_params.beta_att = 0.8;
    uav_params.Gamma_att = diag([1e-3, 1e-3, 1e-3]); uav_params.sigma_att = [0.01; 0.01; 0.01];
    usv_params.gamma_afbc = [0.5; 0.05; 0.02]; usv_params.l_hat_max = [10.0; 1.0; 0.5];

    % --- 状态初始化 (保持不变) ---
    
    
    uav_center_idx = (num_uavs + 1) / 2;
    x_range = [-2, 2]; z_range = [10, 15];
    for k = 1:num_uavs
        y_base = (k - uav_center_idx) * uav_spacing_default;
        uav_state_history(1:3, 1, k) = [x_range(1)+(x_range(2)-x_range(1))*rand(); y_base; z_range(1)+(z_range(2)-z_range(1))*rand()];
        uav_state_history(4:12, 1, k) = 0;
    end
    usv_state.x = (rand(1, num_usvs) - 0.5) * 20;
    usv_state.psi = (rand(1, num_usvs) - 0.5) * (pi/2);
    usv_center_idx = (num_usvs + 1) / 2;
    for j = 1:num_usvs; usv_state.y(j) = (j - usv_center_idx) * usv_spacing_default + (rand - 0.5) * 10; end
    usv_state.u = ones(1, num_usvs) * maneuver.v_c; usv_state.r = zeros(1, num_usvs);
    usv_state.l_hat = repmat([0.1; 0.01; 0.01], 1, num_usvs);
    
    usv_integral_error = zeros(2, num_usvs);

    % --- 航线与通信参数 ---
    plan = build_trajectory_plan(T_sim, maneuver.v_c, maneuver.R_c, [10; 0; 20], 0);
    comm_params.rho_max = 0.2; comm_params.rho_min = 0.1; comm_params.T_relax = 8.0; comm_params.k_accel = 0.1; comm_params.wakeup_time = 0.2;
    smooth_params.Kp = 25.0; smooth_params.Kd = 10.0;
    uav_desired_accel = zeros(3, num_uavs);

    % --- 关键记录数组 (UAV 1 与 USV 8) ---
    hist_uav1_U = zeros(4, N);
    hist_uav1_ang = zeros(3, N);
    hist_usv8_tau = zeros(2, N);
    hist_usv8_ang = zeros(1, N);

    % --- 主循环 ---
    for i = 1:N-1
        t = time_vec(i);
        wind_force_x = 120 + 70*sin(0.4*t) + 30*cos(1.5*t);
        wind_force_y = 40 + 50*cos(0.5*t) + 5*sin(1.8*t);
        
        [pos_d_leader, vel_d_leader, p_ratio] = get_trajectory_state(t, plan, maneuver.v_c, maneuver.R_c);
        desired_state_leader.pos = pos_d_leader; desired_state_leader.vel = vel_d_leader; desired_state_leader.acc = [0;0;0];
        target_spacing_uav = uav_spacing_default + (s_uav_new - uav_spacing_default) * p_ratio;
        target_spacing_usv = usv_spacing_default + (s_usv_new - usv_spacing_default) * p_ratio;
        psi_form = atan2(vel_d_leader(2), vel_d_leader(1));
        R_z_form = [cos(psi_form), -sin(psi_form); sin(psi_form), cos(psi_form)];

        % UAV 更新
        for k = 1:num_uavs
            curr_uav = uav_state_history(:, i, k);
            offset_i = R_z_form * [0; (k - uav_center_idx) * target_spacing_uav];
            des_uav = desired_state_leader; des_uav.pos(1:2) = des_uav.pos(1:2) + offset_i;
            
            % == 核心调用逻辑 ==
            if algo_idx == 1 % 未平滑
                [U1, des_a, l_dot] = UAV_position_controller_AFBC_sm(curr_uav, des_uav, uav_l_hat(:, k), uav_params);
                [U2, U3, U4, latt_dot, att_out] = UAV_attitude_controller_AFBC_s(curr_uav, des_a, uav_l_hat_att(:,k), uav_att_state{k}, dt, uav_params);
            else % 平滑后
                [U1, des_a, l_dot] = UAV_position_controller_AFBC(curr_uav, des_uav, uav_l_hat(:, k), uav_params);
                [U2, U3, U4, latt_dot, att_out] = UAV_attitude_controller_AFBC(curr_uav, des_a, uav_l_hat_att(:,k), uav_att_state{k}, dt, uav_params);
            end
            
            uav_desired_accel(:, k) = des_a;
            uav_l_hat_att(:, k) = max(0, uav_l_hat_att(:, k) + latt_dot * dt);
            uav_att_state{k} = att_out;
            uav_l_hat(:, k) = max(0, uav_l_hat(:, k) + l_dot * dt);
            
            U = [U1; max(-5,min(U2,5)); max(-5,min(U3,5)); max(-2,min(U4,2))];
            U(1) = max(0, min(U(1), uav_params.max_thrust));
            
            if k == 1
                hist_uav1_U(:, i) = U;
                hist_uav1_ang(:, i) = curr_uav(4:6);
            end
            
            s_dot = UAV_dynamic_function(curr_uav, U);
            s_dot(7:9) = s_dot(7:9) + ([wind_force_x; wind_force_y; 0]*0.005) / uav_params.m;
            uav_state_history(:, i+1, k) = curr_uav + s_dot * dt;
        end

        % USV 更新 (包含DETC逻辑)
        is_wakeup = any(abs(t - ([plan.t_start] - comm_params.wakeup_time)) < dt/2);
        for j = 1:num_usvs
            leader_uav_idx = 1; if j>3; leader_uav_idx=2; end; if j>5; leader_uav_idx=3; end
            curr_uav = uav_state_history(:, i, leader_uav_idx);
            p_real = curr_uav(1:2); psi_real = curr_uav(6); r_real = curr_uav(12);
            v_real = [curr_uav(7)*cos(psi_real)-curr_uav(8)*sin(psi_real); curr_uav(7)*sin(psi_real)+curr_uav(8)*cos(psi_real)];
            
            p_pred = usv_comm_state{j}.p_uav_last + usv_comm_state{j}.v_uav_last*(t - usv_comm_state{j}.t_last);
            rho_dot = (comm_params.rho_max - usv_comm_state{j}.rho)/comm_params.T_relax - comm_params.k_accel*norm(uav_desired_accel(1:2, leader_uav_idx));
            rho_new = max(comm_params.rho_min, min(comm_params.rho_max, usv_comm_state{j}.rho + rho_dot*dt));
            usv_comm_state{j}.rho = rho_new;
            
            if is_wakeup || norm(p_real - p_pred) >= rho_new
                usv_comm_state{j}.p_uav_last = p_real; usv_comm_state{j}.v_uav_last = v_real;
                usv_comm_state{j}.psi_uav_last = psi_real; usv_comm_state{j}.r_uav_last = r_real;
                usv_comm_state{j}.t_last = t; p_pred = p_real;
            end
            
            v_pred = usv_comm_state{j}.v_uav_last;
            a_sm = smooth_params.Kp*(p_pred - usv_comm_state{j}.p_smooth) + smooth_params.Kd*(v_pred - usv_comm_state{j}.v_smooth);
            usv_comm_state{j}.v_smooth = usv_comm_state{j}.v_smooth + a_sm*dt;
            usv_comm_state{j}.p_smooth = usv_comm_state{j}.p_smooth + usv_comm_state{j}.v_smooth*dt;
            p_leader_for_usv = usv_comm_state{j}.p_smooth;
            
            [~, vd_curr, ~] = get_trajectory_state(t, plan, maneuver.v_c, maneuver.R_c);
            [~, vd_next, next_p] = get_trajectory_state(t+dt, plan, maneuver.v_c, maneuver.R_c);
            psi_f_curr = atan2(vd_curr(2), vd_curr(1));
            r_form = wrapToPi(atan2(vd_next(2), vd_next(1)) - psi_f_curr) / dt;
            
            if j<=3; sub_idx=j-2; elseif j<=5; sub_idx=j-4.5; else; sub_idx=j-7; end
            pos_d_usv = p_leader_for_usv + [cos(psi_f_curr), -sin(psi_f_curr); sin(psi_f_curr), cos(psi_f_curr)] * [0; sub_idx * target_spacing_usv];
            
            xb=usv_state.x(j); yb=usv_state.y(j); psib=usv_state.psi(j);
            err_ix = pos_d_usv(1)-xb; err_iy = pos_d_usv(2)-yb;
            ye = -err_ix*sin(psib) + err_iy*cos(psib);
            z1 = [err_ix*cos(psib) + err_iy*sin(psib); ye; wrapToPi(psi_f_curr - psib)];
            if abs(ye) < 2.0; usv_integral_error(1, j) = usv_integral_error(1, j) + ye*dt; end
            dyn_s = [usv_state.u(j); 0; 0; 0; usv_state.r(j)];
            speed_d = maneuver.v_c - r_form*(sub_idx * target_spacing_usv);
            
            % == 核心调用逻辑 ==
            if algo_idx == 1 % 未平滑
                [tau, z2] = control(dyn_s, z1, usv_state.l_hat(:, j), usv_integral_error(:, j), speed_d, r_form);
            else % 平滑后
                [tau, z2] = control_s(dyn_s, z1, usv_state.l_hat(:, j), usv_integral_error(:, j), speed_d, r_form);
            end
            
            if j == 8
                hist_usv8_tau(:, i) = tau([1,3]);
                hist_usv8_ang(i) = psib;
            end
            
            dis_usv = [120+70*sin(0.4*t)+30*cos(1.5*t); 40+50*cos(0.5*t)+5*sin(1.8*t); 0; 0; 50+50*sin(0.3*t)+15*cos(1.2*t)];
            acc = dynamic_function(tau, dyn_s, dis_usv, [xb; yb; 0; 0]);
            usv_state.u(j) = usv_state.u(j) + acc(1)*dt; usv_state.r(j) = usv_state.r(j) + acc(5)*dt;
            usv_state.psi(j) = wrapToPi(psib + usv_state.r(j)*dt);
            usv_state.x(j) = xb + usv_state.u(j)*cos(usv_state.psi(j))*dt; usv_state.y(j) = yb + usv_state.u(j)*sin(usv_state.psi(j))*dt;
            
            nv = norm([usv_state.u(j); usv_state.r(j)]);
            l_hat_dot = [usv_params.gamma_afbc(1)*norm(z2); usv_params.gamma_afbc(2)*norm(z2)*nv; usv_params.gamma_afbc(3)*norm(z2)*(nv^2)];
            usv_state.l_hat(:, j) = min(max(usv_state.l_hat(:, j) + l_hat_dot*dt, 1e-6), usv_params.l_hat_max);
        end
    end
    Results.(['A' num2str(algo_idx)]).U = hist_uav1_U;
    Results.(['A' num2str(algo_idx)]).ang_uav = hist_uav1_ang;
    Results.(['A' num2str(algo_idx)]).tau = hist_usv8_tau;
    Results.(['A' num2str(algo_idx)]).ang_usv = hist_usv8_ang;
end

% ========================= 抖振总变差 (TV) 分析 =========================
% 剔除初始大收敛(前20s)，统计稳态及转弯抗扰阶段的真实总变差
eval_idx = (20/dt):N-2; 

% 计算公式：TV = sum(abs(diff(signal)))
tv_U2_ori = sum(abs(diff(Results.A1.U(2, eval_idx)))); tv_U2_sm = sum(abs(diff(Results.A2.U(2, eval_idx))));
tv_U3_ori = sum(abs(diff(Results.A1.U(3, eval_idx)))); tv_U3_sm = sum(abs(diff(Results.A2.U(3, eval_idx))));
tv_U4_ori = sum(abs(diff(Results.A1.U(4, eval_idx)))); tv_U4_sm = sum(abs(diff(Results.A2.U(4, eval_idx))));

tv_phi_ori = sum(abs(diff(Results.A1.ang_uav(1, eval_idx)))); tv_phi_sm = sum(abs(diff(Results.A2.ang_uav(1, eval_idx))));
tv_the_ori = sum(abs(diff(Results.A1.ang_uav(2, eval_idx)))); tv_the_sm = sum(abs(diff(Results.A2.ang_uav(2, eval_idx))));
tv_psi_uav_ori = sum(abs(diff(Results.A1.ang_uav(3, eval_idx)))); tv_psi_uav_sm = sum(abs(diff(Results.A2.ang_uav(3, eval_idx)))); % 修复：补充了UAV偏航角计算

tv_taur_ori = sum(abs(diff(Results.A1.tau(2, eval_idx)))); tv_taur_sm = sum(abs(diff(Results.A2.tau(2, eval_idx))));
tv_psi_usv_ori = sum(abs(diff(Results.A1.ang_usv(eval_idx)))); tv_psi_usv_sm = sum(abs(diff(Results.A2.ang_usv(eval_idx))));

fprintf('\n=== 平滑指标量化对比 (总变差 TV) ===\n');
fprintf('指标 \t\t\t Ori-AFBC(未平滑) \t Sm-AFBC(平滑后) \n');
fprintf('UAV U2 TV: \t\t %.4f \t\t\t %.4f\n', tv_U2_ori, tv_U2_sm);
fprintf('UAV U3 TV: \t\t %.4f \t\t\t %.4f\n', tv_U3_ori, tv_U3_sm);
fprintf('UAV U4 TV: \t\t %.4f \t\t\t %.4f\n', tv_U4_ori, tv_U4_sm);
fprintf('UAV Phi角度 TV: \t %.4f \t\t\t %.4f\n', tv_phi_ori, tv_phi_sm);
fprintf('UAV Psi角度 TV: \t %.4f \t\t\t %.4f\n', tv_psi_uav_ori, tv_psi_uav_sm);
fprintf('USV Tau_r TV: \t\t %.4f \t\t\t %.4f\n', tv_taur_ori, tv_taur_sm);
fprintf('USV Psi角度 TV: \t %.4f \t\t\t %.4f\n', tv_psi_usv_ori, tv_psi_usv_sm);

% ========================= 对比绘图 =========================
% 重新布局：4行2列，移除TV描述，并在姿态角(右侧列)追加局部放大图
figure('Name', '力矩控制抖振与角度平滑度对比', 'Color', 'w', 'Position', [100, 50, 1200, 1000]);
t_plt = time_vec(1:end-1); 
x_lim = [50 200]; % 主图显示区间
zoom_lim = [159 160]; % 局部放大图显示区间

% --- 第1行：UAV 滚转通道 ---
% 1. UAV 滚转力矩 U_2 (左)
ax1 = subplot(4, 2, 1); hold on; grid on; box on;
plot(t_plt, Results.A1.U(2, 1:end-1), 'Color', colors{1}, 'LineWidth'  , 1.5); 
plot(t_plt, Results.A2.U(2, 1:end-1), 'Color', colors{2}, 'LineWidth', 1.5);
title('UAV 1 滚转力矩 U_2'); ylabel('U_2 (N·m)'); xlim(x_lim); 
legend(algo_list, 'Location', 'best');

% 2. UAV 滚转角 \phi (右)
ax2 = subplot(4, 2, 2); hold on; grid on; box on;
plot(t_plt, rad2deg(Results.A1.ang_uav(1, 1:end-1)), 'Color', colors{1}, 'LineWidth', 1.5); 
plot(t_plt, rad2deg(Results.A2.ang_uav(1, 1:end-1)), 'Color', colors{2}, 'LineWidth', 1.5);
title('UAV 1 滚转角 \phi'); ylabel('\phi (deg)'); xlim(x_lim);
% -- 滚转角局部放大 --
axes('Position', [0.6, 0.83, 0.1, 0.08]); hold on; grid on; box on;
plot(t_plt, rad2deg(Results.A1.ang_uav(1, 1:end-1)), 'Color', colors{1}, 'LineWidth', 1.5); 
plot(t_plt, rad2deg(Results.A2.ang_uav(1, 1:end-1)), 'Color', colors{2}, 'LineWidth', 1.5);
xlim(zoom_lim); title(' ', 'FontSize', 8); set(gca, 'FontSize', 8);

% --- 第2行：UAV 俯仰通道 ---
% 3. UAV 俯仰力矩 U_3 (左)
ax3 = subplot(4, 2, 3); hold on; grid on; box on;
plot(t_plt, Results.A1.U(3, 1:end-1), 'Color', colors{1}, 'LineWidth', 1.5); 
plot(t_plt, Results.A2.U(3, 1:end-1), 'Color', colors{2}, 'LineWidth', 1.5);
title('UAV 1 俯仰力矩 U_3'); ylabel('U_3 (N·m)'); xlim(x_lim);

% 4. UAV 俯仰角 \theta (右)
ax4 = subplot(4, 2, 4); hold on; grid on; box on;
plot(t_plt, rad2deg(Results.A1.ang_uav(2, 1:end-1)), 'Color', colors{1}, 'LineWidth', 1.5); 
plot(t_plt, rad2deg(Results.A2.ang_uav(2, 1:end-1)), 'Color', colors{2}, 'LineWidth', 1.5);
title('UAV 1 俯仰角 \theta'); ylabel('\theta (deg)'); xlim(x_lim);
% -- 俯仰角局部放大 --
axes('Position', [0.6, 0.62, 0.1, 0.08]); hold on; grid on; box on;
plot(t_plt, rad2deg(Results.A1.ang_uav(2, 1:end-1)), 'Color', colors{1}, 'LineWidth', 1.5); 
plot(t_plt, rad2deg(Results.A2.ang_uav(2, 1:end-1)), 'Color', colors{2}, 'LineWidth', 1.5);
xlim(zoom_lim); title(' ', 'FontSize', 8); set(gca, 'FontSize', 8);

% --- 第3行：UAV 偏航通道 ---
% 5. UAV 偏航力矩 U_4 (左)
ax5 = subplot(4, 2, 5); hold on; grid on; box on;
plot(t_plt, Results.A1.U(4, 1:end-1), 'Color', colors{1}, 'LineWidth', 1.5); 
plot(t_plt, Results.A2.U(4, 1:end-1), 'Color', colors{2}, 'LineWidth', 1.5);
title('UAV 1 偏航力矩 U_4'); ylabel('U_4 (N·m)'); xlim(x_lim);

% 6. UAV 偏航角 \psi_u (右)
ax6 = subplot(4, 2, 6); hold on; grid on; box on;
plot(t_plt, rad2deg(unwrap(Results.A1.ang_uav(3, 1:end-1))), 'Color', colors{1}, 'LineWidth', 1.5); 
plot(t_plt, rad2deg(unwrap(Results.A2.ang_uav(3, 1:end-1))), 'Color', colors{2}, 'LineWidth', 1.5);
title('UAV 1 偏航角 \psi_u'); ylabel('\psi_u (deg)'); xlim(x_lim);
% -- UAV 偏航角局部放大 --
axes('Position', [0.6, 0.38, 0.1, 0.08]); hold on; grid on; box on;
plot(t_plt, rad2deg(unwrap(Results.A1.ang_uav(3, 1:end-1))), 'Color', colors{1}, 'LineWidth', 1.5); 
plot(t_plt, rad2deg(unwrap(Results.A2.ang_uav(3, 1:end-1))), 'Color', colors{2}, 'LineWidth', 1.5);
xlim(zoom_lim); title(' ', 'FontSize', 8); set(gca, 'FontSize', 8);

% --- 第4行：USV 偏航通道 ---
% 7. USV 偏航力矩 \tau_r (左)
ax7 = subplot(4, 2, 7); hold on; grid on; box on;
plot(t_plt, Results.A1.tau(2, 1:end-1), 'Color', colors{1}, 'LineWidth', 1.5); 
plot(t_plt, Results.A2.tau(2, 1:end-1), 'Color', colors{2}, 'LineWidth', 1.5);
title('USV 8 偏航力矩 \tau_r'); xlabel('Time (s)'); ylabel('\tau_r (N·m)'); xlim(x_lim);

% 8. USV 艏向角 \psi (右)
ax8 = subplot(4, 2, 8); hold on; grid on; box on;
plot(t_plt, rad2deg(unwrap(Results.A1.ang_usv(1:end-1))), 'Color', colors{1}, 'LineWidth', 1.5); 
plot(t_plt, rad2deg(unwrap(Results.A2.ang_usv(1:end-1))), 'Color', colors{2}, 'LineWidth', 1.5);
title('USV 8 艏向角 \psi'); xlabel('Time (s)'); ylabel('\psi (deg)'); xlim(x_lim);
% -- USV 偏航角局部放大 --
axes('Position', [0.6, 0.17, 0.1, 0.08]); hold on; grid on; box on;
plot(t_plt, rad2deg(unwrap(Results.A1.ang_usv(1:end-1))), 'Color', colors{1}, 'LineWidth', 1.5); 
plot(t_plt, rad2deg(unwrap(Results.A2.ang_usv(1:end-1))), 'Color', colors{2}, 'LineWidth', 1.5);
xlim(zoom_lim); title(' ', 'FontSize', 8); set(gca, 'FontSize', 8);
% 依赖函数
function plan = build_trajectory_plan(T_max, v_c, R_c, start_pos, start_psi)
    T_turn = (R_c * pi / 2) / v_c;  
    block = [T_turn, 1, 1, 1; 5, 0, 1, 1; T_turn, 1, 1, 1; 5, 0, 1, 1; 27, 0, 1, 0; 100, 0, 0, 0; 27, 0, 0, 1; 
             T_turn, -1, 1, 1; 5, 0, 1, 1; T_turn, -1, 1, 1; 5, 0, 1, 1; 27, 0, 1, 0; 100, 0, 0, 0; 27, 0, 0, 1];
    seq = [132, 0, 0, 0; 27, 0, 0, 1];
    while sum(seq(:,1)) < T_max; seq = [seq; block]; end
    num_segments = size(seq, 1); plan(num_segments).t_start = 0; current_t = 0; current_pos = start_pos; current_psi = start_psi;
    for idx = 1:num_segments
        dur = seq(idx, 1); type = seq(idx, 2); plan(idx).t_start = current_t; plan(idx).t_end = current_t + dur;
        plan(idx).duration = dur; plan(idx).type = type; plan(idx).p_start = seq(idx, 3); plan(idx).p_end = seq(idx, 4);
        plan(idx).x0 = current_pos(1); plan(idx).y0 = current_pos(2); plan(idx).z0 = current_pos(3); plan(idx).psi0 = current_psi;
        if type == 0; current_pos(1) = current_pos(1) + v_c*dur*cos(current_psi); current_pos(2) = current_pos(2) + v_c*dur*sin(current_psi);
        elseif type == 1; cx = current_pos(1)-R_c*sin(current_psi); cy = current_pos(2)+R_c*cos(current_psi); current_psi = current_psi+pi/2; current_pos(1) = cx+R_c*sin(current_psi); current_pos(2) = cy-R_c*cos(current_psi);
        elseif type == -1; cx = current_pos(1)+R_c*sin(current_psi); cy = current_pos(2)-R_c*cos(current_psi); current_psi = current_psi-pi/2; current_pos(1) = cx-R_c*sin(current_psi); current_pos(2) = cy+R_c*cos(current_psi); end
        current_t = current_t + dur;
    end
end

function [pos_d, vel_d, p] = get_trajectory_state(t, plan, v_c, R_c)
    idx = find(t >= [plan.t_start] & t < [plan.t_end], 1, 'last');
    if isempty(idx); idx = length(plan); t = plan(idx).t_end; end
    seg = plan(idx); dt = t - seg.t_start;
    if seg.p_start == seg.p_end; p = seg.p_start; else; tau = dt / seg.duration; p = seg.p_start + (seg.p_end - seg.p_start) * (6*tau^5 - 15*tau^4 + 10*tau^3); end
    if seg.type == 0
        pos_d = [seg.x0 + v_c*dt*cos(seg.psi0); seg.y0 + v_c*dt*sin(seg.psi0); seg.z0]; vel_d = [v_c*cos(seg.psi0); v_c*sin(seg.psi0); 0];
    elseif seg.type == 1
        cx = seg.x0 - R_c*sin(seg.psi0); cy = seg.y0 + R_c*cos(seg.psi0); psi_t = seg.psi0 + (v_c/R_c)*dt;
        pos_d = [cx + R_c*sin(psi_t); cy - R_c*cos(psi_t); seg.z0]; vel_d = [v_c*cos(psi_t); v_c*sin(psi_t); 0];
    elseif seg.type == -1
        cx = seg.x0 + R_c*sin(seg.psi0); cy = seg.y0 - R_c*cos(seg.psi0); psi_t = seg.psi0 - (v_c/R_c)*dt;
        pos_d = [cx - R_c*sin(psi_t); cy + R_c*cos(psi_t); seg.z0]; vel_d = [v_c*cos(psi_t); v_c*sin(psi_t); 0];
    end
end