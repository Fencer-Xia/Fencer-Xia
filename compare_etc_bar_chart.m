% FILE: compare_etc_bar_chart.m
% =========================================================================
% 独立模块：对比 连续通信、静态事件触发(SETC)、动态事件触发(NDETC)
% 依赖：主程序运行完毕后工作区中的 uav_state_history, time_vec 等变量
% =========================================================================
clc;
fprintf('正在后台分析三种通信机制的触发次数...\n');

target_usvs = [1, 4, 8];
num_targets = length(target_usvs);

% 初始化统计数组 [连续, 静态, 动态]
trigger_counts = zeros(num_targets, 3);

% 静态触发的固定阈值 (取一个合理的经验值，例如 0.15m)
static_threshold = 0.1; 

for idx = 1:num_targets
    j = target_usvs(idx);
    if j <= 3; leader_uav_idx = 1; elseif j <= 5; leader_uav_idx = 2; else; leader_uav_idx = 3; end
    
    % 1. 连续通信 (Continuous)
    % 只要系统在运行，每个 dt 都算一次通信
    trigger_counts(idx, 1) = length(time_vec) - 1;
    
    % 2. 动态事件触发 (Dynamic ETC - 本文算法)
    % 直接读取主程序已经算好的触发次数
    trigger_counts(idx, 3) = length(usv_comm_state{j}.comm_events);
    
    % 3. 静态事件触发 (Static ETC)
    % 我们在后台快速“回放”一遍静态逻辑
    static_trigger_count = 0;
    p_last = uav_state_history(1:2, 1, leader_uav_idx);
    v_last = [2.5; 0];
    t_last = 0;
    
    for i = 1:(length(time_vec)-1)
        t = time_vec(i);
        p_real = uav_state_history(1:2, i, leader_uav_idx);
        
        % 获取真实速度用于更新预测器
        psi_uav = uav_state_history(6, i, leader_uav_idx);
        u_uav = uav_state_history(7, i, leader_uav_idx);
        v_uav = uav_state_history(8, i, leader_uav_idx);
        v_real = [u_uav * cos(psi_uav) - v_uav * sin(psi_uav);
                  u_uav * sin(psi_uav) + v_uav * cos(psi_uav)];
              
        % 零阶预测
        p_pred = p_last + v_last * (t - t_last);
        error_norm = norm(p_real - p_pred);
        
        % 静态触发判断
        if error_norm >= static_threshold
            static_trigger_count = static_trigger_count + 1;
            p_last = p_real;
            v_last = v_real;
            t_last = t;
        end
    end
    trigger_counts(idx, 2) = static_trigger_count;
end

% ========================= 绘制精美柱状图 =========================
figure('Name', '不同通信机制触发次数对比', 'Color', 'w', 'Position', [400, 200, 800, 500]);
bar_handle = bar(trigger_counts, 'grouped');

% 设置颜色：连续(灰), 静态(橙), 动态(蓝)
bar_handle(1).FaceColor = [0.7 0.7 0.7]; 
bar_handle(2).FaceColor = [0.9 0.5 0.2]; 
bar_handle(3).FaceColor = [0.2 0.5 0.8]; 

set(gca, 'XTickLabel', {'USV 1 (左翼)', 'USV 4 (中心)', 'USV 8 (右翼)'}, 'FontSize', 12);
ylabel('累积通信触发次数 (次)', 'FontSize', 12);
title('通信机制效能对比：时间触发 vs 静态事件触发 vs 动态事件触发', 'FontSize', 14);
legend('时间触发 (Time-Triggered)', '静态事件触发 (SETC)', '动态事件触发 (NDETC, 本文)', 'Location', 'northwest');
grid on;
set(gca, 'YScale', 'log'); % 使用对数坐标轴，因为连续通信(几万次)远大于触发通信(几百次)

% 在柱子上标注具体数值
for i = 1:3
    xtips1 = bar_handle(i).XEndPoints;
    ytips1 = bar_handle(i).YEndPoints;
    labels1 = string(bar_handle(i).YData);
    text(xtips1, ytips1, labels1, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 10, 'FontWeight', 'bold');
end

fprintf('分析完毕！请查看柱状图对比结果。\n');