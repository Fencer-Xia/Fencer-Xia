% FILE: control.m
% 增加输入参数 r_form
function [output, z2] = control(dynamic, z1, l_hat, integral_err_kin, desired_speed, r_form)
% =========================================================================
%               纯粹的 USV 自适应固定时间控制器 (Pure AFBC)
%
% 【V9 最终性能调优版】:
%   - 核心目标: 在V8版本稳定无震荡的基础上，解决瞬时误差过大和
%     稳态误差消除过慢的问题。
%   - 解决方案: 对运动学PID控制器的三项增益进行平衡性重调。
%     1. 提高 P(kp_ye) 和 D(kd_psi) 增益，加快瞬态响应速度，减小误差峰值。
%     2. 提高 I(ki_ye) 增益，增强消除稳态误差的能力。
%   - 这是为实现快速、精确、稳定控制效果的最终版本。
% =========================================================================

% --- 步骤 0: 解析输入参数 ---

xe = z1(1); ye = z1(2); psie = z1(3);
u = dynamic(1); r = dynamic(5);

% % =========================================================================
% %             >>>>>  位置环/运动学环 (PID升级) <<<<<
% % =========================================================================
% % --- 期望线速度 ud 的设计 (保持不变) ---
% k_xe1 = 0.8; k_xe2 = 0.6; beta_kin = 1.2;
% smooth_sign_xe = @(x, delta) x ./ sqrt(x.^2 + delta.^2);
% ud_feedback = k_xe1 * tanh(xe) + k_xe2 * abs(xe)^(beta_kin-1) * smooth_sign_xe(xe, 0.1);
% ud = desired_speed * cos(psie) + ud_feedback;
% 
% % --- 【重要升级】期望角速度 rd 的设计 (最终性能PID参数) ---
% kp_ye = 1.2;   % 比例增益 (P) - 提高以加快响应 (原为0.5)
% ki_ye = 0.2;   % 积分增益 (I) - 提高以消除稳态误差 (原为0.1)
% kd_psi = 1.5;  % 微分增益 (D) - 提高以提供更好阻尼 (原为1.0)
% 
% rd = +kp_ye * tanh(ye) + ki_ye * integral_err_kin(1) + kd_psi * tanh(psie);
% 固定时间算子 sig(x)^a = |x|^a * sign(x)
sig = @(x, a) sign(x) * (abs(x)^a);
% sig = @(x, a) tanh(x/0.1) * (abs(x)^a);
% epsilon = 0.05;  % 边界层厚度
% sig = @(x, a) (abs(x) >= epsilon) * (sign(x) * abs(x)^a) + ...
%               (abs(x) < epsilon)  * (x * epsilon^(a - 1));
% --- 步骤 1: 外环虚拟控制律设计 (严格符合 Fixed-Time 结构) ---
m_kin = 1.2; n_kin = 0.6;

% 【分离增益，打破内耗】
% 1. 纵向位置增益 (保持 1.0，确保不掉队)
kx_1 = 1.0; 
kx_2 = 1.0; 
ud = desired_speed * cos(psie) + kx_1 * sig(xe, m_kin) + kx_2 * sig(xe, n_kin);

% 2. 横向与航向增益 (破局核心：重赏横向纠偏，轻罚航向误差)
ky_gain = 1.0;  % 大幅提高横向权重，让船对偏离航线极其敏感
kpsi_1 = 1.2;   % 大幅降低航向权重，允许船在转弯时灵活偏转船头
kpsi_2 = 0.4;   

% 完美补齐反步法推导中 \dot{\psi}_d 对应的前馈项 r_form
rd = r_form + ky_gain * sig(ye, n_kin) + kpsi_1 * sig(psie, m_kin) + kpsi_2 * sig(psie, n_kin);

% --- 步骤 2: 动力学误差定义 ---
u_e = u - ud; 
r_e = r - rd; 
z2 = [u_e; r_e];
% =========================================================================
%          >>>>>    动力学环 (肌肉 - 平顺执行)   <<<<<
% =========================================================================
u_e = u - ud; r_e = r - rd; z2 = [u_e; r_e];
K1 = diag([1.0, 0.8]); K2 = diag([0.5, 0.4]);
beta = 1.2; m_matrix = diag([198, 181]);
current_v_norm = norm([u; r]);
sigma = l_hat(1) + l_hat(2) * current_v_norm + l_hat(3) * (current_v_norm^2);
smooth_sign_dyn = @(x, delta) x ./ sqrt(x.^2 + delta.^2);
s_z2 = smooth_sign_dyn(z2, 0.05);
tanh_z2 = tanh(z2);
feedback_term = - sigma * tanh_z2 - K1 * tanh_z2 - K2 * diag(abs(z2).^(beta - 1)) * s_z2;
tau = m_matrix * feedback_term;
tau_max = [500; 300]; tau = max(min(tau, tau_max), -tau_max);
output = [tau(1); 0; tau(2)];
end
% END OF FILE: control.m