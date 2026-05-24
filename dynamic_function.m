function output = dynamic_function(tau_a, VV, dis, kinematic)
    % 动力学参数
    m = 185; g = 9.8; W = m * g; B = W;
    zG = 0.01; zB = -0.01;
    GG = zG * W - zB * B;
    
    % 当前姿态
    theta = kinematic(4);  % pitch angle

    % 状态变量
    u = VV(1); v = VV(2); w = VV(3);
    q = VV(4); r = VV(5);

    % 干扰项
    dis1 = dis(1); dis2 = dis(2); dis3 = dis(3);
    dis4 = dis(4); dis5 = dis(5);

    % 水动力参数
    Iyy = 95; Izz = 95;
    Xu_dot = -13; Yv_dot = -257; Zw_dot = -257;
    Mq_dot = -86; Nr_dot = -86;
    Xu = -49; Xuu = -16;
    Yv = -243; Yvv = -542;
    Zw = -230; Zww = -422;
    Mq = -140; Mqq = -62;
    Nr = -86; Nrr = -78;

    % 质量矩阵 M
    M = diag([m - Xu_dot, m - Yv_dot, m - Zw_dot, Iyy - Mq_dot, Izz - Nr_dot]);

    % 科氏矩阵 C(v)
    C = [  0                0                0               (m - Zw_dot) * w    -(m - Yv_dot) * v;
           0                0                0                0                   (m - Xu_dot) * u;
           0                0                0               -(m - Xu_dot) * u    0;
          -(m - Zw_dot) * w 0                (m - Xu_dot) * u 0                   0;
           (m - Yv_dot) * v -(m - Xu_dot) * u 0               0                   0 ];

    % 阻尼矩阵 D(v)
    D = -diag([ Xu + Xuu * abs(u), ...
                Yv + Yvv * abs(v), ...
                Zw + Zww * abs(w), ...
                Mq + Mqq * abs(q), ...
                Nr + Nrr * abs(r)]);

    % 重力与浮力矩阵
    G = diag([0, 0, 0, GG * sin(theta), 0]);

    % 合成变量
    vv = [u; v; w; q; r];
    tau = [tau_a(1); 0; 0; tau_a(2); tau_a(3)];
    taud = [dis1; dis2; dis3; dis4; dis5];

    % 加速度计算
    av = pinv(M) * (tau - C * vv - D * vv - taud - G * ones(5, 1));

    % 输出
    output = av;
end
