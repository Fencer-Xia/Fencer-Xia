function wrapped_angle = wrapToPi(angle)
    % 将角度归一化到 [-pi, pi] 范围
    wrapped_angle = angle;
    while wrapped_angle > pi
        wrapped_angle = wrapped_angle - 2*pi;
    end
    while wrapped_angle < -pi
        wrapped_angle = wrapped_angle + 2*pi;
    end
end