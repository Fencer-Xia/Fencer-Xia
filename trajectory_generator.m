function [pos_d, vel_d] = trajectory_generator(t, mode, params)
switch lower(mode)
    case 'circle'
        R = params.R;
        w = params.w;
        h = params.h;
        pos_d = [R * cos(w * t); R * sin(w * t); h];
        vel_d = [-R * w * sin(w * t); R * w * cos(w * t); 0];
    case 'line'
        v = params.v;
        start = params.start_point;
        dir = params.direction;
        pos_d = start + v * t * dir;
        vel_d = v * dir;
    otherwise
        error('Unsupported mode');
end
end
