clc; clear; close all;
a = 55 / (2 * pi);  
v = 100;  
total_time = 300;  
dt = 1; 
initial_theta = 16 * 2 * pi;  
% 龙的结构
head_length = 341 - 2 * 27.5;  
body_tail_length = 220 - 2 * 27.5;  
total_sections = 224;  

x_coords_all = zeros(total_sections, total_time+1);  
y_coords_all = zeros(total_sections, total_time+1);  
v0 = zeros(total_sections-1, total_time+1);  

fig = figure('Position', [100, 100, 800, 800]);
axis equal
hold on
grid on
title('板凳龙运动模拟');
xlabel('X (cm)');
ylabel('Y (cm)');

theta_bg = linspace(0, initial_theta, 500000);
x_bg = a * theta_bg .* cos(theta_bg);
y_bg = a * theta_bg .* sin(theta_bg);
plot(x_bg, y_bg, 'k--', 'LineWidth', 0.5);

%GIF
gif_filename = 'dragon_dance_precise.gif';
frame_delay = 0.1;  

for t = 0:dt:total_time
    s = v * t; 
    current_theta = calculate_theta(initial_theta, s, a);  
    
    x_coords = zeros(1, total_sections);
    y_coords = zeros(1, total_sections);
    

    [x, y] = get_position(current_theta, a);
    x_coords(1) = x;
    y_coords(1) = y;
    

    prev_theta = current_theta;  
    v_front = v; 


    dx_dtheta = a * (cos(current_theta) - current_theta * sin(current_theta));
    dy_dtheta = a * (sin(current_theta) + current_theta * cos(current_theta));
    tangent_vector_front = [dx_dtheta, dy_dtheta];

    for i = 2:total_sections
        if i == 2
            L = head_length;  % 第一节与龙头的长度
        else
            L = body_tail_length;  % 其他节长度
        end
        
        % 使用牛顿法
        delta_theta = newton_method_delta_theta(prev_theta, L, a, x_coords(i-1), y_coords(i-1), 1e-9);
        prev_theta = prev_theta + delta_theta;
        [new_x, new_y] = get_position(prev_theta, a);
        x_coords(i) = new_x;
        y_coords(i) = new_y;


        flag=0;
        if prev_theta>16*2*pi
            flag=1;
        end
        dx_dtheta = a * (cos(prev_theta) - prev_theta * sin(prev_theta));
        dy_dtheta = a * (sin(prev_theta) + prev_theta * cos(prev_theta));
        tangent_vector_later = [dx_dtheta, dy_dtheta];
    
        % 计算速度
        board_vector = [x_coords(i-1) - x_coords(i), y_coords(i-1) - y_coords(i)];
        cos_alpha = board_vector * tangent_vector_front' / (norm(board_vector) * norm(tangent_vector_front));
        cos_beta = board_vector * tangent_vector_later' / (norm(board_vector) * norm(tangent_vector_later));
        v_later = v_front * cos_alpha / cos_beta;
        v0(i-1, t+1) = v_later; 

    if flag==1
    x_coords(i) = NaN;
        y_coords(i) = NaN;
        v0(i-1, t+1) = NaN;
    end
    %剔除未进入点

        tangent_vector_front = tangent_vector_later;
        v_front = v_later;
    end
    
    % 存储所有节的坐标
    x_coords_all(:, t+1) = x_coords';
    y_coords_all(:, t+1) = y_coords';
    
    if t > 0
        delete(dragon_plot);
        delete(head_plot);
        delete(tail_plot);
    end
    dragon_plot = plot(x_coords, y_coords, 'b-', 'LineWidth', 2);
    head_plot = plot(x_coords(1), y_coords(1), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    tail_plot = plot(x_coords(end), y_coords(end), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    
    % 更新图形
    drawnow
    
    frame = getframe(fig);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if t == 0
        imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', frame_delay);
    else
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', frame_delay);
    end
end
v0=v0/100;
x_coords_all=x_coords_all/100;
y_coords_all=y_coords_all/100;
% 牛顿法
function delta_theta = newton_method_delta_theta(prev_theta, L, a, x_prev, y_prev, tol)
    delta_theta = 0.001;  % 初始猜测
    max_iter = 10000;  % 最大迭代次数
    for iter = 1:max_iter
        f = distance_constraint(prev_theta, delta_theta, a, L, x_prev, y_prev);
        % 计算导数（使用有限差分法）
        h = 1e-9; 
        f_prime = (distance_constraint(prev_theta, delta_theta + h, a, L, x_prev, y_prev) - f) / h;
        
        delta_theta_new = delta_theta - f / f_prime;
        if abs(delta_theta_new - delta_theta) < tol
            break;
        end
        
        delta_theta = delta_theta_new;
    end
end

% 距离约束
function dist_error = distance_constraint(prev_theta, delta_theta, a, L, x_prev, y_prev)
    theta_new = prev_theta + delta_theta;
    [new_x, new_y] = get_position(theta_new, a);
    dist_error = sqrt((new_x - x_prev)^2 + (new_y - y_prev)^2) - L;
end

% 根据theta计算位置
function [x, y] = get_position(theta, a)
    x = a * theta * cos(theta);
    y = a * theta * sin(theta);
end

% 根据弧长计算 theta
function theta = calculate_theta(theta0, s, a)
    arclength = @(t) (a / 2) * (t * sqrt(1 + t^2) + log(t + sqrt(1 + t^2)));
    arclength_derivative = @(t) (a / 2) * (sqrt(1 + t^2) + (t^2 / sqrt(1 + t^2)));
    
    theta = theta0;  
    tol = 1e-12;  % 精度要求
    max_iter = 10000;  % 最大迭代次数
    for iter = 1:max_iter
        % 牛顿法更新theta
        f = arclength(theta) - arclength(theta0) + s;
        f_prime = arclength_derivative(theta);
        theta_new = theta - f / f_prime;
        if abs(theta_new - theta) < tol
            break;
        end
        theta = theta_new;
    end
end
