clc; clear; close all;
% 参数设置
canshu = 55;  
a = canshu / (2 * pi);
v = 100;  
total_time_start = 0;  % 开始时间
total_time_end = 425;  % 结束时间
dt = 0.001;  % 时间步长 0.1 秒
initial_theta = 16 * 2 * pi;  
found_collision = false; 
collision_times = [];  

% 龙的结构
head_length = 341 - 2 * 27.5;  % 龙头实际连接长度，单位 cm
body_tail_length = 220 - 2 * 27.5;  % 龙身和龙尾实际连接长度，单位 cm
total_sections = 224;  % 总节数（包括龙头、龙身和龙尾共224节）

% 创建保存结果的矩阵，维度是 (total_sections * 时间步数)，4列分别存储时间、x坐标、y坐标和节点编号
result_matrix = zeros(total_sections * (total_time_end - total_time_start) / dt + 1, 4);
index = 1;  % 矩阵索引，用于记录每一帧数据

% 主循环，遍历时间
for t = total_time_start:dt:total_time_end
    % 计算当前龙头的角度
    s = v * t;  % 当前时刻的弧长（龙头前进的距离）
    current_theta = calculate_theta(initial_theta, s, a);  % 根据弧长计算龙头的角度
    
    % 初始化当前帧的坐标数组
    x_coords = zeros(1, total_sections);
    y_coords = zeros(1, total_sections);
    
    % 计算龙头的位置
    [x, y] = get_position(current_theta, a);  % 获取龙头的 x 和 y 坐标
    x_coords(1) = x;
    y_coords(1) = y;
    
    % 计算龙身和龙尾的其他节点位置
    prev_theta = current_theta;  % 从龙头开始，计算每节龙的角度和位置
    for i = 2:total_sections
        % 根据节数选择长度，龙头与第一节的距离是 head_length，其他部分是 body_tail_length
        if i == 2
            L = head_length;  % 龙头和第一节之间的长度
        else
            L = body_tail_length;  % 龙身和龙尾之间的长度
        end
        
        % 使用牛顿法找到合适的 delta_theta 使得相邻两节龙的距离正确
        delta_theta = newton_method_delta_theta(prev_theta, L, a, x_coords(i-1), y_coords(i-1), 1e-6);
        
        % 更新节的角度和位置
        prev_theta = prev_theta + delta_theta;
        [new_x, new_y] = get_position(prev_theta, a);  % 获取新节点的位置
        x_coords(i) = new_x;
        y_coords(i) = new_y;
    end
    
    % 将当前时刻所有节点的位置保存到结果矩阵
    for o = 1:total_sections
        result_matrix(index, :) = [t, x_coords(o), y_coords(o), o];  % 记录时间、位置和节点编号
        index = index + 1;
    end
end

% 碰撞检测部分
qsqsqs = false;
timess = (total_time_end - total_time_start)/dt + 1;
k_neighbors = 20;  % 只检查距离最近的 20 个节点

% 遍历每个时间步长，检测舞龙是否发生碰撞
for i = 0 : timess-1
    B = [];  % 用于存储每个节点的矩形信息
    
    % 计算每个节点的矩形信息
    for f = 1:total_sections-1
        x_cords = [result_matrix(i*224+f, 2), result_matrix(i*224+f+1, 2)];
        y_cords = [result_matrix(i*224+f, 3), result_matrix(i*224+f+1, 3)];
        
        % 假设龙的宽度和长度
        width = 30;  
        length = 27.5;  
        
        B(f, :) = output_dot(x_cords, y_cords, width, length);
    end
    
    % 针对龙头节点进行最近点查找
    current_node = 1;  % 龙头
    current_x = result_matrix(i*224 + current_node, 2);
    current_y = result_matrix(i*224 + current_node, 3);

    % 将所有节点的坐标（除去龙头自身）存入矩阵
    all_coords = result_matrix(i*224 + 2:i*224 + total_sections, 2:3);  
    head_coord = [current_x, current_y];  % 龙头的坐标
    
    % 使用 knnsearch 查找最近的 20 个节点
    [idx, distances] = knnsearch(all_coords, head_coord, 'K', 20);
    
    % 获取这些最近节点的矩形信息
    closest_indices = idx;  % 最近的节点索引
    
    rect1 = [B(current_node,1), B(current_node,2);
             B(current_node,3), B(current_node,4);
             B(current_node,5), B(current_node,6);
             B(current_node,7), B(current_node,8)];
    
    foundd = false;
    for idx = 3:20
        rect2 = [B(idx,1), B(idx,2);
                 B(idx,3), B(idx,4);
                 B(idx,5), B(idx,6);
                 B(idx,7), B(idx,8)];
        
        % 判断是否相交
        B(idx,9) = check_rectangles_intersect(rect1, rect2);
        if B(idx,9) == 1
            foundd = true;
            collision_times = [collision_times; i * dt];  % 记录发生碰撞的时间
            break;
        end
    end
    
    if foundd
        break;  % 如果检测到相交，跳出循环
    end
end

% 输出发生碰撞的时间
if ~isempty(collision_times)
    disp('碰撞发生的时间（秒）：');
    disp(collision_times);
else
    disp('未发生碰撞。');
end

% 辅助函数：计算theta
function theta = calculate_theta(theta0, s, a)
    % 定义弧长函数
    arclength = @(t) (a / 2) * (t * sqrt(1 + t^2) + log(t + sqrt(1 + t^2)));
    % 使用fzero求解
    equation = @(t) arclength(t) - arclength(theta0) + s;
    theta = fzero(equation, theta0 - s / (a * theta0));
end

% 获取位置
function [x, y] = get_position(theta, a)
    x = a * theta * cos(theta);
    y = a * theta * sin(theta);
end

% 计算矩形的角点
function result = output_dot(x_coords, y_coords, width, length)
    direction_vector = [x_coords(1) - x_coords(2), y_coords(1) - y_coords(2)];
    direction_length = sqrt(sum(direction_vector.^2));
    direction_unit_vector = direction_vector / direction_length;
    
    normal_vector = [-direction_unit_vector(2), direction_unit_vector(1)];
    
    p1 = [x_coords(1), y_coords(1)] + (width / 2) * normal_vector + length * direction_unit_vector;
    p2 = [x_coords(1), y_coords(1)] - (width / 2) * normal_vector + length * direction_unit_vector;
    p3 = [x_coords(2), y_coords(2)] - (width / 2) * normal_vector - length * direction_unit_vector;
    p4 = [x_coords(2), y_coords(2)] + (width / 2) * normal_vector - length * direction_unit_vector;
    
    result = [p1(1), p1(2), p2(1), p2(2), p3(1), p3(2), p4(1), p4(2)];
end

% 检查矩形是否相交
function isIntersect = check_rectangles_intersect(rect1, rect2)
    isIntersect = false; % 初始化是否相交标志
    % 获取两个长方形的 4 条边
    edges1 = [rect1; rect1(1,:)]; 
    edges2 = [rect2; rect2(1,:)]; 
    
    % 遍历所有的边
    for i = 1:4
        for j = 1:4
            % 获取边 (i, i+1) 和 (j, j+1) 的两个点
            A1 = edges1(i, :);
            A2 = edges1(i+1, :);
            B1 = edges2(j, :);
            B2 = edges2(j+1, :);
            

            [a, b] = check_segments_intersect(A1, A2, B1, B2);
            % 如果 a 和 b 都在 [0,1] 之间，则说明这两条边相交
            if a >= 0 && a <= 1 && b >= 0 && b <= 1
                isIntersect = true;
                return;
            end
        end
    end
end

% 检查线段是否相交的辅助函数
function [a, b] = check_segments_intersect(A1, A2, B1, B2)
    % 计算线段向量
    dA = A2 - A1;
    dB = B2 - B1;
    
    % 设置方程组 Ax = b，求解 x = [a; b]
    A = [dA(1), -dB(1); dA(2), -dB(2)];
    bVec = B1' - A1';

    if abs(det(A)) < 1e-10
        a = inf; 
        b = inf;
        return;
    end
    
    % 求解 a 和 b
    solution = A\bVec;
    a = solution(1);
    b = solution(2);
end

% 牛顿法求解 delta_theta
function delta_theta = newton_method_delta_theta(prev_theta, L, a, x_prev, y_prev, tol)
    delta_theta = 0.01;  % 初始猜测
    max_iter = 1000;  % 最大迭代次数
    for iter = 1:max_iter
        % 计算当前 delta_theta 下的距离误差
        f = distance_constraint(prev_theta, delta_theta, a, L, x_prev, y_prev);
        h = 1e-8; 
        f_prime = (distance_constraint(prev_theta, delta_theta + h, a, L, x_prev, y_prev) - f) / h;
        
        % 更新 delta_theta
        delta_theta_new = delta_theta - f / f_prime;
        

        if abs(delta_theta_new - delta_theta) < tol
            break;
        end
        
        % 更新 delta_theta
        delta_theta = delta_theta_new;
    end
end


function dist_error = distance_constraint(prev_theta, delta_theta, a, L, x_prev, y_prev)
    theta_new = prev_theta + delta_theta;
    [new_x, new_y] = get_position(theta_new, a);
    dist_error = sqrt((new_x - x_prev)^2 + (new_y - y_prev)^2) - L;
end
