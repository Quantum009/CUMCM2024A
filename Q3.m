clc; clear; close all;
% 参数设置
for canshu = 40:47
a = canshu / (2 * pi);  % 螺线参数
v = 100;  % 速度 1 m/s = 100 cm/s
total_time_start = 0;  % 开始时间
total_time_end = 400;    % 结束时间
qsqsqs = false;
dt = 1;  % 时间步长1秒
initial_theta = 16 * 2 * pi;  % 初始位置在第16圈结束点
foundd = false;%%引入控制变量
% 龙的结构
head_length = 341 - 2 * 27.5;  % 龙头实际连接长度 cm
body_tail_length = 220 - 2 * 27.5;  % 龙身和龙尾实际连接长度 cm
total_sections = 224;  % 总节数（1龙头 + 221龙身 + 1龙尾）
% 创建保存结果的矩阵，维度是223*100行，4列
% 第一列是时间，第二列是x，第三列是y，第四列是节点编号
result_matrix = zeros(total_sections * (total_time_end - total_time_start ) / dt+1, 4);
index = 1;  % 矩阵的索引
% 主循环，遍历时间
for t = total_time_start:dt:total_time_end
    % 计算当前theta
    s = v * t;  % 行进的弧长
    current_theta = calculate_theta(initial_theta, s, a); % 龙头当前的theta
    
    % 初始化当前帧的坐标数组
    x_coords = zeros(1, total_sections);
    y_coords = zeros(1, total_sections);
    
    % 计算龙头位置
    [x, y] = get_position(current_theta, a);
    x_coords(1) = x;
    y_coords(1) = y;
    
    % 计算其他部分的位置
    prev_theta = current_theta; % 从龙头开始
    for i = 2:total_sections
        % 根据节数选择合适的长度
        if i == 2
            L = head_length; % 第一节与龙头之间的长度
        else
            L = body_tail_length; % 其他节之间的长度
        end
        
        % 初始化角度增量
        delta_theta = 0.01; % 初始增量，越小越精确
        current_distance = 0; % 初始化距离
        % 使用牛顿法找到合适的 delta_theta
        delta_theta = newton_method_delta_theta(prev_theta, L, a, x_coords(i-1), y_coords(i-1), 1e-6);
        % 更新位置
        prev_theta = prev_theta + delta_theta;
        [new_x, new_y] = get_position(prev_theta, a);
        x_coords(i) = new_x;
        y_coords(i) = new_y;
    end
    
    % 将当前时刻所有节点的位置保存到矩阵
    for o = 1:total_sections
        result_matrix(index, :) = [t, x_coords(o), y_coords(o), o];
        index = index + 1;
    end
end
tic
timess = (total_time_end - total_time_start)/dt + 1;
for i = 0 : timess-1
    % 初始化 B 数组为空
    B = [];
    
    % 第一个循环，计算矩形信息
    for f = 1:total_sections-1
        % 获取当前和下一个节点的坐标
        x_cords = [result_matrix(i*224+f, 2), result_matrix(i*224+f+1, 2)];
        y_cords = [result_matrix(i*224+f, 3), result_matrix(i*224+f+1, 3)];
        
        % 假设的宽度和长度
        width = 30;  % 假设的宽度
        length = 27.5;  % 假设的长度
        
        % 调用 output_dot 函数并将结果存入 B
        B(f, :) = output_dot(x_cords, y_cords, width, length);
    end
    
    % 进行矩形相交的判断
    rect1 = [B(1,1),B(1,2); B(1,3),B(1,4);
             B(1,5),B(1,6); B(1,7),B(1,8)];
    
    foundd = false;  % 初始化 foundd 为 false
    for w = 3:223
        rect2 = [B(w,1),B(w,2); B(w,3),B(w,4);
                 B(w,5),B(w,6); B(w,7),B(w,8)];
        
        % 判断矩形是否相交
        B(w,9) = check_rectangles_intersect(rect1, rect2);
        if B(w,9) == 1 && result_matrix(i*224+1,2)^2+result_matrix(i*224+1,3)^2<=450^2 ...
                && result_matrix(i*224+w,2)^2+result_matrix(i*224+w,3)^2<=(16*canshu)^2 ...
                % && result_matrix(i*224+224,2)^2+result_matrix(i*224+224,3)^2<=(16*canshu)^2
            foundd = true;  % 如果相交，设置 foundd 为 true
            break;          % 跳出内层循环
        end
         if B(w,9) == 1 &&  result_matrix(i*224+w,2)^2+result_matrix(i*224+w,3)^2<=(16*canshu)^2
            qsqsqs = true;  % 如果相交，设置 foundd 为 true
            break;          % 跳出内层循环
         end
    end
    
    if foundd
        result_matrix(i*total_sections+w,:)
        B(w,:)
        i% 输出相交的位置
        canshu
        break;  % 跳出外层循环
    end
    if qsqsqs == true
        break;
    end
end
toc
end
% row = size(B(:,1));
% for i = 1:row
%     if(B(i,9)==1)
%         disp(i)
%         disp(result_matrix(i,:))
%     end
% end
% 辅助函数：计算theta
function theta = calculate_theta(theta0, s, a)
    % 定义弧长函数
    arclength = @(t) (a / 2) * (t * sqrt(1 + t^2) + log(t + sqrt(1 + t^2)));
    % 定义要求解的方程
    equation = @(t) arclength(t) - arclength(theta0) + s;
    % 使用fzero求解
    theta = fzero(equation, theta0 - s / (a * theta0));
end
% 辅助函数：获取位置
function [x, y] = get_position(theta, a)
    x = a * theta * cos(theta);
    y = a * theta * sin(theta);
end
function result = output_dot(x_coords, y_coords, width, length)
    % 计算方向向量（从第一个点到第二个点）
    direction_vector = [x_coords(1) - x_coords(2), y_coords(1) - y_coords(2)];
    direction_length = sqrt(sum(direction_vector.^2));  % 计算方向向量的长度
    direction_unit_vector = direction_vector / direction_length;  % 单位方向向量
    % 计算法向量（垂直于方向向量，用来计算宽度）
    normal_vector = [-direction_unit_vector(2), direction_unit_vector(1)];  % 顺时针旋转90度
    % 计算长方形的四个角
    p1 = [x_coords(1), y_coords(1)] + (width / 2) * normal_vector + length * direction_unit_vector; % 右上
    p2 = [x_coords(1), y_coords(1)] - (width / 2) * normal_vector + length * direction_unit_vector; % 左上
    p3 = [x_coords(2), y_coords(2)] - (width / 2) * normal_vector - length * direction_unit_vector;                                 % 左下
    p4 = [x_coords(2), y_coords(2)] + (width / 2) * normal_vector - length * direction_unit_vector;                                 % 右下
    % 将四个点的坐标按顺序平铺成 1 行 8 列
    result = [p1(1), p1(2), p2(1), p2(2), p3(1), p3(2), p4(1), p4(2)];
end
function isIntersect = check_rectangles_intersect(rect1, rect2)
    % rect1 和 rect2 分别是两个长方形的 4 个顶点坐标，按照顺序给出
    % rect1 和 rect2 都是 4x2 的矩阵，表示 (x, y) 坐标
    isIntersect = false; % 初始化是否相交标志
    % 获取两个长方形的 4 条边
    edges1 = [rect1; rect1(1,:)]; % 将第一个点复制到最后，形成闭环
    edges2 = [rect2; rect2(1,:)]; % 将第一个点复制到最后，形成闭环
    
    % 遍历所有的边
    for i = 1:4
        for j = 1:4
            % 获取边 (i, i+1) 和 (j, j+1) 的两个点
            A1 = edges1(i, :);
            A2 = edges1(i+1, :);
            B1 = edges2(j, :);
            B2 = edges2(j+1, :);
            
            % 判断这两条边是否相交
            [a, b] = check_segments_intersect(A1, A2, B1, B2);
            % 如果 a 和 b 都在 [0,1] 之间，则说明这两条边相交
            if a > 0 && a < 1 && b > 0 && b < 1
                isIntersect = true;
                return; % 如果找到相交的边，直接返回
            end
        end
    end
end
function [a, b] = check_segments_intersect(A1, A2, B1, B2)
    % 计算参数 a 和 b 的值
    % 两条线段的方程 (1-a)A1 + aA2 = (1-b)B1 + bB2
    % 返回 a 和 b，判断是否有解，并且 a 和 b 是否在 [0,1] 范围内
    
    % 线段向量
    dA = A2 - A1;
    dB = B2 - B1;
    
    % 设置方程组 Ax = b，求解 x = [a; b]
    A = [dA(1), -dB(1); dA(2), -dB(2)];
    bVec = B1' - A1';
    
    % 检查行列式是否为 0，防止求解无效
    if abs(det(A)) < 1e-10
        a = inf; % 如果行列式为 0，说明线段平行或重叠
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
        % 计算导数（使用有限差分法）
        h = 1e-8;  % 很小的变化，确保导数计算精确
        f_prime = (distance_constraint(prev_theta, delta_theta + h, a, L, x_prev, y_prev) - f) / h;
        
        % 更新 delta_theta
        delta_theta_new = delta_theta - f / f_prime;
        
        % 检查是否收敛到 10^-6 精度
        if abs(delta_theta_new - delta_theta) < tol
            break;
        end
        
        % 更新 delta_theta
        delta_theta = delta_theta_new;
    end
end
% 辅助函数：距离约束，用于牛顿法求解
function dist_error = distance_constraint(prev_theta, delta_theta, a, L, x_prev, y_prev)
    theta_new = prev_theta + delta_theta;
    [new_x, new_y] = get_position(theta_new, a);
    dist_error = sqrt((new_x - x_prev)^2 + (new_y - y_prev)^2) - L;
end
 