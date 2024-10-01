clc; clear;
% 参数设置
a = 170/2/pi; % 螺线参数
theta = 90/17*pi; % 入螺线的角度参数
big_circle_radius = 450; % 大圆区域的半径
for v = 110:0.1:200 % 速度 (cm/s)
    bool = false;
    dt = 0.1;  % 时间步长
    % total_time = 2000;  % 总时间
    initial_theta = 16 * 2 * pi;  % 初始位置在第16圈结束点
    rect_width = 30;  % 长方形的宽度
    rect_length = 27.5;  % 长方形的长度
    n_test = 2;
    % 龙的结构
    head_length = 341 - 2 * 27.5;  % 龙头实际连接长度 cm
    body_tail_length = 220 - 2 * 27.5;  % 龙身和龙尾实际连接长度 cm
    total_sections = 224;  % 总节数（1龙头 + 221龙身 + 1龙尾） % 可以理解为224个板，都取每个板的第一个人
    % 计算入螺线的盘入点坐标和切线
    coordinate = [a*theta*cos(theta), a*theta*sin(theta)]; % 入螺线的坐标
    % 入螺线起点的角度
    theta0 = atan(coordinate(2)/coordinate(1)); % 使用 atan 函数计算起点角度
    % 设置 delta_theta，取值区间为 [0, pi)
    delta_theta = pi / 43.84; % 最佳理想值
    % 计算出螺线的终点坐标
    gamma = theta0 + delta_theta;
    coordinate_out = [big_circle_radius*cos(gamma), big_circle_radius*sin(gamma)];
    out_center = coordinate_out+coordinate;
    % 计算切线方向的单位向量
    vector_in = [-a*(cos(theta)-theta*sin(theta)), -a*(sin(theta)+theta*cos(theta))];
    vector_in = vector_in/norm(vector_in);
    % 垂直方向单位向量
    vector_perp = [vector_in(2), -vector_in(1)];
    % 计算 a 和 b
    a_value = norm(cross([vector_in 0], [coordinate_out - coordinate 0]));
    b_value = dot((coordinate_out - coordinate), vector_in);
    % 计算对应的 R
    n_test = 2;
    R = (a_value^2 + b_value^2) / (2 * a_value) / (n_test + 1);
    % 计算圆心 M 和 N 的位置
    M = coordinate + n_test * R * vector_perp;
    N = coordinate_out - R * vector_perp;
    % 计算切点位置
    vector_MN = (N - M) / norm(N - M);
    P_M = M + R * n_test * vector_MN;  % 圆 M 的切点
    P_N = N - R * vector_MN;           % 圆 N 的切点
    P = P_N;
    % 计算圆弧的角度范围和坐标
    theta_M = linspace(atan2(coordinate(2)-M(2), coordinate(1)-M(1)), atan2(P_M(2)-M(2), P_M(1)-M(1)) - 2*pi, 100);
    x_M = M(1) + n_test * R * cos(theta_M);
    y_M = M(2) + n_test * R * sin(theta_M);
    theta_N = linspace(atan2(P_N(2)-N(2), P_N(1)-N(1)), atan2(coordinate_out(2)-N(2), coordinate_out(1)-N(1)), 100);
    x_N = N(1) + R * cos(theta_N);
    y_N = N(2) + R * sin(theta_N);
    % 初始化 龙头刚到coordinate点所有点的速度
    v0 = zeros(224, 1);  
    v0(1) = v;
    % 初始化当前帧的坐标数组
    x_coords = zeros(1, total_sections);
    y_coords = zeros(1, total_sections);
    % 得到龙头位置
    x_coords(1) = coordinate(1);
    y_coords(1) = coordinate(2);
    current_theta = norm(coordinate)/a;
    % 计算其他部分的位置
    prev_theta = current_theta; % 从龙头开始
    v_front = v;
    % 初始化前者的切线向量
    dx_dtheta = a * (cos(current_theta) - current_theta * sin(current_theta));
    dy_dtheta = a * (sin(current_theta) + current_theta * cos(current_theta));
    tangent_vector_front = [dx_dtheta, dy_dtheta];
    for i = 2:total_sections
        % 根据节数选择合适的长度
        if i == 2
            L = head_length; % 第一节与龙头之间的长度
        else
            L = body_tail_length; % 其他节之间的长度
        end
        % 使用牛顿法找到合适的 delta_theta
        delta_theta = newton_method_delta_theta(prev_theta, L, a, x_coords(i-1), y_coords(i-1), 1e-9);
        % 更新位置
        prev_theta = prev_theta + delta_theta;
        [new_x, new_y] = get_position(prev_theta, a);
        x_coords(i) = new_x;
        y_coords(i) = new_y;
        
        % 计算当前节的切线向量
        dx_dtheta = a * (cos(prev_theta) - prev_theta * sin(prev_theta));
        dy_dtheta = a * (sin(prev_theta) + prev_theta * cos(prev_theta));
        tangent_vector_later = [dx_dtheta, dy_dtheta];
        board_vector = [x_coords(i-1) - x_coords(i), y_coords(i-1) - y_coords(i)];
        cos_aphla = board_vector * tangent_vector_front'/norm(board_vector)/norm(tangent_vector_front);
        cos_beta = board_vector * tangent_vector_later'/norm(board_vector)/norm(tangent_vector_later);
        v_later = v_front * cos_aphla / cos_beta;
        v0(i) = v_later;  % 使用round确保数组索引正确
        % 更新前后者向量
        tangent_vector_front = tangent_vector_later;
        v_front = v_later;
    end
    cs_zuobiao = [x_coords ;y_coords]';
    all_matrix = [cs_zuobiao ,v0];
    s_to_coordinate = zeros(1, total_sections);
    for i = 2:total_sections
        theta_i = get_theta(x_coords(i), y_coords(i), a);
        s_to_coordinate(i) = calculate_s(theta_i, current_theta, a);
    end
    all_matrix = [all_matrix ,s_to_coordinate'];
    all_matrix = [all_matrix ones(224,1)];
    all_matrix(1,5) = 2;
    
    % % 初始化图形
    % fig = figure('Position', [100, 100, 800, 800]);
    % axis equal
    % hold on
    % grid on
    % title('舞龙队运动模拟');
    % xlabel('X (cm)');
    % ylabel('Y (cm)');
    %
    % % 绘制螺线背景
    % theta_bg = linspace(0, initial_theta, 1000);
    % x_bg = a * theta_bg .* cos(theta_bg);
    % y_bg = a * theta_bg .* sin(theta_bg);
    % plot(x_bg, y_bg, 'k--', 'LineWidth', 0.5);
    %
    % % 准备GIF
    % gif_filename = 'dragon_dance_with_rectangles.gif';
    % frame_delay = 0.1;  % GIF帧间延迟（秒）
    
    for ppp = 1:1600 % ppp是次数
        cla; % 清除当前坐标区的所有对象
        s = v * ppp * dt; % 当前行进的弧长
    
        % 确定龙头的当前位置和方向
        if s > 0 && s  <= n_test * R * abs(theta_M(1) - theta_M(end))
            all_matrix(1,5) = 2;
            % 沿蓝色弧线运动
            alpha_M = s/n_test/R;
            x_current = M(1) + n_test * R * cos(theta_M(1) - alpha_M); % 当前圆弧上的 x 坐标
            y_current = M(2) + n_test * R * sin(theta_M(1) - alpha_M); % 当前圆弧上的 y 坐标
            all_matrix(1,1:2)= [x_current, y_current];
            direct_vector(1,:) = get_direct_vector_2(x_current, y_current, M);
        elseif s > n_test * R * abs(theta_M(1) - theta_M(end)) && s - n_test * R * (theta_M(1) - theta_M(end)) <= R * abs(theta_N(1) - theta_N(end))
            all_matrix(1,5) = 3;
            % 沿绿色弧线运动
            alpha_N = (s - n_test * R * (theta_M(1) - theta_M(end)))/R;
            x_current = N(1) + R * cos(theta_N(1) + alpha_N); % 当前圆弧上的 x 坐标
            y_current = N(2) + R * sin(theta_N(1) + alpha_N); % 当前圆弧上的 y 坐标
            all_matrix(1,1:2)= [x_current, y_current];
            direct_vector(1,:) = get_direct_vector_3(x_current, y_current, N);
        elseif s - n_test * R * (theta_M(1) - theta_M(end)) > R * abs(theta_N(1) - theta_N(end))
            all_matrix(1,5) = 4;
            % 沿出螺线运动
            ds = s - n_test * R * (theta_M(1) - theta_M(end)) - R * abs(theta_N(1) - theta_N(end));
            theta_out = calculate_theta_ni(norm(coordinate_out)/a, ds, a); % 出螺线的起始点
            x_current = - a * theta_out * cos(theta_out) + out_center(1); % 当前出螺线上的 x 坐标
            y_current = - a * theta_out * sin(theta_out) + out_center(2); % 当前出螺线上的 y 坐标
            all_matrix(1,1:2)= [x_current, y_current];
            direct_vector(1,:) = get_direct_vector_4(x_current, y_current, a, coordinate, coordinate_out);
        end
    
        % 更新其他点的坐标
        for qqq = 2:224
            % 根据状态更新点的坐标
            if all_matrix(qqq,5) == 1
                ds = all_matrix(qqq, 3) * dt;
                [x_current, y_current, x_over] = in_hx(all_matrix(qqq,1), all_matrix(qqq,2), ds, coordinate, a);
                if x_over == 0
                    all_matrix(qqq,1:2)= [x_current, y_current];
                    direct_vector(qqq,:) = get_direct_vector_1(all_matrix(qqq,1), all_matrix(qqq,2), a);
                else
                    [x_current, y_current, x_over] = in_circle(coordinate(1),coordinate(2),x_over,M,R,P);
                    all_matrix(qqq,1:2) = [x_current,y_current];
                    all_matrix(qqq,5) = 2;
                    direct_vector(qqq,:) = get_direct_vector_2(all_matrix(qqq,1), all_matrix(qqq,2), M);
                end
            elseif all_matrix(qqq,5) == 2
                % 更新在圆弧上的位置
                x = all_matrix(qqq,1);
                y = all_matrix(qqq,2);
                ds = all_matrix(qqq, 3) * dt;
                [x_current, y_current, x_over] = in_circle(x, y, ds, M, R, P);
                if x_over == 0
                    all_matrix(qqq,1:2)= [x_current, y_current];
                    direct_vector(qqq,:) = get_direct_vector_2(all_matrix(qqq,1), all_matrix(qqq,2), M);
                else
                    [x_current, y_current, ~] = out_circle(P(1), P(2), x_over, N, R, coordinate_out);
                    all_matrix(qqq,1:2) = [x_current,y_current];
                    all_matrix(qqq,5) = 3;
                    direct_vector(qqq,:) = get_direct_vector_3(all_matrix(qqq,1), all_matrix(qqq,2), N);
                end
            elseif all_matrix(qqq,5) == 3
                % 更新在第二个圆弧上的位置
                x = all_matrix(qqq,1);
                y = all_matrix(qqq,2);
                ds = all_matrix(qqq, 3) * dt;
                [x_current, y_current, x_over] = out_circle(x, y, ds, N, R, coordinate_out);
                if x_over == 0
                    all_matrix(qqq,1:2)= [x_current,y_current];
                    direct_vector(qqq,:) = get_direct_vector_3(all_matrix(qqq,1), all_matrix(qqq,2), N);
                else
                    [x_current, y_current] = out_hx(coordinate_out(1),coordinate_out(2),x_over,a,coordinate,coordinate_out);
                    all_matrix(qqq,1:2) = [x_current,y_current];
                    all_matrix(qqq,5) = 4;
                    direct_vector(qqq,:) = get_direct_vector_4(all_matrix(qqq,1), all_matrix(qqq,2), a, coordinate, coordinate_out);
                end
            elseif all_matrix(qqq,5) == 4
                % 更新在出螺线上的位置
                x = all_matrix(qqq,1);
                y = all_matrix(qqq,2);
                ds = all_matrix(qqq, 3) * dt;
                [x_current, y_current] = out_hx(x,y,ds,a,coordinate,coordinate_out);
                all_matrix(qqq,1:2) = [x_current,y_current];
                direct_vector(qqq,:) = get_direct_vector_4(all_matrix(qqq,1), all_matrix(qqq,2), a, coordinate, coordinate_out);
            end
        end
    
        % 更新速度和方向
        v_front = v;
        tangent_vector_front = direct_vector(1,:);
    
        for i = 2:total_sections
            % 计算每节的速度
            tangent_vector_later = direct_vector(i,:);
            board_vector = [all_matrix(i-1,1) - all_matrix(i,1), all_matrix(i-1,2) - all_matrix(i,2)];
            cos_aphla = board_vector * tangent_vector_front' / norm(board_vector) / norm(tangent_vector_front);
            cos_beta = board_vector * tangent_vector_later' / norm(board_vector) / norm(tangent_vector_later);
            v_later = v_front * cos_aphla / cos_beta;
            all_matrix(i,3) = v_later;
            tangent_vector_front = tangent_vector_later;
            v_front = v_later;
        end
        % % 绘制所有点的路径
        % plot(all_matrix(:,1), all_matrix(:,2), 'o-', 'MarkerSize', 3, 'LineWidth', 1);
        % drawnow; % 实时更新图形
        % pause(0.001); % 增加延迟以观察变化
    
        if max(all_matrix(:,3)) > 200
            bool = true;
            disp(v)
            break
        end
    end
    if bool == true
        break
    end
end
%% 第一个圆弧部分
function [x_current, y_current, ds_over] = in_circle(x, y, ds, M, R, P)
    ds_over = 0;
    d_theta = ds/2/R;
    pre_theta = atan2( y-M(2), x-M(1));%atan2先输入y的值再输入x
    now_theta = pre_theta - d_theta;
    x_current = M(1) + 2 * R * cos(now_theta);
    y_current = M(2) + 2 * R * sin(now_theta);
    next_point = [x_current, y_current];
    pre_point = [x, y];
    judge_vector = cross([next_point - M, 0], [P - M, 0]);
    if judge_vector(3) >= 0 % 说明跑过头了
        theta0 = real(acos((P-M)*(pre_point-M)'/norm(P-M)/norm(pre_point-M)));
        ds_over = ds - theta0*2*R;
        x_current = P(1);
        y_current = P(2);
    end
end
%% 第二个圆弧部分
function [x_current, y_current, ds_over] = out_circle(x, y, ds, N, R, coordinate_out)
    ds_over = 0;
    d_theta = ds/R;
    pre_theta = atan2( y-N(2),x-N(1));
    now_theta = pre_theta + d_theta;
    x_current = N(1) + R * cos(now_theta);
    y_current = N(2) + R * sin(now_theta);
    next_point = [x_current, y_current];
    pre_point = [x, y];
    judge_vector = cross([next_point - N, 0], [coordinate_out - N, 0]);
    if judge_vector(3) <= 0 % 说明跑过头了
        theta0 = real(acos(dot(coordinate_out-N, pre_point-N)/norm(coordinate_out-N)/norm(pre_point-N)));
        ds_over = ds - theta0*R;
        x_current = coordinate_out(1);
        y_current = coordinate_out(2);
    end
end
%% 入螺线部分
function [x, y, x_over] = in_hx(x_position, y_position, ds, coordinate, a)%输入x,y的原有坐标以及运行的长度 前面那一段的螺旋
    x_over = 0;
    theta = get_theta(x_position, y_position, a);
    theta_0 = calculate_theta(theta, ds, a);%最后算出来的theta
    theta_std = get_theta(coordinate(1), coordinate(2), a);%%标准的 不能小于的theta
    if theta_0 < theta_std %如果小于了说明要变化！
        x_over = ds - calculate_s(theta, theta_std, a);
        theta_0 = theta_std;
        % disp('fuck you');
    end
    [x,y] = get_position(theta_0,a);
end
%% 出螺线部分
function [x,y] = out_hx(x_position,y_position,ds,a,coordinate,coordinate_out)
    theta = get_theta_ni(x_position,y_position,a,coordinate,coordinate_out);
    theta_0 = calculate_theta_ni(theta,ds,a); %最后算出来的theta
    [x,y] = get_position_ni(theta_0, a, coordinate, coordinate_out);
end
function [x, y] = get_position_ni(theta, a, coordinate, coordinate_out)
    center_new = coordinate+coordinate_out;
    % 计算逆时针方向的 x, y 坐标
    x = -a * theta * cos(theta) + center_new(1);  % 保持与顺时针相同
    y = -a * theta * sin(theta) + center_new(2);  % 保持与顺时针相同
end
function delta_theta = newton_method_delta_theta(prev_theta, L, a, x_prev, y_prev, tol)
    delta_theta = 0.01;  % 初始猜测
    max_iter = 10000;  % 最大迭代次数
    for iter = 1:max_iter
        % 计算当前 delta_theta 下的距离误差
        f = distance_constraint(prev_theta, delta_theta, a, L, x_prev, y_prev);
        % 计算导数（使用有限差分法）
        h = 1e-8;  % 很小的变化，确保导数计算精确
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
% 辅助函数：距离约束，用于牛顿法求解
function dist_error = distance_constraint(prev_theta, delta_theta, a, L, x_prev, y_prev)
    theta_new = prev_theta + delta_theta;
    [new_x, new_y] = get_position(theta_new, a);
    dist_error = sqrt((new_x - x_prev)^2 + (new_y - y_prev)^2) - L;
end
function theta = calculate_theta(theta0, s, a)
    % 定义弧长函数
    arclength = @(t) (a / 2) * (t * sqrt(1 + t^2) + log(t + sqrt(1 + t^2)));
    % 定义弧长函数的导数（使用解析导数）
    arclength_derivative = @(t) (a / 2) * (sqrt(1 + t^2) + (t^2 / sqrt(1 + t^2)));
    % 初始theta
    theta = theta0;
    % 牛顿法参数
    tol = 1e-9;  % 精度要求
    max_iter = 10000;  % 最大迭代次数
    for iter = 1:max_iter
        % 计算方程的当前值
        f = arclength(theta) - arclength(theta0) + s;
        % 计算导数的当前值
        f_prime = arclength_derivative(theta);
        % 使用牛顿法更新theta
        theta_new = theta - f / f_prime;
        % 检查是否满足精度要求
        if abs(theta_new - theta) < tol
            break;
        end
        % 更新theta
        theta = theta_new;
    end
    if iter == max_iter
        warning('牛顿法未在最大迭代次数内收敛');
    end
end
% function theta = calculate_theta(theta0, s, a)
%     % 定义弧长函数
%     arclength = @(t) (a / 2) * (t * sqrt(1 + t^2) + log(t + sqrt(1 + t^2)));
%
%     % 定义要求解的方程
%     equation = @(t) arclength(t) - arclength(theta0) + s;
%
%     % 定义二分法的初始边界
%     lower_bound = theta0 - 10;  % 可以根据具体情况调整
%     upper_bound = theta0 + 10;  % 可以根据具体情况调整
%     tol = 1e-3;  % 精度要求
%     max_iter = 1000;  % 最大迭代次数
%
%     % 检查初始边界是否适合
%     if equation(lower_bound) * equation(upper_bound) > 0
%         error('初始边界未包含根，请调整初始区间。');
%     end
%
%     % 二分法求解
%     for iter = 1:max_iter
%         mid_point = (lower_bound + upper_bound) / 2;
%         f_mid = equation(mid_point);
%
%         if abs(f_mid) < tol
%             theta = mid_point;
%             return;
%         elseif equation(lower_bound) * f_mid < 0
%             upper_bound = mid_point;
%         else
%             lower_bound = mid_point;
%         end
%     end
%
%     error('二分法未在最大迭代次数内收敛');
% end
function [x, y] = get_position(theta, a)
    x = a * theta * cos(theta);
    y = a * theta * sin(theta);
end
function s = calculate_s(theta0, theta1, a) % theta0是上限，theta1是下限
    % 定义弧长函数
    arclength = @(t) (a / 2) * (t * sqrt(1 + t^2) + log(t + sqrt(1 + t^2)));
    % 定义要求解的方程
    s = arclength(theta0) - arclength(theta1);
end
function theta = get_theta(x, y, a)
    theta = sqrt(x^2 + y^2) / a;
end
% function theta = calculate_theta_ni(theta0, s, a)
%     % 逆时针旋转角度计算，theta 角度增大
%     % 定义弧长函数
%     arclength = @(t) (a / 2) * (t * sqrt(1 + t^2) + log(t + sqrt(1 + t^2)));
%
%     % 定义弧长函数的导数（使用解析导数）
%     arclength_derivative = @(t) (a / 2) * (sqrt(1 + t^2) + (t^2 / sqrt(1 + t^2)));
%
%     % 初始theta
%     theta = theta0;
%
%     % 牛顿法参数
%     tol = 1e-9;  % 精度要求
%     max_iter = 1000;  % 最大迭代次数
%
%     for iter = 1:max_iter
%         % 计算方程的当前值 (逆时针方向，弧长增加)
%         f = arclength(theta) - arclength(theta0) - s;  % 注意这里减去 s，使 theta 增大
%
%         % 计算导数的当前值
%         f_prime = arclength_derivative(theta);
%
%         % 使用牛顿法更新theta
%         theta_new = theta - f / f_prime;
%
%         % 检查是否满足精度要求
%         if abs(theta_new - theta) < tol
%             break;
%         end
%
%         % 更新theta
%         theta = theta_new;
%     end
%
%     if iter == max_iter
%         warning('牛顿法未在最大迭代次数内收敛');
%     end
% end
function theta = calculate_theta_ni(theta0, s, a)
    % 逆时针旋转角度计算，theta 角度增大
    % 定义弧长函数
    arclength = @(t) (a / 2) * (t * sqrt(1 + t^2) + log(t + sqrt(1 + t^2)));
    
    % 定义弧长函数的导数（使用解析导数）
    arclength_derivative = @(t) (a / 2) * (sqrt(1 + t^2) + (t^2 / sqrt(1 + t^2)));
    
    % 初始theta
    theta = theta0 + 17/2/pi;
    
    % 牛顿法参数
    tol = 1e-9;  % 精度要求
    max_iter = 1000;  % 最大迭代次数
    
    for iter = 1:max_iter
        % 计算方程的当前值 (逆时针方向，弧长增加)
        f = arclength(theta) - arclength(theta0) - s;  % 注意这里减去 s，使 theta 增大
        
        % 惩罚条件：确保 theta 大于 theta0
        if theta == theta0
            f = f + abs(theta - theta0); % 增大误差以驱动 theta 增加
        end
        
        % 计算导数的当前值
        f_prime = arclength_derivative(theta);
        
        % 使用牛顿法更新theta
        theta_new = theta - f / f_prime;
        
        % 检查是否满足精度要求
        if abs(theta_new - theta) < tol
            break;
        end
        
        % 更新theta
        theta = theta_new;
    end
    
    if iter == max_iter
        warning('牛顿法未在最大迭代次数内收敛');
    end
end
function theta = get_theta_ni(x,y,a,coordinate,coordinate_out)
    new_yuandian = coordinate+coordinate_out;
    theta = sqrt((new_yuandian(1)-x)^2+(new_yuandian(2)-y)^2)/a ;
end
function direct_vector = get_direct_vector_1(x, y, a)
    theta = sqrt(x^2+y^2)/a;
    d_x = a*(cos(theta)-theta*sin(theta));
    d_y = a*(sin(theta)+theta*cos(theta));
    direct_vector = [-d_x, -d_y];
    direct_vector = direct_vector/norm(direct_vector);
end
function direct_vector = get_direct_vector_4(x, y, a, coordinate, coordinate_out)
    x1 = coordinate(1)+coordinate_out(1)-x;
    y1 = coordinate(2)+coordinate_out(2)-y;
    direct_vector = get_direct_vector_1(x1, y1, a);
end
function direct_vector = get_direct_vector_2(x, y, M)
    current_coordinate = [x,y];
    TM = M - current_coordinate;
    direct_vector = [-TM(2), TM(1)];
    direct_vector = direct_vector/norm(direct_vector);
end
function direct_vector = get_direct_vector_3(x, y, N)
    current_coordinate = [x,y];
    TN = N - current_coordinate;
    direct_vector = [TN(2), -TN(1)];
    direct_vector = direct_vector/norm(direct_vector);
end