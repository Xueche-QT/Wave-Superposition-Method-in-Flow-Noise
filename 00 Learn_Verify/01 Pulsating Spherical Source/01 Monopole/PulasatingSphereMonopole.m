%{
 *=======================================================================================
 *========================================【M FILE】=====================================
 * Copyright 流体力学与声学技术实验室
 * ALL right reserved.See COPYRIGHT for detailed Information.
 *
 * @File:       PulasatingSphereMonopole.m
 * @Brief:      1. 以【球体】为研究对象，导入[节点坐标]
 *              2. 基于二分法获取[最佳缩放系数]，以[均方根误差]为判断准则
 *              3. 求解场点频域声压
 *
 * @Author:     Haiger
 * @date:       2023.06.12
 *=======================================================================================
%}

clc;
clear;
close all;

%% ------------------------------【1 导入数据 / 输入参数】------------------------------
%{
    导入结构离散节点坐标文件
%}
OriginPoint_Coord_Table = readtable("Sphere_R1m_58.txt", 'VariableNamingRule', 'preserve');             % 读取节点坐标
OriginPoint_Num_Array = OriginPoint_Coord_Table{:, 1};                                                  % 节点索引
OriginPoint_X_Array = OriginPoint_Coord_Table{:, 2};                                                    % 节点X轴坐标Array
OriginPoint_Y_Array = OriginPoint_Coord_Table{:, 3};                                                    % 节点Y轴坐标Array
OriginPoint_Z_Array = OriginPoint_Coord_Table{:, 4};                                                    % 节点Z轴坐标Array
OriginPoint_Coord_Array = [OriginPoint_X_Array, OriginPoint_Y_Array, OriginPoint_Z_Array];              % 节点三轴坐标Array

% 结构信息
Sphere_Radius = 1;                                                                                      % 【脉动球源】直径
Velocity_Amplitude = 1;                                                                                 % 表面振速幅值，默认沿着法向
OriginPoint_NodeVelocity_Array = ones(length(OriginPoint_Num_Array), 1) * Velocity_Amplitude;           % 【脉动球源】离散节点 [速度幅值]

% 介质信息
Density_Medium = 1.225;                                                                                 % 介质密度
AcousticVelocity_Medium = 340;                                                                          % 介质声速
Vibra_Fre_Origin = 100;                                                                                 % 振动频率[初始]，后续考虑【全波数】范围内【振动频率】会随波数的变化而变化
Omega_Origin = 2 * pi * Vibra_Fre_Origin;                                                               % 圆频率[初始]，后续考虑【全波数】范围内【圆频率】会随波数的变化而变化
k_Origin = Omega_Origin / AcousticVelocity_Medium;                                                      % 波数[初始]，后续考虑【全波数】范围内【波数】会随波数的变化而变化

% 辐射声压检测场点的位置和数量
FieldPoints_Num = 50;                                                                                   % 自定义场点数量
FieldPoints(:, 1) = linspace(10 * Sphere_Radius, 20 * Sphere_Radius, FieldPoints_Num)';                 % 建立场点距离向量 X 轴坐标，范围为10倍至20倍脉动球源半径，转置为列向量
FieldPoints(:, 2) = 0;
FieldPoints(:, 3) = 0;

% 【等效简单源】参数
Scale_Factor = 0.88;                                                                                     % 缩比系数，即【等效简单源】所在球面半径 / 【脉动球源】半径，默认为0.5

% 波数范围
k_Range = linspace(0, 20, 50);                                                                          % 定义波数 k 的范围
Monitor_Point = [0, 0, 1];                                                                              % 【波数变化】监测场点

%% ------------------------------【2 结构表面节点可视化】------------------------------
%{
    结构表面节点可视化
%}
figure;                                                                                                 % 图(1.1)
subplot(2, 2, 1);
scatter3(OriginPoint_X_Array, OriginPoint_Y_Array, OriginPoint_Z_Array, 'r', 'filled');
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
xlim([-Sphere_Radius Sphere_Radius]);                                                                   % 以脉动球源半径为横纵坐标幅值
ylim([-Sphere_Radius Sphere_Radius]);
zlim([-Sphere_Radius Sphere_Radius]);
title('图1.1 结构表面离散节点图');

subplot(2, 2, 2);                                                                                       % 图(1.2)
[x_OriginSphere,y_OriginSphere,z_OriginSphere] = sphere;                                                % 生成半径为1的球体模型，便于显示离散点
x_OriginSphere = x_OriginSphere * Sphere_Radius;                                                        % 按【脉动球源半径】比例缩放
y_OriginSphere = y_OriginSphere * Sphere_Radius;
z_OriginSphere = z_OriginSphere * Sphere_Radius;
surf(x_OriginSphere, y_OriginSphere, z_OriginSphere);                                                   % 绘制球体
hold on;
scatter3(OriginPoint_X_Array, OriginPoint_Y_Array, OriginPoint_Z_Array, 'r', 'filled');
hold off;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
xlim([-Sphere_Radius Sphere_Radius]);                                                                   % 以脉动球源半径为横纵坐标幅值
ylim([-Sphere_Radius Sphere_Radius]);
zlim([-Sphere_Radius Sphere_Radius]);
title('图1.2 结构表面离散节点图');

%% ------------------------------【3 求解结构表面离散节点外法向量】------------------------------
%{
    基于点云求解结构表面离散节点外法向
%}
ptCloud = pointCloud(OriginPoint_Coord_Array);                                                          % 读取点云数据
OriginPoint_Normals_Array = pcnormals(ptCloud);                                                         % 计算法向量(部分点为内法向，部分点为外法向)
centroid = mean(ptCloud.Location, 1);                                                                   % 计算点云的中心
vectors = centroid - ptCloud.Location;                                                                  % 计算从每个点到中心的向量
dotProduct = sum(vectors .* OriginPoint_Normals_Array, 2);                                              % 计算向量和法向量的点积
OriginPoint_Normals_Array(dotProduct > 0, :) = -OriginPoint_Normals_Array(dotProduct > 0, :);           % 取反指向中心的法向量

scale = 0.1;                                                                                            % 法向量(箭头)缩放系数

ptCloud_X = ptCloud.Location(:, 1);                                                                     % 点云的三轴坐标
ptCloud_Y = ptCloud.Location(:, 2);
ptCloud_Z = ptCloud.Location(:, 3);

Normal_X = OriginPoint_Normals_Array(:, 1);                                                             % 获取法向量的三轴坐标
Normal_Y = OriginPoint_Normals_Array(:, 2);
Normal_Z = OriginPoint_Normals_Array(:, 3);

subplot(2, 2, 3);                                                                                       % 图(1.3)
pcshow(ptCloud, 'MarkerSize', 200);
hold on;
quiver3(ptCloud_X, ptCloud_Y, ptCloud_Z, scale * Normal_X, scale * Normal_Y, scale * Normal_Z);
hold off;
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
xlim([-1.5 * Sphere_Radius 1.5 * Sphere_Radius]);                                                       % 以脉动球源半径为横纵坐标幅值
ylim([-1.5 * Sphere_Radius 1.5 * Sphere_Radius]);
zlim([-1.5 * Sphere_Radius 1.5 * Sphere_Radius]);
title('图1.3 结构表面离散节点及外法线图');
whitebg('white');
axis on;

subplot(2, 2, 4);                                                                                       % 图(1.4)
pcshow(ptCloud, 'MarkerSize', 200);
hold on;
quiver3(ptCloud_X, ptCloud_Y, ptCloud_Z, scale * Normal_X, scale * Normal_Y, scale * Normal_Z);
surf(x_OriginSphere, y_OriginSphere, z_OriginSphere);                                                   % 绘制球体
hold off;
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
xlim([-1.5 * Sphere_Radius 1.5 * Sphere_Radius]);                                                       % 以脉动球源半径为横纵坐标幅值
ylim([-1.5 * Sphere_Radius 1.5 * Sphere_Radius]);
zlim([-1.5 * Sphere_Radius 1.5 * Sphere_Radius]);
title('图1.4 结构表面离散节点及外法线图');
whitebg('white');
axis on;

%% ------------------------------【4 结构辐射声压解析解】------------------------------
%{
    基于解析式求解【脉动球源】辐射声压
    1. 以[距离]为步长
    2. 以[波数]为步长
%}

% 1. 以[距离]为步长
P_Distance = zeros(FieldPoints_Num, 1);                                                                 % 初始化【解析解】辐射声压矢量(距离)

% 计算【解析解】辐射声压(距离)
for i = 1 : FieldPoints_Num                                                                             % 遍历各个测点
    r_i = norm(FieldPoints(i));                                                                         % 测点与【脉动球源】球心的距离
    P_Distance(i) = (k_Origin * Density_Medium * AcousticVelocity_Medium * Velocity_Amplitude * Sphere_Radius^2) / (1i + k_Origin * Sphere_Radius) * exp(1i * k_Origin * (r_i - Sphere_Radius)) / r_i;       % 解析解
end

% 绘制【脉动球源】辐射声压图
figure;
subplot(2, 2, 1);                                                                                       % 图(2.1)
scatter3(FieldPoints(:, 1), FieldPoints(:, 2), FieldPoints(:, 3), 'g', 'filled');
hold on;
scatter3(OriginPoint_X_Array, OriginPoint_Y_Array, OriginPoint_Z_Array, 'r', 'filled');
surf(x_OriginSphere, y_OriginSphere, z_OriginSphere);                                                   % 绘制球体
hold off;
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('图2.1 【解析解】结构及测点视图(距离)');

subplot(2, 2, 2);                                                                                       % 图(2.2)
plot(FieldPoints(:, 1), abs(P_Distance), 'LineWidth', 1.5);
xlabel('离球心的距离 (m)');
ylabel('声压幅值 (Pa)');
title('图2.2 【解析解】结构辐射声压幅值(距离)');
xlim([10 * Sphere_Radius 20 * Sphere_Radius]); 

subplot(2, 2, 3);                                                                                       % 图(2.3)
plot(FieldPoints(:, 1), real(P_Distance), 'Color', '#8b7d6b', 'LineWidth', 1.5);
xlabel('离球心的距离 (m)');
ylabel('实部声压 (Pa)');
title('图2.3 【解析解】结构辐射声压实部(距离)');
xlim([10 * Sphere_Radius 20 * Sphere_Radius]); 

subplot(2, 2, 4);                                                                                       % 图(2.4)
plot(FieldPoints(:, 1), imag(P_Distance), 'Color', '#5f9ea0', 'LineWidth', 1.5);
xlabel('离球心的距离 (m)');
ylabel('虚部声压 (Pa)');
title('图2.4 【解析解】结构辐射声压虚部(距离)');
xlim([10 * Sphere_Radius 20 * Sphere_Radius]); 

% 2. 以[波数]为步长
P_kRange = zeros(size(k_Range));                                                                        % 初始化【解析解】辐射声压矢量(波数)

% 计算【解析解】辐射声压(波数)
r_i = norm(Monitor_Point);
for i = 1 : length(k_Range)
    k = k_Range(i);
    P_kRange(i) = (k * Density_Medium * AcousticVelocity_Medium * Velocity_Amplitude * Sphere_Radius^2) / (1i + k * Sphere_Radius) * exp(1i * k * (r_i - Sphere_Radius)) / r_i;       % 解析解
end

% 绘制【脉动球源】辐射声压图
figure;
subplot(2, 2, 1);                                                                                       % 图(3.1)
scatter3(Monitor_Point(:, 1), Monitor_Point(:, 2), Monitor_Point(:, 3), 'g', 'filled');
hold on;
scatter3(OriginPoint_X_Array, OriginPoint_Y_Array, OriginPoint_Z_Array, 'r', 'filled');
surf(x_OriginSphere, y_OriginSphere, z_OriginSphere);                                                   % 绘制球体
hold off;
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('图3.1 【解析解】结构及测点视图(波数)');

subplot(2, 2, 2);                                                                                       % 图(2.2)
plot(k_Range, abs(P_kRange), 'LineWidth', 1.5);
xlabel('波数 k');
ylabel('声压幅值 (Pa)');
title('图3.2 【解析解】结构辐射声压幅值(波数)');
xlim([min(k_Range) max(k_Range)]); 

subplot(2, 2, 3);                                                                                       % 图(2.3)
plot(k_Range, real(P_kRange), 'Color', '#8b7d6b', 'LineWidth', 1.5);
xlabel('波数 k');
ylabel('实部声压 (Pa)');
title('图3.3 【解析解】结构辐射声压实部(波数)');
xlim([min(k_Range) max(k_Range)]); 

subplot(2, 2, 4);                                                                                       % 图(2.4)
plot(k_Range, imag(P_kRange), 'Color', '#5f9ea0', 'LineWidth', 1.5);
xlabel('波数 k');
ylabel('虚部声压 (Pa)');
title('图3.4 【解析解】结构辐射声压虚部(波数)');
xlim([min(k_Range) max(k_Range)]); 


%% ------------------------------【5 等效简单源辐射声压】------------------------------
%{
    基于波叠加法，求解偶极矩阵、源强度矩阵，进而求解等效压力场
%}

EquivSource_Index_Array = OriginPoint_Num_Array;                                                        % 【等效简单源】索引
EquivSource_X_Array = Scale_Factor .* OriginPoint_X_Array;                                              % 【等效简单源】X坐标
EquivSource_Y_Array = Scale_Factor .* OriginPoint_Y_Array;                                              % 【等效简单源】X坐标
EquivSource_Z_Array = Scale_Factor .* OriginPoint_Z_Array;                                              % 【等效简单源】X坐标
EquivSource_Coord_Array = [EquivSource_X_Array, EquivSource_Y_Array, EquivSource_Z_Array];              % 【等效简单源】三轴坐标

figure;
subplot(2, 2, 1);                                                                                       % 图(4.1)
scatter3(OriginPoint_X_Array, OriginPoint_Y_Array, OriginPoint_Z_Array, 'r', 'filled');
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
xlim([-Sphere_Radius Sphere_Radius]);                                                                   % 以脉动球源半径为横纵坐标幅值
ylim([-Sphere_Radius Sphere_Radius]);
zlim([-Sphere_Radius Sphere_Radius]);
title('图4.1 结构表面离散节点图');

subplot(2, 2, 2);                                                                                       % 图(4.2)
surf(x_OriginSphere, y_OriginSphere, z_OriginSphere);                                                   % 绘制【原尺寸】球体
hold on;
scatter3(OriginPoint_X_Array, OriginPoint_Y_Array, OriginPoint_Z_Array, 'r', 'filled');
hold off;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
xlim([-Sphere_Radius Sphere_Radius]);                                                                   % 以脉动球源半径为横纵坐标幅值
ylim([-Sphere_Radius Sphere_Radius]);
zlim([-Sphere_Radius Sphere_Radius]);
title('图4.2 结构表面离散节点图');

subplot(2, 2, 3);                                                                                       % 图(4.3)
scatter3(EquivSource_X_Array, EquivSource_Y_Array, EquivSource_Z_Array, 'b', 'filled');
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
xlim([-Sphere_Radius Sphere_Radius]);                                                                   % 以脉动球源半径为横纵坐标幅值
ylim([-Sphere_Radius Sphere_Radius]);
zlim([-Sphere_Radius Sphere_Radius]);
title('图4.3 等效简单源图');

subplot(2, 2, 4);                                                                                       % 图(4.4)
[x_EquivSphere,y_EquivSphere,z_EquivSphere] = sphere;                                                   % 生成半径为1的球体模型，便于显示离散点
x_EquivSphere = x_EquivSphere * Sphere_Radius;                                                          % 按【脉动球源半径】比例缩放
y_EquivSphere = y_EquivSphere * Sphere_Radius;
z_EquivSphere = z_EquivSphere * Sphere_Radius;
surf(x_EquivSphere, y_EquivSphere, z_EquivSphere);                                                      % 绘制【缩比】球体
hold on;
scatter3(EquivSource_X_Array, EquivSource_Y_Array, EquivSource_Z_Array, 'b', 'filled');
hold off;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
xlim([-Sphere_Radius Sphere_Radius]);                                                                   % 以脉动球源半径为横纵坐标幅值
ylim([-Sphere_Radius Sphere_Radius]);
zlim([-Sphere_Radius Sphere_Radius]);
title('图4.4 等效简单源图');

% 基于[法向振速]正求[源强度矢量]
[Equiv_D_Radiation, Equiv_Q_Radiation] = Fun_BuildEquiv_D_Q(OriginPoint_Coord_Array, EquivSource_Coord_Array, k_Origin, OriginPoint_NodeVelocity_Array, OriginPoint_Normals_Array);
% 基于[源强度矢量]正求节点[辐射声压]
[Equiv_P_Radiation, Equiv_M] = Fun_BuildEquiv_M_P(FieldPoints, EquivSource_Coord_Array, k_Origin, AcousticVelocity_Medium, Density_Medium, Equiv_Q_Radiation);

% 绘制【脉动球源】辐射声压图
figure;
subplot(2, 2, 1);                                                                                       % 图(5.1)
scatter3(FieldPoints(:, 1), FieldPoints(:, 2), FieldPoints(:, 3), 'g', 'filled');
hold on;
scatter3(OriginPoint_X_Array, OriginPoint_Y_Array, OriginPoint_Z_Array, 'r', 'filled');
surf(x_OriginSphere, y_OriginSphere, z_OriginSphere);                                                   % 绘制【原尺寸】球体
hold off;
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('图2.1 【解析解】结构及测点视图(距离)');

subplot(2, 2, 2);                                                                                       % 图(5.2)
plot(FieldPoints(:, 1), abs(P_Distance), 'LineWidth', 1.5);
hold on;
plot(FieldPoints(:, 1), abs(Equiv_P_Radiation), '--r');
hold off;
xlabel('离球心的距离 (m)');
ylabel('声压幅值 (Pa)');
title('图5.2 结构辐射声压幅值对比(距离)');
xlim([10 * Sphere_Radius 20 * Sphere_Radius]); 

subplot(2, 2, 3);                                                                                       % 图(5.3)
plot(FieldPoints(:, 1), real(P_Distance), 'Color', '#8b7d6b', 'LineWidth', 1.5);
hold on;
plot(FieldPoints(:, 1), real(Equiv_P_Radiation), '--r');
hold off;
xlabel('离球心的距离 (m)');
ylabel('实部声压 (Pa)');
title('图5.3 结构辐射声压实部对比(距离)');
xlim([10 * Sphere_Radius 20 * Sphere_Radius]); 

subplot(2, 2, 4);                                                                                       % 图(5.4)
plot(FieldPoints(:, 1), imag(P_Distance), 'Color', '#5f9ea0', 'LineWidth', 1.5);
hold on;
plot(FieldPoints(:, 1), imag(Equiv_P_Radiation), '--r');
hold off;
xlabel('离球心的距离 (m)');
ylabel('虚部声压 (Pa)');
title('图5.4 结构辐射声压虚部对比(距离)');
xlim([10 * Sphere_Radius 20 * Sphere_Radius]); 

Equiv_P_kRange = zeros(size(k_Range));

for i = 1 : length(k_Range)
    k = k_Range(i);
    [Equiv_D_Radiation, Equiv_Q_Radiation] = Fun_BuildEquiv_D_Q(OriginPoint_Coord_Array, EquivSource_Coord_Array, k, OriginPoint_NodeVelocity_Array, OriginPoint_Normals_Array);
    [Equiv_P_kRange(i), ~] = Fun_BuildEquiv_M_P(Monitor_Point, EquivSource_Coord_Array, k, AcousticVelocity_Medium, Density_Medium, Equiv_Q_Radiation);
end

figure;
subplot(2, 2, 1);                                                                                       % 图(6.1)
scatter3(Monitor_Point(:, 1), Monitor_Point(:, 2), Monitor_Point(:, 3), 'g', 'filled');
hold on;
scatter3(OriginPoint_X_Array, OriginPoint_Y_Array, OriginPoint_Z_Array, 'r', 'filled');
surf(x_OriginSphere, y_OriginSphere, z_OriginSphere);                                                   % 绘制【原尺寸】球体
hold off;
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('图6.1 【解析解】结构及测点视图(波数)');

subplot(2, 2, 2);                                                                                       % 图(6.2)
plot(k_Range, abs(P_kRange), 'LineWidth', 1.5);
hold on;
plot(k_Range, abs(Equiv_P_kRange), '--r');
hold off;
xlabel('波数 k');
ylabel('声压幅值 (Pa)');
title('图6.2 【解析解】结构辐射声压幅值(波数)');
xlim([min(k_Range) max(k_Range)]); 

subplot(2, 2, 3);                                                                                       % 图(6.3)
plot(k_Range, real(P_kRange), 'Color', '#8b7d6b', 'LineWidth', 1.5);
hold on;
plot(k_Range, real(Equiv_P_kRange), '--r');
hold off;
xlabel('波数 k');
ylabel('实部声压 (Pa)');
title('图6.3 【解析解】结构辐射声压实部(波数)');
xlim([min(k_Range) max(k_Range)]); 

subplot(2, 2, 4);                                                                                       % 图(6.4)
plot(k_Range, imag(P_kRange), 'Color', '#5f9ea0', 'LineWidth', 1.5);
hold on;
plot(k_Range, imag(Equiv_P_kRange), '--r');
hold off;
xlabel('波数 k');
ylabel('虚部声压 (Pa)');
title('图6.4 【解析解】结构辐射声压虚部(波数)');
xlim([min(k_Range) max(k_Range)]); 