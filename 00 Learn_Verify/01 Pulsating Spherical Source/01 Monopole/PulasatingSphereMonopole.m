%{
 *=======================================================================================
 *========================================【M FILE】=====================================
 * Copyright 流体力学与声学技术实验室
 * ALL right reserved.See COPYRIGHT for detailed Information.
 *
 * @File:       PulasatingSphereMonopole.m
 * @Brief:      1. 以【球体】为研究对象，导入[节点坐标](分三种节点密度：58、87、123)
 *              2. 基于【点云处理】计算得到节点[外法向量]，并进行可视化
 *              3. 基于【解析式】求解【脉动球源】辐射声压(距离、波数)
 *              4. 基于【波叠加法】求解【脉动球源】辐射声压(距离、波数)
 *              5. 对比【解析式】和【波叠加法】计算结果
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
    导入结构离散节点坐标文件(分三种节点密度：58、87、123)
%}
OriginPoint_Coord_Table = readtable("Sphere_R1m_58.txt", 'VariableNamingRule', 'preserve');           % 读取节点坐标
% OriginPoint_Coord_Table = readtable("Sphere_R1m_87.txt", 'VariableNamingRule', 'preserve');             % 读取节点坐标
% OriginPoint_Coord_Table = readtable("Sphere_R1m_123.txt", 'VariableNamingRule', 'preserve');          % 读取节点坐标
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
Scale_Factor = 0.5;                                                                                     % 缩比系数，即【等效简单源】所在球面半径 / 【脉动球源】半径，默认为0.5

% 波数范围
k_Range = linspace(0, 20, 200);                                                                         % 定义波数 k 的范围
Monitor_Point = [0, 0, 1];                                                                              % 【波数变化】监测场点

%% ------------------------------【2 结构表面节点可视化】------------------------------
%{
    结构表面节点可视化
%}
figure;                                                                                                 % 图(1.1) 结构表面离散节点图
subplot(2, 2, 1);
Fun_MultiPlot3(1, OriginPoint_X_Array, OriginPoint_Y_Array, OriginPoint_Z_Array, 'X', 'Y', 'Z', '图1.1 结构表面离散节点图', true, -Sphere_Radius, Sphere_Radius, -Sphere_Radius, Sphere_Radius, -Sphere_Radius, Sphere_Radius, 'n', 'r', 'n', true);

subplot(2, 2, 2);                                                                                       % 图(1.2) 结构表面离散节点图
[x_OriginSphere,y_OriginSphere,z_OriginSphere] = sphere;                                                % 生成半径为1的球体模型，便于显示离散点
x_OriginSphere = x_OriginSphere * Sphere_Radius;                                                        % 按【脉动球源半径】比例缩放
y_OriginSphere = y_OriginSphere * Sphere_Radius;
z_OriginSphere = z_OriginSphere * Sphere_Radius;
surf(x_OriginSphere, y_OriginSphere, z_OriginSphere);                                                   % 绘制球体
hold on;
Fun_MultiPlot3(1, OriginPoint_X_Array, OriginPoint_Y_Array, OriginPoint_Z_Array, 'X', 'Y', 'Z', '图1.2 结构表面离散节点图', true, -Sphere_Radius, Sphere_Radius, -Sphere_Radius, Sphere_Radius, -Sphere_Radius, Sphere_Radius, 'n', 'r', 'n', true);
hold off;

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

subplot(2, 2, 3);                                                                                       % 图(1.3) 结构表面离散节点及外法线图
pcshow(ptCloud, 'MarkerSize', 200);
hold on;
quiver3(ptCloud_X, ptCloud_Y, ptCloud_Z, scale * Normal_X, scale * Normal_Y, scale * Normal_Z);
hold off;
Fun_GraphSet('X', 'Y', 'Z', '图1.3 结构表面离散节点及外法线图', true, -1.5 * Sphere_Radius, 1.5 * Sphere_Radius, -1.5 * Sphere_Radius, 1.5 * Sphere_Radius, -1.5 * Sphere_Radius, 1.5 * Sphere_Radius);
whitebg('white');
axis on;

subplot(2, 2, 4);                                                                                       % 图(1.4) 结构表面离散节点及外法线图
pcshow(ptCloud, 'MarkerSize', 200);
hold on;
quiver3(ptCloud_X, ptCloud_Y, ptCloud_Z, scale * Normal_X, scale * Normal_Y, scale * Normal_Z);
Fun_MultiPlot3(2, x_OriginSphere, y_OriginSphere, z_OriginSphere, 'X', 'Y', 'Z', '图1.4 结构表面离散节点及外法线图', true, -1.5 * Sphere_Radius, 1.5 * Sphere_Radius, -1.5 * Sphere_Radius, 1.5 * Sphere_Radius, -1.5 * Sphere_Radius, 1.5 * Sphere_Radius, 'n', 'n', 'n', 'n');
hold off;
whitebg('white');
axis on;

% 法向量【数值】
% Normal_X = 2 .* OriginPoint_X_Array;
% Normal_Y = 2 .* OriginPoint_Y_Array;
% Normal_Z = 2 .* OriginPoint_Z_Array;
% OriginPoint_Normals_Array = [Normal_X, Normal_Y, Normal_Z];

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
subplot(2, 2, 1);                                                                                       % 图(2.1) 【解析解】结构及测点视图(距离)
scatter3(FieldPoints(:, 1), FieldPoints(:, 2), FieldPoints(:, 3), 'g', 'filled');
hold on;
surf(x_OriginSphere, y_OriginSphere, z_OriginSphere);                                                   % 绘制球体
Fun_MultiPlot3(1, OriginPoint_X_Array, OriginPoint_Y_Array, OriginPoint_Z_Array, 'X', 'Y', 'Z', '图2.1 【解析解】结构及测点视图(距离)', true, 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'r', 'n', true);
hold off;

subplot(2, 2, 2);                                                                                       % 图(2.2) 【解析解】结构辐射声压幅值(距离)
Fun_MultiPlot2(1, FieldPoints(:, 1), abs(P_Distance), '离球心的距离 (m)', '声压幅值 (Pa)', '图2.2 【解析解】结构辐射声压幅值(距离)', false, 10 * Sphere_Radius, 20 * Sphere_Radius);
subplot(2, 2, 3);                                                                                       % 图(2.3) 【解析解】结构辐射声压实部(距离)
Fun_MultiPlot2(1, FieldPoints(:, 1), real(P_Distance), '离球心的距离 (m)', '实部声压 (Pa)', '图2.3 【解析解】结构辐射声压实部(距离)', false, 10 * Sphere_Radius, 20 * Sphere_Radius, 'n', 'n', 'n', '#8b7d6b', 'n');
subplot(2, 2, 4);                                                                                       % 图(2.4) 【解析解】结构辐射声压虚部(距离)
Fun_MultiPlot2(1, FieldPoints(:, 1), imag(P_Distance), '离球心的距离 (m)', '虚部声压 (Pa)', '图2.4 【解析解】结构辐射声压虚部(距离)', false, 10 * Sphere_Radius, 20 * Sphere_Radius, 'n', 'n', 'n', '#5f9ea0', 'n'); 

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
subplot(2, 2, 1);                                                                                       % 图(3.1) 【解析解】结构及测点视图(波数)
scatter3(Monitor_Point(:, 1), Monitor_Point(:, 2), Monitor_Point(:, 3), 'g', 'filled');
hold on;
surf(x_OriginSphere, y_OriginSphere, z_OriginSphere);                                                   % 绘制球体
Fun_MultiPlot3(1, OriginPoint_X_Array, OriginPoint_Y_Array, OriginPoint_Z_Array, 'X', 'Y', 'Z', '图3.1 【解析解】结构及测点视图(波数)', true, 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'r', 'n', true);
hold off;

subplot(2, 2, 2);                                                                                       % 图(3.2) 【解析解】结构辐射声压幅值(波数)
Fun_MultiPlot2(1, k_Range, abs(P_kRange), '波数 k', '声压幅值 (Pa)', '图3.2 【解析解】结构辐射声压幅值(波数)', false, min(k_Range), max(k_Range)); 
subplot(2, 2, 3);                                                                                       % 图(3.3) 【解析解】结构辐射声压实部(波数)
Fun_MultiPlot2(1, k_Range, real(P_kRange), '波数 k', '实部声压 (Pa)', '图3.3 【解析解】结构辐射声压实部(波数)', false, min(k_Range), max(k_Range), 'n', 'n', 'n', '#8b7d6b', 'n');
subplot(2, 2, 4);                                                                                       % 图(3.4) 【解析解】结构辐射声压虚部(波数)
Fun_MultiPlot2(1, k_Range, imag(P_kRange), '波数 k', '虚部声压 (Pa)', '图3.4 【解析解】结构辐射声压虚部(波数)', false, min(k_Range), max(k_Range), 'n', 'n', 'n', '#5f9ea0', 'n');


%% ------------------------------【5 等效简单源辐射声压】------------------------------
%{
    基于波叠加法，求解偶极矩阵、源强度矩阵，进而求解等效压力场
%}

EquivSource_Index_Array = OriginPoint_Num_Array;                                                        % 【等效简单源】索引
EquivSource_X_Array = Scale_Factor .* OriginPoint_X_Array;                                              % 【等效简单源】X坐标
EquivSource_Y_Array = Scale_Factor .* OriginPoint_Y_Array;                                              % 【等效简单源】X坐标
EquivSource_Z_Array = Scale_Factor .* OriginPoint_Z_Array;                                              % 【等效简单源】X坐标
EquivSource_Coord_Array = [EquivSource_X_Array, EquivSource_Y_Array, EquivSource_Z_Array];              % 【等效简单源】三轴坐标

[x_EquivSphere,y_EquivSphere,z_EquivSphere] = sphere;                                                   % 生成半径为1的球体模型，便于显示离散点
x_EquivSphere = x_EquivSphere * Sphere_Radius * Scale_Factor;                                           % 按【脉动球源半径】比例缩放
y_EquivSphere = y_EquivSphere * Sphere_Radius * Scale_Factor;
z_EquivSphere = z_EquivSphere * Sphere_Radius * Scale_Factor;

figure;
subplot(2, 2, 1);                                                                                       % 图(4.1) 结构表面离散节点图
Fun_MultiPlot3(1, OriginPoint_X_Array, OriginPoint_Y_Array, OriginPoint_Z_Array, 'X', 'Y', 'Z', '图4.1 结构表面离散节点图', true, 'n', 'n', 'n', 'n', 'n', 'n', 'n', 'r', 'n', true);

subplot(2, 2, 2);                                                                                       % 图(4.2) 结构表面离散节点图
surf(x_OriginSphere, y_OriginSphere, z_OriginSphere);                                                   % 绘制【原尺寸】球体
hold on;
Fun_MultiPlot3(1, OriginPoint_X_Array, OriginPoint_Y_Array, OriginPoint_Z_Array, 'X', 'Y', 'Z', '图4.2 结构表面离散节点图', true, -Sphere_Radius, Sphere_Radius, -Sphere_Radius, Sphere_Radius, -Sphere_Radius, Sphere_Radius, 'n', 'r', 'n', true);
hold off;

subplot(2, 2, 3);                                                                                       % 图(4.3) 等效简单源图
Fun_MultiPlot3(1, EquivSource_X_Array, EquivSource_Y_Array, EquivSource_Z_Array, 'X', 'Y', 'Z', '图4.3 等效简单源图', true, -Sphere_Radius, Sphere_Radius, -Sphere_Radius, Sphere_Radius, -Sphere_Radius, Sphere_Radius, 'n', 'b', 'n', true);

subplot(2, 2, 4);                                                                                       % 图(4.4) 等效简单源图
surf(x_EquivSphere, y_EquivSphere, z_EquivSphere);                                                      % 绘制【缩比】球体
hold on;
scatter3(EquivSource_X_Array, EquivSource_Y_Array, EquivSource_Z_Array, 'b', 'filled');
Fun_MultiPlot3(1, EquivSource_X_Array, EquivSource_Y_Array, EquivSource_Z_Array, 'X', 'Y', 'Z', '图4.4 等效简单源图', true, -Sphere_Radius, Sphere_Radius, -Sphere_Radius, Sphere_Radius, -Sphere_Radius, Sphere_Radius, 'n', 'b', 'n', true);
hold off;

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
title('图5.1 【解析解】结构及测点视图(距离)');

subplot(2, 2, 2);                                                                                       % 图(5.2) 结构辐射声压幅值对比(距离)
plot(FieldPoints(:, 1), abs(P_Distance), 'LineWidth', 1.5);
hold on;
Fun_MultiPlot2(1, FieldPoints(:, 1), abs(Equiv_P_Radiation), '离球心的距离 (m)', '声压幅值 (Pa)', '图5.2 结构辐射声压幅值对比(距离)', false, 10 * Sphere_Radius, 20 * Sphere_Radius, 'n', 'n', '--', 'r', 1.5);
hold off;

subplot(2, 2, 3);                                                                                       % 图(5.3) 结构辐射声压实部对比(距离)
plot(FieldPoints(:, 1), real(P_Distance), 'Color', '#8b7d6b', 'LineWidth', 1.5);
hold on;
Fun_MultiPlot2(1, FieldPoints(:, 1), real(Equiv_P_Radiation), '离球心的距离 (m)', '实部声压 (Pa)', '图5.3 结构辐射声压实部对比(距离)', false, 10 * Sphere_Radius, 20 * Sphere_Radius, 'n', 'n', '--', 'r', 1.5);
hold off;
 
subplot(2, 2, 4);                                                                                       % 图(5.4) 结构辐射声压虚部对比(距离)
plot(FieldPoints(:, 1), imag(P_Distance), 'Color', '#5f9ea0', 'LineWidth', 1.5);
hold on;
Fun_MultiPlot2(1, FieldPoints(:, 1), imag(Equiv_P_Radiation), '离球心的距离 (m)', '虚部声压 (Pa)', '图5.4 结构辐射声压虚部对比(距离)', false, 10 * Sphere_Radius, 20 * Sphere_Radius, 'n', 'n', '--', 'r', 1.5);
hold off;

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

subplot(2, 2, 2);                                                                                       % 图(6.2) 结构辐射声压幅值对比(波数)
plot(k_Range, abs(P_kRange), 'LineWidth', 1.5);
hold on;
Fun_MultiPlot2(1, k_Range, abs(Equiv_P_kRange), '波数 k', '声压幅值 (Pa)', '图6.2 结构辐射声压幅值对比(波数)', false, min(k_Range), max(k_Range), 'n', 'n', '--', 'r', 1.5);
hold off;

subplot(2, 2, 3);                                                                                       % 图(6.3) 结构辐射声压实部对比(波数)
plot(k_Range, real(P_kRange), 'Color', '#8b7d6b', 'LineWidth', 1.5);
hold on;
Fun_MultiPlot2(1, k_Range, real(Equiv_P_kRange), '波数 k', '实部声压 (Pa)', '图6.3 结构辐射声压实部对比(波数)', false, min(k_Range), max(k_Range), 'n', 'n', '--', 'r', 1.5);
hold off;

subplot(2, 2, 4);                                                                                       % 图(6.4) 结构辐射声压虚部对比(波数)
plot(k_Range, imag(P_kRange), 'Color', '#5f9ea0', 'LineWidth', 1.5);
hold on;
Fun_MultiPlot2(1, k_Range, imag(Equiv_P_kRange), '波数 k', '虚部声压 (Pa)', '图6.4 结构辐射声压虚部对比(波数)', false, min(k_Range), max(k_Range), 'n', 'n', '--', 'r', 1.5);
hold off;
