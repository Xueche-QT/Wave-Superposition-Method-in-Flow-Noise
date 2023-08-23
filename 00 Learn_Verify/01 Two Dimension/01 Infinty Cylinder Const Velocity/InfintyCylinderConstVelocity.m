%{
 *=======================================================================================
 *========================================【M FILE】=====================================
 * Copyright 流体力学与声学技术实验室
 * ALL right reserved.See COPYRIGHT for detailed Information.
 *
 * @File:       InfintyCylinderConstVelocity.m
 * @Brief:      1. 以【脉动球源】为研究对象，导入[节点坐标](分三种节点密度：58、87、123)
 *              2. 基于【点云处理】计算得到节点[外法向量]，并进行可视化
 *              3. 基于【解析式】求解【脉动球源】辐射声压(距离、波数)
 *              4. 基于【波叠加法】求解【脉动球源】辐射声压(距离、波数)
 *              5. 对比【解析式】和【波叠加法】计算结果
 *
 * @Author:     Haiger
 * @date:       2023.08.22
 *=======================================================================================
%}

clc;
clear;
close all;

%% ------------------------------【1 导入数据 / 输入参数】------------------------------
%{
    结构及介质参数
%}

% 结构信息
Circle_Radius = 1;                                                                                      % 【柱面】直径
Discrete_NodeNum = 30;                                                                                  % 【柱面】离散节点数目
Velocity_Amplitude = 0.001;                                                                                 % 表面振速幅值，默认沿着法向（径向）
OriginPoint_NodeVelocity_Array = ones(Discrete_NodeNum, 1) * Velocity_Amplitude;                        % 【柱面】离散节点 [速度幅值]

% 介质信息
Density_Medium = 1.225;                                                                                 % 介质密度
AcousticVelocity_Medium = 340;                                                                          % 介质声速
Vibra_Fre_Origin = 100;                                                                                 % 振动频率[初始]，后续考虑【全波数】范围内【振动频率】会随波数的变化而变化
Omega_Origin = 2 * pi * Vibra_Fre_Origin;                                                               % 圆频率[初始]，后续考虑【全波数】范围内【圆频率】会随波数的变化而变化
k_Origin = Omega_Origin / AcousticVelocity_Medium;                                                      % 波数[初始]，后续考虑【全波数】范围内【波数】会随波数的变化而变化

% 辐射声压检测场点的位置和数量
FieldPoints_Num = 50;                                                                                   % 自定义场点数量
FieldPoints_Cartesian(:, 1) = linspace(10 * Circle_Radius, 20 * Circle_Radius, FieldPoints_Num)';       % 【笛卡尔】建立场点距离向量 X 轴坐标，范围为10倍至20倍半径，转置为列向量
FieldPoints_Cartesian(:, 2) = 0;
FieldPoints_Polar = Cartesian_to_Polar_2D(FieldPoints_Cartesian);                                       % 【极坐标】场点数组[角坐标，半径坐标]

% 波数范围
k_Range = linspace(0, 20, 200);                                                                         % 定义波数 k 的范围
Monitor_Point = [0, 1];                                                                                 % 【波数变化】监测场点

% 【等效简单源】参数
Scale_Factor = 0.1;                                                                                     % 缩比系数，即【等效简单源】所在圆面半径 / 【柱面】半径，默认为0.5

%% ------------------------------【2 结构表面及简单源节点可视化】------------------------------
%{
    结构表面节点可视化
%}
OriginPoint_Theta_Array = linspace(0, 2 * pi, Discrete_NodeNum)';                                       % 【柱坐标】节点角坐标
OriginPoint_R_Array = Circle_Radius .* ones(1, Discrete_NodeNum)';                                      % 【柱坐标】节点半径坐标

OriginPoint_X_Array = Circle_Radius .* cos(OriginPoint_Theta_Array);                                    % 【直角坐标】节点X轴坐标Array
OriginPoint_Y_Array = Circle_Radius .* sin(OriginPoint_Theta_Array);                                    % 【直角坐标】节点Y轴坐标Array
OriginPoint_Cartesian_Array = [OriginPoint_X_Array, OriginPoint_Y_Array];                               % 【直角坐标】节点双轴坐标Array

figure;
subplot(2, 2, 1);                                                                                       % 图(1.1) 结构表面离散节点
Fun_MultiPlot2(1, OriginPoint_X_Array, OriginPoint_Y_Array, "X", "Y", "图1.1 结构表面离散节点", true, true, round(min(OriginPoint_X_Array) - 1), round(max(OriginPoint_X_Array) + 1), round(min(OriginPoint_Y_Array) - 1), round(max(OriginPoint_Y_Array) + 1), '-');
hold on;
plot(OriginPoint_X_Array, OriginPoint_Y_Array, 'r.', 'MarkerSize', 20);
hold off;

EquivSource_X_Array = Scale_Factor .* OriginPoint_X_Array;                                              % 【直角坐标】[等效简单源]X坐标Array
EquivSource_Y_Array = Scale_Factor .* OriginPoint_Y_Array;                                              % 【直角坐标】[等效简单源]Y坐标Array
EquivSource_Cartesian_Array = [EquivSource_X_Array, EquivSource_Y_Array];                               % 【直角坐标】[等效简单源]双轴坐标Array

subplot(2, 2, 3);                                                                                       % 图(1.3) 结构表面离散节点及内部简单源
Fun_MultiPlot2(1, OriginPoint_X_Array, OriginPoint_Y_Array, "X", "Y", "图1.3 结构表面离散节点及内部简单源", true, true, round(min(OriginPoint_X_Array) - 1), round(max(OriginPoint_X_Array) + 1), round(min(OriginPoint_Y_Array) - 1), round(max(OriginPoint_Y_Array) + 1), '-');
hold on;
plot(OriginPoint_X_Array, OriginPoint_Y_Array, 'r.', 'MarkerSize', 20);
plot(EquivSource_X_Array, EquivSource_Y_Array, 'b.', 'MarkerSize', 20);
plot(EquivSource_X_Array, EquivSource_Y_Array, 'm-', 'LineWidth', 1.5);
hold off;

%% ------------------------------【3 求解结构表面离散节点外法向量】------------------------------
%{
    求解结构表面离散节点外法向
%}

% [法一] 解析公式求解法向量
OriginPoint_Normal_Vector = zeros(Discrete_NodeNum, 2);
for ii = 1 : Discrete_NodeNum
    OriginPoint_Normal_Vector(ii, 1) = OriginPoint_X_Array(ii) / sqrt(OriginPoint_X_Array(ii)^2 + OriginPoint_Y_Array(ii)^2);
    OriginPoint_Normal_Vector(ii, 2) = OriginPoint_Y_Array(ii) / sqrt(OriginPoint_X_Array(ii)^2 + OriginPoint_Y_Array(ii)^2);
end

% [法二] 相邻段的法向量的和

% 可视化
Scale_Normal = 0.1;
subplot(2, 2, 2);                                                                                       % 图(1.2) 结构表面离散节点及外法线图
Fun_MultiPlot2(1, OriginPoint_X_Array, OriginPoint_Y_Array, "X", "Y", "图1.2 结构表面离散节点及外法线图", true, true, round(min(OriginPoint_X_Array) - 1), round(max(OriginPoint_X_Array) + 1), round(min(OriginPoint_Y_Array) - 1), round(max(OriginPoint_Y_Array) + 1), '-');
hold on;
plot(OriginPoint_X_Array, OriginPoint_Y_Array, 'r.', 'MarkerSize', 20);
quiver(OriginPoint_X_Array, OriginPoint_Y_Array, Scale_Normal * OriginPoint_Normal_Vector(:, 1), Scale_Normal * OriginPoint_Normal_Vector(:, 2), 'c');
hold off;

subplot(2, 2, 4);                                                                                       % 图(1.4) 结构表面离散节点外法线及内部简单源
Fun_MultiPlot2(1, OriginPoint_X_Array, OriginPoint_Y_Array, "X", "Y", "图1.4 结构表面离散节点外法线及内部简单源", true, true, round(min(OriginPoint_X_Array) - 1), round(max(OriginPoint_X_Array) + 1), round(min(OriginPoint_Y_Array) - 1), round(max(OriginPoint_Y_Array) + 1), '-');
hold on;
plot(OriginPoint_X_Array, OriginPoint_Y_Array, 'r.', 'MarkerSize', 20);
plot(EquivSource_X_Array, EquivSource_Y_Array, 'b.', 'MarkerSize', 20);
plot(EquivSource_X_Array, EquivSource_Y_Array, 'm-', 'LineWidth', 1.5);
quiver(OriginPoint_X_Array, OriginPoint_Y_Array, Scale_Normal * OriginPoint_Normal_Vector(:, 1), Scale_Normal * OriginPoint_Normal_Vector(:, 2), 'c');
hold off;


%% ------------------------------【4 结构辐射声压解析解】------------------------------
%{
    基于解析式求解【脉动球源】辐射声压
    1. 以[距离]为步长
    2. 以[波数]为步长
%}
% 原方程为【柱坐标】系下，优先以柱坐标节点输入

% 1. 以[距离]为步长
P_Distance = zeros(FieldPoints_Num, 1);                                                                 % 初始化【解析解】辐射声压矢量(距离)
for ii = 1 : FieldPoints_Num
    P_Distance(ii, 1) = 1i * Density_Medium * AcousticVelocity_Medium * Velocity_Amplitude * besselh(0, k_Origin * FieldPoints_Polar(ii, 2)) / ((besselj(-1, k_Origin * Circle_Radius) - besselj(1, k_Origin * Circle_Radius)) / 2 + 1i * ((bessely(-1, k_Origin * Circle_Radius) - bessely(1, k_Origin * Circle_Radius)) / 2));
end

figure;
subplot(2, 2, 1);                                                                                       % 图(2.1) 【解析解】结构及测点视图(距离)
Fun_MultiPlot2(1, OriginPoint_X_Array, OriginPoint_Y_Array, "X", "Y", "图2.1 【解析解】结构及测点视图(距离)", true, false, -2 * Circle_Radius, 20 * Circle_Radius, -2 * Circle_Radius, 2 * Circle_Radius, '-');
hold on;
plot(OriginPoint_X_Array, OriginPoint_Y_Array, 'r.', 'MarkerSize', 20);
plot(FieldPoints_Cartesian(:, 1), FieldPoints_Cartesian(:, 2), 'g.', 'MarkerSize', 20);
hold off;

subplot(2, 2, 2);                                                                                       % 图(2.2) 【解析解】结构辐射声压幅值(距离)
Fun_MultiPlot2(1, FieldPoints_Cartesian(:, 1), abs(P_Distance), "离圆心的距离(m)", "声压幅值(pa)", "图2.2 【解析解】结构辐射声压幅值(距离)", false, false, 10 * Circle_Radius, 20 * Circle_Radius);
subplot(2, 2, 3);                                                                                       % 图(2.3) 【解析解】结构辐射声压实部(距离)
Fun_MultiPlot2(1, FieldPoints_Cartesian(:, 1), real(P_Distance), '离圆心的距离 (m)', '实部声压 (Pa)', '图2.3 【解析解】结构辐射声压实部(距离)', false, false, 10 * Circle_Radius, 20 * Circle_Radius, 'n', 'n', 'n', '#8b7d6b', 'n');
subplot(2, 2, 4);                                                                                       % 图(2.4) 【解析解】结构辐射声压虚部(距离)
Fun_MultiPlot2(1, FieldPoints_Cartesian(:, 1), imag(P_Distance), '离圆心的距离 (m)', '虚部声压 (Pa)', '图2.4 【解析解】结构辐射声压虚部(距离)', false, false, 10 * Circle_Radius, 20 * Circle_Radius, 'n', 'n', 'n', '#5f9ea0', 'n'); 

% 2. 以[波数]为步长
P_kRange = zeros(size(k_Range));                                                                        % 初始化【解析解】辐射声压矢量(波数)

% 计算【解析解】辐射声压(波数)
r_i = norm(Monitor_Point);
for ii = 1 : length(k_Range)
    k = k_Range(ii);
    P_kRange(ii) = 1i * Density_Medium * AcousticVelocity_Medium * Velocity_Amplitude * besselh(0, k * r_i) / ((besselj(-1, k * Circle_Radius) - besselj(1, k * Circle_Radius)) / 2 + 1i * ((bessely(-1, k * Circle_Radius) - bessely(1, k * Circle_Radius)) / 2));
end

figure;
subplot(2, 2, 1);                                                                                       % 图(3.1) 【解析解】结构及测点视图(波数)
Fun_MultiPlot2(1, OriginPoint_X_Array, OriginPoint_Y_Array, "X", "Y", "图3.1 【解析解】结构及测点视图(波数)", true, true, -2 * Circle_Radius, 2 * Circle_Radius, -2 * Circle_Radius, 2 * Circle_Radius, '-');
hold on;
plot(OriginPoint_X_Array, OriginPoint_Y_Array, 'r.', 'MarkerSize', 20);
plot(Monitor_Point(:, 1), Monitor_Point(:, 2), 'g.', 'MarkerSize', 20);
hold off;

subplot(2, 2, 2);                                                                                       % 图(3.2) 【解析解】结构辐射声压幅值(波数)
Fun_MultiPlot2(1, k_Range, abs(P_kRange), "波数", "声压幅值(pa)", "图3.2 【解析解】结构辐射声压幅值(波数)", false, false, min(k_Range), max(k_Range));
subplot(2, 2, 3);                                                                                       % 图(3.3) 【解析解】结构辐射声压实部(波数)
Fun_MultiPlot2(1, k_Range, real(P_kRange), '波数', '实部声压 (Pa)', '图3.3 【解析解】结构辐射声压实部(波数)', false, false, min(k_Range), max(k_Range), 'n', 'n', 'n', '#8b7d6b', 'n');
subplot(2, 2, 4);                                                                                       % 图(3.4) 【解析解】结构辐射声压虚部(波数)
Fun_MultiPlot2(1, k_Range, imag(P_kRange), '波数', '虚部声压 (Pa)', '图3.4 【解析解】结构辐射声压虚部(波数)', false, false, min(k_Range), max(k_Range), 'n', 'n', 'n', '#5f9ea0', 'n'); 

%% ------------------------------【5 等效简单源辐射声压】------------------------------
%{
    基于波叠加法，求解偶极矩阵、源强度矩阵，进而求解等效压力场
%}
% 1. 以[距离]为步长
% 基于[法向振速]正求[源强度矢量]
[Equiv_D_Radiation, Equiv_Q_Radiation] = Fun_BuildEquiv_D_Q_2D(OriginPoint_Cartesian_Array, EquivSource_Cartesian_Array, k_Origin, OriginPoint_NodeVelocity_Array, OriginPoint_Normal_Vector);
% 基于[源强度矢量]正求节点[辐射声压]
[Equiv_P_Radiation, Equiv_M] = Fun_BuildEquiv_M_P_2D(FieldPoints_Cartesian, EquivSource_Cartesian_Array, k_Origin, AcousticVelocity_Medium, Density_Medium, Equiv_Q_Radiation);

figure;
subplot(2, 2, 1);                                                                                       % 图(4.1) 结构及测点视图(距离)
Fun_MultiPlot2(1, OriginPoint_X_Array, OriginPoint_Y_Array, "X", "Y", "图4.1 结构及测点视图(距离)", true, false, -2 * Circle_Radius, 20 * Circle_Radius, -2 * Circle_Radius, 2 * Circle_Radius, '-');
hold on;
plot(OriginPoint_X_Array, OriginPoint_Y_Array, 'r.', 'MarkerSize', 20);
plot(FieldPoints_Cartesian(:, 1), FieldPoints_Cartesian(:, 2), 'g.', 'MarkerSize', 20);
hold off;

subplot(2, 2, 2);                                                                                       % 图(4.2) 结构辐射声压幅值对比(距离)
plot(FieldPoints_Cartesian(:, 1), abs(P_Distance), 'LineWidth', 1.5);
hold on;
Fun_MultiPlot2(1, FieldPoints_Cartesian(:, 1), abs(Equiv_P_Radiation), "离圆心的距离(m)", "声压幅值(pa)", "图4.2 结构辐射声压幅值对比(距离)", false, false, 10 * Circle_Radius, 20 * Circle_Radius, 'n', 'n', '--', 'r', 1.5);
hold off;
legend('解析解','波叠加法');

subplot(2, 2, 3);                                                                                       % 图(4.3) 结构辐射声压实部对比(距离)
plot(FieldPoints_Cartesian(:, 1), real(P_Distance), 'Color', '#8b7d6b', 'LineWidth', 1.5);
hold on;
Fun_MultiPlot2(1, FieldPoints_Cartesian(:, 1), real(Equiv_P_Radiation), '离圆心的距离 (m)', '实部声压 (Pa)', '图4.3 结构辐射声压实部对比(距离)', false, false, 10 * Circle_Radius, 20 * Circle_Radius, 'n', 'n', '--', 'r', 1.5);
hold off;
legend('解析解','波叠加法');


subplot(2, 2, 4);                                                                                       % 图(4.4) 结构辐射声压虚部对比(距离)
plot(FieldPoints_Cartesian(:, 1), imag(P_Distance), 'Color', '#5f9ea0', 'LineWidth', 1.5);
hold on;
Fun_MultiPlot2(1, FieldPoints_Cartesian(:, 1), imag(Equiv_P_Radiation), '离圆心的距离 (m)', '虚部声压 (Pa)', '图4.4 结构辐射声压虚部对比(距离)', false, false, 10 * Circle_Radius, 20 * Circle_Radius, 'n', 'n', '--', 'r', 1.5); 
hold off;
legend('解析解','波叠加法');

% 2. 以[波数]为步长
Equiv_P_kRange = zeros(size(k_Range));

for i = 1 : length(k_Range)
    k = k_Range(i);
    [Equiv_D_Radiation, Equiv_Q_Radiation] = Fun_BuildEquiv_D_Q_2D(OriginPoint_Cartesian_Array, EquivSource_Cartesian_Array, k, OriginPoint_NodeVelocity_Array, OriginPoint_Normal_Vector);
    [Equiv_P_kRange(i), ~] = Fun_BuildEquiv_M_P_2D(Monitor_Point, EquivSource_Cartesian_Array, k, AcousticVelocity_Medium, Density_Medium, Equiv_Q_Radiation);
end

figure;
subplot(2, 2, 1);                                                                                       % 图(5.1) 结构及测点视图(波数)
Fun_MultiPlot2(1, OriginPoint_X_Array, OriginPoint_Y_Array, "X", "Y", "图5.1 结构及测点视图(波数)", true, true, -2 * Circle_Radius, 2 * Circle_Radius, -2 * Circle_Radius, 2 * Circle_Radius, '-');
hold on;
plot(OriginPoint_X_Array, OriginPoint_Y_Array, 'r.', 'MarkerSize', 20);
plot(Monitor_Point(:, 1), Monitor_Point(:, 2), 'g.', 'MarkerSize', 20);
hold off;

subplot(2, 2, 2);                                                                                       % 图(5.2) 结构辐射声压幅值对比(波数)
plot(k_Range, abs(P_kRange), 'LineWidth', 1.5);
hold on;
Fun_MultiPlot2(1, k_Range, abs(Equiv_P_kRange), "波数", "声压幅值(pa)", "图5.2 结构辐射声压幅值对比(波数)", false, false, min(k_Range), max(k_Range), 'n', 'n', '--', 'r', 1.5);
hold off;
legend('解析解','波叠加法');


subplot(2, 2, 3);                                                                                       % 图(5.3) 结构辐射声压实部对比(波数)
plot(k_Range, real(P_kRange), 'Color', '#8b7d6b', 'LineWidth', 1.5);
hold on;
Fun_MultiPlot2(1, k_Range, real(Equiv_P_kRange), '波数', '实部声压 (Pa)', '图5.3 结构辐射声压实部对比(波数)', false, false, min(k_Range), max(k_Range), 'n', 'n', '--', 'r', 1.5);
hold off;
legend('解析解','波叠加法');

subplot(2, 2, 4);                                                                                       % 图(5.4) 结构辐射声压虚部对比(波数)
plot(k_Range, imag(P_kRange), 'Color', '#5f9ea0', 'LineWidth', 1.5);
hold on;
Fun_MultiPlot2(1, k_Range, imag(Equiv_P_kRange), '波数', '虚部声压 (Pa)', '图5.4 结构辐射声压虚部对比(波数)', false, false, min(k_Range), max(k_Range), 'n', 'n', '--', 'r', 1.5);
hold off;
legend('解析解','波叠加法');