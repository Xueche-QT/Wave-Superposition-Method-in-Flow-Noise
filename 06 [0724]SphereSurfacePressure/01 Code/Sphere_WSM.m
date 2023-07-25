%{
 *=======================================================================================
 *========================================【M FILE】=====================================
 * Copyright 流体力学与声学技术实验室
 * ALL right reserved.See COPYRIGHT for detailed Information.
 *
 * @File:       Sphere_WSM.m
 * @Brief:      1. 以【Sphere】为研究对象，导入[节点坐标]和[节点某时刻压力(稳态)]
 *              2. 基于二分法获取[最佳缩放系数]，以[均方根误差]为判断准则
 *              3. 求解场点频域声压
 *
 * @Author:     Haiger
 * @date:       2023.07.24
 *=======================================================================================
%}

clc;
close all;
clear;

%% ------------------------------【1 导入数据 / 输入参数】------------------------------
%{
    输入几何参数及声学参数
%}
OriginPoint_Data_Table = readtable("..\02 Files\Wall_Pressure_1100", 'VariableNamingRule', 'preserve');         % 读取节点索引、坐标和静压数据Table
OriginPoint_Coord_Table = OriginPoint_Data_Table(:, 2:4);                                                       % 节点三轴坐标Table

OriginPoint_X_Array = OriginPoint_Coord_Table{:, 1};                                                            % 节点X轴坐标Array
OriginPoint_Y_Array = OriginPoint_Coord_Table{:, 2};                                                            % 节点Y轴坐标Array
OriginPoint_Z_Array = OriginPoint_Coord_Table{:, 3};                                                            % 节点Z轴坐标Array

OriginPoint_Coord_Array = [1000 .* OriginPoint_X_Array, 1000 .* OriginPoint_Y_Array, 1000 .* OriginPoint_Z_Array];  % 节点三轴坐标Array

OriginPoint_X_Array = OriginPoint_Coord_Array(:, 1);                                                            % 节点X轴坐标Array
OriginPoint_Y_Array = OriginPoint_Coord_Array(:, 2);                                                            % 节点Y轴坐标Array
OriginPoint_Z_Array = OriginPoint_Coord_Array(:, 3);                                                            % 节点Z轴坐标Array

OriginPoint_Pressure_Array = OriginPoint_Data_Table{:, 5};                                                      % 节点压力Array

% 介质参数
Rho_Medium = 1;                                                                                                 % 介质密度
SoundVelocity_Medium = 340;                                                                                     % 介质声速

%% ------------------------------【2 结构表面节点可视化】------------------------------
%{
    结构表面节点可视化
%}
figure;
subplot(2, 2, 1);
scatter3(OriginPoint_X_Array, OriginPoint_Y_Array, OriginPoint_Z_Array, 'r', 'filled');
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('结构表面离散节点');

%% ------------------------------【3 求解结构表面离散节点外法向量】------------------------------
%{
    基于点云求解结构表面离散节点外法向
%}
ptCloud = pointCloud(OriginPoint_Coord_Array);                          % 读取点云数据
normals = pcnormals(ptCloud);                                           % 计算法向量（部分点为内法向，部分点为外法向）
centroid = mean(ptCloud.Location, 1);                                   % 计算点云的中心
vectors = centroid - ptCloud.Location;                                  % 计算从每个点到中心的向量
dotProduct = sum(vectors .* normals, 2);                                % 计算向量和法向量的点积
normals(dotProduct > 0, :) = -normals(dotProduct > 0, :);               % 取反指向中心的法向量

subplot(2, 2, 3);
pcshow(ptCloud, 'MarkerSize', 200);                                     % 绘制点云图
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('结构表面节点点云图');
whitebg('white');

scale = 0.1;                                                            % 法向量缩放系数

% 点云的三轴坐标
x = ptCloud.Location(:, 1);
y = ptCloud.Location(:, 2);
z = ptCloud.Location(:, 3);

% 获取法向量的三轴坐标
u = normals(:, 1);
v = normals(:, 2);
w = normals(:, 3);

subplot(2, 2, 4);
pcshow(ptCloud, 'MarkerSize', 200);
hold on;
quiver3(x, y, z, scale * u, scale * v, scale * w);
hold off;
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('SUBOFF表面节点及外法向');
whitebg('white');

%% ------------------------------【4 方法1：迭代计算均方根误差最小值】------------------------------
%{
    通过缩比的方式，获得等效简单源的坐标，等效简单源数目和结构面离散节点数目一致
    以[缩放系数]Scale_Factor为自变量，使用【最小二乘法】迭代求解均方根误差，得到误差最小的[缩放系数]
%}

Calculate_Fre = 100;                                            % [计算频率]
Calculate_Omega = 2 * pi * Calculate_Fre;                       % 该[计算频率]下[圆频率]
Calculate_k = Calculate_Omega / SoundVelocity_Medium;           % 该[计算频率]下[波数]

UpperBound = 1.0;               % 上界
LowerBound = 0.0;               % 下界
iterations = 10;                 % 迭代次数
Scale_Factor_Opt = 0.0;         % 初始化最佳缩放系数
Error_RMS_Opt = Inf;            % 初始化最优均方根误差

disp("------------------------------------------------------------");
fprintf('Frequency %.2f\n', Calculate_Fre);
fprintf('\n');

for j = 1 : 1 : iterations                                                      % 遍历[迭代次数]
    Scale_Factor = (UpperBound + LowerBound) / 2;                           % 取上界和下界的中点作为缩放系数
%     Scale_Factor = 0.58;                           % 取上界和下界的中点作为缩放系数
    Equiv_Simple_Source_X_Array = Scale_Factor .* OriginPoint_X_Array;      % 【等效简单源】X坐标
    Equiv_Simple_Source_Y_Array = Scale_Factor .* OriginPoint_Y_Array;      % 【等效简单源】X坐标
    Equiv_Simple_Source_Z_Array = Scale_Factor .* OriginPoint_Z_Array;      % 【等效简单源】X坐标
    Equiv_Simple_Source_Coord_Array = [Equiv_Simple_Source_X_Array, Equiv_Simple_Source_Y_Array, Equiv_Simple_Source_Z_Array];  % 【等效简单源】X、Y、Z坐标

    % 基于[节点压力]反求[源强度矢量]
    EquivReverse_Q = Fun_ReverseBuild_Q(OriginPoint_Coord_Array, Equiv_Simple_Source_Coord_Array, Calculate_k, SoundVelocity_Medium, Rho_Medium, OriginPoint_Pressure_Array);
    % 基于[源强度矢量]反求[法向振速]
    [EquivReverse_D, EquivReverse_U_normal] = Fun_ReverseBuild_D_U(OriginPoint_Coord_Array, Equiv_Simple_Source_Coord_Array, Calculate_k, EquivReverse_Q, normals);
    % 基于[法向振速]正求[源强度矢量]
    [Equiv_D_Radiation, Equiv_Q_Radiation] = Fun_BuildEquiv_D_Q(OriginPoint_Coord_Array, Equiv_Simple_Source_Coord_Array, Calculate_k, EquivReverse_U_normal, normals);
    % 基于[源强度矢量]正求节点[辐射声压]
    [Equiv_P_Radiation, Equiv_M] = Fun_BuildEquiv_M_P(OriginPoint_Coord_Array, Equiv_Simple_Source_Coord_Array, Calculate_k, SoundVelocity_Medium, Rho_Medium, Equiv_Q_Radiation);

    Equiv_P_Radiation_ABS = abs(Equiv_P_Radiation);                         % 节点[辐射声压]幅值
    Error_RMS = sqrt(mean((Equiv_P_Radiation_ABS - OriginPoint_Pressure_Array).^2));                     % [均方根误差]

    disp("----------------------------");
    fprintf('iteration %d\n', j);
    fprintf('Scale_Factor %f\n', Scale_Factor);
    fprintf('Error_RMS %e\n', Error_RMS);
    fprintf('Min %e\n', min(Equiv_P_Radiation_ABS));
    fprintf('Max %e\n', max(Equiv_P_Radiation_ABS));
    fprintf('Mean %e\n', mean(Equiv_P_Radiation_ABS));
    disp("----------------------------");

    if (Error_RMS < Error_RMS_Opt) && (max(Equiv_P_Radiation_ABS) > max(OriginPoint_Pressure_Array))     % 判断标准
        LowerBound = Scale_Factor;                              % 更换界限
        if min(Equiv_P_Radiation_ABS) > min(OriginPoint_Pressure_Array)
            Error_RMS_Opt = Error_RMS;
            Scale_Factor_Opt = Scale_Factor;
        end
    else
        UpperBound = Scale_Factor;
        if min(Equiv_P_Radiation_ABS) > min(OriginPoint_Pressure_Array)
            Error_RMS_Opt = Error_RMS;
            Scale_Factor_Opt = Scale_Factor;
        end
    end
end

fprintf('\n');
fprintf('Scale_Factor_Opt %f\n', Scale_Factor_Opt);
fprintf('Error_RMS_Opt %e\n', Error_RMS_Opt);
disp("------------------------------------------------------------");
