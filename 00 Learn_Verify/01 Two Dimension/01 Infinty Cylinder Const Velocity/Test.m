%{
 *=======================================================================================
 *========================================【M FILE】=====================================
 * Copyright 流体力学与声学技术实验室
 * ALL right reserved.See COPYRIGHT for detailed Information.
 *
 * @File:       Test.m
 * @Brief:      子代码编写测试
 *
 * @Author:     Haiger
 * @date:       2023.08.22
 *=======================================================================================
%}

clc;
clear;
close all;

OriginPoint_R_Array = 1 .* ones(1, 30)';
OriginPoint_Theta_Array = linspace(0, 2 * pi, 30)';
OriginPoint_X_Array = 1 .* cos(OriginPoint_Theta_Array);                                    % 【直角坐标】节点X轴坐标Array
OriginPoint_Y_Array = 1 .* sin(OriginPoint_Theta_Array);                                    % 【直角坐标】节点Y轴坐标Array

Cartesian_2D = [OriginPoint_X_Array, OriginPoint_Y_Array];

Point_Num = length(Cartesian_2D);           % 输入点的数目
Polar_2D = zeros(Point_Num, 2);          % 二维极坐标数组初始化

Polar_2D(:, 1) = atan2(Cartesian_2D(:, 2), Cartesian_2D(:, 1));
Polar_2D(:, 2) = sqrt(Cartesian_2D(:, 1).^2 + Cartesian_2D(:, 2).^2);



polarplot(Polar_2D(:, 1), Polar_2D(:, 2));
