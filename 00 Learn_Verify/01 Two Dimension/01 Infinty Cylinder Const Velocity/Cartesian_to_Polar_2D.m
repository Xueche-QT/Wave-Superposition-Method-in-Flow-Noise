%{
 *------------------------------------------------------------------------------------------
 *---------------------------------------【Fun FILE】---------------------------------------
 * Copyright 流体力学与声学技术实验室
 * ALL right reserved.See COPYRIGHT for detailed Information.
 *
 * @File:       Cartesian_to_Polar_2D.m
 * @Brief:      [二维]笛卡尔坐标转换为极坐标
 * 
 * @Input:      Cartesian_2D                    二维笛卡尔坐标                         N×2 double
 *
 * @Output:     Polar_2D                        二维极坐标                            N×2 double 
 * 
 * @Author:     Haiger
 * @date:       2023.08.22
 *------------------------------------------------------------------------------------------
%}

function Polar_2D = Cartesian_to_Polar_2D(Cartesian_2D)

Point_Num = length(Cartesian_2D);                                           % 输入点的数目
Polar_2D = zeros(Point_Num, 2);                                             % 二维【极坐标】数组初始化

Polar_2D(:, 1) = atan2(Cartesian_2D(:, 2), Cartesian_2D(:, 1));             % 【极坐标】角坐标
Polar_2D(:, 2) = sqrt(Cartesian_2D(:, 1).^2 + Cartesian_2D(:, 2).^2);       % 【极坐标】半径坐标

end