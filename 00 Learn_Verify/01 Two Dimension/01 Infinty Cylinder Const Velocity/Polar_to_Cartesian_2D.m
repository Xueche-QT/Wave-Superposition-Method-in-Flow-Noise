%{
 *------------------------------------------------------------------------------------------
 *---------------------------------------【Fun FILE】---------------------------------------
 * Copyright 流体力学与声学技术实验室
 * ALL right reserved.See COPYRIGHT for detailed Information.
 *
 * @File:       Polar_to_Cartesian_2D.m
 * @Brief:      [二维]极坐标转换为笛卡尔坐标
 * 
 * @Input:      Polar_2D                        二维极坐标                            N×2 double
 *
 * @Output:     Cartesian_2D                    二维笛卡尔坐标                         N×2 double
 * 
 * @Author:     Haiger
 * @date:       2023.08.22
 *------------------------------------------------------------------------------------------
%}

function Cartesian_2D = Polar_to_Cartesian_2D(Polar_2D)

Point_Num = length(Polar_2D);                                       % 输入点的数目
Cartesian_2D = zeros(Point_Num, 2);                                 % 二维【直角】坐标数组初始化

Cartesian_2D(:, 1) = Polar_2D(:, 1) .* cos(Polar_2D(:, 2));         % 二维【直角】X坐标
Cartesian_2D(:, 2) = Polar_2D(:, 1) .* sin(Polar_2D(:, 2));         % 二维【直角】Y坐标

end