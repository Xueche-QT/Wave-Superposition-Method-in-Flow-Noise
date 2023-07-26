%{
 *------------------------------------------------------------------------------------------
 *---------------------------------------【Fun FILE】---------------------------------------
 * Copyright 流体力学与声学技术实验室
 * ALL right reserved.See COPYRIGHT for detailed Information.
 *
 * @File:       Fun_BuildEquiv_M_P.m
 * @Brief:      
 *              [文献：]
 *              1. 
 * 
 * @Input:      FieldPoints                     X轴声压检测场点坐标                  FieldPoints_Num × 1数组
 *              Equiv_Simple_Source_Coord       【等效简单源】X、Y、Z坐标            M × 3数组
 *              k                               波数
 *              SoundVelocity                   声速
 *              Rho                             密度
 * 
 * @Output:     Equiv_M                         系数矩阵
 *              Q                               源强度矩阵
 * 
 * @Author:     Haiger
 * @date:       2023.05.01
 *------------------------------------------------------------------------------------------
%}

function [Equiv_P, Equiv_M] = Fun_BuildEquiv_M_P(FieldPoints, Equiv_Simple_Source_Coord, k, SoundVelocity, Rho, Equiv_Q)

[FieldPoints_Num, ~] = size(FieldPoints);                                   % 声压监测场点数目
Equiv_Simple_Source_Num = length(Equiv_Simple_Source_Coord);                % 获取【等效简单源】三轴坐标数组行数，即【等效简单源】数目N
Omega = SoundVelocity * k;                                                  % 初始化测点声压矩阵
Equiv_P = zeros(FieldPoints_Num, 1);
Equiv_M = zeros(Equiv_Simple_Source_Num, FieldPoints_Num);                  % 初始化【系数矩阵】
for i = 1 : FieldPoints_Num
    r_i = FieldPoints(i, :);                                                % 声压监测点坐标
    for j = 1 : Equiv_Simple_Source_Num
        r_vector = r_i - [Equiv_Simple_Source_Coord(j, 1), Equiv_Simple_Source_Coord(j, 2), Equiv_Simple_Source_Coord(j, 3)];     % 观测点矢量与等效源矢量差(新的向量)
        r_distance = norm(r_vector);                                        % 观测点与等效源的距离(新的向量)的模
        Equiv_M(j, i) = - 1i * Omega * Rho * exp(1i * k * r_distance) / (4 * pi * r_distance);                                               % 系数矩阵M
    end
    Equiv_P(i) = sum(Equiv_M(:, i) .* Equiv_Q);
end

end