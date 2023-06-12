%{
 *------------------------------------------------------------------------------------------
 *---------------------------------------【Fun FILE】---------------------------------------
 * Copyright 流体力学与声学技术实验室
 * ALL right reserved.See COPYRIGHT for detailed Information.
 *
 * @File:       Fun_ReverseBuild_Q.m
 * @Brief:      
 *              [文献：]
 *              1. 
 * 
 * @Input:      FieldPoints                     声压检测场点坐标                  FieldPoints_Num × 1数组
 *              Equiv_Simple_Source_Coord       【等效简单源】X、Y、Z坐标            M × 3数组
 *              k                               波数
 *              SoundVelocity                   声速
 *              Rho                             密度
 *              FieldPoint_Pressure             场点声压
 * 
 * @Output:     Q                               源强度矩阵
 * 
 * @Author:     Haiger
 * @date:       2023.05.19
 *------------------------------------------------------------------------------------------
%}

function EquivReverse_Q = Fun_ReverseBuild_Q(FieldPoints, Equiv_Simple_Source_Coord, k, SoundVelocity, Rho, FieldPoint_Pressure)
[FieldPoints_Num, ~] = size(FieldPoints);                                   % 声压监测场点数目
Equiv_Simple_Source_Num = length(Equiv_Simple_Source_Coord);                % 获取【等效简单源】三轴坐标数组行数，即【等效简单源】数目N
Omega = SoundVelocity * k;
% 构建系数矩阵
Equiv_M_Matrix = zeros(FieldPoints_Num, Equiv_Simple_Source_Num);
for i = 1 : FieldPoints_Num
    r_i = FieldPoints(i, :);
    for j = 1 : Equiv_Simple_Source_Num
        r_vector = r_i - [Equiv_Simple_Source_Coord(j, 1), Equiv_Simple_Source_Coord(j, 2), Equiv_Simple_Source_Coord(j, 3)];
        r_distance = norm(r_vector);
        Equiv_M_Matrix(i, j) = - 1i * Omega * Rho * exp(1i * k * r_distance) / (4 * pi * r_distance);
    end
end

% 求解线性系统
EquivReverse_Q = Equiv_M_Matrix \ FieldPoint_Pressure;

end