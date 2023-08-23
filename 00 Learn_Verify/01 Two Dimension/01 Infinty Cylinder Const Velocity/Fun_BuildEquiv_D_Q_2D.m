%{
 *------------------------------------------------------------------------------------------
 *---------------------------------------【Fun FILE】---------------------------------------
 * Copyright 流体力学与声学技术实验室
 * ALL right reserved.See COPYRIGHT for detailed Information.
 *
 * @File:       Fun_BuildEquiv_D_Q_2D.m
 * @Brief:      【二维问题】求解[声偶极矩阵]和[源强度矩阵]
 * 
 * @Input:      Origin_Discrete_Node_Coord      【结构源】表面离散X、Y坐标           N×2数组
 *              Equiv_Simple_Source_Coord       【等效简单源】X、Y坐标               M×2数组
 *              k                               波数
 *              U_normal                        【结构源】表面离散节点法向振速幅值    N×1数组
 *              Origin_Discrete_Node_Normal     【结构源】表面离散节点外法向向量      N×2数组
 * 
 * @Output:     Equiv_D                         声偶极矩阵
 *              Equiv_Q                         源强度矩阵
 * 
 * @Author:     Haiger
 * @date:       2023.08.23
 *------------------------------------------------------------------------------------------
%}

function [Equiv_D, Equiv_Q] = Fun_BuildEquiv_D_Q_2D(Origin_Discrete_Node_Coord, Equiv_Simple_Source_Coord, k, U_normal, Origin_Discrete_Node_Normal)

Origin_Discrete_Node_Num = length(Origin_Discrete_Node_Coord);              % 获取【结构源】表面离散节点双轴坐标数组行数，即【结构源】表面离散节点数目M
Equiv_Simple_Source_Num = length(Equiv_Simple_Source_Coord);                % 获取【等效简单源】双轴坐标数组行数，即【等效简单源】数目N
Equiv_D = zeros(Equiv_Simple_Source_Num, Origin_Discrete_Node_Num);         % 初始化【声偶极矩阵】

% for i = 1 : Equiv_Simple_Source_Num
%     r_i = [Equiv_Simple_Source_Coord(i, 1), Equiv_Simple_Source_Coord(i, 2)];                                           % 【等效简单源】上第 i 点的位置矢量
%     for j = 1 : Origin_Discrete_Node_Num
%         r_j = [Origin_Discrete_Node_Coord(j, 1), Origin_Discrete_Node_Coord(j, 2)];                                     % 【结构源】表面上第 j 个离散节点的位置矢量
%         r_difference = r_j - r_i;                                           % 位置矢量差(新的向量)
%         r_distance = norm(r_difference);                                    % 两点间距离
% %         normal_j = r_j / norm(r_j);                                         
%         normal_j = Origin_Discrete_Node_Normal(j, :) / norm(Origin_Discrete_Node_Normal(j, :));                         % 【结构源】表面第 j 个离散节点的法向向量
% %         theta_ij = acos(dot(r_difference, normal_j) / r_distance);          % 位置矢量差(新的向量)与【结构源】表面第 j 个离散节点法向向量的夹角
%         
%         theta_ij = acos(min(1, max(-1, dot(r_difference, normal_j) / r_distance)));
%         Equiv_D(i,j) = - 1 / (4 * pi) * (1j * k * r_distance - 1) / r_distance^2 * exp(1j * k * r_distance) * cos(theta_ij);    % 声偶极矩阵
%     end
% end

for ii = 1 : Origin_Discrete_Node_Num
    r_ii = [Origin_Discrete_Node_Coord(ii, 1), Origin_Discrete_Node_Coord(ii, 2)];                                          % 【结构源】表面上第 ii 个离散节点的位置矢量
    normal_ii = Origin_Discrete_Node_Normal(ii, :) / norm(Origin_Discrete_Node_Normal(ii, :));                              % 【结构源】表面第 j 个离散节点的法向向量
    for jj = 1 : Equiv_Simple_Source_Num
        r_jj = [Equiv_Simple_Source_Coord(jj, 1), Equiv_Simple_Source_Coord(jj, 2)];                                        % 【等效简单源】上第 jj 点的位置矢量
        r_difference = r_ii - r_jj;                                         % 位置矢量差(新的向量)
        r_distance = norm(r_difference);                                    % 两点间距离
        Equiv_D(ii,jj) = (1i * pi * k) / 2 * (besselh(-1, k * r_distance) - besselh(1, k * r_distance)) / r_distance * dot(r_difference, normal_ii);
    end
end

Equiv_Q = Equiv_D \ U_normal;
% Equiv_Q = pinv(Equiv_D) * U_normal;                                                           % 求解源强度矩阵，D 的广义逆矩阵 * U_normal
end