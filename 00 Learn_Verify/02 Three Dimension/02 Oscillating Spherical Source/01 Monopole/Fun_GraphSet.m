%{
 *------------------------------------------------------------------------------------------
 *---------------------------------------【Fun FILE】---------------------------------------
 * Copyright 流体力学与声学技术实验室
 * ALL right reserved.See COPYRIGHT for detailed Information.
 *
 * @File:       Fun_GraphSet.m
 * @Brief:      二维或三维图像属性设置
 * 
 * @Input:      XLabel                          横轴标签                              字符串
 *              YLabel                          纵轴标签                              字符串
 *              ZLabel                          竖轴标签                              字符串
 *              Title                           图谱标题                              字符串
 *              Axis_Option                     判断是否加载网格线                    bool
 *                  true    axis equal          三轴缩放系数相同
 *                  flase                       三轴缩放系数默认
 *              X_Min                           横轴最小值【可选】                    double
 *              X_Max                           横轴最大值【可选】                    double
 *              Y_Min                           纵轴最小值【可选】                    double
 *              Y_Max                           纵轴最大值【可选】                    double
 *              Z_Min                           竖轴最小值【可选】                    double
 *              Z_Max                           竖轴最大值【可选】                    double
 *
 * @Output:     
 * 
 * @Author:     Haiger
 * @date:       2023.07.26
 *------------------------------------------------------------------------------------------
%}

function [] = Fun_GraphSet(XLabel, YLabel, ZLabel, Title, Axis_Option, X_Min, X_Max, Y_Min, Y_Max, Z_Min, Z_Max)

    xlabel(XLabel);
    ylabel(YLabel);

    if ~exist('ZLabel', 'var') || strcmp(ZLabel, 'n')                     % 未指定【横轴最小值】
        
    else
        zlabel(ZLabel);
    end

    title(Title);
    if Axis_Option
        axis equal;
    end
    % 判断是否输入了X、Y、Z轴的限定范围

    if ~exist('X_Min', 'var') || strcmp(X_Min, 'n')                     % 未指定【横轴最小值】
        X_Min_Flag = 0;
    else
        X_Min_Flag = 1;                                                 % 指定了【横轴最小值】
    end
    if ~exist('X_Max', 'var') || strcmp(X_Max, 'n')                     % 未指定【横轴最大值】
        X_Max_Flag = 0;
    else                                                                % 指定了【横轴最大值】
        X_Max_Flag = 1;
    end
    if ~exist('Y_Min', 'var') || strcmp(Y_Min, 'n')
        Y_Min_Flag = 0;
    else
        Y_Min_Flag = 1;
    end
    if ~exist('Y_Max', 'var') || strcmp(Y_Max, 'n')
        Y_Max_Flag = 0;
    else
        Y_Max_Flag = 1;
    end
    if ~exist('Z_Min', 'var') || strcmp(Z_Min, 'n')
        Z_Min_Flag = 0;
    else
        Z_Min_Flag = 1;
    end
    if ~exist('Z_Max', 'var') || strcmp(Z_Max, 'n')
        Z_Max_Flag = 0;
    else
        Z_Max_Flag = 1;
    end

    if X_Min_Flag && X_Max_Flag                                         % 当【横轴】最大值和最小值均指定时，才会进行坐标轴限制
        xlim([X_Min X_Max]);
    end
    if Y_Min_Flag && Y_Max_Flag
        ylim([Y_Min Y_Max]);
    end
    if Z_Min_Flag && Z_Max_Flag
        zlim([Z_Min Z_Max]);
    end

end