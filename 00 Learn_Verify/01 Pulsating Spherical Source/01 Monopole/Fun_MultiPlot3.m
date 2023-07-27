%{
 *------------------------------------------------------------------------------------------
 *---------------------------------------【Fun FILE】---------------------------------------
 * Copyright 流体力学与声学技术实验室
 * ALL right reserved.See COPYRIGHT for detailed Information.
 *
 * @File:       Fun_MultiPlot3.m
 * @Brief:      选择不同类型的三维图像绘制
 * 
 * @Input:      Plot_Selection                  绘图函数选择                          int
 *                  1   scatter3                散点图
 *                  2   surf                    曲面图
 *              X_Series                        横轴数据                              N×1double
 *              Y_Series                        纵轴数据                              N×1double
 *              Z_Series                        竖轴数据                              N×1double
 *              XLabel                          横轴标签                              字符串
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
 *              Size                            散点大小【可选】                      double
 *              Color                           散点或曲面颜色【可选】                字符串
 *              MarketType                      散点形状【可选】                      字符串
 *                  '+'                         加号符
 *                  'o'                         空心圆
 *                  '*'                         星号
 *                  '.'                         点
 *                  'x'                         加号符
 *                  '_'                         水平线条
 *                  '|'                         垂直线条
 *                  'square'                    正方形
 *                  'diamond'                   菱形
 *                  '^'                         上三角形
 *                  'v'                         下三角形
 *                  '>'                         右三角形
 *                  '<'                         左三角形
 *                  'pentagram'                 五角星
 *                  'hexagram'                  六边形
 *                  'none'                      无标记
 *              Filled_Option                   是否填充散点【可选】                  bool
 *                  true    'filled'            填充
 *                  false                       不填充
 *
 * @Output:     
 * 
 * @Author:     Haiger
 * @date:       2023.07.26
 *------------------------------------------------------------------------------------------
%}

function [] = Fun_MultiPlot3(Plot_Selection, X_Series, Y_Series, Z_Series, XLabel, YLabel, ZLabel, Title, Axis_Option, X_Min, X_Max, Y_Min, Y_Max, Z_Min, Z_Max, Size, Color, MarketType, Filled_Option)

if ~exist('Size', 'var') || strcmp(Size, 'n')
    Size = 36;
end
if ~exist('Color', 'var') || strcmp(Color, 'n')
    Color = '#0d79c0';
end
if ~exist('MarketType', 'var') || strcmp(MarketType, 'n')
    MarketType = 'o';
end
if ~exist('Filled_Option', 'var') || strcmp(Filled_Option, 'n') 
    Filled_Option = false;
end

switch Plot_Selection
    case 1
        if Filled_Option == true
            scatter3(X_Series, Y_Series, Z_Series, Size, Color, MarketType, 'filled');
        else
            scatter3(X_Series, Y_Series, Z_Series, Size, Color, MarketType);
        end
    case 2
        surf(X_Series, Y_Series, Z_Series);
end

    xlabel(XLabel);
    ylabel(YLabel);
    zlabel(ZLabel);
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
