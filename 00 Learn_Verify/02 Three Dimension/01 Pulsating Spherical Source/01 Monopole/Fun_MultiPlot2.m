%{
 *------------------------------------------------------------------------------------------
 *---------------------------------------【Fun FILE】---------------------------------------
 * Copyright 流体力学与声学技术实验室
 * ALL right reserved.See COPYRIGHT for detailed Information.
 *
 * @File:       Fun_MultiPlot2.m
 * @Brief:      选择不同类型的二维图像绘制
 * 
 * @Input:      Plot_Selection                  绘图函数选择                          int
 *                  1   plot                    坐标轴为线性
 *                  2   semilogx                横坐标为对数
 *              X_Series                        横轴数据                              N×1double
 *              Y_Series                        纵轴数据                              N×1double
 *              XLabel                          横轴标签                              字符串
 *              YLabel                          纵轴标签                              字符串
 *              Title                           图谱标题                              字符串
 *              Grid_Option                     判断是否加载网格线                     bool
 *                  true    grid on             加载网格线
 *                  flase   grid off            不加载网格线
 *              X_Min                           横轴最小值【可选】                     double
 *              X_Max                           横轴最大值【可选】                     double
 *              Y_Min                           纵轴最小值【可选】                     double
 *              Y_Max                           纵轴最大值【可选】                     double
 *              Linetype                        线型                                  字符串
 *                  '-'                         实线
 *                  '--'                        双划线
 *                  ':'                         虚线
 *                  ':.'                        点划线
 *
 *                  '+'                         加号符
 *                  'o'                         空心圆
 *                  '*'                         星号
 *                  '.'                         实心圆
 *                  'x'                         加号符
 *                  's'                         正方形
 *                  'd'                         菱形
 *                  '^'                         上三角形
 *                  'v'                         下三角形
 *                  '>'                         右三角形
 *                  '<'                         左三角形
 *                  'p'                         五角星
 *                  'h'                         六边形
 *              LineColor                       线条颜色                              字符串
 *                  '#8b7d6b'
 *              LineWidth                       线条粗细                              double
 *                  1.5
 * 
 * @Output:     
 * 
 * @Author:     Haiger
 * @date:       2023.07.26
 *------------------------------------------------------------------------------------------
%}

function [] = Fun_MultiPlot2(Plot_Selection, X_Series, Y_Series, XLabel, YLabel, Title, Grid_Option, X_Min, X_Max, Y_Min, Y_Max, Linetype, LineColor, LineWidth)

if ~exist('Linetype', 'var') || strcmp(Linetype, 'n')
    Linetype = '-';
end
if ~exist('LineColor', 'var') || strcmp(LineColor, 'n')
    LineColor = '#0d79c0';
end
if ~exist('LineWidth', 'var') || strcmp(LineWidth, 'n')
    LineWidth = 1.5;
end

switch Plot_Selection
    case 1
        plot(X_Series, Y_Series, Linetype, 'Color', LineColor, 'LineWidth', LineWidth);
    case 2
        semilogx(X_Series, Y_Series, Linetype, 'Color', LineColor, 'LineWidth', LineWidth);
end

    xlabel(XLabel);
    ylabel(YLabel);
    title(Title);
    if Grid_Option
        grid on;
    end
    % 判断是否输入了X、Y轴的限定范围
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
    
    if X_Min_Flag && X_Max_Flag                                         % 当【横轴】最大值和最小值均指定时，才会进行坐标轴限制
        xlim([X_Min X_Max]);
    end
    if Y_Min_Flag && Y_Max_Flag
        ylim([Y_Min Y_Max]);
    end
   
end
