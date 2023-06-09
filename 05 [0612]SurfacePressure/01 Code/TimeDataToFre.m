%{
 *------------------------------------------------------------------------------------------
 *---------------------------------------【Fun FILE】---------------------------------------
 * Copyright 流体力学与声学技术实验室
 * ALL right reserved.See COPYRIGHT for detailed Information.
 *
 * @File:       TimeDataToFre.m
 * @Brief:      对[节点压力]进行傅里叶变换和修正，得到频域下的压力(复数、幅值)
 * 
 * @Input:      Fs                              采样频率
 *              Point_Pressure                  [节点压力]            M × N表格table[double]
 * 
 * @Output:     Data_Spectrum_Complex           频域压力复数          M × N表格table[complex]
 *              Data_Spectrum_ABS               频域压力幅值          M × N表格table[double]
 * 
 * @Author:     Haiger
 * @date:       2023.05.19
 *------------------------------------------------------------------------------------------
%}

function [Data_Spectrum_Complex, Data_Spectrum_ABS] = TimeDataToFre(Fs, Point_Pressure)

OriginPoint_Num = width(Point_Pressure) - 1;                                % 节点数目，第一列为时间列
Point_Pressure_CutMean = Point_Pressure;                                    % Point_Pressure_CutMean表格初始化
Mean_values = mean(Point_Pressure{:, 2:end});                               % 求各个[节点压力]的时均值

for i = 1:length(Mean_values)                                               % 使用每一列的均值来修正数据
    Point_Pressure_CutMean{:, 1+i} = Point_Pressure{:, 1+i} - Mean_values(i);
end

% 初始化一个新的表，用来存储频谱数据
% 列数保持和Point_Pressure_CutMean一样，行数为(Point_Pressure_CutMean行数除以2)
Fre_ColsNum = OriginPoint_Num + 1;
Fre_RowsNum = floor(height(Point_Pressure_CutMean) / 2);
Data_Spectrum_ABS = array2table(nan(Fre_RowsNum, Fre_ColsNum));             % 初始化Data_Spectrum_ABS表格，全设为NaN
Data_Spectrum_Complex = array2table(nan(Fre_RowsNum, Fre_ColsNum));         % 初始化Data_Spectrum_Complex表格，全设为NaN

% 设定列名
Data_Spectrum_ABS.Properties.VariableNames = ["频率", cellstr(strcat("幅值Point_", string(1:OriginPoint_Num)))];
Data_Spectrum_Complex.Properties.VariableNames = ["频率", cellstr(strcat("复值Point_", string(1:OriginPoint_Num)))];

for i = 2 : Fre_ColsNum                                                     % 遍历所有数据行
    spectrum_complex = fft(Point_Pressure_CutMean{:, i});                   % 计算傅里叶变换
    spectrum_complex = spectrum_complex / height(Point_Pressure_CutMean);   % 除以样本点数
    spectrum_complex = spectrum_complex(1:Fre_RowsNum);                     % 由于频谱是对称的，只需要前一半的数据
    spectrum_complex(2 : end - 1 ) = 2 * spectrum_complex(2 : end - 1 );    % 除首尾两个元素，其余元素幅值乘以2，完成修正
    Data_Spectrum_Complex{:, i} = spectrum_complex;
    
    spectrum_abs = abs(spectrum_complex);                                   % 修正后的幅值
    Data_Spectrum_ABS{:, i} = spectrum_abs;
end

f = Fs * (0:Fre_RowsNum-1) / (2 * Fre_RowsNum);                             % 计算频率值得到一个包含频率值的数组
f_ColVector = f.';

Data_Spectrum_ABS{:, 1} = f_ColVector;
Data_Spectrum_Complex{:, 1} = f_ColVector;
end