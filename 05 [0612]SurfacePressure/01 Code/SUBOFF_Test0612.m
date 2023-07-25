%{
 *=======================================================================================
 *========================================【M FILE】=====================================
 * Copyright 流体力学与声学技术实验室
 * ALL right reserved.See COPYRIGHT for detailed Information.
 *
 * @File:       SUBOFF_Main0612.m
 * @Brief:      1. 以【SUBOFF】为研究对象，导入[节点]和[节点时域压力]
 *              2. 基于二分法获取[最佳缩放系数]，以[均方根误差]为判断准则
 *              3. 求解场点频域声压
 *
 * @Author:     Haiger
 * @date:       2023.06.12
 *=======================================================================================
%}

clc;
clear;

%% ------------------------------【1 导入数据 / 输入参数】------------------------------
%{
    输入几何参数及声学参数
%}
OriginPoint_Coord = readtable("Point_Coord.txt");                                       % 读取节点坐标
OriginPoint_Pressure = readtable("p.txt");                                              % 读取节点时域压力

%------------------------------------------------------------------------------------------
% 预处理 1——在[节点坐标]中删除(0,0,0)节点
%------------------------------------------------------------------------------------------
OriginPoint_Coord.Properties.VariableNames = ["节点X坐标", "节点Y坐标", "节点Z坐标"];   % 设置列属性
% 在[节点坐标]中找到三列全为0的行
OriginPoint_Coord_RowsRemoveFor0 = OriginPoint_Coord.("节点X坐标") == 0 & OriginPoint_Coord.("节点Y坐标") == 0 & OriginPoint_Coord.("节点Z坐标") == 0;
OriginPoint_Coord_RowsRemoveFor0_Index = find(OriginPoint_Coord_RowsRemoveFor0);        % 因坐标包含(0,0,0)而需要删除的行索引
OriginPoint_Coord(OriginPoint_Coord_RowsRemoveFor0_Index, :) = [];                      % 删除该行，即某个[节点]

% 对应在[节点时域压力]中删除该[节点]对应的列
OriginPoint_Pressure_ColsRemoveFor0_Index = OriginPoint_Coord_RowsRemoveFor0_Index + 1; % 计算需要删除的列索引，为上述行索引+1，因第1列为时间列
OriginPoint_Pressure(:,OriginPoint_Pressure_ColsRemoveFor0_Index) = [];                 % 删除列

%------------------------------------------------------------------------------------------
% 预处理 2，在[节点时域压力]中找到并删除含有绝对值大于500的节点列
%------------------------------------------------------------------------------------------
OriginPoint_Pressure_BoolForBig500 = abs(OriginPoint_Pressure{:, 2:end}) > 500;         % 创建绝对值大于 500 的布尔矩阵
OriginPoint_Pressure_ColsRemoveForBig500 = any(OriginPoint_Pressure_BoolForBig500, 1);  % 找出至少有一个元素大于 500 的列
OriginPoint_Pressure_ColsRemoveForBig500_Index = find(OriginPoint_Pressure_ColsRemoveForBig500);        % 需要删除的列索引，注意此时的列索引没有加上时间列

OriginPoint_Coord_RowsRemoveForBig500_Index = OriginPoint_Pressure_ColsRemoveForBig500_Index;           % [节点坐标]需删除的行索引
OriginPoint_Coord(OriginPoint_Coord_RowsRemoveForBig500_Index, :) = [];                 % 删除该行，即某个[节点]

OriginPoint_Coord_ColsRemoveForBig500_Index = OriginPoint_Coord_RowsRemoveForBig500_Index + 1;          % [节点时域压力]需删除的列索引，加上了时间列
OriginPoint_Pressure(:, OriginPoint_Coord_ColsRemoveForBig500_Index) = [];              % 删除该列，即某个[节点时域压力]

%------------------------------------------------------------------------------------------
% 预处理 3——在[节点时域压力]中找到并删除含有绝对值大于均值绝对值20%的节点列
%------------------------------------------------------------------------------------------
OriginPoint_Pressure_Mean = mean(OriginPoint_Pressure{:, 2:end});                       % 计算第2列到最后1列的平均值
OriginPoint_Pressure_MeanAbsMult = abs(OriginPoint_Pressure_Mean) * 1.4;                  % 均值的倍数，作为判断标准
OriginPoint_Pressure_BoolForMean = abs(OriginPoint_Pressure{:, 2:end}) > OriginPoint_Pressure_MeanAbsMult; % 建立绝对值大于 均值倍数 的布尔矩阵
% 找出大于 均值倍数 的坐标索引
[OriginPoint_Pressure_RowsRemoveForMean_Index, OriginPoint_Pressure_ColsRemoveForMean_Index] = find(OriginPoint_Pressure_BoolForMean);
OriginPoint_Pressure_ColsRemoveForMean_Index = OriginPoint_Pressure_ColsRemoveForMean_Index + 1;           % 列索引 + 1 ，因第1列是时间列

OriginPoint_Pressure_Delete = OriginPoint_Pressure;                                     % 复制[节点时域压力]表格，用于剔除绝对值大于 均值倍数 的坐标
% 将这些值设为NaN
for i = 1 : length(OriginPoint_Pressure_RowsRemoveForMean_Index)
    OriginPoint_Pressure_Delete{OriginPoint_Pressure_RowsRemoveForMean_Index(i), OriginPoint_Pressure_ColsRemoveForMean_Index(i)} = NaN;
end
OriginPoint_Pressure_MeanDelete = mean(OriginPoint_Pressure_Delete{:, 2:end}, 'omitnan');                % 计算新的均值

OriginPoint_Pressure_New = OriginPoint_Pressure;
% 将这些值设为刚刚所求得的均值
for i = 1 : length(OriginPoint_Pressure_RowsRemoveForMean_Index)
    OriginPoint_Pressure_New{OriginPoint_Pressure_RowsRemoveForMean_Index(i), OriginPoint_Pressure_ColsRemoveForMean_Index(i)} = OriginPoint_Pressure_MeanDelete(OriginPoint_Pressure_ColsRemoveForMean_Index(i) - 1);
end

%------------------------------------------------------------------------------------------
% 预处理 4——在[节点时域压力]中找到并删除全部为0的列
%------------------------------------------------------------------------------------------
OriginPoint_Pressure_NewColsRemoveFor0 = all(OriginPoint_Pressure_New{:, 2:end} == 0, 1);               % 建立[节点时域压力]中全为0的列的布尔矩阵
OriginPoint_Pressure_NewColsRemoveFor0_Index = find(OriginPoint_Pressure_NewColsRemoveFor0);            % 在[节点时域压力]中找到全为0的列索引，注意此时并未考虑时间列

OriginPoint_Coord_NewRowsRemoveFor0_Index = OriginPoint_Pressure_NewColsRemoveFor0_Index;               % 转化为[节点坐标]需删除的行索引
OriginPoint_Coord(OriginPoint_Coord_NewRowsRemoveFor0_Index, :) = [];                                   % [节点坐标]删除对应行，即对应节点坐标

OriginPoint_Pressure_NewColsRemoveFor0_Index = OriginPoint_Pressure_NewColsRemoveFor0_Index + 1;        % 考虑时间列，[节点时域压力]需删除的列索引
OriginPoint_Pressure_New(:, OriginPoint_Pressure_NewColsRemoveFor0_Index) = [];                         % [节点时域压力]删除对应列，即对应节点时域压力

% 将table类型的[节点坐标]转化为数组
OriginPoint_Num = height(OriginPoint_Coord);                                            
OriginPoint_X = OriginPoint_Coord{:, 1};
OriginPoint_Y = OriginPoint_Coord{:, 2};
OriginPoint_Z = OriginPoint_Coord{:, 3};
OriginPoint_Coord_Array = [OriginPoint_X, OriginPoint_Y, OriginPoint_Z];

OriginPoint_Pressure_New.Properties.VariableNames = ["Time", cellstr(strcat("Point_", string(1:OriginPoint_Num)))];
Time_Interval = OriginPoint_Pressure_New{2, 1} - OriginPoint_Pressure_New{1, 1};
Fs = 1 / Time_Interval;                                                                                 % 采样频率

[Data_Spectrum_Complex, Data_Spectrum_ABS] = TimeDataToFre(Fs, OriginPoint_Pressure_New);               % 转换得到各个节点频谱压力，其中傅里叶变换是经过修正的

% 介质参数
Rho_Medium = 1000;                                                                                      % 介质密度，默认为水
SoundVelocity_Medium = 1500;                                                                            % 介质声速，默认为水中的声速

figure;
plot(OriginPoint_Pressure_New{:, 1}, OriginPoint_Pressure_New{:, 2}, 'LineWidth', 1.5);
hold on;
plot(OriginPoint_Pressure_New{:, 1}, OriginPoint_Pressure_New{:, 3}, 'LineWidth', 1.5);
plot(OriginPoint_Pressure_New{:, 1}, OriginPoint_Pressure_New{:, 4}, 'LineWidth', 1.5);
hold off;

%% ------------------------------【2 SUBOFF表面节点可视化】------------------------------
%{
    
%}
figure;
subplot(2, 2, 1);
scatter3(OriginPoint_X, OriginPoint_Y, OriginPoint_Z, 'r', 'filled');
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('SUBOFF表面节点');

%% ------------------------------【3 求解结构表面离散节点外法向量】------------------------------
%{
    
%}
ptCloud = pointCloud(OriginPoint_Coord_Array);                          % 读取点云数据
normals = pcnormals(ptCloud);                                           % 计算法向量（部分点为内法向，部分点为外法向）
centroid = mean(ptCloud.Location, 1);                                   % 计算点云的中心
vectors = centroid - ptCloud.Location;                                  % 计算从每个点到中心的向量
dotProduct = sum(vectors .* normals, 2);                                % 计算向量和法向量的点积
normals(dotProduct > 0, :) = -normals(dotProduct > 0, :);               % 取反指向中心的法向量

subplot(2, 2, 3);
pcshow(ptCloud, 'MarkerSize', 200);                                     % 绘制点云图
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('SUBOFF表面节点点云图');
whitebg('white');

scale = 0.1;                                                            % 法向量缩放系数

% 点云的三轴坐标
x = ptCloud.Location(:, 1);
y = ptCloud.Location(:, 2);
z = ptCloud.Location(:, 3);

% 获取法向量的三轴坐标
u = normals(:, 1);
v = normals(:, 2);
w = normals(:, 3);

subplot(2, 2, 4);
pcshow(ptCloud, 'MarkerSize', 200);
hold on;
quiver3(x, y, z, scale * u, scale * v, scale * w);
hold off;
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('SUBOFF表面节点及外法向');
whitebg('white');

% %% ------------------------------【4 迭代计算均方根误差最小值】------------------------------
% %{
%     通过缩比的方式，获得等效简单源的坐标，等效简单源数目和结构面离散节点数目一致
%     以[缩放系数]Scale_Factor为自变量，使用【最小二乘法】迭代求解均方根误差，得到误差最小的[缩放系数]
% %}
% 
% %{
%     取Data_Spectrum_Complex前11列，即0~1000Hz，压力值设为每行的平均值
% %}
% Data_Spectrum_Complex_Partial = Data_Spectrum_Complex(1 : 11, :);
% Data_Spectrum_Complex_PartialMean = Data_Spectrum_Complex_Partial;
% for i = 2 : 1 : height(Data_Spectrum_Complex_PartialMean)
%     for j = 2 : 1 : width(Data_Spectrum_Complex_PartialMean)
%         Data_Spectrum_Complex_PartialMean{i, j} = mean(Data_Spectrum_Complex_Partial{i, 2 : end});
%     end
% end
% Data_Spectrum_Abs_PartialMean = Data_Spectrum_Complex_PartialMean;
% Data_Spectrum_Abs_PartialMean{2 : end, 2 : end} = abs(Data_Spectrum_Abs_PartialMean{2 : end, 2 : end});
% 
% OriginPoint_Fre_ScaleError = zeros(3, height(Data_Spectrum_Complex_PartialMean) - 1);                                   % 各个频率下(删去0Hz)的最佳缩放系数、最优均方根误差
% 
% % for i = 2 : 1 : height(Data_Spectrum_Complex_PartialMean)                       % 遍历[频率范围]
% for i = 2 : 1 : 2
%     UpperBound = 1.0;               % 上界
%     LowerBound = 0.9;               % 下界
%     iterations = 10;                % 迭代次数
%     Scale_Factor_Opt = 0.0;         % 最佳缩放系数
%     Error_RMS_Opt = Inf;            % 最优均方根误差0
% 
%     Calculate_Fre = Data_Spectrum_Complex_PartialMean{i, 1};                    % 计算频率，对应Data_Spectrum_Complex的第i行第1列
%     Calculate_Row = i;                                              % 计算行数，针对Data_Spectrum_Complex数组
%     Calculate_OriginPoint_Complex = Data_Spectrum_Complex_PartialMean{Calculate_Row, 2 : end}.';                        % 该[计算频率]下各个[节点压力]复数
%     Calculate_OriginPoint_ABS = Data_Spectrum_Abs_PartialMean{Calculate_Row, 2 : end}.';                                % 该[计算频率]下各个[节点压力]幅值
%     Calculate_Omega = 2 * pi * Calculate_Fre;                       % 该[计算频率]下[圆频率]
%     Calculate_k = Calculate_Omega / SoundVelocity_Medium;           % 该[计算频率]下[波数]
%     
%     disp("------------------------------------------------------------");
%     fprintf('Frequency %.2f\n', Calculate_Fre);
%     fprintf('\n');
% 
%     for j = 1 : 1 : iterations                                      % 遍历[迭代次数]
%         Scale_Factor = (UpperBound + LowerBound) / 2;               % 取上界和下界的中点作为缩放系数
%         
%         Equiv_Simple_Source_X = Scale_Factor .* OriginPoint_X;      % 【等效简单源】X坐标
%         Equiv_Simple_Source_Y = Scale_Factor .* OriginPoint_Y;      % 【等效简单源】X坐标
%         Equiv_Simple_Source_Z = Scale_Factor .* OriginPoint_Z;      % 【等效简单源】X坐标
%         Equiv_Simple_Source_Coord = [Equiv_Simple_Source_X, Equiv_Simple_Source_Y, Equiv_Simple_Source_Z];  % 【等效简单源】X、Y、Z坐标
%         
%         % 基于[节点压力]反求[源强度矢量]
%         EquivReverse_Q = Fun_ReverseBuild_Q(OriginPoint_Coord_Array, Equiv_Simple_Source_Coord, Calculate_k, SoundVelocity_Medium, Rho_Medium, Calculate_OriginPoint_Complex);
%         % 基于[源强度矢量]反求[法向振速]
%         [EquivReverse_D, EquivReverse_U_normal] = Fun_ReverseBuild_D_U(OriginPoint_Coord_Array, Equiv_Simple_Source_Coord, Calculate_k, EquivReverse_Q, normals);
%         % 基于[法向振速]正求[源强度矢量]
%         [Equiv_D_Radiation, Equiv_Q_Radiation] = Fun_BuildEquiv_D_Q(OriginPoint_Coord_Array, Equiv_Simple_Source_Coord, Calculate_k, EquivReverse_U_normal, normals);
%         % 基于[源强度矢量]正求节点[辐射声压]
%         Equiv_P_Radiation = Fun_BuildEquiv_M_P(OriginPoint_Coord_Array, Equiv_Simple_Source_Coord, Calculate_k, SoundVelocity_Medium, Rho_Medium, Equiv_Q_Radiation);
%     
%         Equiv_P_Radiation_ABS = abs(Equiv_P_Radiation);             % 节点[辐射声压]幅值
%         Error_RMS = sqrt(mean((Equiv_P_Radiation_ABS - Calculate_OriginPoint_ABS).^2));                     % [均方根误差]
%         
%         disp("----------------------------");
%         fprintf('iteration %d\n', j);
%         fprintf('Scale_Factor %f\n', Scale_Factor);
%         fprintf('Error_RMS %e\n', Error_RMS);
%         fprintf('Min %e\n', min(Equiv_P_Radiation_ABS));
%         fprintf('Max %e\n', max(Equiv_P_Radiation_ABS));
%         fprintf('Mean %e\n', mean(Equiv_P_Radiation_ABS));
%         disp("----------------------------");
%     
%         if (Error_RMS < Error_RMS_Opt) && (max(Equiv_P_Radiation_ABS) > max(Calculate_OriginPoint_ABS))     % 判断标准
%             LowerBound = Scale_Factor;                              % 更换界限
%             Error_RMS_Opt = Error_RMS;
%             Scale_Factor_Opt = Scale_Factor;
%         else
%             UpperBound = Scale_Factor;
%             Error_RMS_Opt = Error_RMS;
%             Scale_Factor_Opt = Scale_Factor;
%         end
%     end
%     
%     fprintf('\n');
%     fprintf('Scale_Factor_Opt %f\n', Scale_Factor_Opt);
%     fprintf('Error_RMS_Opt %e\n', Error_RMS_Opt);
%     disp("------------------------------------------------------------");
% 
%     % 存储各个频率下(删去0Hz)的最佳缩放系数、最优均方根误差
%     OriginPoint_Fre_ScaleError(1, i - 1) = Calculate_Fre;           % 计算频率
%     OriginPoint_Fre_ScaleError(2, i - 1) = Scale_Factor_Opt;        % 最佳缩放系数
%     OriginPoint_Fre_ScaleError(3, i - 1) = Error_RMS_Opt;           % 最优均方根误差
% 
% end

%% ------------------------------【4 计算当前缩放系数下表面压力】------------------------------
% 给定【缩放系数】，即距离
Scale_Factor = 0.995313;

% 辐射声压检测场点的位置
% FieldPoints = [0.5004, 1.07, 0];
FieldPoints = [0.206, 0, 0];

FieldPoints_Fre_P = zeros(3, height(Data_Spectrum_Complex));

Equiv_Simple_Source_X = Scale_Factor .* OriginPoint_X;      % 【等效简单源】X坐标
Equiv_Simple_Source_Y = Scale_Factor .* OriginPoint_Y;      % 【等效简单源】X坐标
Equiv_Simple_Source_Z = Scale_Factor .* OriginPoint_Z;      % 【等效简单源】X坐标
Equiv_Simple_Source_Coord = [Equiv_Simple_Source_X, Equiv_Simple_Source_Y, Equiv_Simple_Source_Z];  % 【等效简单源】X、Y、Z坐标

% 给定简单源【等效源强度矢量矩阵】
EquivReverse_Q_Scale = 5*10^(-4);
for i = 2 : 1 : height(Data_Spectrum_Complex)
% for i = 2 : 1 : 10
    Calculate_Fre_Current = Data_Spectrum_Complex{i, 1};

    disp("------------------------------------------------------------");
    fprintf('Frequency %.2f\n', Calculate_Fre_Current);
    disp("------------------------------------------------------------");
    
    Calculate_OriginPoint_ComplexCurrent = Data_Spectrum_Complex{i, 2 : end}.';                        % 该[计算频率]下各个[节点压力]复数
    Calculate_OriginPoint_ABSCurrent = Data_Spectrum_ABS{i, 2 : end}.';                                % 该[计算频率]下各个[节点压力]幅值
    Calculate_Omega_Current = 2 * pi * Calculate_Fre_Current;                       % 该[计算频率]下[圆频率]
    Calculate_k_Current = Calculate_Omega_Current / SoundVelocity_Medium;           % 该[计算频率]下[波数]

    EquivReverse_Q = EquivReverse_Q_Scale .* Data_Spectrum_Complex{i , 2:end}';
    
    %------------------------------------------------------------------------------------------
    % 预处理 ——    在[等效源强度矩阵]中找到幅值大于幅值均值20%的索引
    %              大于1.2的取均值的1.2倍
    %              小于0.8的取均值的0.8倍
    %------------------------------------------------------------------------------------------
    
    Mean_Q_ABS = mean(abs(EquivReverse_Q));
    Mean_Q_Com = mean((EquivReverse_Q));
    for j = 1 : 1 : height(EquivReverse_Q)
        if (abs(EquivReverse_Q(j)) / Mean_Q_ABS >= 1.5)
            EquivReverse_Q(j) = 1.5 * Mean_Q_Com;
        elseif (abs(EquivReverse_Q(j)) / Mean_Q_ABS <= 0.5)
            EquivReverse_Q(j) = 0.5 * Mean_Q_Com;
        else
            EquivReverse_Q(j) = EquivReverse_Q(j);
        end
    end

    % 基于[源强度矢量]反求[法向振速]
    [EquivReverse_D, EquivReverse_U_normal] = Fun_ReverseBuild_D_U(OriginPoint_Coord_Array, Equiv_Simple_Source_Coord, Calculate_k_Current, EquivReverse_Q, normals);
    % 基于[法向振速]正求[源强度矢量]
    [Equiv_D_Radiation, Equiv_Q_Radiation] = Fun_BuildEquiv_D_Q(OriginPoint_Coord_Array, Equiv_Simple_Source_Coord, Calculate_k_Current, EquivReverse_U_normal, normals);
    % 基于[源强度矢量]正求节点[辐射声压]
    Equiv_P_Radiation = Fun_BuildEquiv_M_P(FieldPoints, Equiv_Simple_Source_Coord, Calculate_k_Current, SoundVelocity_Medium, Rho_Medium, Equiv_Q_Radiation);

    Equiv_P_Radiation_ABS = abs(Equiv_P_Radiation);             % 节点[辐射声压]幅值

    FieldPoints_Fre_P(1, i) = Calculate_Fre_Current;% 频率，排除了0Hz
    FieldPoints_Fre_P(2, i) = Equiv_P_Radiation;% 复数声压
    FieldPoints_Fre_P(3, i) = 20 * log10(Equiv_P_Radiation_ABS / (10^(-6)));% 幅值声压


end

figure;
plot(FieldPoints_Fre_P(1, :), FieldPoints_Fre_P(3, :), 'LineWidth', 1.5);
xlabel('频率');
ylabel('声压幅值(Pa)');
title('测点声压幅值');

figure;
subplot(1, 2, 1);
plot(FieldPoints_Fre_P(1, :), real(FieldPoints_Fre_P(2, :)), 'LineWidth', 1.5);
xlabel('频率');
ylabel('实部声压(Pa)');
title('测点实部声压');

subplot(1, 2, 2);
plot(FieldPoints_Fre_P(1, :), imag(FieldPoints_Fre_P(2, :)), 'LineWidth', 1.5);
xlabel('频率');
ylabel('虚部声压(Pa)');
title('测点虚部声压');
