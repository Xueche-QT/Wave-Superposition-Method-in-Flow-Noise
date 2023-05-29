clc;
clear;

%% ------------------------------【1 导入数据 / 输入参数】------------------------------
%{
    输入几何参数及声学参数
%}
OriginPoint_Coord = readtable("Point_Coord.txt");
OriginPoint_Pressure = readtable("Point_Pressure.txt");
% OriginPoint_Pressure_Array = OriginPoint_Pressure{:, 1 : 1767};
% OriginPoint_Pressure = OriginPoint_Pressure(:, 1 : 1767);

% dlmwrite("Point_Pressure_1767.txt", OriginPoint_Pressure_Array, 'delimiter', '\t', 'precision', '%.6f');

OriginPoint_Coord.Properties.VariableNames = ["节点X坐标", "节点Y坐标", "节点Z坐标"];
OriginPoint_Coord_RowsToRemove = OriginPoint_Coord.("节点X坐标") == 0 & OriginPoint_Coord.("节点Y坐标") == 0 & OriginPoint_Coord.("节点Z坐标") == 0;
OriginPoint_Coord(OriginPoint_Coord_RowsToRemove, :) = [];
OriginPoint_Num = height(OriginPoint_Coord);
OriginPoint_X = OriginPoint_Coord{:, 1};
OriginPoint_Y = OriginPoint_Coord{:, 2};
OriginPoint_Z = OriginPoint_Coord{:, 3};
OriginPoint_Coord_Array = [OriginPoint_X, OriginPoint_Y, OriginPoint_Z];

% % 打开文件进行写入
% fileID = fopen('output.txt','w');
% 
% % 遍历数据并写入文件
% for i = 1:size(OriginPoint_Coord_Array,1)
%     fprintf(fileID, '(%f %f %f)\n', OriginPoint_Coord_Array(i,1), OriginPoint_Coord_Array(i,2), OriginPoint_Coord_Array(i,3));
% end
% 
% % 关闭文件
% fclose(fileID);    

OriginPoint_Pressure_ColsToRemove = find(OriginPoint_Coord_RowsToRemove) + 1;  % 计算列索引
OriginPoint_Pressure(:,OriginPoint_Pressure_ColsToRemove) = [];  % 删除列
OriginPoint_Pressure.Properties.VariableNames = ["Time", cellstr(strcat("Point_", string(1:OriginPoint_Num)))];
Time_Interval = OriginPoint_Pressure{2, 1} - OriginPoint_Pressure{1, 1};
Fs = 1 / Time_Interval;


[Data_Spectrum_Complex, Data_Spectrum_ABS] = TimeDataToFre(Fs, OriginPoint_Pressure);

% 介质参数
Rho_Medium = 1000;                                                  % 介质密度，默认为空气
SoundVelocity_Medium = 1500;                                         % 介质声速，默认为空气中的声速

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
ptCloud = pointCloud(OriginPoint_Coord_Array);                       % 读取点云数据
normals = pcnormals(ptCloud);                                  % 计算法向量（部分点为内法向，部分点为外法向）
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

%% ------------------------------【5 迭代计算均方根误差最小值】------------------------------
%{
    通过缩比的方式，获得等效简单源的坐标，等效简单源数目和椭球面离散节点数目一致
    以[缩放系数]Scale_Factor为自变量，迭代求解均方根误差，得到误差最小的[缩放系数]
%}

Error_ABS_Opt = Inf;
for Scale_Factor = 0.976 : 0.00001 : 0.976
    Equiv_Simple_Source_X = Scale_Factor .* OriginPoint_X;      % 【等效简单源】X坐标
    Equiv_Simple_Source_Y = Scale_Factor .* OriginPoint_Y;      % 【等效简单源】X坐标
    Equiv_Simple_Source_Z = Scale_Factor .* OriginPoint_Z;      % 【等效简单源】X坐标
    Equiv_Simple_Source_Coord = [Equiv_Simple_Source_X, Equiv_Simple_Source_Y, Equiv_Simple_Source_Z];            % 【等效简单源】X、Y、Z坐标

    Calculate_Fre = 100;
    Calculate_Row = 2;

    Calculate_OriginPoint_Complex = Data_Spectrum_Complex{Calculate_Row, 2 : end}.';
    Calculate_OriginPoint_ABS = Data_Spectrum_ABS{Calculate_Row, 2 : end}.';
    Calculate_Omega = 2 * pi * Calculate_Fre;
    Calculate_k = Calculate_Omega / SoundVelocity_Medium;

    EquivReverse_Q = Fun_ReverseBuild_Q(OriginPoint_Coord_Array, Equiv_Simple_Source_Coord, Calculate_k, SoundVelocity_Medium, Rho_Medium, Calculate_OriginPoint_Complex);
    
    [EquivReverse_D, EquivReverse_U_normal] = Fun_ReverseBuild_D_U(OriginPoint_Coord_Array, Equiv_Simple_Source_Coord, Calculate_k, EquivReverse_Q, normals);

    [Equiv_D_Radiation, Equiv_Q_Radiation] = Fun_BuildEquiv_D_Q(OriginPoint_Coord_Array, Equiv_Simple_Source_Coord, Calculate_k, EquivReverse_U_normal, normals);

    Equiv_P_Radiation = Fun_BuildEquiv_M_P(OriginPoint_Coord_Array, Equiv_Simple_Source_Coord, Calculate_k, SoundVelocity_Medium, Rho_Medium, Equiv_Q_Radiation);
    Equiv_P_Radiation_ABS = abs(Equiv_P_Radiation);
    
    Error_ABS_Array = Equiv_P_Radiation_ABS - Calculate_OriginPoint_ABS;

    if max(Error_ABS_Array) < 0
        Error_ABS_Curr = min(Error_ABS_Array);
        disp("全为负数");
        
    elseif min(Error_ABS_Array) > 0
        disp("全为正数");
        Error_ABS_Curr = max(Error_ABS_Array);
    else
        disp("正负数均有");
        Error_ABS_Curr = max(abs(Error_ABS_Array));
    end
    disp(Scale_Factor)
    disp(Error_ABS_Curr)
%     if Error_ABS_Curr > 0 && Error_ABS_Curr < Error_ABS_Opt
%         Error_ABS_Opt = Error_ABS_Curr;
%         Scale_Factor_Opt = Scale_Factor;
%     end
end

% disp("最小绝对值误差：");
% disp(Error_ABS_Opt);
% 
% disp("最佳缩放系数为：");
% disp(Scale_Factor_Opt);

% Scale_Factor_Opt = 0.999;
% Equiv_Simple_Source_X = Scale_Factor_Opt .* OriginPoint_X;      % 【等效简单源】X坐标
% Equiv_Simple_Source_Y = Scale_Factor_Opt .* OriginPoint_Y;      % 【等效简单源】X坐标
% Equiv_Simple_Source_Z = Scale_Factor_Opt .* OriginPoint_Z;      % 【等效简单源】X坐标
% Equiv_Simple_Source_Coord = [Equiv_Simple_Source_X, Equiv_Simple_Source_Y, Equiv_Simple_Source_Z];                      % 【等效简单源】X、Y、Z坐标
% 
% Calculate_Fre = 100;
% Calculate_Row = 2;
% 
% Calculate_OriginPoint_Complex = Data_Spectrum_Complex{Calculate_Row, 2 : end}.';
% Calculate_OriginPoint_ABS = Data_Spectrum_ABS{Calculate_Row, 2 : end}.';
% Calculate_Omega = 2 * pi * Calculate_Fre;
% Calculate_k = Calculate_Omega / SoundVelocity_Medium; 
% 
% EquivReverse_Q = Fun_ReverseBuild_Q(OriginPoint_Coord_Array, Equiv_Simple_Source_Coord, Calculate_k, SoundVelocity_Medium, Rho_Medium, Calculate_OriginPoint_Complex);
% 
% [EquivReverse_D, EquivReverse_U_normal] = Fun_ReverseBuild_D_U(OriginPoint_Coord_Array, Equiv_Simple_Source_Coord, Calculate_k, EquivReverse_Q, normals);

%% ------------------------------【7 正向求解结构辐射声压】------------------------------
%{
    基于结构表面法向振速，求解结构辐射声压
%}
% 辐射声压检测场点的位置和数量
% FieldPoints_Num = 50;                                               % 自定义场点数量
% FieldPoints = zeros(FieldPoints_Num, 3);
% FieldPoints(:, 2) = linspace(10 * (max(OriginPoint_Y)), 20 * (max(OriginPoint_Y)), FieldPoints_Num)';          % 建立场点距离向量 X 轴坐标，范围为10倍至20倍脉动球源半径，转置为列向量
% FieldPoints(:, 1) = centroid(1);
% FieldPoints(:, 3) = centroid(3);
% 
% [Equiv_D_Radiation, Equiv_Q_Radiation] = Fun_BuildEquiv_D_Q(OriginPoint_Coord_Array, Equiv_Simple_Source_Coord, Calculate_k, EquivReverse_U_normal, normals);
% Equiv_P_Radiation = Fun_BuildEquiv_M_P(OriginPoint_Coord_Array, Equiv_Simple_Source_Coord, Calculate_k, SoundVelocity_Medium, Rho_Medium, Equiv_Q_Radiation);
% 
% Equiv_P_Radiation_ABS = abs(Equiv_P_Radiation);
% Error_RMS_Calculate = sqrt(mean((Equiv_P_Radiation_ABS - Calculate_OriginPoint_ABS).^2));
% 
% Equiv_P_Radiation_FieldPoints = Fun_BuildEquiv_M_P(FieldPoints, Equiv_Simple_Source_Coord, Calculate_k, SoundVelocity_Medium, Rho_Medium, Equiv_Q_Radiation);
% 
% figure;
% plot(FieldPoints(:, 2), abs(Equiv_P_Radiation_FieldPoints), 'LineWidth', 1.5);
% xlabel('测点(Y轴)');
% ylabel('声压幅值(Pa)');
% title('等效简单源辐射声压实部');