
function [Data_Spectrum_Complex, Data_Spectrum_ABS] = TimeDataToFre(Fs, Point_Pressure)
OriginPoint_Num = width(Point_Pressure) - 1;
Point_Pressure_CutMean = Point_Pressure;
Mean_values = mean(Point_Pressure{:, 2:end});

for i = 1:length(Mean_values)                                       % 使用每一列的均值来修正数据
    Point_Pressure_CutMean{:, 1+i} = Point_Pressure{:, 1+i} - Mean_values(i);
end

% 初始化一个新的表，用来存储频谱数据
% 列数保持和Point_Pressure_CutMean一样，行数为(Point_Pressure_CutMean行数除以2)
Fre_ColsNum = OriginPoint_Num + 1;
Fre_RowsNum = floor(height(Point_Pressure_CutMean) / 2);
Data_Spectrum_ABS = array2table(nan(Fre_RowsNum, Fre_ColsNum));
Data_Spectrum_Complex = array2table(nan(Fre_RowsNum, Fre_ColsNum));

% 设定列名
Data_Spectrum_ABS.Properties.VariableNames = ["频率", cellstr(strcat("幅值Point_", string(1:OriginPoint_Num)))];
Data_Spectrum_Complex.Properties.VariableNames = ["频率", cellstr(strcat("复值Point_", string(1:OriginPoint_Num)))];

for i = 2 : Fre_ColsNum                                                  % 遍历所有数据行
% for i = 3 : 3
    spectrum_complex = fft(Point_Pressure_CutMean{:, i});                % 计算傅里叶变换
    spectrum_complex = spectrum_complex(1:Fre_RowsNum);            % 由于频谱是对称的，只需要前一半的数据
    Data_Spectrum_Complex{:, i} = spectrum_complex;
    
    spectrum_abs = abs(spectrum_complex);
    Data_Spectrum_ABS{:, i} = spectrum_abs;
end

f = Fs * (0:Fre_RowsNum-1) / (2 * Fre_RowsNum);                        % 计算频率值得到一个包含频率值的数组
f_ColVector = f.';

Data_Spectrum_ABS{:, 1} = f_ColVector;
Data_Spectrum_Complex{:, 1} = f_ColVector;
end