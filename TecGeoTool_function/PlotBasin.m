function PlotBasin(sbd_dir,sbd_fld)
%
% 用法：
%   绘制指定流域的细分盆地结果
%    sbd_dir：大流域文件夹位置 
%    sbd_fld: 细分结果文件夹名称
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数作者：Yarong Zhang - 更新日期：2024年12月18日 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off

 % 弹出对话框，输入流域编号
    prompt = {'请输入流域编号 (例如: 1, 2, 3):'};
    dlgTitle = '流域编号输入';
    numLines = 1;  defaultAnswer = {'1'};  % 默认值
    answer = inputdlg(prompt, dlgTitle, numLines, defaultAnswer);
    
    if isempty(answer)  % 如果用户取消对话框
        disp('用户取消了输入');
        return;
    end
    
    % 获取输入的流域编号
    basinNumber = str2double(answer{1});  % 将输入的流域编号转换为数字
    
    if isnan(basinNumber)  % 检查输入是否有效
        disp('请输入有效的数字');
        return;
    end

 asciiFileName = sprintf('Basin_%d_DEM.txt', basinNumber);  % 生成文件名
    filePath = fullfile(sbd_dir, asciiFileName);  % 完整的文件路径
    
    % 检查文件是否存在
    if ~isfile(filePath)
        disp('指定的 DEM 文件不存在');
        return;
    end
    
 % 设置文件夹路径
asciiFileName = sprintf('Basin_%d_DEM.txt', basinNumber);  % 生成文件名
filePath = fullfile(sbd_dir, asciiFileName);  % 完整的文件路径

% 1. 读取ASCII文件
fid = fopen(filePath, 'r');
headerInfo = textscan(fid, '%s %f', 6); % 读取前6行的头部信息
fclose(fid);

% 提取头部信息
ncols = headerInfo{2}(1);       % 栅格列数
nrows = headerInfo{2}(2);       % 栅格行数
xllcorner = headerInfo{2}(3);   % 左下角X坐标
yllcorner = headerInfo{2}(4);   % 左下角Y坐标
cellsize = headerInfo{2}(5);    % 栅格单元大小
nodata_value = headerInfo{2}(6);% 无数据值

% 3. 读取DEM数据
demData = dlmread(filePath, ' ', 6, 0); % 从第7行开始读取DEM数据
demData(demData == nodata_value) = NaN; % 将无数据值替换为 NaN

% 5. 生成地理坐标
xCoords = xllcorner + (0:ncols-1) * cellsize; % X坐标
yCoords = yllcorner + (0:nrows-1) * cellsize; % Y坐标
DEM = GRIDobj(xCoords, yCoords, demData);

figure;
[RGB] = imageschs(DEM, DEM, 'colormap', 'gray');
[~, R] = GRIDobj2im(DEM);
imshow(RGB, R); hold on;
%imageschs(DEM,DEM,'colormap','gray');hold on;

%%
MatFileName = sprintf('Basin_%d_Data.mat', basinNumber);
filePath = fullfile(sbd_dir, MatFileName);  
fid = load(filePath);
SL=fid.SLc;
plot(SL,'-w'); hold on;
%%
folderPath= fullfile(sbd_dir, sbd_fld);
matFiles = dir(fullfile(folderPath, '*.mat'));  % 获取所有.mat文件
% 1. 初始化变量存储 KSN_stats 和颜色
KSN_stats_all = zeros(length(matFiles), 1);  % 存储所有的 KSN_stats 值

% 2. 遍历所有 .mat 文件，提取 KSN_stats
for i = 1:length(matFiles)
    matFile = fullfile(matFiles(i).folder, matFiles(i).name);
    data = load(matFile);  % 加载.mat文件
    
    % 提取 KSN_stats 值
    KSN_stats = data.KSNc_stats(1);
    KSN_stats_all(i) = KSN_stats;
end
colorMap = ksncolor(length(matFiles));  % 获取颜色映射
assignin('base', 'colorMap', colorMap);  % 将变量赋值到基础工作区
% 3. 按 KSN_stats 排序
[sortedValues] = sort(KSN_stats_all, 'descend');
assignin('base', 'sortedValues', sortedValues);  % 将变量赋值到基础工作区

% 4. 遍历并绘制 DEM 区域
for j = 1:length(sortedValues)
    value = KSN_stats_all(j);
    idx = find(sortedValues == value);
    % 构建对应的ASCII文件路径
    asciiFileName = sprintf('Basin_%d_DataSubset_%d_DEM.txt',basinNumber,j);
    asciiFilePath = fullfile(folderPath, asciiFileName);
    
    % 读取ASCII文件
    fid = fopen(asciiFilePath, 'r');
    headerInfo = textscan(fid, '%s %f', 6); % 读取头部信息
    fclose(fid);
    
    % 提取头部信息
    ncols = headerInfo{2}(1);
    nrows = headerInfo{2}(2);
    xllcorner = headerInfo{2}(3);
    yllcorner = headerInfo{2}(4);
    cellsize = headerInfo{2}(5);
    nodata_value = headerInfo{2}(6);
    
    % 读取DEM数据
    demData = dlmread(asciiFilePath, ' ', 6, 0);  % 跳过前6行头部信息
    demData(demData == nodata_value) = NaN;       % 将nodata值替换为NaN
    
    % 生成网格坐标
    xCoords = xllcorner + (0:ncols-1) * cellsize; % X坐标
    yCoords = yllcorner + (0:nrows-1) * cellsize; % Y坐标
    yCoords = flip(yCoords);  % 翻转 Y 坐标

    % 提取有效区域边界
    validMask = ~isnan(demData);                 % 有效区域掩膜
    boundaries = bwboundaries(validMask);        % 提取有效区域边界

    % 绘制DEM边界
    for k = 1:length(boundaries)
        boundary = boundaries{k};
        rows = boundary(:, 1);  % 行索引
        cols = boundary(:, 2);  % 列索引

        % 行列索引转换为地理坐标
        x = xCoords(cols);     y = yCoords(rows);

        % 创建 polyshape 并绘制
        %y = max(yCoords) - (y - min(yCoords));
        polygon = polyshape(x, y, 'Simplify', true);
        if ~isempty(polygon.Vertices)
            hold on
            plot(polygon, 'FaceColor', colorMap(idx,:), 'EdgeColor', 'none', 'FaceAlpha', 0.8);        
        end
    end
end

% 设置颜色映射和色条
set(gca, 'YDir', 'normal'); % 保证 Y 轴从下到上
colormap(ksncolor(20));
caxis([min(KSN_stats_all), max(KSN_stats_all)]);  % 设置颜色条范围
c1 = colorbar;
ylabel(c1, 'Normalized Channel Steepness');
hold off;
%%
while true
    % 弹出对话框，询问是否继续
    choice = questdlg('操作完成，是否继续处理其他流域？', ...
                      '继续操作', ...
                      '继续', '退出', '继续'); % 三个选项，默认"继续"

    % 根据用户的选择执行操作
    switch choice
        case '继续'
            disp('继续新的流域操作...');
            PlotBasin(sbd_dir, sbd_fld); % 调用自身函数，重复操作
            return; % 结束当前操作并重新开始
        case '退出'
            disp('操作已结束。');
            break; % 退出 while 循环，程序结束
    end
end
