function RiverProInfo(Elev, FlowDir, FlowAcc, cellsize, headX, headY, LeftX, DownY, eleWin, mVSn, str1)
% 使用方法:
% RiverProInfo(Elev, FlowDir, FlowAcc, cellsize, headX, headY, LeftX, DownY, eleWin, mVSn, str1)
%
% 描述:
% 本函数用于从DEM数据中提取河流纵剖面信息，支持多起点追踪。通过流向矩阵进行河道溯源分析，
% 计算各节点的溯源距离、Chi值等形态参数，输出包含完整纵剖面信息的文本文件。
%
% 必要输入:
% Elev - 高程矩阵(MxN)，数据类型为double，负值表示无效区域
% FlowDir - D8流向矩阵(MxN)，取值范围{1,2,4,8,16,32,64,128}，0表示无效值
% FlowAcc - 汇水面积矩阵(MxN)，单位为栅格数，需为整数
% cellsize - 栅格单元物理尺寸（米），正实数
% headX - 河流起点X坐标数组（地图坐标系，与LeftX同基准）
% headY - 河流起点Y坐标数组（地图坐标系，与DownY同基准）
% LeftX - 栅格左下角X坐标（地图坐标系）
% DownY - 栅格左下角Y坐标（地图坐标系）
% eleWin - 高程采样窗口（米），控制纵剖面点间距的垂直阈值
% mVSn - 河床侵蚀模型m/n比率值，用于计算Chi值
% str1 - 输出文件前缀字符串，建议包含路径信息
%
% 输出文件:
% [str1]_RiverPro[N].txt - 纵剖面数据文件，每行包含：
% Y坐标 X坐标 溯源距离(m) 高程(m) 流向 汇水面积(m²) Chi值
%
% 关键处理流程:
% 1. 坐标转换：将地理坐标转换为矩阵行列号
% 2. 河道溯源：根据D8流向矩阵进行河道追踪
% 3. 数据采样：按高程窗口eleWin控制剖面点密度
% 4. Chi计算：基于积分公式 Chi = ∫(A0/A)^(m/n) dx
%
% 注意事项:
% 1. 输入矩阵需满足：
% - Elev与FlowDir/FlowAcc同尺寸
% - FlowDir严格遵循D8编码规范
% - FlowAcc应为累积栅格数
% 2. 坐标系要求：
% - headX/headY应与LeftX/DownY使用相同地图投影
% - Y轴需满足北向递增（通常为UTM或高斯投影）
% 3. 参数建议：
% - eleWin一般设置为0.000001
% - m/n值参考范围0.3-0.6（构造活跃区取较高值）
% 4. 异常处理：
% - 自动跳过高程为负的栅格
% - 流向为0或非法值时终止追踪
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数编写者: 王一舟 - 更新日期: 2025年2月4日
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


A0 = 1; % 无维度化汇水面积1m�0�5
[m, n] = size(Elev); % 获取高程矩阵尺寸
lenheadX = length(headX); % 起点数量
upY = DownY + m * cellsize; % 栅格的顶部 Y 坐标


for j = 1:lenheadX
    % 初始化起点行列号
    midRow = floor((abs(upY - headY(j))) / cellsize) + 1;
    midCol = floor((abs(headX(j) - LeftX)) / cellsize) + 1;
    RiverLen = 0; 
    RiverPro1 = []; 
    kk = 1;

    % 首先填充顺流矩阵 RiverPro1
    RiverPro1 = [RiverPro1; midRow, midCol, RiverLen, Elev(midRow, midCol), FlowDir(midRow, midCol), FlowAcc(midRow, midCol)];

    % 检查背景值并替换为 NaN
FlowDir(FlowDir == 0) = NaN;

% 初始化
RiverLen = 0;
visited = false(size(FlowDir)); % 记录访问过的坐标
validDirections = [1, 2, 4, 8, 16, 32, 64, 128];

while midRow > 1 && midRow < size(Elev, 1) && midCol > 1 && midCol < size(Elev, 2)
    % 检查是否越界或重复访问
    if midRow < 1 || midRow > size(Elev, 1) || midCol < 1 || midCol > size(Elev, 2) || visited(midRow, midCol)
        break;
    end
    visited(midRow, midCol) = true;

    % 检查 Elev 和 FlowDir 是否有效
    if Elev(midRow, midCol) < 0 || isnan(FlowDir(midRow, midCol)) || ~ismember(FlowDir(midRow, midCol), validDirections)
        break;
    end

    % 根据流向更新坐标和河流长度
    switch FlowDir(midRow, midCol)
        case 1 % 东
            midCol = midCol + 1; 
            RiverLen = RiverLen + 1;
        case 2 % 东南
            midRow = midRow + 1; 
            midCol = midCol + 1; 
            RiverLen = RiverLen + 1.414;
        case 4 % 南
            midRow = midRow + 1; 
            RiverLen = RiverLen + 1;
        case 8 % 西南
            midRow = midRow + 1; 
            midCol = midCol - 1; 
            RiverLen = RiverLen + 1.414;
        case 16 % 西
            midCol = midCol - 1; 
            RiverLen = RiverLen + 1;
        case 32 % 西北
            midRow = midRow - 1; 
            midCol = midCol - 1; 
            RiverLen = RiverLen + 1.414;
        case 64 % 北
            midRow = midRow - 1; 
            RiverLen = RiverLen + 1;
        case 128 % 东北
            midRow = midRow - 1; 
            midCol = midCol + 1; 
            RiverLen = RiverLen + 1.414;
        otherwise
            break; % 异常流向值
    end



        % 添加数据点
        if Elev(midRow, midCol) < 0
            break; % 高程为负，退出
        end
        if (RiverPro1(kk, 4) - Elev(midRow, midCol) >= eleWin) % 超过高程窗口
            RiverPro1 = [RiverPro1; midRow, midCol, RiverLen, Elev(midRow, midCol), FlowDir(midRow, midCol), FlowAcc(midRow, midCol)];
            kk = kk + 1;
        end
    end

    % 初始化溯流矩阵
    jRiverPro = zeros(kk, 7); 
    minEle = RiverPro1(kk, 4); 
    MaxRiverLen = RiverPro1(kk, 3);
    for i = 1:kk
        jRiverPro(i, 1) = RiverPro1(kk - i + 1, 1); % 行号
        jRiverPro(i, 2) = RiverPro1(kk - i + 1, 2); % 列号
        jRiverPro(i, 3) = MaxRiverLen - RiverPro1(kk - i + 1, 3); % 溯源距离
        jRiverPro(i, 4) = RiverPro1(kk - i + 1, 4); % 高程
        jRiverPro(i, 5) = RiverPro1(kk - i + 1, 5); % 流向
        jRiverPro(i, 6) = RiverPro1(kk - i + 1, 6); % 汇水面积
    end

    % 转换为地理坐标
    jRiverPro(:, 1) = upY - cellsize .* jRiverPro(:, 1); % Y 坐标
    jRiverPro(:, 2) = LeftX + cellsize .* jRiverPro(:, 2); % X 坐标
    jRiverPro(:, 3) = jRiverPro(:, 3) .* cellsize; % 溯源距离
    jRiverPro(:, 6) = jRiverPro(:, 6) .* (cellsize^2); % 汇水面积

    % 计算 Chi 距离
    for i = 2:kk
        jRiverPro(i, 7) = jRiverPro(i - 1, 7) + ((A0 ./ jRiverPro(i - 1, 6)).^mVSn) .* (jRiverPro(i, 3) - jRiverPro(i - 1, 3));
    end

    % 保存到文件
    str2 = '_RiverPro'; 
    str3 = num2str(j);
    str4 = '.txt';
    str_j = strcat(str1, str2, str3, str4);
    fid = fopen(str_j, 'w');
    for i = 1:kk
        fprintf(fid, '%f %f %f %f %f %f %f\n', jRiverPro(i, :));
    end
    fclose(fid);
end
end
