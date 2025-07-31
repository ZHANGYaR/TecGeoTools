function head = Search_WaterHead(FlowAcc,FlowDir,LeftX,DownY,cellsize,ouletX,ouletY,Ac,str1)
% 使用方法:
% head = Search_WaterHead(FlowAcc, FlowDir, LeftX, DownY, cellsize, ouletX, ouletY, Ac, str1)
%
% 描述:
% 该函数通过逆向追踪法自动识别满足汇水面积阈值的河道源头（水头），基于D8算法反向搜索，
% 结合地形参数生成水头点坐标文件。适用于河网自动提取及地形演化研究。
%
% 必要输入:
% FlowAcc - 汇水面积矩阵(m×n)，数值为各栅格点汇水面积栅格数（需预先乘以cellsize^2转换为真实面积）
% FlowDir - D8流向矩阵(m×n)，取值范围{1,2,4,8,16,32,64,128}，遵循ArcGIS流向标准
% LeftX - 栅格左下角X坐标（地图坐标系）
% DownY - 栅格左下角Y坐标（地图坐标系）
% cellsize - 栅格单元物理尺寸（米）
% ouletX - 目标河段出口点X坐标数组（地图坐标系）
% ouletY - 目标河段出口点Y坐标数组（地图坐标系）
% Ac - 汇水面积阈值(m²)，用于识别起始水头点的临界值
% str1 - 输出文件路径前缀，示例：'D:/data/output_'
%
% 输出结果:
% head - 水头点坐标矩阵(N×2)，格式为：[X坐标, Y坐标]
% 生成文件 [str1]_Head.txt - ASCII文本文件，每行存储水头点坐标（X,Y浮点格式）
%
% 算法流程:
% 1. 将输入坐标转换为矩阵行列号
% 2. 初始化栈结构进行逆向追踪：
% - 从每个出口点出发，沿反向D8流向搜索上游支流
% - 动态记录各节点的汇入支流数和有效支流数
% 3. 水头判定条件：
% - 当某点所有汇入支路的汇水面积 ≤ Ac 时，标记为水头
% 4. 坐标转换：
% - 将行列号转换为地理坐标（带栅格中心点校正）
%
% 示例调用:
% % 基础参数设置
% cellsize = 30; % DEM分辨率30米
% oulet = [120.5, 25.3]; % 关注河口坐标
% Ac = 1e6; % 汇水面积阈值1km²
% [heads] = Search_WaterHeader(FlowAcc, FlowDir, 120000, 2500000,...
% cellsize, oulet(1), oulet(2), Ac, 'output/watershed1');
%
% 关键参数建议:
% 1. Ac选择：
% - 研究大流域：1e6 ~ 1e7 m²（1~10 km²）
% - 小流域分析：1e4 ~ 1e5 m²（0.01~0.1 km²）
% 2. 流向处理：
% - 确保FlowDir中背景区域值为0，已自动过滤
% - 输入前应通过填洼预处理DEM
%
% 注意事项:
% 1. 地理坐标系统：
% - 要求Y轴北向递增（适合UTM坐标系）
% - 坐标转换公式：Y = UpY - row×cellsize
% 2. 输入矩阵要求：
% - FlowAcc需为整型栅格数，先转换为真实面积累加结果
% - 流向矩阵需严格遵循D8编码规范
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数编写者: 王一舟 - 更新日期: 2025年2月4日
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long;
[m,n]=size(FlowAcc);
UpY=DownY+cellsize*m;
Ac=floor(Ac/(cellsize^2)); % Ac位汇水面积临界值
FlowDir = double(FlowDir);
len_ouletX=length(ouletX);
stepPoint=[]; WaterHead=[];

num_ToPoint=0; % 汇入该点的栅格数
num_slope=0; % 汇入该点的栅格中，<Ac的栅格数目
% 当二者相等，说明该点为水头

for ilen=1:len_ouletX
    % 起始搜索点的行、列号
    startRow=floor((UpY-ouletY(ilen))/cellsize)+1; 
    startCol=floor((ouletX(ilen)-LeftX)/cellsize)+1; 

    stepPoint=[stepPoint;[startRow,startCol,FlowAcc(startRow,startCol),0,0]];
    empty_stepPoint=1;
    % 搜索点行、列号；该点的汇水面积；该点汇入区域的行、列号（0,0表示无汇入）
    while empty_stepPoint
        startRow=stepPoint(1,1);
        startCol=stepPoint(1,2);
        
        if FlowAcc(startRow,startCol)>=Ac
            % 确保该点是河流点
        if startCol<n && FlowDir(startRow,startCol+1)==2^4 % && FlowAcc(startRow,startCol)>Ac % && FlowAcc(startRow,startCol+1)>Ac
            num_ToPoint=num_ToPoint+1; % 汇入该点
            % 该点右侧(>Ac)，入栈
            if FlowAcc(startRow,startCol+1)>Ac; stepPoint=[stepPoint;[startRow,startCol+1,FlowAcc(startRow,startCol+1),startRow,startCol]]; end
            % 该点右侧(<Ac)，该点有作为水头的可能
            if FlowAcc(startRow,startCol+1)<=Ac; num_slope=num_slope+1; end
        end  
        if startRow<m && startCol<n && FlowDir(startRow+1,startCol+1)==2^5 % && FlowAcc(startRow,startCol)>Ac % && FlowAcc(startRow+1,startCol+1)>Ac
             num_ToPoint=num_ToPoint+1; % 汇入该点
             if FlowAcc(startRow+1,startCol+1)>Ac;stepPoint=[stepPoint;[startRow+1,startCol+1,FlowAcc(startRow+1,startCol+1),startRow,startCol]];end
             if FlowAcc(startRow+1,startCol+1)<=Ac;num_slope=num_slope+1;end
        end
        if startRow<m && FlowDir(startRow+1,startCol)==2^6 % && FlowAcc(startRow,startCol)>Ac % && FlowAcc(startRow+1,startCol)>Ac
            num_ToPoint=num_ToPoint+1; % 汇入该点
            if FlowAcc(startRow+1,startCol)>Ac;stepPoint=[stepPoint;[startRow+1,startCol,FlowAcc(startRow+1,startCol),startRow,startCol]];end
            if FlowAcc(startRow+1,startCol)<=Ac;num_slope=num_slope+1;end
        end
        if startRow<m && startCol>1 && FlowDir(startRow+1,startCol-1)==2^7 % && FlowAcc(startRow,startCol)>Ac  % && FlowAcc(startRow+1,startCol-1)>Ac
            num_ToPoint=num_ToPoint+1; % 汇入该点
            if FlowAcc(startRow+1,startCol-1)>Ac;stepPoint=[stepPoint;[startRow+1,startCol-1,FlowAcc(startRow+1,startCol-1),startRow,startCol]];end
            if FlowAcc(startRow+1,startCol-1)<=Ac;num_slope=num_slope+1;end
        end
        if startCol>1 && FlowDir(startRow,startCol-1)==2^0 % && FlowAcc(startRow,startCol)>Ac % && FlowAcc(startRow,startCol-1)>Ac
            num_ToPoint=num_ToPoint+1; % 汇入该点
            if FlowAcc(startRow,startCol-1)>Ac;stepPoint=[stepPoint;[startRow,startCol-1,FlowAcc(startRow,startCol-1),startRow,startCol]];end
            if FlowAcc(startRow,startCol-1)<=Ac;num_slope=num_slope+1;end
        end
        if startRow>1 && startCol>1 && FlowDir(startRow-1,startCol-1)==2^1 % && FlowAcc(startRow,startCol)>Ac % && FlowAcc(startRow-1,startCol-1)>Ac
            num_ToPoint=num_ToPoint+1; % 汇入该点
            if FlowAcc(startRow-1,startCol-1)>Ac;stepPoint=[stepPoint;[startRow-1,startCol-1,FlowAcc(startRow-1,startCol-1),startRow,startCol]];end
            if FlowAcc(startRow-1,startCol-1)<=Ac;num_slope=num_slope+1;end
        end
        if startRow>1 && FlowDir(startRow-1,startCol)==2^2 % && FlowAcc(startRow,startCol)>Ac % && FlowAcc(startRow-1,startCol)>Ac
            num_ToPoint=num_ToPoint+1; % 汇入该点
            if FlowAcc(startRow-1,startCol)>Ac;stepPoint=[stepPoint;[startRow-1,startCol,FlowAcc(startRow-1,startCol),startRow,startCol]];end
            if FlowAcc(startRow-1,startCol)<=Ac;num_slope=num_slope+1;end
        end
        if startRow>1 && startCol<n && FlowDir(startRow-1,startCol+1)==2^3 % && FlowAcc(startRow,startCol)>Ac  % && FlowAcc(startRow-1,startCol+1)>Ac
            num_ToPoint=num_ToPoint+1; % 汇入该点
            if FlowAcc(startRow-1,startCol+1)>Ac;stepPoint=[stepPoint;[startRow-1,startCol+1,FlowAcc(startRow-1,startCol+1),startRow,startCol]];end
            if FlowAcc(startRow-1,startCol+1)<=Ac;num_slope=num_slope+1;end
        end
        % 汇入点数目 与 汇入点<Ac数目 相等，则该点为水头
        if num_slope==num_ToPoint; WaterHead=[WaterHead;[startRow,startCol]]; end
        end
        num_slope=0; num_ToPoint=0;
        stepPoint(1,:)=[]; % % 每一轮搜索完成后，原来的起始点 出栈, 第一个点为上一步的搜索点
        empty_stepPoint=~isempty(stepPoint);
    end
end
WaterHead(:,1)=UpY-WaterHead(:,1).*cellsize+cellsize/2; % row，Y坐标
WaterHead(:,2)=LeftX+WaterHead(:,2).*cellsize-cellsize/2; % Col，X坐标
length_WaterHead = length(WaterHead);

% 文件路径
str2 = '_Head';    str3 = '.txt';     
Pathi_Out = strcat(str1, str2, str3);
% 确保文件夹存在
folderPath = fileparts(Pathi_Out);  % 获取文件夹路径
if exist(folderPath, 'dir') == 0
    mkdir(folderPath);  % 如果文件夹不存在，则创建
end

% 打开文件
fid = fopen(Pathi_Out, 'w');  
if fid == -1
    error('无法打开文件: %s', Pathi_Out);  % 如果文件打开失败，输出错误信息
end
for j = 1:length_WaterHead
    fprintf(fid, '%f\t', WaterHead(j, 2));  % 写入 X 坐标
    fprintf(fid, '%f\n', WaterHead(j, 1));  % 写入 Y 坐标
end
fclose(fid);
head = WaterHead(:, [2, 1]);
end























