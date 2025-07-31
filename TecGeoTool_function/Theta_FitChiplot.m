function Theta_FitChiplot(Path1,startRiver,endRiver,str1)
% 使用方法:
% Theta_FitChiplot(Path1, startRiver, endRiver, str1)
%
% 描述:
% 本函数基于Mudd et al (2018)方法，通过分析多河流纵剖面的Chi-高程关系，
% 计算使高程散点最小化的最优河道凹度θ值。支持批量处理河流数据，生成θ-MisFit
% 关系曲线及不同θ值的Chi-高程对比图。
%
% 必要输入:
% Path1 - 河流纵剖面数据存储路径，字符串格式。要求路径下包含[N].txt格式文件
% startRiver - 起始河流编号，正整数，文件命名起始序号
% endRiver - 终止河流编号，正整数，文件命名结束序号
% str1 - 输出文件前缀，字符串，建议包含完整存储路径
%
% 输出文件:
% [str1]_theta-min_Scatter_rather.txt - 最优θ值记录文件，包含：
% 最小MisFit值(m) 最佳θ值
%
% 图形输出:
% 1) θ-MisFit关系曲线图
% 2) 最佳θ值对应的Chi-高程分布图
% 3) θ-0.2时的Chi-高程对比图
% 4) θ+0.2时的Chi-高程对比图
%
% 关键处理流程:
% 1. 数据读取：批量加载指定编号的河流剖面文件，文件需包含：
% 上游距离(m) | 高程(m) | 汇水面积(m²) 等关键参数
% 2. Chi计算：根据公式 χ = ∫(A0/A)^θ dx 生成不同θ值的Chi序列
% 3. 主干识别：选择Chi最长的河流作为基准剖面
% 4. 误差计算：通过线性插值比较各支流与主干的高程偏差，计算均方根误差(MisFit)
% 5. 参数优化：寻找使总MisFit最小的θ值作为最优河道凹度
%
% 注意事项:
% 1. 输入文件要求：
% - 文本格式，每列包含：X坐标 Y坐标 上游距离 高程 流向 汇水面积
% - 文件命名需连续编号，如1.txt, 2.txt,...n.txt
% 2. 参数设置：
% - θ搜索范围固定为0.1-1.0，步长0.05，需修改需调整代码
% - 高程标准化采用最小值归零处理
% 3. 图形特性：
% - 子图1显示θ-MisFit曲线，自动标定y轴范围
% - 子图2-4使用灰色绘制各河流Chi-高程曲线
% 4. 异常处理：
% - 自动跳过不存在或格式错误的文件
% - θ超出范围时自动取边界值
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数编写者: 王一舟 - 更新日期: 2025年2月4日
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 设置
format long; close all;
theta = 0.1:0.05:1;
len_theta = length(theta);

%% 读取所有河流，并计算不同θ对应的chi-z
Path3 = '.txt';
sumRiver = endRiver - startRiver + 1; %河流个数
AllRiverPro = cell(sumRiver,1); % 利用cell存储所有的河流纵剖面矩阵
lengthChi = [];

% 读取每条河流的数据
for i = startRiver:endRiver
    Path2 = num2str(i);
    Path_i = strcat(Path1, Path2, Path3);
    RiverPro_i = textread(Path_i); % 读取数据
    Ele_i = RiverPro_i(:,4); 
    Area_i = RiverPro_i(:,6); 
    upLen_i = RiverPro_i(:,3);

    chiLen_i = cal_chiLen_i(Area_i, upLen_i, theta); % 计算chiLen矩阵
    
    minEle_i = min(Ele_i); 
    Ele_i = Ele_i - minEle_i; % 标准化高度
    AllRiverPro{i - startRiver + 1, 1} = [chiLen_i, Ele_i]; % 存储所有河流数据
end

%% 计算每个theta的MisFit值，并绘制图形
MisFit_i = zeros(1, len_theta);
figure;

for i = 1:len_theta
    RiverPro_Chosen = cell(sumRiver, 1);
    for j = 1:sumRiver
        ChiZ_ij = AllRiverPro{j, 1};
        RiverPro_Chosen{j, 1} = [ChiZ_ij(:, i), ChiZ_ij(:, len_theta + 1)]; % 提取当前θ下的Chi-z
    end
    MisFit_i(i) = co_stem(RiverPro_Chosen); % 计算MisFit
end

% 绘制θ与MisFit的关系图
subplot(2, 2, 1);
plot(theta, MisFit_i); 
xlabel('theta');
ylabel('Elevation scatter (m)');
[min_MisFit_i, num_MisFit_i] = min(MisFit_i);
min_theta = theta(num_MisFit_i);
title(strcat('Min MisFit: ', num2str(min_MisFit_i), ', theta = ', num2str(min_theta)));
max_MisFit_i = max(MisFit_i);
ylim_down = floor(min_MisFit_i / 10) * 10;
ylim_up = ceil(max_MisFit_i / 10) * 10;
ylim([ylim_down, ylim_up]);

%% 绘制最佳凹度下所有的chi与高程图
subplot(2, 2, 2);
for j = 1:sumRiver
    RiverPro_ij = AllRiverPro{j, 1};
    plot(RiverPro_ij(:, num_MisFit_i), RiverPro_ij(:, len_theta + 1),'Color', [0.5 0.5 0.5]);
    hold on;
end
xlabel('Chi Length');
ylabel('Elevation');
title(strcat('Chi-Z for best theta  = ', num2str(min_theta)));

%% 绘制最佳凹度 - 0.2 时的chi与高程图
theta_minus_02 = min_theta - 0.2;
subplot(2, 2, 3);
for j = 1:sumRiver
    RiverPro_ij = AllRiverPro{j, 1};
    [~, idx] = min(abs(theta - theta_minus_02)); % 找到最接近min_theta-0.2的theta索引
    plot(RiverPro_ij(:, idx), RiverPro_ij(:, len_theta + 1), 'Color', [0.5 0.5 0.5]);
    hold on;
end
xlabel('Chi Length');
ylabel('Elevation');
title(strcat('Chi-Z for theta  = ', num2str(theta_minus_02)));

%% 绘制最佳凹度 + 0.2 时的chi与高程图
theta_plus_02 = min_theta + 0.2;
subplot(2, 2, 4);
for j = 1:sumRiver
    RiverPro_ij = AllRiverPro{j, 1};
    [~, idx] = min(abs(theta - theta_plus_02)); % 找到最接近min_theta+0.2的theta索引
    plot(RiverPro_ij(:, idx), RiverPro_ij(:, len_theta + 1), 'Color', [0.5 0.5 0.5]);
    hold on;
end
xlabel('Chi Length');
ylabel('Elevation');
title(strcat('Chi-Z for theta  = ', num2str(theta_plus_02)+0.2));

%% 保存最小误差值
str2 = 'theta-min_Scatter_rather.txt';
str = strcat(str1, str2);
fid = fopen(str, 'a+');
fprintf(fid, '%f\t', min_MisFit_i);       
fprintf(fid, '%f\n', min_theta);
fclose(fid);

end

%%  计算不同θ值对应的Chi
function chiLen_i=cal_chiLen_i(Area_new,upLen_new,concavity)

chiLen_i=[];

len_con=length(concavity);      mLen_new=length(Area_new);
for i=1:len_con
    Chi=[]; sum_i=0; Chi=[Chi;sum_i];
    for j=2:mLen_new
        sum_i=sum_i+(1/Area_new(j-1)).^(concavity(i))*(upLen_new(j)-upLen_new(j-1));
        Chi=[Chi;sum_i]; % Chi是列向量
    end
    chiLen_i=[chiLen_i,Chi];
end
end
%% 计算Chi-z的misFit,最接近主干道
function MisFit_i=co_stem(RiverPro_Chosen)
%RiverPro_Chosen是cell，每个元素都是chi-z

sumRiver=length(RiverPro_Chosen);

% 寻找stem，即Chi最长的河段
max_Chi=0;
for i=1:sumRiver
    RiverPro_i=RiverPro_Chosen{i,1};    Chi_i=RiverPro_i(:,1);      max_Chi_i=max(Chi_i);
    if max_Chi_i>max_Chi
        max_Chi=max_Chi_i;      num_max_Chi=i;
    end
end
RiverPro_Stem=RiverPro_Chosen{num_max_Chi,1};    Stem_Chi=RiverPro_Stem(:,1);   Stem_z=RiverPro_Stem(:,2);

% 计算MisFit_i
num_node=0;     MisFit_i=0;
for j=1:sumRiver
    if j ~=  num_max_Chi
        RiverProj=RiverPro_Chosen{j,1};    Chi_j=RiverProj(:,1);    z_j=RiverProj(:,2); % 都是列向量
        z_j_recal=interp1(Stem_Chi,Stem_z,Chi_j);
        MisFit_i=MisFit_i+sum((z_j-z_j_recal).^2);      num_node=num_node+length(Chi_j);
    end
end
MisFit_i=sqrt(MisFit_i/num_node);
end

%%


