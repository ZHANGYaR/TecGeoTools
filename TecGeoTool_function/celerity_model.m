function celerity_model(inputFile,folderPath)
    % celerity_model函数用于根据特定的数据读取规则（从文本文件中读取相关信息），并结合用户交互输入的参数，进行类似于原函数的分析、计算及结果展示与保存操作。
    % 函数主要操作流程包括：
    %   1. 从指定文本文件读取河道ID、裂点x坐标等基础信息，按照新规则处理相同ID的情况（只记录一个ID，同ID有多行时只保存第一行的KnickDATA_x）。
    %   2. 根据河道ID读取对应的面积和河道深度数据文件，获取相应列的数据，并添加根据x坐标定位提取对应数据的功能，将面积和深度数据分别整理存储。
    %   3. 通过用户交互输入获取K、m、Age等参数，然后进行相应的计算分析（类似原函数中的波速相关分析），调整代码以适配Age为多个值的情况，先完成所有年龄运算再统一绘图，将相关结果绘制在一幅包含三个子图的图中。
    %   4. 绘制相关图表展示结果，并保存分析结果数据以及重要图形。

    % 注意事项：
    %   - 代码依赖于文本文件的正确格式以及存在性，需确保文件名为"DEM_Pnt_analysis.txt"以及形如"DEM_RiverPro+ID.txt"（其中ID为具体数值）的文件能正确被读取，且文件中各列数据顺序和类型符合代码中的使用要求。
    %   - 交互输入的参数需符合相应的数学和物理意义，以保证计算的正确性。
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 函数作者：张亚荣 - 更新日期：2024年12月5日 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 初始化存储ID、裂点x坐标的数组
    ID = [];
    KnickDATA_x = [];
    % 用于记录已经出现过的ID，方便后续判断是否重复
    recordedIDs = [];
    % 打开输入文件（inputFile），获取文件标识符
    fid = fopen(inputFile, 'r');
    if fid == -1
        error('无法打开文件，请检查文件是否存在及权限。');
    end
    % 逐行读取文件内容，解析并存储ID和裂点x坐标，按照新规则处理相同ID情况
    while ~feof(fid)
        line = fgetl(fid);
        if ~isempty(line)
            data = textscan(line, '%f%f%f%f%f%f', 'Delimiter', '\t'); % 假设文件中每行数据是6列，以制表符分隔，根据实际情况调整格式说明符
            currentID = data{1};
            % 判断当前ID是否已经记录过，如果没有记录过，则添加到ID和KnickDATA_x数组中，并记录该ID已出现过
            if ~ismember(currentID, recordedIDs)
                ID(end + 1) = currentID;
                KnickDATA_x(end + 1) = data{2};
                recordedIDs(end + 1) = currentID;
            end
        end
    end
    fclose(fid);

    % 获取河道数量，即ID的数量
    TotalNumChannel = length(ID);

    % 初始化用于存储所有河道相关数据的单元数组，第一列存面积数据，第二列存河道深度数据
    chandataSelectedAll = cell(TotalNumChannel, 2); 

     [~, fileName, ~] = fileparts(inputFile);
    prefix = strsplit(fileName, '_');
    prefix = prefix{1};

    % 循环读取每个河道对应的面积和河道深度数据文件
    for i = 1:TotalNumChannel
        % 根据当前河道ID生成对应的文件名
        riverProFileName = fullfile(folderPath, sprintf('%s_RiverPro%d.txt', prefix, ID(i)));
        % 尝试打开文件，若失败给出提示信息
        riverProFid = fopen(riverProFileName, 'r');
        if riverProFid == -1
            warning(['无法打开文件 ', riverProFileName,'，跳过该河道对应数据读取。']);
            continue;
        end
        % 初始化用于存储当前河道的面积和河道深度数据的数组
        areaData = [];
        dfmData = [];
        % 逐行读取文件，提取第6列（面积）和第3列（河道深度）的数据，并进行定位提取裂点相关数据操作
        while ~feof(riverProFid)
            line = fgetl(riverProFid);
            if ~isempty(line)
                data = textscan(line, '%f%f%f%f%f%f', 'Delimiter', '\t'); % 假设每行数据有6列，以制表符分隔，根据实际调整格式说明符
                areaData(end + 1) = data{6};
                dfmData(end + 1) = data{3};
                % 检查是否满足定位条件
                if isempty(KnickDATA_x)
                    KnickDATA_x = data{2};
                elseif any(data{2} == KnickDATA_x)
                    % 如果x的数据与riverProFid第二列相同，则提取riverProFid该行第3列的数据到KnickDATA中
                    % 这里暂时不做额外处理，因为KnickDATA的赋值逻辑在之前已确定主要通过读取inputFile来初始化
                end
            end
        end
        fclose(riverProFid);
        % 将当前河道的面积和河道深度数据分别存入chandataSelectedAll的对应列中
        chandataSelectedAll{i, 1} = areaData;
        chandataSelectedAll{i, 2} = dfmData;
    end

    K_gr = 10.^[linspace(-10, -4, 500)];
    m_gr = linspace(0.2, 1.2, 40);

    [K,m] = meshgrid(K_gr, m_gr);

    % 修改此处，提示用户输入多个年龄值，以向量形式输入，例如 [1 2 3]*10^6
     % 1. 定义输入对话框的提示信息、标题、输入框尺寸以及默认输入值
    prompt = {'请输入起始年龄：',...
              '请输入终止年龄：',...
              '请输入时间间隔：'};
    dlgtitle = '输入年龄参数相关信息';
    dims = [1 50];
    definput = {'10^6', '10^7', '10^6'}; % 设置默认输入值，可根据实际需求调整
    answer = inputdlg(prompt, dlgtitle, dims, definput);

    % 2. 将用户输入的字符串转换为数值
    startAge = str2num(answer{1});
    endAge = str2num(answer{2});
    interval = str2num(answer{3});
    if ~(isnumeric(startAge) && isnumeric(endAge) && isnumeric(interval))
        error('输入的年龄参数需要是数值，请检查输入格式。');
    end

    % 3. 根据输入的起始年龄、终止年龄和时间间隔生成年龄向量
    Age = startAge:interval:endAge;
    orgSize = size(K);

    % 用于存储每个年龄对应的所有计算结果，每个元素对应一个年龄的结果，结构与原来单个年龄计算结果类似
    AllResults = cell(length(Age), 1);

    %% 进行相关计算分析（类似原函数中的celerity_analysis部分逻辑，并添加定位提取数据功能），先完成所有年龄的运算并存储结果
    for ageIndex = 1:length(Age)
        ObsKP = KnickDATA_x; % 假设裂点x坐标对应原函数中的观测尼克点数据，这里明确赋值，可根据实际情况调整含义

        K_col = K(:);  % 转换为列矩阵
        m_col = m(:);

        FinalMatrix = cell(TotalNumChannel, 1);
        for NumChannel = 1:TotalNumChannel
            % 分别获取当前河道的面积和河道深度数据
            Area = chandataSelectedAll{NumChannel, 1}; 
            dfm = chandataSelectedAll{NumChannel, 2}; 

            dX = diff(dfm);
            Area_Avg = (Area(1:end - 1) + Area(2:end)) / 2;

            EstKP = zeros(size(K_col));
            parfor k = 1:numel(K_col)
                EstKP(k) = Estimate_KnickPoint(Area_Avg, dX, K_col(k), m_col(k), Age(ageIndex));
            end

            FinalMatrix{NumChannel} = [K_col, m_col, (EstKP - ObsKP(NumChannel)).^2, EstKP];
        end

        Misfit = zeros(size(K_col));
        for j = 1:TotalNumChannel
            Misfit = Misfit + FinalMatrix{j}(:, 3);
        end
        [Min_Misfit, IndexMin] = min(Misfit);
        C_est = K_col(IndexMin);
        p_est = m_col(IndexMin);

        V = zeros(TotalNumChannel, 1);
        for i = 1:TotalNumChannel
            EstKP_F(i) = FinalMatrix{i}(IndexMin, 4);
            V(i) = (EstKP_F(i) * 1000) / Age(ageIndex);
        end
        V_mean = mean(V);
        V_max = max(V);
        V_min = min(V);
        V_std = std(V);

        Res{1} = log10(reshape(K_col, orgSize));
        Res{2} = reshape(m_col, orgSize);
        Res{3} = reshape(log10(Misfit), orgSize);

        % 将当前年龄对应的所有计算结果存储到AllResults中
        AllResults{ageIndex} = {Min_Misfit, p_est, C_est, V_mean, V_max, V_min, V_std, Res};
    end

    %% 
    % 图1：包含三个子图展示不同参数随年龄变化情况
    figure(1)
    % 子图1：提取最小拟合误差数据
    minMisfitData = zeros(length(Age), 1);
    for ageIndex = 1:length(Age)
        minMisfitData(ageIndex) = AllResults{ageIndex}{1};
    end
    % 最小拟合误差随年龄变化
    subplot(3, 1, 1)
    hold on
    plot(Age, minMisfitData, '*-', 'LineWidth', 1.5)
    title('最小拟合误差随年龄变化')
    xlabel('年龄')
    ylabel('最小拟合误差')

    % 子图2：提取p_est数据
    pEstData = zeros(length(Age), 1);
    for ageIndex = 1:length(Age)
        pEstData(ageIndex) = AllResults{ageIndex}{2};
    end
    % p_est随年龄变化
    subplot(3, 1, 2)
    hold on
    plot(Age, pEstData, '+-', 'LineWidth', 1.5)
    title('p_{est}随年龄变化') % 添加 'Interpreter', 'latex' 确保LaTeX语法正确解析显示下标
    xlabel('年龄')
    ylabel('p_{est}')

    % 子图3：提取C_est数据
    CEstData = zeros(length(Age), 1);
    for ageIndex = 1:length(Age)
        CEstData(ageIndex) = AllResults{ageIndex}{3};
    end
    % C_est随年龄变化
    subplot(3, 1, 3)
    hold on
    plot(Age, CEstData, 'o-', 'LineWidth', 1.5)
    title('C_{est}随年龄变化')
    xlabel('年龄')
    ylabel('C_{est}')

    %% 绘制等高线图展示相关参数关系及最优参数标记
   figure(2)
    numSubplots = length(Age); % 获取年龄参数的数量，即要绘制的子图数量
    if numSubplots > 20 % 限制最多绘制20个子图（4x5布局最多容纳20个），可根据实际需求调整
        error('输入的年龄数量过多，超出当前预设的4x5子图布局容纳范围，请减少年龄数量。');
    end
    for i = 1:numSubplots
        subplot(4, 5, i) % 创建4x5布局的子图
        hold on
        contourf(AllResults{i}{8}{1}, AllResults{i}{8}{2}, AllResults{i}{8}{3}, 'LineStyle', 'none')
        plot(log10(AllResults{i}{3}), AllResults{i}{2}, '+', 'MarkerSize', 2, 'LineWidth', 10, 'MarkerFaceColor','red')
        plot(log10(AllResults{i}{3}), AllResults{i}{2}, '*r', 'MarkerSize', 4, 'LineWidth', 6, 'MarkerFaceColor','red')
        if i <= 15 % 根据子图位置合理设置坐标轴标签等，避免重复标注导致图形混乱，这里仅在每行最后一个子图显示x轴标签
            xlabel('log_{10}(C_{est})')
        end
        if i >= 16 % 仅在最后一行子图显示y轴标签
            ylabel('p_{est}')
        end
    end


    %% 绘制平均速度随年龄变化图
    % 提取平均速度数据
    VMeanData = zeros(length(Age), 1);
    for ageIndex = 1:length(Age)
        VMeanData(ageIndex) = AllResults{ageIndex}{4};
    end
    figure(3)
    hold on
    plot(Age, VMeanData, '-*', 'LineWidth', 1.5)
    title('平均速度随年龄变化')
    xlabel('年龄')
    ylabel('平均速度')

    save(fullfile(folderPath, 'All_Analysis_Results.mat'), 'Age', 'AllResults');
end

function [EstimatedKP_dist] = Estimate_KnickPoint(Area, dX, K, m, Age)
    % Estimate_KnickPoint函数用于估计裂点距离，基于给定的面积、差分距离、系数、指数以及年龄参数进行计算。
    dT = dX./(K * Area.^m);
    Cumulative_Time = cumsum(dT);
    Index = (find(Cumulative_Time <= Age));
    if isempty(Index)
        EstimatedKP_dist = 0;
    else
        EstimatedKP_dist = sum(dX(Index));
    end
end