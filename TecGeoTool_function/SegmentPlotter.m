function SegmentPlotter(basin_nums, varargin)
    %
    % 用法:
    %   SegmentPlotter(basin_nums);
    %   SegmentPlotter(basin_nums,'name',value,...);
    %
    % 描述:
    %   此函数用于绘制由 'SegmentPicker' 函数挑选出的河网段的 chi-Z 关系、纵向剖面图和坡度-面积图。
    %
    % 必须输入:
    %   basin_nums - 需要一起绘制的流域编号的行向量或列向量。假定 .mat 文件位于 MATLAB 路径的活动目录下。
    % 
    % 可选输入:
    %   separate [false] - 是否将所有河段作为单独的图形绘制。
    %   subset [] - 指定的河流编号列表，仅适用于单一流域编号。
    %   label [false] - 是否标注单独的河流编号。
    %   names [] - 当 'label' 为 true 时，为流域中的河流添加标识名称。
    %   in_dir [pwd] - 输入数据所在目录。
    %   out_dir [pwd] - 图形保存目录。
    %   save_fig [false] - 是否保存图形为 PDF 文件。
    %
    % 示例:
    %   SegmentPlotter([4,5],'label',true,'names',{'big','little'});
    %   SegmentPlotter(4,'label',true,'names','big');
    %   SegmentPlotter(1,'subset',[3,5,7,8]);
    
    % 解析输入参数
    p = inputParser;
    p.FunctionName = 'SegmentPlotter';
    addRequired(p, 'basin_nums', @(x) isnumeric(x));
    addParameter(p, 'separate', false, @(x) islogical(x));
    addParameter(p, 'subset', [], @(x) isnumeric(x) || isempty(x));
    addParameter(p, 'label', false, @(x) islogical(x));
    addParameter(p, 'names', [], @(x) ischar(x) || iscell(x) || isempty(x));
    addParameter(p, 'in_dir', pwd, @(x) ischar(x) && isfolder(x));
    addParameter(p, 'out_dir', pwd, @(x) ischar(x) && isfolder(x));
    addParameter(p, 'save_fig', false, @(x) islogical(x));

    parse(p, basin_nums, varargin{:});
    params = p.Results;

    % 检查subset参数使用条件
    if ~isempty(params.subset) && numel(basin_nums) > 1
        warning('参数 subset 仅在提供单一流域编号时适用。忽略 subset 输入。');
        params.subset = [];
    end

    % 确定颜色映射
    num_basins = numel(basin_nums);
    cMap = jet(num_basins);

    % 处理每个流域
    for ii = 1:num_basins
        basin_num = basin_nums(ii);
        fileName = fullfile(params.in_dir, ['PickedSegments_' num2str(basin_num) '.mat']);
        if ~isfile(fileName)
            warning('文件 %s 不存在。跳过流域 %d。', fileName, basin_num);
            continue;
        end
        data = load(fileName);
        
        % 确定河段编号
        if isfield(data, 'Heads')
            seg_num = data.Heads(:, 3);
        elseif isfield(data, 'Outlets')
            seg_num = data.Outlets(:, 3);
        else
            warning('在 %s 中未找到河段信息。', fileName);
            continue;
        end

        % 根据subset筛选河段
        if ~isempty(params.subset)
            idx = ismember(seg_num, params.subset);
            data.ChiSgmnts = data.ChiSgmnts(idx);
            data.SlpAreaSgmnts = data.SlpAreaSgmnts(idx, :);
            seg_num = seg_num(idx);
        end
        
        % 绘制河段
        if ~params.separate
            plotSegmentsTogether(data, seg_num, params, basin_num, cMap(ii, :));
        else
            plotSegmentsSeparately(data, seg_num, params, basin_num, ii);
        end
    end
end

function plotSegmentsTogether(data, seg_num, params, basin_num, color)
    % 创建图形窗口
    f = figure;
    set(f, 'Units', 'inches', 'Position', [1.0, 1, 10, 10], 'renderer', 'painters', 'PaperSize', [10, 10]);

    % 绘制每个河段
    for jj = 1:numel(data.ChiSgmnts)
        C = data.ChiSgmnts{jj};
        bSA = data.SlpAreaSgmnts{jj, 1};
        aSA = data.SlpAreaSgmnts{jj, 2};

        % 绘制 chi-Z 关系图
        subplot(3, 1, 1);
        hold on;
        plot(C.chi, C.elev, 'Color', color);
        if params.label
            [mc, ix] = max(C.chi);
            text(mc, C.elev(ix), num2str(seg_num(jj)), 'Color', color);
        end
        xlabel('\chi');
        ylabel('高程 (米)');
        title('\chi - Z 关系图');
        hold off;

        % 绘制纵向剖面图
        subplot(3, 1, 2);
        hold on;
        plot(C.distance, C.elev, 'Color', color);
        xlabel('距河口距离 (米)');
        ylabel('高程 (米)');
        title('纵向剖面图');
        hold off;

        % 绘制坡度-面积图
        subplot(3, 1, 3);
        hold on;
        scatter(aSA(:, 2), aSA(:, 1), 5, color, '+');
        scatter(bSA(:, 2), bSA(:, 1), 20, 'k', 'filled');
        set(gca, 'XScale', 'log', 'YScale', 'log', 'XDir', 'reverse');
        xlabel('对数流域面积');
        ylabel('对数坡度');
        title('坡度-面积图');
        hold off;
    end
    
    % 如果指定了保存图形
    if params.save_fig
        figFile = fullfile(params.out_dir, ['StreamSegments_' num2str(basin_num) '.pdf']);
        print(f, '-dpdf', figFile, '-fillpage');
    end
end

function plotSegmentsSeparately(data, seg_num, params, basin_num, fig_num)
    % 单独绘制每个河段
    for jj = 1:numel(data.ChiSgmnts)
        f = figure;
        set(f, 'Units', 'inches', 'Position', [1.0, 1, 10, 10], 'renderer', 'painters', 'PaperSize', [10, 10]);
        C = data.ChiSgmnts{jj};
        bSA = data.SlpAreaSgmnts{jj, 1};
        aSA = data.SlpAreaSgmnts{jj, 2};

        % 绘制 chi-Z 关系图
        subplot(3, 1, 1);
        plot(C.chi, C.elev, '-k');
        xlabel('\chi');
        ylabel('高程 (米)');
        title(sprintf('流域 %d - 河段 %d', basin_num, seg_num(jj)));

        % 绘制纵向剖面图
        subplot(3, 1, 2);
        plot(C.distance, C.elev, '-k');
        xlabel('距河口距离 (米)');
        ylabel('高程 (米)');

        % 绘制坡度-面积图
        subplot(3, 1, 3);
        scatter(aSA(:, 2), aSA(:, 1), 5, [0.5, 0.5, 0.5], '+');
        scatter(bSA(:, 2), bSA(:, 1), 20, 'k', 'filled');
        set(gca, 'XScale', 'log', 'YScale', 'log', 'XDir', 'reverse');
        xlabel('对数流域面积');
        ylabel('对数坡度');

        % 如果指定了保存图形
        if params.save_fig
            figFile = fullfile(params.out_dir, sprintf('StreamSegments_%d_%d.pdf', basin_num, fig_num));
            print(f, '-dpdf', figFile, '-fillpage');
        end
        fig_num = fig_num + 1;
    end
end