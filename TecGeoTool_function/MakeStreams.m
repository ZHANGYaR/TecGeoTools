function [DEM,FD,A,S]=MakeStreams(dem,threshold_area,varargin)
% 使用方法:
%   [DEM,FD,A,S]=MakeStreams(dem,threshold_area);
%   [DEM,FD,A,S]=MakeStreams(dem,threshold_area,'name',value,...);
%
% 说明:
%  此函数接收一个数字高程模型（DEM），并输出其他TopoToolbox函数所需的基本数据集。
%  输入的DEM如果其网格分辨率（即像元大小）不是整数，有时会在伴随函数中引发问题。
%  如果提供的DEM的像元大小不是整数，代码将警告用户（但不会采取任何操作）。
%  如果你想修复像元大小的问题，可以在GIS程序中重新投影，或者使用此代码（将'resample_grid'设置为true）自动完成。
%
% 必需输入:
%   dem - DEM文件的完整路径，可以是ASCII文本文件（推荐）或GeoTIFF文件，或DEM的GRIDobj对象
%   threshold_area - 定义溪流的最小汇水面积，单位为平方米
%
% 可选输入:
%   file_name [] - 包含DEM、FD、A和S以及河流网络的shapefile的mat文件名称。
%       如果未提供file_name，则函数假定用户不希望将结果保存到mat文件（结果仍将出现在工作区）或shapefile中。
%   precip_grid [] - 可选输入降水的GRIDobj对象。如果提供此参数，代码将使用该参数生成加权流量累积网格。
%   rr_grid [] - 可选输入径流系数的GRIDobj对象。如果提供此参数，代码将使用该参数和'precip_grid'的输入生成加权流量累积网格。
%   no_data_exp [] - 定义无数据条件的输入。期望输入为一个使用变量DEM的有效等式字符串或'auto'。例如，如果你希望将小于或等于0的任何高程设置为无数据，
%   可以提供'DEM<=0'，或者如果你希望将小于500且大于1000的高程设置为无数据，
%   可以提供'DEM<500 | DEM>1000'。如果表达式无效，用户将收到警告，但代码将继续并忽略此操作。如果你提供'auto'，
%   代码将使用梯度的对数来识别真正连通的平坦区域，并将这些区域设置为NaN。如果你想更精确地控制多高程的平坦区域（例如内部排水盆地），请考虑使用'RemoveFlats'。
%   min_flat_area [1e8] - 当'no_data_exp'设置为'auto'时，用于识别平坦区域并将其设置为NaN的最小面积（单位为m^2）。如果没有调用'no_data_exp'或提供了有效的逻辑表达式，则忽略'min_flat_area'的输入。
%   resample_grid [false] - 网格重采样标志。如果未提供new_cellsize的输入，则网格将被重采样为最接近的整数大小。
%   new_cellsize [] - 新像元大小的值（单位为地图单位）。
%   mex [false] - 是否使用编译的mex程序进行流向计算的逻辑标志，仅在已为目标机器编译了TopoToolbox的mex文件时有效，请参阅TopoToolbox中的'compilemexfiles.m'了解更多信息。
%
% 输出:
%   DEM - DEM的GRIDobj对象
%   FD - 由提供的DEM生成的FLOWobj对象
%   A - 流量累积网格（GRIDobj）
%   S - 由DEM生成的STREAMobj对象
%
% 示例: 
%   [DEM,FD,A,S]=MakeStreams('/Users/forte/topo/dem.tif',1e6);
%   [DEM,FD,A,S]=MakeStreams('/Users/forte/topo/dem.tif',1e6,'file_name','AreaFiles');
%   [DEM,FD,A,S]=MakeStreams(DEMgrid,1e6,'resample_grid',true); % DEMgrid为一个GRIDobj对象	
%   [DEM,FD,A,S]=MakeStreams('/Users/forte/topo/dem.tif',1e6,'no_data_exp','DEM<=-100 | DEM>10000'); %将低于-100m或高于10000m的高程设置为NaN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数由Adam M. Forte编写 - 更新日期: 2018年6月18日
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 新增自动填洼功能 by ZhangYR  
% 更新日期：2024年9月24日
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 解析输入参数
    p = inputParser;         
    p.FunctionName = 'MakeStreams';
    addRequired(p,'dem',@(x) isa(x,'GRIDobj') | ischar(x));
    addRequired(p,'threshold_area', @(x) isscalar(x));

    addParameter(p,'file_name',[],@(x) ischar(x));
    addParameter(p,'no_data_exp',[],@(x) ischar(x) || isempty(x));
    addParameter(p,'min_flat_area',1e8,@(x) isnumeric(x) && isscalar(x));
    addParameter(p,'resample_grid',false,@(x) isscalar(x) && islogical(x));
    addParameter(p,'new_cellsize',[],@(x) isscalar(x) && isnumeric(x) || isempty(x));
    addParameter(p,'precip_grid',[],@(x) isa(x,'GRIDobj') || isempty(x));
    addParameter(p,'rr_grid',[],@(x) isa(x,'GRIDobj') || isempty(x));
    addParameter(p,'mex',false,@(x) isscalar(x) && islogical(x));
    addParameter(p,'silent',false,@(x) isscalar(x) && islogical(x));

    parse(p,dem,threshold_area,varargin{:});
    dem=p.Results.dem;
    threshold_area=p.Results.threshold_area;

    file_name=p.Results.file_name;
    no_data_exp=p.Results.no_data_exp;
    min_flat_area=p.Results.min_flat_area;
    resample_grid=p.Results.resample_grid;
    new_cellsize=p.Results.new_cellsize;
    precip_grid=p.Results.precip_grid;
    rr_grid=p.Results.rr_grid;
    mexP=p.Results.mex;
    silent=p.Results.silent;

    % 检查保存文件名
    if isempty(file_name)
        save_output=false;
    else
        save_output=true;
    end

    % 处理输入类型
    if isa(dem,'GRIDobj')
        DEM=dem;
    elseif ischar(dem)
        if ~silent
            disp('正在加载和处理DEM数据')
        end
        DEM=GRIDobj(dem);
    else
        if isdeployed
            errordlg('输入的dem未被识别为GRIDobj对象或字符路径')
        end
        error('输入的dem未被识别为GRIDobj对象或字符路径')
    end

    % 执行栅格重采样
    if resample_grid && isempty(new_cellsize)
        if ~silent
            disp('正在重采样DEM - 可能需要较长时间，请稍候')
        end
        DEM=resample(DEM,ceil(DEM.cellsize),'bicubic');
    elseif resample_grid && ~isempty(new_cellsize)
        if ~silent
            disp('正在重采样DEM - 可能需要较长时间，请稍候')
        end
        DEM=resample(DEM,new_cellsize,'bicubic');        
    end

    % 检查像元大小
    if mod(DEM.cellsize,1)~=0 && ~resample_grid
        if isdeployed
            warndlg('网格像元大小不是整数，可能引发后续处理问题，建议启用resample_grid选项')
        else
            warning('网格像元大小不是整数，可能引发后续处理问题，建议启用resample_grid选项')
        end
    end

    % 数据清洗处理
    if ~silent
        disp('正在进行DEM数据清洗')
    end
    if ~isempty(no_data_exp) && ~strcmp(no_data_exp,'auto')
        try 
            IDX=eval(no_data_exp);
            DEM.Z(IDX.Z)=nan;
            DEM=crop(DEM); % 裁剪NaN边界
        catch
            if isdeployed
                warndlg('提供的无效数据表达式无效，已跳过该过滤条件')
            else
                warning('提供的无效数据表达式无效，已跳过该过滤条件')
            end
        end
    elseif strcmp(no_data_exp,'auto')
        [DEM]=AutoFlat(DEM,min_flat_area);
    else
        DEM=crop(DEM); % 裁剪NaN边界
    end

    % 新增填洼处理
    if ~silent
        disp('正在进行填洼处理')
    end
    DEM = fillsinks(DEM);

    % 保存中间结果
    if save_output
        fileNameBase=file_name;
        MatFileName=[fileNameBase '.mat'];
        ShpFileName=[fileNameBase '.shp'];
        save(MatFileName,'DEM','-v7.3');
    end

    % 计算流向
    if ~silent
        disp('正在计算水流方向')
    end

    if mexP
        try 
            if silent
                FD=FLOWobj(DEM,'preprocess','carve','mex',true,'verbose',false);
            else
                FD=FLOWobj(DEM,'preprocess','carve','verbose',true,'mex',true);
            end
        catch
            if silent
                FD=FLOWobj(DEM,'preprocess','carve','verbose',false);
            else
                FD=FLOWobj(DEM,'preprocess','carve','verbose',true);
                disp('Mex加速不可用，已切换至普通模式')
            end
        end
    else
        if silent
            FD=FLOWobj(DEM,'preprocess','carve');
        else
            FD=FLOWobj(DEM,'preprocess','carve','verbose',true);
        end
    end

    if save_output
        save(MatFileName,'FD','-append');
    end

    % 计算流量累积
    if ~silent
        disp('正在计算流量累积')
    end

    if isempty(precip_grid)
        A=flowacc(FD);
    elseif ~isempty(precip_grid) && isempty(rr_grid)
        if ~validatealignment(precip_grid,DEM)
            precip_grid=resample(precip_grid,DEM,'nearest');
        end
        A=flowacc(FD,precip_grid);
    elseif isempty(precip_grid) && ~isempty(rr_grid)
        precip_grid=GRIDobj(DEM);
        precip_grid.Z=ones(DEM.size);
        if ~validatealignment(rr_grid,DEM)
            rr_grid=resample(rr_grid,DEM,'nearest');
        end
        A=flowacc(FD,precip_grid,rr_grid);            
    elseif ~isempty(precip_grid) && ~isempty(rr_grid)
        if ~validatealignment(precip_grid,DEM)
            precip_grid=resample(precip_grid,DEM,'nearest');
        end
        if ~validatealignment(rr_grid,DEM)
            rr_grid=resample(rr_grid,DEM,'nearest');
        end
        A=flowacc(FD,precip_grid,rr_grid);
    end

    % 提取河网
    if ~silent
        disp('正在提取河网数据')
    end
    S=STREAMobj(FD,'unit','mapunits','minarea',threshold_area);

    % 保存最终结果
    if save_output
        if ~silent
            disp('正在保存输出结果')
        end
        save(MatFileName,'A','S','-append');
        MS=STREAMobj2mapstruct(S);
        shapewrite(MS,ShpFileName);
    end
end

% 平坦区域自动处理子函数
function [DEMn] = AutoFlat(DEM,min_area)
    num_pix=round(min_area/(DEM.cellsize^2));
    LG=log10(gradient8(DEM));
    BW=isnan(LG.Z) | isinf(LG.Z);
    CC=bwconncomp(BW);
    FLATS=GRIDobj(DEM,'logical');

    for ii=1:numel(CC.PixelIdxList)
        if numel(CC.PixelIdxList{ii})>=num_pix
            idx=CC.PixelIdxList{ii};
            FLATS.Z(idx)=true;
        end
    end

    DEMn=DEM;
    DEMn.Z(FLATS.Z)=nan;
end