function [MS]=Basin2Shape(DEM,location_of_data_files,varargin)
    %
    % 使用方法：
    %   [MapStruct]=Basin2Shape(DEM, location_of_data_files);
    %   [MapStruct]=Basin2Shape(DEM, location_of_data_files, 'name', value, ...);
    %
    % 描述：
    %   该函数用于将 'ProcessRiverBasins' 和 'SubDivideBigBasins' 的输出转换为一个显示多边形轮廓的单一 shapefile，并包含常见的属性字段。
    %   如果在 'ProcessRiverBasins' 中提供了附加的网格，那么这些网格的均值和标准误差值将自动填充到 shapefile 中，字段名将使用附加网格输入的第二列字符数组。
    %   此函数还允许您输入希望包含的附加字段（请参见下文的可选输入部分）。如果您更希望创建一个带有指定值的 GRIDobj，请使用 'Basin2Raster'。
    %
    % 必要输入：
    %   DEM - 用于 'ProcessRiverBasins' 的输入 DEM 的 GRIDobj，用于感兴趣的流域。
    %   location_of_data_files - 包含 'ProcessRiverBasins' 结果的 .mat 文件的文件夹的完整路径。
    %
    % 可选输入：
    %   location_of_subbasins ['SubBasins'] - 包含感兴趣子流域的文件夹的名称（如果您使用 "SubDivideBigBasins" 创建了子流域），
    %     该文件夹应位于提供的 "location_of_data_files" 文件夹中。如果没有提供正确的子流域目录名，子流域的值将不会被包含在输出中。
    %   shape_name ['basins'] - 要导出的 shapefile 的名称，必须没有空格以符合 ArcGIS 的有效名称，并且不要包括 '.shp' 后缀；
    %   include ['all'] - 用于指定要包含哪些流域构建 shapefile 的参数。默认值 'all' 会包括指定文件夹中的所有流域 .mat 文件。
    %     提供 'subdivided' 会检查给定的主流域是否通过 'SubdivideBigBasins' 被细分，只有该流域的细分版本会被包含在 shapefile 中（即原始主流域将不会包含在 shapefile 中）。
    %     提供 'bigonly' 会仅包含通过 'ProcessRiverBasins' 生成的原始流域，即使 'SubDivideBigBasins' 已经运行。如果未运行 'SubDivideBigBasins'，则 'all' 和 'bigonly' 的结果相同。
    %   extra_field_values [] - 您希望包含的额外字段值的单元格数组。该数组的第一列必须是流域编号（即 'ProcessRiverBasins' 中 RiverMouth 输入的第三列所提供的标识编号，或在 'SubDivideBigBasins' 中为流域生成的编号）。
    %     每个流域编号只能有一行数据，且在处理的所有流域中，必须为每个流域编号关联一个值。附加列将被视为您希望填充的额外字段的值。这些值可以是字符数组或数字，其他值将导致错误。
    %   extra_field_names [] - 一个 1 x m 的单元格数组，包含字段名称（字符型，不能有空格），与 extra_field_values 的字段值对应。
    %     这些名称必须与 extra_field_values 中的值的顺序一致。例如，如果您的 extra_field_values 数组有三列，包括流域编号、样本名称和侵蚀率，则您的 extra_field_names 数组应按顺序包括 'sample_name' 和 'erosion_rate'。
    %   new_concavity [] - 一个 1 x m 的数组，包含新的凹度值，用于重新计算标准化的河道陡度统计（均值、标准误差和/或标准差）。
    %     默认方法非常快速，但仅为近似。如果您希望新的 ksn 统计量在新凹度下是准确的，可以将 'new_ksn_method' 设置为 'exact'。
    %   new_ksn_method ['approximate'] - 用于控制如果提供了 'new_concavity'，如何计算新的凹度。选项有 'approximate'（默认）和 'exact'。
    %     将此选项设置为 'exact' 将大大降低计算速度。
    %   segment_length [1000] - 如果提供了新的凹度并且 'new_ksn_method' 设置为 'exact'，则此参数为 ksn 平滑距离，否则被忽略。
    %   uncertainty ['se'] - 控制包含的误差度量，期望 'se' 为标准误差（默认值），'std' 为标准差，或 'both' 同时包括标准误差和标准差。
    %   populate_categories [false] - 逻辑标志，指示是否添加条目，表示每个类别在流域内所占的百分比。如果您在运行 'ProcessRiverBasins' 时提供了 'add_cat_grids' 参数，该参数是一个地质图，包含三种单元：'Q'、'Mz' 和 'Pz'，并且将 'populate_categories' 设置为 true，
    %     那么结果的 shapefile 中将包含名为 'Q'、'Mz' 和 'Pz' 的字段，且这些字段对应每个流域中每个单元所占的百分比。如果没有在运行 'ProcessRiverBasins' 时提供 'add_cat_grids' 参数，则设置 'populate_categories' 为 true 不会产生任何影响。
    %
    % 输出：
    %   输出一个 mapstructure（MS），并保存一个 shapefile，包含以下默认字段：
    %     river_mouth - 提供给 ProcessRiverBasins 的河口编号
    %     drainage_area - 流域面积（单位：km^2）
    %     center_x - 流域中心点的 x 坐标（投影坐标）
    %     center_y - 流域中心点的 y 坐标（投影坐标）
    %     outlet_elevation - 排水点的海拔（单位：m）
    %     mean_el - 流域的平均海拔（单位：米）
    %     max_el - 流域的最高海拔（单位：米）
    %     mean_ksn - 平均河道陡度
    %     mean_gradient - 平均梯度
    %   根据 'uncertainty' 参数的值，填充海拔、ksn 和梯度的标准误差、标准差或两者。
    %   还将为任何附加网格填充均值和标准误差/标准差/两者。
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 函数编写：Adam M. Forte - 更新：06/18/18 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 解析输入
    p = inputParser;
    p.FunctionName = 'Basin2Shape';
    addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
    addRequired(p,'location_of_data_files',@(x) isdir(x));

    addParameter(p,'location_of_subbasins','SubBasins',@(x) ischar(x) || isempty(x));
    addParameter(p,'shape_name','basins',@(x) ischar(x));
    addParameter(p,'include','all',@(x) ischar(validatestring(x,{'all','subdivided','bigonly'})));
    addParameter(p,'extra_field_values',[],@(x) isa(x,'cell'));
    addParameter(p,'extra_field_names',[],@(x) isa(x,'cell') & size(x,1)==1);
    addParameter(p,'new_concavity',[],@(x) isnumeric(x));
    addParameter(p,'new_ksn_method','approximate',@(x) ischar(validatestring(x,{'approximate','exact'})));
    addParameter(p,'segment_length',1000,@(x) isscalar(x) && isnumeric(x));
    addParameter(p,'uncertainty','se',@(x) ischar(validatestring(x,{'se','std','both'})));
    addParameter(p,'populate_categories',false,@(x) isscalar(x) && islogical(x))
    addParameter(p,'suppress_shape_write',false,@(x) isscalar(x) && islogical(x))

    parse(p,DEM,location_of_data_files,varargin{:});
    DEM=p.Results.DEM;
    location_of_data_files=p.Results.location_of_data_files;

    location_of_subbasins=p.Results.location_of_subbasins;
    shape_name=p.Results.shape_name;
    include=p.Results.include;
    efv=p.Results.extra_field_values;
    efn=p.Results.extra_field_names;
    new_concavity=p.Results.new_concavity;
    new_ksn_method=p.Results.new_ksn_method;
    segment_length=p.Results.segment_length;
    uncertainty=p.Results.uncertainty;
    pc=p.Results.populate_categories;
    ssw=p.Results.suppress_shape_write;

    % 处理位置格式的可变性
    [sub_head,~,~]=fileparts(location_of_subbasins);
    if isempty(sub_head)
        location_of_subbasins=[location_of_data_files filesep location_of_subbasins];
    end

    % 根据要包含的流域进行切换
    switch include
        case 'all'
            FileList1=dir([location_of_data_files filesep '*_Data.mat']);
            FileList2=dir([location_of_subbasins filesep '*_DataSubset*.mat']);
            FileList=vertcat(FileList1,FileList2);
            num_files=numel(FileList);
        case 'bigonly'
            FileList=dir([location_of_data_files filesep '*_Data.mat']);
            num_files=numel(FileList);
        case 'subdivided'
            AllFullFiles=dir([location_of_data_files filesep '*_Data.mat']);
            num_basins=numel(AllFullFiles);
            basin_nums=zeros(num_basins,1);
            for jj=1:num_basins
                fileName=AllFullFiles(jj,1).name;
                basin_nums(jj)=sscanf(fileName,'%*6s %i'); %%%
            end

            FileCell=cell(num_basins,1);
            for kk=1:num_basins
                basin_num=basin_nums(kk);
                SearchAllString=[location_of_data_files filesep '*_' num2str(basin_num) '_Data.mat'];
                SearchSubString=[location_of_subbasins filesep '*_' num2str(basin_num) '_DataSubset*.mat'];

                if numel(dir(SearchSubString))>0
                    Files=dir(SearchSubString);
                else
                    Files=dir(SearchAllString);
                end

                FileCell{kk}=Files;
            end
            FileList=vertcat(FileCell{:});
            num_files=numel(FileList);
    end

    % 按排水面积排序以确保正确的绘制顺序
    dalist=zeros(num_files,1);
    for ii=1:num_files
        FileName=[FileList(ii,1).folder filesep FileList(ii,1).name];
        load(FileName,'drainage_area');
        dalist(ii,1)=drainage_area;
    end
    [~,six]=sort(dalist,'descend');
    FileList=FileList(six);

    % 初始化地图结构体
    MS=struct;

    % 开始主循环
    w1=waitbar(0,'构建多边形');
    for ii=1:num_files
        FileName=[FileList(ii,1).folder filesep FileList(ii,1).name];
        DB=GRIDobj(DEM);

        load(FileName,'DEMoc','RiverMouth','drainage_area','out_el','KSNc_stats','Zc_stats','Gc_stats','Centroid','hyps','Chic','DEMcc','Sc','Ac','theta_ref');

        I=~isnan(DEMoc.Z);
        [X,Y]=getcoordinates(DEMoc);
        xmat=repmat(X,numel(Y),1);
        ymat=repmat(Y,1,numel(X));

        xix=xmat(I);
        yix=ymat(I);

        ix=coord2ind(DEM,xix,yix);

        DB.Z(ix)=RiverMouth(:,3);

        [ms_temp,~,~]=GRIDobj2polygon(DB);

        % 填充输出地图结构体中的默认字段
        MS(ii,1).Geometry='Polygon';
        MS(ii,1).X=ms_temp.X;
        MS(ii,1).Y=ms_temp.Y;
        MS(ii,1).ID=ii;
        MS(ii,1).river_mouth=ms_temp.gridval;
        MS(ii,1).center_x=Centroid(1);
        MS(ii,1).center_y=Centroid(2);
        MS(ii,1).drainage_area=drainage_area;
        MS(ii,1).outlet_elevation=out_el;
        MS(ii,1).mean_el=Zc_stats(1);
        MS(ii,1).max_el=Zc_stats(5);
        MS(ii,1).mean_ksn=KSNc_stats(1);
        MS(ii,1).mean_gradient=Gc_stats(1);

        switch uncertainty
            case 'se'
                MS(ii,1).se_el=Zc_stats(2);
                MS(ii,1).se_ksn=KSNc_stats(2);
                MS(ii,1).se_gradient=Gc_stats(2);
            case 'std'
                MS(ii,1).std_el=Zc_stats(3);
                MS(ii,1).std_ksn=KSNc_stats(3);
                MS(ii,1).std_gradient=Gc_stats(3);
            case 'both'
                MS(ii,1).se_el=Zc_stats(2);
                MS(ii,1).se_ksn=KSNc_stats(2);
                MS(ii,1).se_gradient=Gc_stats(2);
                MS(ii,1).std_el=Zc_stats(3);
                MS(ii,1).std_ksn=KSNc_stats(3);
                MS(ii,1).std_gradient=Gc_stats(3);
        end

		if ~isempty(new_concavity)
			load(FileName,'MSNc');
			for jj=1:numel(new_concavity)
				switch new_ksn_method
				case 'approximate'
					[mean_ksn,std_ksn,se_ksn]=ksn_convert_approx(MSNc,new_concavity(jj));
				case 'exact'
					[mean_ksn,std_ksn,se_ksn]=ksn_convert_exact(FileName,segment_length,new_concavity(jj));
				end
				ksn_cat_name=matlab.lang.makeValidName(['mnKSN_' num2str(new_concavity(jj))]);
				MS(ii,1).(ksn_cat_name)=mean_ksn;
				switch uncertainty
				case 'se'
					ksn_cat_name_se=matlab.lang.makeValidName(['seKSN_' num2str(new_concavity(jj))]);
					MS(ii,1).(ksn_cat_name_se)=se_ksn;
				case 'std'
					ksn_cat_name_std=matlab.lang.makeValidName(['stdKSN_' num2str(new_concavity(jj))]);
					MS(ii,1).(ksn_cat_name_std)=std_ksn;
				case 'both'
					ksn_cat_name_se=matlab.lang.makeValidName(['seKSN_' num2str(new_concavity(jj))]);
					MS(ii,1).(ksn_cat_name_se)=se_ksn;
					ksn_cat_name_std=matlab.lang.makeValidName(['stdKSN_' num2str(new_concavity(jj))]);
					MS(ii,1).(ksn_cat_name_std)=std_ksn;					
				end
			end
		end		

MS(ii,1).hyp_int=double(abs(trapz((hyps(:,2)-min(hyps(:,2)))/(max(hyps(:,2))-min(hyps(:,2))),hyps(:,1)/100)));
MS(ii,1).theta=Chic.mn;

c=chiplot(Sc,DEMcc,Ac,'a0',1,'mn',theta_ref,'plot',false);
MS(ii,1).chi_r2=c.R2;

% 确定是否存在地理参考结构，如果存在则生成采样点的经纬度位置
if ~isempty(DEM.georef)
    try
        % 检查投影信息的存储方式（兼容旧版 TopoToolbox 方法）
        if isfield(DEM.georef,'mstruct')
            proj=DEM.georef.mstruct;
        else
            proj=DEM.georef;
        end
        % 进行投影反变换
        [s_lat,s_lon]=projinv(proj,RiverMouth(:,1),RiverMouth(:,2));
        % 存储河口的经纬度坐标
        MS(ii,1).outlet_lat=s_lat;
        MS(ii,1).outlet_lon=s_lon;
    catch
        str=sprintf('投影不被支持或你没有安装 Mapping Toolbox，无法将河口坐标转换为经纬度');
        disp(str);
    end
else
    str=sprintf('GRIDobj 没有投影信息，无法将河口坐标转换为经纬度');
    disp(str);
end

% 检查 ProcessRiverBasins 输出中是否存在额外的网格
VarList=whos('-file',FileName);

VarInd=find(strcmp(cellstr(char(VarList.name)),'KSNQc_stats'));
if ~isempty(VarInd)
    load(FileName,'KSNQc_stats');

    MS(ii,1).mean_ksn_q=KSNQc_stats(:,1);

    switch uncertainty
        case 'se'
            MS(ii,1).se_ksn_q=KSNQc_stats(:,2);
        case 'std'
            MS(ii,1).std_ksn_q=KSNQc_stats(:,3);
        case 'both'
            MS(ii,1).se_ksn_q=KSNQc_stats(:,2);
            MS(ii,1).std_ksn_q=KSNQc_stats(:,3);
    end
end

VarInd=find(strcmp(cellstr(char(VarList.name)),'AGc'));
if ~isempty(VarInd)
    load(FileName,'AGc','AGc_stats');
    num_grids=size(AGc,1);

    for kk=1:num_grids
        mean_prop_name=['mean_' AGc{kk,2}];
        MS(ii,1).(mean_prop_name)=double(AGc_stats(kk,1));

        switch uncertainty
            case 'se'
                se_prop_name=['se_' AGc{kk,2}];
                MS(ii,1).(se_prop_name)=double(AGc_stats(kk,2));
            case 'std'
                std_prop_name=['std_' AGc{kk,2}];
                MS(ii,1).(std_prop_name)=double(AGc_stats(kk,3));
            case 'both'
                se_prop_name=['se_' AGc{kk,2}];
                MS(ii,1).(se_prop_name)=double(AGc_stats(kk,2));
                std_prop_name=['std_' AGc{kk,2}];
                MS(ii,1).(std_prop_name)=double(AGc_stats(kk,3));
        end
    end
end

VarInd=find(strcmp(cellstr(char(VarList.name)),'rlf'));
if ~isempty(VarInd)
    load(FileName,'rlf','rlf_stats');
    num_grids=size(rlf,1);

    for kk=1:num_grids
        mean_prop_name=['mean_rlf' num2str(rlf{kk,2})];
        se_prop_name=['se_rlf' num2str(rlf{kk,2})];
        MS(ii,1).(mean_prop_name)=double(rlf_stats(kk,1));
        MS(ii,1).(se_prop_name)=double(rlf_stats(kk,2));

        switch uncertainty
            case 'se'
                se_prop_name=['se_rlf' num2str(rlf{kk,2})];
                MS(ii,1).(se_prop_name)=double(rlf_stats(kk,2));
            case 'std'
                std_prop_name=['std_rlf' num2str(rlf{kk,2})];
                MS(ii,1).(std_prop_name)=double(rlf_stats(kk,3));
            case 'both'
                se_prop_name=['se_rlf' num2str(rlf{kk,2})];
                MS(ii,1).(se_prop_name)=double(rlf_stats(kk,2));
                std_prop_name=['std_rlf' num2str(rlf{kk,2})];
                MS(ii,1).(std_prop_name)=double(rlf_stats(kk,3));
        end
    end
end

% 检查输入时是否提供了额外的字段
if ~isempty(efv)
    bnl=cell2mat(efv(:,1));

    ix=find(bnl==RiverMouth(:,3));
    % 检查是否为每个流域编号提供了唯一的条目
    if ~isempty(ix) && numel(ix)==1
        efvOI=efv(ix,2:end); % 去除流域编号列
        num_efv=size(efvOI,2);

        for kk=1:num_efv
            field_name=efn{kk};
            field_value=efvOI{kk};
            % 检查字段值是数字还是字符串
            if ischar(field_value)
                MS(ii,1).(field_name)=field_value;
            elseif isnumeric(field_value)
                MS(ii,1).(field_name)=double(field_value);
            else
                if isdeployed
                    errordlg(['为 ' field_name ' 提供的额外字段值既不是数字也不是字符串'])
                end
                error(['为 ' field_name ' 提供的额外字段值既不是数字也不是字符串']);
            end
        end
    elseif numel(ix)>1
        if isdeployed
            errordlg(['为流域 ' num2str(RiverMouth(:,3)) ' 的额外字段提供了多个条目'])
        end
        error(['为流域 ' num2str(RiverMouth(:,3)) ' 的额外字段提供了多个条目']);
    elseif isempty(ix)
        if isdeployed
            errordlg(['未为流域 ' num2str(RiverMouth(:,3)) ' 的额外字段值提供任何条目'])
        end
        error(['未为流域 ' num2str(RiverMouth(:,3)) ' 的额外字段值提供任何条目']);
    end
end

VarInd=find(strcmp(cellstr(char(VarList.name)),'ACGc'));
if ~isempty(VarInd)
    load(FileName,'ACGc','ACGc_stats');
    num_grids=size(ACGc,1);

    for kk=1:num_grids
        mode_prop_name=['mode_' ACGc{kk,3}];
        perc_prop_name=['mode_' ACGc{kk,3} '_pct'];
        ix=find(ACGc{kk,2}.Numbers==ACGc_stats(kk,1),1);
        MS(ii,1).(mode_prop_name)=ACGc{kk,2}.Categories{ix};
        total_nodes=sum(ACGc{kk,2}.Counts);
        MS(ii,1).(perc_prop_name)=double((ACGc{kk,2}.Counts(ix)/total_nodes)*100);

        if pc
            ACG_T=ACGc{kk,2};
            total_nodes=sum(ACG_T.Counts);
            for ll=1:numel(ACG_T.Categories)
                cat_name=matlab.lang.makeValidName(ACG_T.Categories{ll});
                MS(ii,1).(cat_name)=double((ACG_T.Counts(ll)/total_nodes)*100);
            end
        end
    end
end

waitbar(ii/num_files);
end
close(w1);

[head_dir,~,~]=fileparts(location_of_data_files);
if ~isempty(head_dir)
    out_shape_name=[head_dir filesep shape_name '.shp'];
else
    out_shape_name=[shape_name '.shp'];
end
shapewrite(MS,out_shape_name);
end

function [mean_ksn,std_ksn,se_ksn]=ksn_convert_approx(okm,new_ref_concavity)

    g=[okm.gradient];
    a=[okm.uparea];

    ksn_calc=g./a.^-new_ref_concavity;

    mean_ksn=mean(ksn_calc,'omitnan');
    std_ksn=std(ksn_calc,'omitnan');
    se_ksn=std_ksn/sqrt(numel(ksn_calc));

end

function [mean_ksn,std_ksn,se_ksn]=ksn_convert_exact(FN,segment_length,new_ref_concavity)
    % 确定ksn计算方法
    load(FN,'DEMoc','DEMcc','FDc','Ac','Sc','ksn_method');

    % 计算ksn
    switch ksn_method
        case 'quick'
            [MSNc]=KSN_Quick(DEMoc,DEMcc,Ac,Sc,new_ref_concavity,segment_length);
        case 'trunk'
            [MSNc]=KSN_Trunk(DEMoc,DEMcc,Ac,Sc,new_ref_concavity,segment_length,min_order);
        case 'trib'
            % 若流域面积非常小，则更改选择，因为KSN_Trib对于小流域会失效
            if drainage_area>2.5
                [MSNc]=KSN_Trib(DEMoc,DEMcc,FDc,Ac,Sc,new_ref_concavity,segment_length);
            else
                [MSNc]=KSN_Quick(DEMoc,DEMcc,Ac,Sc,new_ref_concavity,segment_length);
            end
    end

    % 计算全流域的ksn统计量
    mean_ksn=mean([MSNc.ksn],'omitnan');
    std_ksn=std([MSNc.ksn],'omitnan');
    se_ksn=std_ksn/sqrt(numel(MSNc)); % 标准误差
end

function [ksn_ms]=KSN_Quick(DEM,DEMc,A,S,theta_ref,segment_length)
    g=gradient(S,DEMc);
    G=GRIDobj(DEM);
    G.Z(S.IXgrid)=g;

    Z_RES=DEMc-DEM;

    ksn=G./(A.*(A.cellsize^2)).^(-theta_ref);

    SD=GRIDobj(DEM);
    SD.Z(S.IXgrid)=S.distance;

    ksn_ms=STREAMobj2mapstruct(S,'seglength',segment_length,'attributes',...
        {'ksn' ksn @mean 'uparea' (A.*(A.cellsize^2)) @mean 'gradient' G @mean 'cut_fill' Z_RES @mean...
        'min_dist' SD @min 'max_dist' SD @max});

    seg_dist=[ksn_ms.max_dist]-[ksn_ms.min_dist];
    distcell=num2cell(seg_dist');
    [ksn_ms(1:end).seg_dist]=distcell{:};
    ksn_ms=rmfield(ksn_ms,{'min_dist','max_dist'});
end

function [ksn_ms]=KSN_Trunk(DEM,DEMc,A,S,theta_ref,segment_length,min_order)

    order_exp=['>=' num2str(min_order)];

    Smax=modify(S,'streamorder',order_exp);
    Smin=modify(S,'rmnodes',Smax);

    g=gradient(S,DEMc);
    G=GRIDobj(DEM);
    G.Z(S.IXgrid)=g;

    Z_RES=DEMc-DEM;

    ksn=G./(A.*(A.cellsize^2)).^(-theta_ref);

    SDmax=GRIDobj(DEM);
    SDmin=GRIDobj(DEM);
    SDmax.Z(Smax.IXgrid)=Smax.distance;
    SDmin.Z(Smin.IXgrid)=Smin.distance;

    ksn_ms_min=STREAMobj2mapstruct(Smin,'seglength',segment_length,'attributes',...
        {'ksn' ksn @mean 'uparea' (A.*(A.cellsize^2)) @mean 'gradient' G @mean 'cut_fill' Z_RES @mean...
        'min_dist' SDmin @min 'max_dist' SDmin @max});

    ksn_ms_max=STREAMobj2mapstruct(Smax,'seglength',segment_length,'attributes',...
        {'ksn' ksn @mean 'uparea' (A.*(A.cellsize^2)) @mean 'gradient' G @mean 'cut_fill' Z_RES @mean...
        'min_dist' SDmax @min 'max_dist' SDmax @max});

    ksn_ms=vertcat(ksn_ms_min,ksn_ms_max);
    seg_dist=[ksn_ms.max_dist]-[ksn_ms.min_dist];
    distcell=num2cell(seg_dist');
    [ksn_ms(1:end).seg_dist]=distcell{:};
    ksn_ms=rmfield(ksn_ms,{'min_dist','max_dist'});
end

function [ksn_ms]=KSN_Trib(DEM,DEMc,FD,A,S,theta_ref,segment_length)

    % 定义不相交的河段分段
    [as]=networksegment_slim(DEM,FD,S);
    seg_bnd_ix=as.ix;
    % 预计算或提取后续所需的值
    z=getnal(S,DEMc);
    zu=getnal(S,DEM);
    z_res=z-zu;
    g=gradient(S,DEMc);
    c=chitransform(S,A,'a0',1,'mn',theta_ref);
    d=S.distance;
    da=getnal(S,A.*(A.cellsize^2));
    ixgrid=S.IXgrid;
    % 提取有序的河道节点索引列表，并找出河道间的断点
    s_node_list=S.orderednanlist;
    streams_ix=find(isnan(s_node_list));
    streams_ix=vertcat(1,streams_ix);
    % 为ksn值生成空的节点属性列表
    ksn_nal=zeros(size(d));
    % 开始遍历所有河道
    num_streams=numel(streams_ix)-1;
    seg_count=1;
    for ii=1:num_streams
        % 提取当前关注河道的节点列表
        if ii==1
            snlOI=s_node_list(streams_ix(ii):streams_ix(ii+1)-1);
        else
            snlOI=s_node_list(streams_ix(ii)+1:streams_ix(ii+1)-1);
        end

        % 确定该河道内的分段
        [~,~,dn]=intersect(snlOI,seg_bnd_ix(:,1));
        [~,~,up]=intersect(snlOI,seg_bnd_ix(:,2));
        seg_ix=intersect(up,dn);

        num_segs=numel(seg_ix);
        dn_up=seg_bnd_ix(seg_ix,:);
        for jj=1:num_segs
            % 在节点列表中找到位置
            dnix=find(snlOI==dn_up(jj,1));
            upix=find(snlOI==dn_up(jj,2));
            % 提取所需分段的索引
            seg_ix_oi=snlOI(upix:dnix);
            % 提取流程距离并归一化
            dOI=d(seg_ix_oi);
            dnOI=dOI-min(dOI);
            num_bins=ceil(max(dnOI)/segment_length);
            bin_edges=[0:segment_length:num_bins*segment_length];
            % 遍历分段内的每个子段
            for kk=1:num_bins
                idx=dnOI>bin_edges(kk) & dnOI<=bin_edges(kk+1);
                bin_ix=seg_ix_oi(idx);
                cOI=c(bin_ix);
                zOI=z(bin_ix);
                if numel(cOI)>2
                    [ksn_val,r2]=Chi_Z_Spline(cOI,zOI);
                    ksn_nal(bin_ix)=ksn_val;

                    % 构建地图结构体
                    ksn_ms(seg_count).Geometry='Line';
                    ksn_ms(seg_count).BoundingBox=[min(S.x(bin_ix)),min(S.y(bin_ix));max(S.x(bin_ix)),max(S.y(bin_ix))];
                    ksn_ms(seg_count).X=S.x(bin_ix);
                    ksn_ms(seg_count).Y=S.y(bin_ix);
                    ksn_ms(seg_count).ksn=ksn_val;
                    ksn_ms(seg_count).uparea=mean(da(bin_ix));
                    ksn_ms(seg_count).gradient=mean(g(bin_ix));
                    ksn_ms(seg_count).cut_fill=mean(z_res(bin_ix));
                    ksn_ms(seg_count).seg_dist=max(S.distance(bin_ix))-min(S.distance(bin_ix));
                    ksn_ms(seg_count).chi_r2=r2;

                    seg_count=seg_count+1;
                end
            end
        end
    end
end

function seg = networksegment_slim(DEM,FD,S)
    % TopoToolbox主库中“networksegment”函数的简化版本，同时移除零长度和单节点长度的分段

    %% 识别河道源头、汇合点、b型汇合点和出口
    Vhead = streampoi(S,'channelheads','logical');  ihead=find(Vhead==1);  IXhead=S.IXgrid(ihead);
    Vconf = streampoi(S,'confluences','logical');   iconf=find(Vconf==1);  IXconf=S.IXgrid(iconf);
    Vout = streampoi(S,'outlets','logical');        iout=find(Vout==1);    IXout=S.IXgrid(iout);
    Vbconf = streampoi(S,'bconfluences','logical'); ibconf=find(Vbconf==1);IXbconf=S.IXgrid(ibconf);

    %% 识别与b型汇合点和出口相关的流域
    DB   = drainagebasins(FD,vertcat(IXbconf,IXout));DBhead=DB.Z(IXhead); DBbconf=DB.Z(IXbconf); DBconf=DB.Z(IXconf); DBout=DB.Z(IXout);

    %% 计算流程距离
    D = flowdistance(FD);

    %% 识别河段
    % 河道源头与b型汇合点之间的连接
    [~,ind11,ind12]=intersect(DBbconf,DBhead);
    % 汇合点与b型汇合点之间的连接
    [~,ind21,ind22]=intersect(DBbconf,DBconf);
    % 河道源头与出口之间的连接
    [~,ind31,ind32]=intersect(DBout,DBhead);
    % 汇合点与出口之间的连接
    [~,ind41,ind42]=intersect(DBout,DBconf);
    % 将连接组合成河段
    IX(:,1) = [ IXbconf(ind11)' IXbconf(ind21)' IXout(ind31)'  IXout(ind41)'  ];   ix(:,1)= [ ibconf(ind11)' ibconf(ind21)' iout(ind31)'  iout(ind41)'  ];
    IX(:,2) = [ IXhead(ind12)'  IXconf(ind22)'  IXhead(ind32)' IXconf(ind42)' ];   ix(:,2)= [ ihead(ind12)'  iconf(ind22)'  ihead(ind32)' iconf(ind42)' ];

    % 计算河段的流程长度
    flength=double(abs(D.Z(IX(:,1))-D.Z(IX(:,2))));

    % 移除零长度和单节点长度的元素
    idx=flength>=2*DEM.cellsize;
    seg.IX=IX(idx,:);
    seg.ix=ix(idx,:);
    seg.flength=flength(idx);

    % 河段数量
    seg.n=numel(IX(:,1));
end

function [KSN,R2] = Chi_Z_Spline(c,z)

    % 使用三次样条插值对chi-高程关系进行重采样
    [~,minIX]=min(c);
    zb=z(minIX);
    chiF=c-min(c);
    zabsF=z-min(z);
    chiS=linspace(0,max(chiF),numel(chiF)).';
    zS=spline(chiF,zabsF,chiS);

    % 通过斜率计算ksn
    KSN= chiS\(zS); % 由于a0固定为1，不需要mn

    % 计算R²
    z_pred=chiF.*KSN;
    sstot=sum((zabsF-mean(zabsF)).^2);
    ssres=sum((zabsF-z_pred).^2);
    R2=1-(ssres/sstot);

end