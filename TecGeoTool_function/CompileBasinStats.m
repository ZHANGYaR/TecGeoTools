function [T]=CompileBasinStats(location_of_data_files,varargin)
%
% 用法：
%	[T]=CompileBasinStats(location_of_data_files);
%	[T]=CompileBasinStats(location_of_data_files,'name',value,...);
%
% 描述：
% 	该函数用于从 'ProcessRiverBasins' 和 'SubDivideBigBasins' 的输出中生成总结性统计表格。该表格是 'BasinStatsPlots' 的必要输入。若在 'ProcessRiverBasins' 中提供了额外栅格数据，表格中将包含这些栅格的均值和标准误差。可通过可选参数添加自定义字段（详见下文）。若在 'ProcessRiverBasins' 中提供了分类栅格，则可计算额外的分类相关参数。
%
% 必需输入：
% 		location_of_data_files - 包含 'ProcessRiverBasins' 输出的 .mat 文件的完整路径。
%
% 可选输入：
%		location_of_subbasins ['SubBasins'] - 子流域文件夹名称（当使用 "SubDivideBigBasins" 创建子流域时指定）。该目录应位于主流域目录内。若未正确指定子流域目录，无论 'include' 参数如何设置，表格中将不包含子流域数据。
%		include ['all'] - 控制包含的流域类型。可选值：
%			'all' - 包含指定目录下所有流域 .mat 文件（默认）
%			'subdivided' - 仅包含经 'SubdivideBigBasins' 细分的子流域（不包含原主流域）
%			'bigonly' - 仅包含原始主流域（即使存在细分流域）
%		extra_field_values [] - 自定义字段值的单元格数组。第一列必须为流域编号（对应 'ProcessRiverBasins' 的 RiverMouth 第三列或 'SubDivideBigBasins' 生成的编号），每流域仅允许一行数据。其他列为自定义字段值，支持字符数组或数值类型。
%		extra_field_names [] - 1×m 单元格数组，定义自定义字段名称（需符合 shapefile 属性命名规范，不含空格）。字段顺序需与 extra_field_values 中的列对应。
%		new_concavity [] - 新凹度值数组，用于重新计算标准化河道陡度统计量。默认采用快速近似算法，若需精确计算，请将 'new_ksn_method' 设为 'exact'。
%		new_ksn_method ['approximate'] - 新凹度计算方法，可选：
%			'approximate' - 快速近似计算（默认）
%			'exact' - 精确计算（显著增加计算时间）
%		segment_length [1000] - 精确计算时的河道平滑长度（单位：米），仅当 new_ksn_method 为 'exact' 时生效。
%		uncertainty ['se'] - 不确定性度量方式，可选：
%			'se' - 标准误差（默认）
%			'std' - 标准差
%			'both' - 同时包含标准误差和标准差
%		dist_along_azimuth [] - 指定方位角（0-360度）计算流域沿该方向的延伸距离。
%		filter_by_category [false] - 是否根据分类栅格重新计算均值。需配合 'filter_type', 'cat_grid' 和 'cat_values' 使用。默认对河道陡度、坡度和高程进行过滤计算，同时处理所有附加栅格。
%		filter_type ['exclude'] - 分类过滤方式，可选：
%			'exclude' - 排除指定类别区域
%			'include' - 仅包含指定类别区域
%			'mode' - 按流域主导类别自动过滤（无需指定 cat_values）
%		cat_grid [] - 分类栅格名称，需与 'ProcessRiverBasins' 中 'add_cat_grids' 输入的第三列一致。
%		cat_values [] - 1×m 单元格数组，指定过滤使用的类别值（需匹配分类栅格的类别编码）。
%		populate_categories [false] - 是否添加分类占比字段。例如地质分类栅格将生成各岩性单元的面积百分比字段。
%		means_by_category [] - 按类别计算均值的栅格列表。格式为 1×m 单元格数组，首元素为分类栅格名称，后续为待统计的栅格名称。例：{'geology','ksn','rlf2500'} 表示按地质分类统计河道陡度及2.5km²起伏度的均值。
%
% 输出：
%		T - 结构统计表格，包含以下默认字段：
%			river_mouth : 流域出口编号
%			drainage_area : 集水面积（km²）
%			out_x/out_y : 出口坐标
%			center_x/center_y : 流域质心坐标
%			outlet_elevation : 出口高程（m）
%			mean_el/max_el : 平均/最大高程（m）
%			mean_ksn : 平均河道陡度
%			mean_gradient : 平均坡度
%			根据 uncertainty 设置，包含相应误差统计量。
%			所有附加栅格均添加对应的均值及误差统计字段。
%
% 示例：
%		% 基本用法
%		[T]=CompileBasinStats('/Users/You/basin_files');
%
%		% 按地质分类统计坡度和起伏度
%		[T]=CompileBasinStats('/Users/You/basin_files','means_by_category',{'geology','gradient','rlf2500'});
%
%		% 使用多凹度值重新计算
%		[T]=CompileBasinStats('/Users/You/basin_files','new_concavity',[0.45 0.55]);
%
%		% 排除冲积层区域计算统计量
%		[T]=CompileBasinStats('/Users/You/basin_files','filter_by_category',true,'cat_grid','geology','cat_values',{'Q'},'filter_type','exclude');
%
% 注意：
%		- 分类过滤后的河道陡度值基于插值计算（KsnOBJc），与原始 mean_ksn 的计算方法不同
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数作者：Adam M. Forte - 最后更新于06/18/18 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	

	% 解析输入参数
	p = inputParser;
	p.FunctionName = 'CompileBasinStats';
	addRequired(p,'location_of_data_files',@(x) isdir(x));

	addParameter(p,'location_of_subbasins','SubBasins',@(x) ischar(x) || isempty(x));
	addParameter(p,'include','all',@(x) ischar(validatestring(x,{'all','subdivided','bigonly'})));
	addParameter(p,'extra_field_values',[],@(x) isa(x,'cell') || isempty(x));
	addParameter(p,'extra_field_names',[],@(x) isa(x,'cell') && size(x,1)==1 || isempty(x));
	addParameter(p,'new_concavity',[],@(x) isnumeric(x));
	addParameter(p,'new_ksn_method','approximate',@(x) ischar(validatestring(x,{'approximate','exact'})));
	addParameter(p,'segment_length',1000,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'dist_along_azimuth',[],@(x) isnumeric(x) && isscalar(x) && x>=0 && x<=360 || isempty(x));
	addParameter(p,'uncertainty','se',@(x) ischar(validatestring(x,{'se','std','both'})));
	addParameter(p,'populate_categories',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'means_by_category',[],@(x) isa(x,'cell') && size(x,2)>=2 || isempty(x));
	addParameter(p,'filter_by_category',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'filter_type','exclude',@(x) ischar(validatestring(x,{'exclude','include','mode'})));
	addParameter(p,'cat_grid',[],@(x) ischar(x));
	addParameter(p,'cat_values',[],@(x) isa(x,'cell') && size(x,1)==1);

	parse(p,location_of_data_files,varargin{:});
	location_of_data_files=p.Results.location_of_data_files;

	location_of_subbasins=p.Results.location_of_subbasins;
	include=p.Results.include;
	efv=p.Results.extra_field_values;
	efn=p.Results.extra_field_names;
	new_concavity=p.Results.new_concavity;
	new_ksn_method=p.Results.new_ksn_method;
	segment_length=p.Results.segment_length;
	az=p.Results.dist_along_azimuth;
	uncertainty=p.Results.uncertainty;
	pc=p.Results.populate_categories;
	mbc=p.Results.means_by_category;
	fbc=p.Results.filter_by_category;
	ft=p.Results.filter_type;
	cgn=p.Results.cat_grid;
	cgv=p.Results.cat_values;

	% 检查必要参数完整性
	if fbc & ~strcmp(ft,'mode') && isempty(cgn) |isempty(cgv)
		if isdeployed
			errordlg('使用"include"或"exclude"过滤时，必须同时提供"cat_grid"和"cat_values"')
		end
		error('使用"include"或"exclude"过滤时，必须同时提供"cat_grid"和"cat_values"');
	elseif fbc & strcmp(ft,'mode') & isempty(cgn)
		if isdeployed
			errordlg('使用"mode"过滤时，必须提供"cat_grid"')
		end
		error('使用"mode"过滤时，必须提供"cat_grid"');
	end

	% 处理路径格式
	if ~isempty(location_of_subbasins)
		[sub_head,~,~]=fileparts(location_of_subbasins);
		if isempty(sub_head)
			location_of_subbasins=[location_of_data_files filesep location_of_subbasins];
		end
	end

	% 确定包含的流域类型
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

	% 空文件列表检查
	if num_files==0
		error('输入目录有效但未找到流域文件，请检查路径设置');
	end

	% 初始化表格
	T=table;

	% 初始化进度条
	if ~isempty(mbc)
		w1=waitbar(0,'正在编译表格并按类别计算均值');
	elseif fbc
		w1=waitbar(0,'正在编译表格并计算过滤均值');
	else
		w1=waitbar(0,'正在编译表格');
	end

	warning off
	for ii=1:num_files
		FileName=[FileList(ii,1).folder filesep FileList(ii,1).name];

		load(FileName,'DEMoc','RiverMouth','drainage_area','out_el','KSNc_stats','Zc_stats','Gc_stats','Centroid','hyps','Chic','DEMcc','Sc','Ac','theta_ref');

		% 填充基础字段
		T.ID(ii,1)=ii;
		T.river_mouth(ii,1)=RiverMouth(3);
		T.out_x(ii,1)=RiverMouth(1);
		T.out_y(ii,1)=RiverMouth(2);
		T.center_x(ii,1)=Centroid(1);
		T.center_y(ii,1)=Centroid(2);
		T.drainage_area(ii,1)=drainage_area;
		T.outlet_elevation(ii,1)=out_el;
		T.mean_el(ii,1)=Zc_stats(1);
		T.max_el(ii,1)=Zc_stats(5);
		switch uncertainty
		case 'se'
			T.se_el(ii,1)=Zc_stats(2);
		case 'std'
			T.std_el(ii,1)=Zc_stats(3);
		case 'both'
			T.se_el(ii,1)=Zc_stats(2);
			T.std_el(ii,1)=Zc_stats(3);
		end

		T.mean_ksn(ii,1)=KSNc_stats(1);
		switch uncertainty
		case 'se'
			T.se_ksn(ii,1)=KSNc_stats(2);
		case 'std'
			T.std_ksn(ii,1)=KSNc_stats(3);
		case 'both'
			T.se_ksn(ii,1)=KSNc_stats(2);
			T.std_ksn(ii,1)=KSNc_stats(3);
		end

		% 处理新凹度计算
		if ~isempty(new_concavity)
			load(FileName,'MSNc');
			for jj=1:numel(new_concavity)
				switch new_ksn_method
				case 'approximate'
					[mean_ksn,std_ksn,se_ksn]=ksn_convert_approx(MSNc,new_concavity(jj));
				case 'exact'
					[mean_ksn,std_ksn,se_ksn]=ksn_convert_exact(FileName,segment_length,new_concavity(jj));
				end
				ksn_cat_name=matlab.lang.makeValidName(['mean_ksn_' num2str(new_concavity(jj))]);
				T.(ksn_cat_name)(ii,1)=mean_ksn;
				switch uncertainty
				case 'se'
					ksn_cat_name_se=matlab.lang.makeValidName(['se_ksn_' num2str(new_concavity(jj))]);
					T.(ksn_cat_name_se)(ii,1)=se_ksn;
				case 'std'
					ksn_cat_name_std=matlab.lang.makeValidName(['std_ksn_' num2str(new_concavity(jj))]);
					T.(ksn_cat_name_std)(ii,1)=std_ksn;
				case 'both'
					ksn_cat_name_se=matlab.lang.makeValidName(['se_ksn_' num2str(new_concavity(jj))]);
					T.(ksn_cat_name_se)(ii,1)=se_ksn;
					ksn_cat_name_std=matlab.lang.makeValidName(['std_ksn_' num2str(new_concavity(jj))]);
					T.(ksn_cat_name_std)(ii,1)=std_ksn;					
				end
			end
		end

		% 填充坡度统计
		T.mean_gradient(ii,1)=Gc_stats(1);
		switch uncertainty
		case 'se'
			T.se_gradient(ii,1)=Gc_stats(2);
		case 'std'
			T.std_gradient(ii,1)=Gc_stats(3);
		case 'both'
			T.se_gradient(ii,1)=Gc_stats(2);
			T.std_gradient(ii,1)=Gc_stats(3);
		end

		% 填充地形积分曲线
		T.hypsometry{ii,1}=hyps;
		T.hyp_integral(ii,1)=abs(trapz((hyps(:,2)-min(hyps(:,2)))/(max(hyps(:,2))-min(hyps(:,2))),hyps(:,1)/100));
		T.concavity(ii,1)=Chic.mn;

		% 计算χ方统计量
		c=chiplot(Sc,DEMcc,Ac,'a0',1,'mn',theta_ref,'plot',false);
		T.chi_R_squared(ii,1)=c.R2;
		c_trunk=chiplot(trunk(Sc),DEMcc,Ac,'a0',1,'mn',theta_ref,'plot',false);
		T.chi_R_squared_trunk(ii,1)=c_trunk.R2;

		% 处理附加栅格
		VarList=whos('-file',FileName);
		AgInd=find(strcmp(cellstr(char(VarList.name)),'AGc'));
		RlfInd=find(strcmp(cellstr(char(VarList.name)),'rlf'));
		AcgInd=find(strcmp(cellstr(char(VarList.name)),'ACGc'));
		KsnQInd=find(strcmp(cellstr(char(VarList.name)),'KSNQc_stats'));

		if ~isempty(KsnQInd)
			load(FileName,'KSNQc_stats');

			T.mean_ksn_q(ii,1)=KSNQc_stats(:,1);

			switch uncertainty
			case 'se'
				T.se_ksn_q(ii,1)=KSNQc_stats(:,2);
			case 'std'
				T.std_ksn_q(ii,1)=KSNQc_stats(:,3);
			case 'both'
				T.se_ksn_q(ii,1)=KSNQc_stats(:,2);
				T.std_ksn_q(ii,1)=KSNQc_stats(:,3);
			end
		end

		% 处理附加栅格统计
		if ~isempty(AgInd)
			load(FileName,'AGc','AGc_stats');
			num_grids=size(AGc,1);

			for kk=1:num_grids
				mean_prop_name=['mean_' AGc{kk,2}];		
				T.(mean_prop_name)(ii,1)=double(AGc_stats(kk,1));

				switch uncertainty
				case 'se'
					se_prop_name=['se_' AGc{kk,2}];
					T.(se_prop_name)(ii,1)=double(AGc_stats(kk,2));
				case 'std'
					std_prop_name=['std_' AGc{kk,2}];
					T.(std_prop_name)(ii,1)=double(AGc_stats(kk,3));
				case 'both'
					se_prop_name=['se_' AGc{kk,2}];
					T.(se_prop_name)(ii,1)=double(AGc_stats(kk,2));
					std_prop_name=['std_' AGc{kk,2}];
					T.(std_prop_name)(ii,1)=double(AGc_stats(kk,3));
				end
			end

			ag_flag=true;
		else
			ag_flag=false;
		end		

		% 处理地形起伏统计
		if ~isempty(RlfInd)
			load(FileName,'rlf','rlf_stats');
			num_grids=size(rlf,1);

			for kk=1:num_grids
				mean_prop_name=['mean_rlf' num2str(rlf{kk,2})];
				T.(mean_prop_name)(ii,1)=double(rlf_stats(kk,1));

				switch uncertainty
				case 'se'
					se_prop_name=['se_rlf' num2str(rlf{kk,2})];
					T.(se_prop_name)(ii,1)=double(rlf_stats(kk,2));
				case 'std'
					std_prop_name=['std_rlf' num2str(rlf{kk,2})];
					T.(std_prop_name)(ii,1)=double(rlf_stats(kk,3));
				case 'both'
					se_prop_name=['se_rlf' num2str(rlf{kk,2})];
					T.(se_prop_name)(ii,1)=double(rlf_stats(kk,2));
					std_prop_name=['std_rlf' num2str(rlf{kk,2})];
					T.(std_prop_name)(ii,1)=double(rlf_stats(kk,3));
				end
			end
			rlf_flag=true;
		else
			rlf_flag=false;
		end		

		% 分类过滤计算
		if fbc && ~isempty(AcgInd)
			load(FileName,'ACGc');
			% 提取分类栅格及编码表
			cix=find(strcmp(ACGc(:,3),cgn));
			CG=ACGc{cix,1};
			cgt=ACGc{cix,2};
			% 创建过滤掩膜
			F=GRIDobj(CG,'logical');
			if strcmp(ft,'include')
				% 包含指定类别
				vcix=ismember(cgt.Categories,cgv);
				vnix=cgt.Numbers(vcix);
				F.Z=ismember(CG.Z,vnix);
			elseif strcmp(ft,'exclude')
				% 排除指定类别
				vcix=ismember(cgt.Categories,cgv);
				vnix=cgt.Numbers(vcix);
				F.Z=~ismember(CG.Z,vnix);
			elseif strcmp(ft,'mode')
				% 按主导类别过滤
				load(FileName,'ACGc_stats');
				F.Z=ismember(CG.Z,ACGc_stats(cix,1));
			end
			% 应用过滤计算
			load(FileName,'DEMoc','Goc','MSNc');

			T.mean_el_f(ii,1)=mean(DEMoc.Z(F.Z),'omitnan');
			switch uncertainty
			case 'se'
				T.se_el_f(ii,1)=std(DEMoc.Z(F.Z),'omitnan')/sqrt(sum(~isnan(DEMoc.Z(F.Z))));
			case 'std'
				T.std_el_f(ii,1)=std(DEMoc.Z(F.Z),'omitnan');
			case 'both'
				T.se_el_f(ii,1)=std(DEMoc.Z(F.Z),'omitnan')/sqrt(sum(~isnan(DEMoc.Z(F.Z))));
				T.std_el_f(ii,1)=std(DEMoc.Z(F.Z),'omitnan');
			end

			T.mean_gradient_f(ii,1)=mean(Goc.Z(F.Z),'omitnan');
			switch uncertainty
			case 'se'
				T.se_gradient_f(ii,1)=std(Goc.Z(F.Z),'omitnan')/sqrt(sum(~isnan(Goc.Z(F.Z))));
			case 'std'
				T.std_gradient_f(ii,1)=std(Goc.Z(F.Z),'omitnan');
			case 'both'
				T.se_gradient_f(ii,1)=std(Goc.Z(F.Z),'omitnan')/sqrt(sum(~isnan(Goc.Z(F.Z))));
				T.std_gradient_f(ii,1)=std(Goc.Z(F.Z),'omitnan');
			end

			% 计算过滤后的河道陡度
			KSNG=GRIDobj(CG);
			KSNG.Z(:,:)=NaN;
			for kk=1:numel(MSNc)
				ix=coord2ind(CG,MSNc(kk).X,MSNc(kk).Y);
				KSNG.Z(ix)=MSNc(kk).ksn;
			end

			T.mean_ksn_f(ii,1)=mean(KSNG.Z(F.Z),'omitnan');
			switch uncertainty
			case 'se'
				T.se_ksn_f(ii,1)=std(KSNG.Z(F.Z),'omitnan')/sqrt(sum(~isnan(KSNG.Z(F.Z))));
			case 'std'
				T.std_ksn_f(ii,1)=std(KSNG.Z(F.Z),'omitnan');
			case 'both'
				T.se_ksn_f(ii,1)=std(KSNG.Z(F.Z),'omitnan')/sqrt(sum(~isnan(KSNG.Z(F.Z))));
				T.std_ksn_f(ii,1)=std(KSNG.Z(F.Z),'omitnan');
			end

			% 处理附加栅格的过滤统计
			if ag_flag
				ag_grids=size(AGc,1);
				for kk=1:ag_grids
					agG=AGc{kk,1};
					mean_prop_name=['mean_' AGc{kk,2} '_f'];		
					T.(mean_prop_name)(ii,1)=mean(agG.Z(F.Z),'omitnan');

					switch uncertainty
					case 'se'
						se_prop_name=['se_' AGc{kk,2} '_f'];
						T.(se_prop_name)(ii,1)=std(agG.Z(F.Z),'omitnan')/sqrt(sum(~isnan(agG.Z(F.Z))));
					case 'std'
						std_prop_name=['std_' AGc{kk,2} '_f'];
						T.(std_prop_name)(ii,1)=std(agG.Z(F.Z),'omitnan');
					case 'both'
						se_prop_name=['se_' AGc{kk,2} '_f'];
						T.(se_prop_name)(ii,1)=std(agG.Z(F.Z),'omitnan')/sqrt(sum(~isnan(agG.Z(F.Z))));
						std_prop_name=['std_' AGc{kk,2} '_f'];
						T.(std_prop_name)(ii,1)=std(agG.Z(F.Z),'omitnan');
					end
				end
			end

			if rlf_flag
				rlf_grids=size(rlf,1);
				for kk=1:rlf_grids
					rlfG=rlf{kk,1};
					mean_prop_name=['mean_rlf' num2str(rlf{kk,2}) '_f'];
					T.(mean_prop_name)(ii,1)=mean(rlfG.Z(F.Z),'omitnan');

					switch uncertainty
					case 'se'
						se_prop_name=['se_rlf' num2str(rlf{kk,2}) '_f'];
						T.(se_prop_name)(ii,1)=std(rlfG.Z(F.Z),'omitnan')/sqrt(sum(~isnan(rlfG.Z(F.Z))));
					case 'std'
						std_prop_name=['std_rlf' num2str(rlf{kk,2}) '_f'];
						T.(std_prop_name)(ii,1)=std(rlfG.Z(F.Z),'omitnan');
					case 'both'
						se_prop_name=['se_rlf' num2str(rlf{kk,2}) '_f'];
						T.(se_prop_name)(ii,1)=std(rlfG.Z(F.Z),'omitnan')/sqrt(sum(~isnan(rlfG.Z(F.Z))));
						std_prop_name=['std_rlf' num2str(rlf{kk,2})];
						T.(std_prop_name)(ii,1)=std(rlfG.Z(F.Z),'omitnan');
					end
				end
			end
			% 生成过滤记录列
			if ~strcmp(ft,'mode')
				filt_name=join(cgv);
				filt_name=filt_name{1};
				T.filter{ii,1}=[ft ' ' filt_name];
			else
				T.filter{ii,1}=[cgn ' mode'];
			end

		elseif fbc && isempty(AcgInd)
			if isdeployed
				errordlg('未提供分类栅格数据，无法计算过滤值')
			end
			error('未提供分类栅格数据，无法计算过滤值');
		end

		% 检查自定义字段输入
		if ~isempty(efv)
			bnl=cell2mat(efv(:,1));

			ix=find(bnl==RiverMouth(:,3));
			% 确保每个流域编号仅对应一行数据
			if ~isempty(ix) && numel(ix)==1
				efvOI=efv(ix,2:end); % 移除流域编号列
				num_efv=size(efvOI,2);

				for kk=1:num_efv
					field_name=efn{kk};
					field_value=efvOI{kk};
					% 检查字段值类型
					if ischar(field_value)
						T.(field_name){ii,1}=field_value;
					elseif isnumeric(field_value)
						T.(field_name)(ii,1)=double(field_value);
					else
						if isdeployed
							errordlg(['字段 ' field_name ' 的值必须为数值或字符类型'])
						end
						error(['字段 ' field_name ' 的值必须为数值或字符类型']);
					end
				end
			elseif numel(ix)>1
				if isdeployed
					errordlg(['流域 ' num2str(RiverMouth(:,3)) ' 的自定义字段存在重复条目'])
				end
				error(['流域 ' num2str(RiverMouth(:,3)) ' 的自定义字段存在重复条目']);
			elseif isempty(ix)
				if isdeployed
					errordlg(['未找到流域 ' num2str(RiverMouth(:,3)) ' 的自定义字段'])
				end
				error(['未找到流域 ' num2str(RiverMouth(:,3)) ' 的自定义字段']);
			end
		end

		% 处理分类栅格统计
		if ~isempty(AcgInd)
			load(FileName,'ACGc','ACGc_stats');
			num_grids=size(ACGc,1);

			for kk=1:num_grids
				mode_prop_name=['mode_' ACGc{kk,3}];
				perc_prop_name=['mode_' ACGc{kk,3} '_percent'];
				ix=find(ACGc{kk,2}.Numbers==ACGc_stats(kk,1),1);
				T.(mode_prop_name){ii,1}=ACGc{kk,2}.Categories{ix};
				total_nodes=sum(ACGc{kk,2}.Counts);				
				T.(perc_prop_name)(ii,1)=double((ACGc{kk,2}.Counts(ix)/total_nodes)*100);

				if pc
					ACG_T=ACGc{kk,2};
					total_nodes=sum(ACG_T.Counts);
					for ll=1:numel(ACG_T.Categories)
						cat_name=ACG_T.Categories{ll};
						cat_name=matlab.lang.makeValidName([ACGc{kk,3} '_perc_' cat_name]);
						T.(cat_name)(ii,1)=double((ACG_T.Counts(ll)/total_nodes)*100);
					end
				end

				if ~isempty(mbc)
					warn_flag=false;
					% 解析分类统计参数
					cg=mbc(1);
					dg=mbc(2:end);
					num_dg=numel(dg);
					% 定位目标分类栅格
					cix=find(strcmp(ACGc(:,3),cg));
					ACG=ACGc{cix,1}; % 栅格对象
					ACG_T=ACGc{cix,2}; % 分类编码表
					% 遍历分类类别
					for ll=1:numel(ACG_T.Categories)
						IDX=GRIDobj(ACG,'logical');
						IDX.Z=ismember(ACG.Z,ACG_T.Numbers(ll));
						cat_name=ACG_T.Categories{ll};
						for mm=1:num_dg
							dgOI=dg{mm};
							if strcmp(dgOI,'ksn')
								load(FileName,'MSNc');
								KSNG=GRIDobj(ACG);
								KSNG.Z(:,:)=NaN;
								for oo=1:numel(MSNc)
									ix=coord2ind(ACG,MSNc(oo).X,MSNc(oo).Y);
									KSNG.Z(ix)=MSNc(oo).ksn;
								end
								cat_nameN=matlab.lang.makeValidName(['mksn_' cat_name]);
								T.(cat_nameN)(ii,1)=mean(KSNG.Z(IDX.Z),'omitnan');
							elseif strcmp(dgOI,'gradient')
								load(FileName,'Goc');
								cat_nameN=matlab.lang.makeValidName(['mgrad_' cat_name]);
								T.(cat_nameN)(ii,1)=mean(Goc.Z(IDX.Z),'omitnan');
							elseif regexp(dgOI,regexptranslate('wildcard','rlf*'))
								rlfval=str2num(strrep(dgOI,'rlf',''));
								rlfix=find(cell2mat(rlf(:,2))==rlfval);
								if ~isempty(rlfix)
									Rg=rlf{rlfix,1};
									cat_nameN=matlab.lang.makeValidName(['mr' num2str(rlfval) '_' cat_name]);
									T.(cat_nameN)(ii,1)=mean(Rg.Z(IDX.Z),'omitnan');	
								end								
							else 
								try
									dgix=find(strcmp(AGc(:,2),dgOI));
									AGcOI=AGc{dgix,1};
									cat_nameN=matlab.lang.makeValidName(['m' AGc{dgix,2} '_' cat_name]);
									T.(cat_nameN)(ii,1)=mean(AGcOI.Z(IDX.Z),'omitnan');
								catch
									warn_flag=true;
								end
							end
						end
					end
				end
			end
		end	

		T.file_path{ii,1}=FileName;

		waitbar(ii/num_files);
	end
	warning on

	% 计算沿方位角的距离
	if ~isempty(az)
		% 以数据集中心为原点旋转
		x0=mean(T.center_x);
		y0=mean(T.center_y);
		% 转换方位角
		azn=az-90;
		% 执行坐标旋转
		d=(T.center_x-x0).*cosd(azn)-(T.center_y-y0).*sind(azn);
		% 归一化距离
		d=d-min(d);	
		% 添加至表格
		az_name=['dist_along_' num2str(round(az))];
		T=addvars(T,d,'NewVariableNames',az_name,'After','hyp_integral');
	end

	% 处理分类统计警告
	if ~isempty(mbc)
		if warn_flag==true
			if isdeployed
				warndlg('"means_by_category"中存在无法识别的栅格名称，相关条目未加入表格')
			end
			warning('"means_by_category"中存在无法识别的栅格名称，相关条目未加入表格')
		end
	end

	close(w1);
end

%% 辅助函数

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
		% 如果流域面积过小，则覆盖选择，因为KSN_Trib方法在小流域会失败
		if drainage_area>2.5
			[MSNc]=KSN_Trib(DEMoc,DEMcc,FDc,Ac,Sc,new_ref_concavity,segment_length);
		else
			[MSNc]=KSN_Quick(DEMoc,DEMcc,Ac,Sc,new_ref_concavity,segment_length);
		end
	end

	% 计算流域范围的ksn统计量
	mean_ksn=mean([MSNc.ksn],'omitnan');
	std_ksn=std([MSNc.ksn],'omitnan');
	se_ksn=std_ksn/sqrt(numel(MSNc)); % 标准误
end

function [ksn_ms]=KSN_Quick(DEM,DEMc,A,S,theta_ref,segment_length)
	g=gradient(S,DEMc); % 计算坡度
	G=GRIDobj(DEM);
	G.Z(S.IXgrid)=g;

	Z_RES=DEMc-DEM; % 计算填挖量

	ksn=G./(A.*(A.cellsize^2)).^(-theta_ref); % 计算ksn

	SD=GRIDobj(DEM);
	SD.Z(S.IXgrid)=S.distance; % 存储河道距离
	
	% 将河道对象转换为地理结构体
	ksn_ms=STREAMobj2mapstruct(S,'seglength',segment_length,'attributes',...
		{'ksn' ksn @mean 'uparea' (A.*(A.cellsize^2)) @mean 'gradient' G @mean 'cut_fill' Z_RES @mean...
		'min_dist' SD @min 'max_dist' SD @max});

	% 计算并添加河段距离字段
	seg_dist=[ksn_ms.max_dist]-[ksn_ms.min_dist];
	distcell=num2cell(seg_dist');
	[ksn_ms(1:end).seg_dist]=distcell{:};
	ksn_ms=rmfield(ksn_ms,{'min_dist','max_dist'}); % 移除临时字段
end

function [ksn_ms]=KSN_Trunk(DEM,DEMc,A,S,theta_ref,segment_length,min_order)

	order_exp=['>=' num2str(min_order)]; % 设置最小河道等级阈值

    Smax=modify(S,'streamorder',order_exp); % 提取主干河道
	Smin=modify(S,'rmnodes',Smax); % 获取支流河道

	g=gradient(S,DEMc); % 坡度计算
	G=GRIDobj(DEM);
	G.Z(S.IXgrid)=g;

	Z_RES=DEMc-DEM; % 填挖量计算

	ksn=G./(A.*(A.cellsize^2)).^(-theta_ref); % ksn计算

	% 为不同等级河道创建距离网格
	SDmax=GRIDobj(DEM);
	SDmin=GRIDobj(DEM);
	SDmax.Z(Smax.IXgrid)=Smax.distance;
	SDmin.Z(Smin.IXgrid)=Smin.distance;

	% 分别处理支流和主干河道
	ksn_ms_min=STREAMobj2mapstruct(Smin,'seglength',segment_length,'attributes',...
		{'ksn' ksn @mean 'uparea' (A.*(A.cellsize^2)) @mean 'gradient' G @mean 'cut_fill' Z_RES @mean...
		'min_dist' SDmin @min 'max_dist' SDmin @max});

	ksn_ms_max=STREAMobj2mapstruct(Smax,'seglength',segment_length,'attributes',...
		{'ksn' ksn @mean 'uparea' (A.*(A.cellsize^2)) @mean 'gradient' G @mean 'cut_fill' Z_RES @mean...
		'min_dist' SDmax @min 'max_dist' SDmax @max});

	% 合并结果并处理字段
	ksn_ms=vertcat(ksn_ms_min,ksn_ms_max);
	seg_dist=[ksn_ms.max_dist]-[ksn_ms.min_dist];
	distcell=num2cell(seg_dist');
	[ksn_ms(1:end).seg_dist]=distcell{:};
	ksn_ms=rmfield(ksn_ms,{'min_dist','max_dist'});
end

function [ksn_ms]=KSN_Trib(DEM,DEMc,FD,A,S,theta_ref,segment_length)

	% 定义非交叉河段
	[as]=networksegment_slim(DEM,FD,S);
	seg_bnd_ix=as.ix;
	% 预计算或提取后续需要的值
	z=getnal(S,DEMc);
	zu=getnal(S,DEM);
	z_res=z-zu;
	g=gradient(S,DEMc);
	c=chitransform(S,A,'a0',1,'mn',theta_ref);
	d=S.distance;
	da=getnal(S,A.*(A.cellsize^2));
	ixgrid=S.IXgrid;
	% 提取有序的河道节点列表并查找河道间的断点
	s_node_list=S.orderednanlist;
	streams_ix=find(isnan(s_node_list));
	streams_ix=vertcat(1,streams_ix);
	% 创建空节点属性列表存储ksn值
	ksn_nal=zeros(size(d));
	% 开始主循环处理各河道
	num_streams=numel(streams_ix)-1;
	seg_count=1;
	for ii=1:num_streams
		% 提取当前河道的节点列表
		if ii==1
			snlOI=s_node_list(streams_ix(ii):streams_ix(ii+1)-1);
		else
			snlOI=s_node_list(streams_ix(ii)+1:streams_ix(ii+1)-1);
		end

		% 确定当前河道包含的河段
		[~,~,dn]=intersect(snlOI,seg_bnd_ix(:,1));
		[~,~,up]=intersect(snlOI,seg_bnd_ix(:,2));
		seg_ix=intersect(up,dn);

		num_segs=numel(seg_ix);
		dn_up=seg_bnd_ix(seg_ix,:);
		for jj=1:num_segs
			% 在节点列表中定位起止点
			dnix=find(snlOI==dn_up(jj,1));
			upix=find(snlOI==dn_up(jj,2));
			% 提取目标河段的节点索引
			seg_ix_oi=snlOI(upix:dnix);
			% 提取流动距离并归一化
			dOI=d(seg_ix_oi);
			dnOI=dOI-min(dOI);
			num_bins=ceil(max(dnOI)/segment_length);
			bin_edges=[0:segment_length:num_bins*segment_length];
			% 分箱处理
			for kk=1:num_bins
				idx=dnOI>bin_edges(kk) & dnOI<=bin_edges(kk+1);
				bin_ix=seg_ix_oi(idx);
				cOI=c(bin_ix);
				zOI=z(bin_ix);
					if numel(cOI)>2
						[ksn_val,r2]=Chi_Z_Spline(cOI,zOI);
						ksn_nal(bin_ix)=ksn_val;

						% 构建地理结构体
						ksn_ms(seg_count).Geometry='Line';
						ksm_ms(seg_count).BoundingBox=[min(S.x(bin_ix)),min(S.y(bin_ix));max(S.x(bin_ix)),max(S.y(bin_ix))];
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
	% 精简版网络分段函数，移除了零长度和单节点长度的河段

	%% 识别河道关键节点
	Vhead = streampoi(S,'channelheads','logical');  ihead=find(Vhead==1);  IXhead=S.IXgrid(ihead); % 河源点
	Vconf = streampoi(S,'confluences','logical');   iconf=find(Vconf==1);  IXconf=S.IXgrid(iconf); % 汇合点
	Vout = streampoi(S,'outlets','logical');        iout=find(Vout==1);    IXout=S.IXgrid(iout);   % 出口点
	Vbconf = streampoi(S,'bconfluences','logical'); ibconf=find(Vbconf==1);IXbconf=S.IXgrid(ibconf); % 支流汇入点

	%% 识别关联流域
	DB   = drainagebasins(FD,vertcat(IXbconf,IXout));DBhead=DB.Z(IXhead); DBbconf=DB.Z(IXbconf); DBconf=DB.Z(IXconf); DBout=DB.Z(IXout);

	%% 计算流动距离
	D = flowdistance(FD);

	%% 识别河道段
	% 河源点到支流汇入点的连接
	[~,ind11,ind12]=intersect(DBbconf,DBhead);
	% 汇合点到支流汇入点的连接
	[~,ind21,ind22]=intersect(DBbconf,DBconf);
	% 河源点到出口点的连接
	[~,ind31,ind32]=intersect(DBout,DBhead);
	% 汇合点到出口点的连接
	[~,ind41,ind42]=intersect(DBout,DBconf);
	% 组合连接段
	IX(:,1) = [ IXbconf(ind11)' IXbconf(ind21)' IXout(ind31)'  IXout(ind41)'  ];   ix(:,1)= [ ibconf(ind11)' ibconf(ind21)' iout(ind31)'  iout(ind41)'  ];
	IX(:,2) = [ IXhead(ind12)'  IXconf(ind22)'  IXhead(ind32)' IXconf(ind42)' ];   ix(:,2)= [ ihead(ind12)'  iconf(ind22)'  ihead(ind32)' iconf(ind42)' ];

	% 计算河段长度
	flength=double(abs(D.Z(IX(:,1))-D.Z(IX(:,2))));

	% 移除无效短河段
	idx=flength>=2*DEM.cellsize;
	seg.IX=IX(idx,:);
	seg.ix=ix(idx,:);
	seg.flength=flength(idx);

	% 河段数量统计
	seg.n=numel(IX(:,1));
end

function [KSN,R2] = Chi_Z_Spline(c,z)

	% 使用三次样条插值重采样χ-高程关系
	[~,minIX]=min(c);
	zb=z(minIX);
	chiF=c-min(c);
	zabsF=z-min(z);
	chiS=linspace(0,max(chiF),numel(chiF)).';
	zS=spline(chiF,zabsF,chiS);

	% 通过斜率计算ksn
	KSN= chiS\(zS); % 因a0固定为1，无需mn参数

	% 计算决定系数R²
	z_pred=chiF.*KSN;
	sstot=sum((zabsF-mean(zabsF)).^2);
	ssres=sum((zabsF-z_pred).^2);
	R2=1-(ssres/sstot);

end