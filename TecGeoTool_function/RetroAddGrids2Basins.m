function RetroAddGrids2Basins(location_of_data_files,DEM,behavior,varargin)
	%
	% 用法:
	%		RetroAddGrids2Basins(数据文件位置,DEM,操作模式,'参数名',参数值,...);
	%
	% 描述:
	%		本函数对'ProcessRiverBasins'和'SubdivideBigBasins'生成的流域/子流域数据集进行操作，
	%		为其添加附加栅格或分类栅格数据。可控制覆盖或追加已有数据。
	%		
	% 必需输入:
	% 		数据文件位置 - 包含流域mat文件的完整路径
	%		DEM - 原始DEM的GRIDobj对象
	%		操作模式 - 数据处理方式：
	%			'overwrite'覆盖现有数据 | 'append'追加新数据
	%
	% 可选参数:
	%		add_grids [] - 附加栅格数据（格式同ProcessRiverBasins的add_grids参数）
	%		add_cat_grids [] - 附加分类栅格（格式同ProcessRiverBasins的add_cat_grids参数）
	%		include ['all'] - 处理范围选择：
	%			'all'所有流域 | 'subdivided'仅子流域 | 'bigonly'仅主流域
	%		location_of_subbasins ['SubBasins'] - 子流域存储目录名
	%		resample_method ['nearest'] - 重采样方法：
	%			'nearest'最近邻 | 'bilinear'双线性 | 'bicubic'双三次
	% 
	% 示例:
	% 		RetroAddGrids2Basins('basins',DEM,'append','add_grids',{PRE,'precip'});
	% 		RetroAddGrids2Basins('basins',DEM,'append','add_grids',{PRE,'precip'},'include','subdivided','location_of_subbasins','my_subs');	
	%	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% 作者：Adam M. Forte - 最后更新：2020/05/16 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% 解析输入参数
	p = inputParser;
	p.FunctionName = 'RetroAddGrids2Basins';
	addRequired(p,'location_of_data_files',@(x) isdir(x));
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'behavior',@(x) ischar(validatestring(x,{'overwrite','append'})));

	addParameter(p,'location_of_subbasins','SubBasins',@(x) ischar(x) || isempty(x));
	addParameter(p,'include','all',@(x) ischar(validatestring(x,{'all','subdivided','bigonly'})));
	addParameter(p,'add_grids',[],@(x) isa(x,'cell') && size(x,2)==2 || isempty(x));
	addParameter(p,'add_cat_grids',[],@(x) isa(x,'cell') && size(x,2)==3 || isempty(x));
	addParameter(p,'resample_method','nearest',@(x) ischar(validatestring(x,{'nearest','bilinear','bicubic'})));

	parse(p,location_of_data_files,DEM,behavior,varargin{:});
	location_of_data_files=p.Results.location_of_data_files;
	DEM=p.Results.DEM;
	behavior=p.Results.behavior;

	location_of_subbasins=p.Results.location_of_subbasins;
	include=p.Results.include;
	AG=p.Results.add_grids;
	ACG=p.Results.add_cat_grids;
	resample_method=p.Results.resample_method;

	% 校验必要输入
	if isempty(AG) & isempty(ACG)
		error('必须提供"add_grids"或"add_cat_grids"参数，否则本函数无执行意义')
	end

	% 确定处理范围
	switch include
	case 'all'
		% 获取主流域和子流域文件列表
		FileList1=dir([location_of_data_files filesep '*_Data.mat']);
		FileList2=dir([location_of_subbasins filesep '*_DataSubset*.mat']);
		FileList=vertcat(FileList1,FileList2);
		num_files=numel(FileList);
	case 'bigonly'
		% 仅处理主流域文件
		FileList=dir([location_of_data_files filesep '*_Data.mat']);
		num_files=numel(FileList);
	case 'subdivided'
		% 智能识别子流域文件
		AllFullFiles=dir([location_of_data_files filesep '*_Data.mat']);
		num_basins=numel(AllFullFiles);
		basin_nums=zeros(num_basins,1);
		for jj=1:num_basins
			fileName=AllFullFiles(jj,1).name;
			basin_nums(jj)=sscanf(fileName,'%*6s %i'); 
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

	% 预处理附加栅格数据
	if ~isempty(AG)
		num_grids=size(AG,1);
		for jj=1:num_grids
			AGoi=AG{jj,1};
			if ~validatealignment(AGoi,DEM);
				disp(['正在将 ' AG{jj,2} ' 栅格重采样至DEM分辨率（方法：' resample_method '）']);
				AG{jj,1}=resample(AGoi,DEM,resample_method);
			end
		end
	end

	% 预处理分类栅格数据（强制使用最近邻采样）
	if ~isempty(ACG)
		num_grids=size(ACG,1);
		for jj=1:num_grids
			ACGoi=ACG{jj,1};
			if ~validatealignment(ACGoi,DEM);
				disp(['正在将 ' ACG{jj,3} ' 分类栅格重采样至DEM分辨率（强制最近邻方法）']);
				ACG{jj,1}=resample(ACGoi,DEM,'nearest');
			end
		end
	end

	% 主处理循环
	w1=waitbar(0,'正在向流域添加栅格...');
	for ii=1:num_files
		% 加载流域DEM
		FileName=[FileList(ii,1).folder filesep FileList(ii,1).name];
		load(FileName,'DEMoc');
		% 构建掩膜
		[x,y]=getcoordinates(DEMoc);
		[X,Y]=meshgrid(x,y);
		ix=isnan(DEMoc.Z(:));
		xx=X(ix); yy=Y(ix);
		gix=coord2ind(DEM,xx,yy);
		I=GRIDobj(DEM,'logical');
		I.Z(gix)=true;

		% 处理附加栅格
		if ~isempty(AG)
			switch behavior
			case 'overwrite'
				% 覆盖模式：完全替换原有数据
				num_grids=size(AG,1);
				AGc=cell(size(AG));
				for jj=1:num_grids
					AGcOI=crop(AG{jj,1},I,nan);
					AGc{jj,1}=AGcOI;
					AGc{jj,2}=AG{jj,2};
					% 计算统计指标
					mean_AGc=mean(AGcOI.Z(:),'omitnan');
					min_AGc=min(AGcOI.Z(:),[],'omitnan');
					max_AGc=max(AGcOI.Z(:),[],'omitnan');
					std_AGc=std(AGcOI.Z(:),'omitnan');
					se_AGc=std_AGc/sqrt(sum(~isnan(AGcOI.Z(:))));
					AGc_stats(jj,:)=[mean_AGc se_AGc std_AGc min_AGc max_AGc];
				end
				save(FileName,'AGc','AGc_stats','-append');
			case 'append'
				% 追加模式：保留原有数据
				vI=who('-file',FileName);
				if any(ismember(vI,'AGc'))
					load(FileName,'AGc','AGc_stats');
					num_existing=size(AGc,1);
					pos_vec=num_existing+1:num_grids+num_existing;
					num_grids=size(AG,1);
					for jj=1:num_grids
						AGcOI=crop(AG{jj,1},I,nan);
						AGc{pos_vec(jj),1}=AGcOI;
						AGc{pos_vec(jj),2}=AG{jj,2};
						mean_AGc=mean(AGcOI.Z(:),'omitnan');
						min_AGc=min(AGcOI.Z(:),[],'omitnan');
						max_AGc=max(AGcOI.Z(:),[],'omitnan');
						std_AGc=std(AGcOI.Z(:),'omitnan');
						se_AGc=std_AGc/sqrt(sum(~isnan(AGcOI.Z(:))));
						AGc_stats(pos_vec(jj),:)=[mean_AGc se_AGc std_AGc min_AGc max_AGc];
					end
					save(FileName,'AGc','AGc_stats','-append');						

				else
					% 首次添加情况
					num_grids=size(AG,1);
					AGc=cell(size(AG));
					for jj=1:num_grids
						AGcOI=crop(AG{jj,1},I,nan);
						AGc{jj,1}=AGcOI;
						AGc{jj,2}=AG{jj,2};
						mean_AGc=mean(AGcOI.Z(:),'omitnan');
						min_AGc=min(AGcOI.Z(:),[],'omitnan');
						max_AGc=max(AGcOI.Z(:),[],'omitnan');
						std_AGc=std(AGcOI.Z(:),'omitnan');
						se_AGc=std_AGc/sqrt(sum(~isnan(AGcOI.Z(:))));
						AGc_stats(jj,:)=[mean_AGc se_AGc std_AGc min_AGc max_AGc];
					end
					save(FileName,'AGc','AGc_stats','-append');					
				end
			end				
		end

		% 处理分类栅格
		if ~isempty(ACG)
			switch behavior
			case 'overwrite'
				% 覆盖模式
				num_grids=size(ACG,1);
				ACGc=cell(size(ACG));
				for jj=1:num_grids
					ACGcOI=crop(ACG{jj,1},I,nan);
					ACGc{jj,1}=ACGcOI;
					ACGc{jj,3}=ACG{jj,3};
					% 构建分类直方图
					edg=ACG{jj,2}.Numbers;
					edg=edg+0.5;
					edg=vertcat(0.5,edg);
					[N,~]=histcounts(ACGcOI.Z(:),edg);
					T=ACG{jj,2};
					T.Counts=N';
					ACGc{jj,2}=T;
					ACGc_stats(jj,1)=[mode(ACGcOI.Z(:))];
				end
				save(FileName,'ACGc','ACGc_stats','-append');
			case 'append'
				% 追加模式
				vI=who('-file',FileName);
				if any(ismember(vI,'ACGc'))
					num_existing=size(ACGc,1);
					pos_vec=num_existing+1:num_grids+num_existing;
					num_grids=size(ACG,1);
					for jj=1:num_grids
						ACGcOI=crop(ACG{jj,1},I,nan);
						ACGc{pos_vec(jj),1}=ACGcOI;
						ACGc{pos_vec(jj),3}=ACG{jj,3};
						edg=ACG{jj,2}.Numbers;
						edg=edg+0.5;
						edg=vertcat(0.5,edg);
						[N,~]=histcounts(ACGcOI.Z(:),edg);
						T=ACG{jj,2};
						T.Counts=N';
						ACGc{pos_vec(jj),2}=T;
						ACGc_stats(pos_vec(jj),1)=[mode(ACGcOI.Z(:))];
					end
					save(FileName,'ACGc','ACGc_stats','-append');
				else
					% 首次添加分类栅格
					num_grids=size(ACG,1);
					ACGc=cell(size(ACG));
					for jj=1:num_grids
						ACGcOI=crop(ACG{jj,1},I,nan);
						ACGc{jj,1}=ACGcOI;
						ACGc{jj,3}=ACG{jj,3};
						edg=ACG{jj,2}.Numbers;
						edg=edg+0.5;
						edg=vertcat(0.5,edg);
						[N,~]=histcounts(ACGcOI.Z(:),edg);
						T=ACG{jj,2};
						T.Counts=N';
						ACGc{jj,2}=T;
						ACGc_stats(jj,1)=[mode(ACGcOI.Z(:))];
					end
					save(FileName,'ACGc','ACGc_stats','-append');					
				end
			end
		end
		waitbar(ii/num_files);
	end
	close(w1);
end