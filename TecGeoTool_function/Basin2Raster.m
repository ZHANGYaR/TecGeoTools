function [OUT]=Basin2Raster(DEM,valueOI,location_of_data_files,varargin)
	%
% 使用方法:
%   [RASTER]=Basin2Raster(DEM,valueOI,location_of_data_files);
%   [RASTER]=Basin2Raster(DEM,valueOI,location_of_data_files,'name',value,...);
%
% 描述:
%   该函数将 'ProcessRiverBasins' 函数的输出作为输入，并生成一个包含各个流域的单个GRIDobj
% （由 'ProcessRiverBasins' 和 'SubDivideBigBasins' 选择）并赋予不同的值
%
% 必要输入:
%   DEM - 包含数据集的完整范围的GRIDobj
%   valueOI - 要赋给流域的值，接受的输入包括：
%     'ksn' - 流域的平均ksn值
%     'gradient' - 流域的平均坡度
%     'elevation' - 流域的平均海拔
%     'relief' - 流域的平均地形起伏（必须使用 'relief_radius' 参数指定兴趣半径）
%     'chir2' - chi-z 拟合的R^2值（用于衡量不平衡）
%     'drainage_area' - 流域的排水面积（以平方千米为单位）
%     'hypsometric_integral' - 流域的高程积分
%     'id' - 流域ID编号（即RiverMouth输出的第三列）
%     'theta' - topo toolbox chiplot函数产生的最佳拟合凹度
%     'NAME' - 额外网格的名称（即 'add_grid' 的第二列或 'add_cat_grid' 的第三列），value输入将为附加网格名称的均值或附加分类网格名称的众数
%   location_of_data_files - 包含来自 'ProcessRiverBasins' 的.mat文件的文件夹的完整路径（作为字符串）
%
% 可选输入:
%   file_name_prefix ['basins'] - 输出文件的前缀，将自动附加输出类型，例如 'ksn'、'elevation' 等
%   location_of_subbasins ['SubBasins'] - 包含感兴趣子流域的文件夹名称（如果使用 "SubDivideBigBasins" 创建了子流域），应位于与 "location_of_data_files" 提供的主流域文件夹中
%   method ['subdivided'] - 用于细分流域的方法。如果使用了 'ProcessRiverBasins' 然后使用了 'SubDivideBigBasins'，或者仅使用了 'ProcessRiverBasins' 但没有选择任何嵌套流域（即 'ProcessRiverBasins' 提供的河口不在其他流域的流域边界内），则应使用 'subdivided'（这是默认值，因此无需为此属性指定值）。如果手动选择了嵌套流域并运行了 'ProcessRiverBasins'，则应使用 'nested'。
%   relief_radius [2500] - 如果 'valueOI' 设置为 'relief'，则使用的地形起伏半径
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数编写者: Adam M. Forte - 更新日期: 2019年1月27日
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'Basin2Raster';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'valueOI',@(x) ischar(x));
	addRequired(p,'location_of_data_files',@(x) isdir(x));

	addParameter(p,'location_of_subbasins','SubBasins',@(x) ischar(x) || isempty(x));
	addParameter(p,'file_name_prefix','basins',@(x) ischar(x));
	addParameter(p,'method','subdivided',@(x) ischar(validatestring(x,{'subdivided','nested'})));
	addParameter(p,'relief_radius',2500,@(x) isscalar(x) && isnumeric(x));

	parse(p,DEM,valueOI,location_of_data_files,varargin{:});
	DEM=p.Results.DEM;
	valueOI=p.Results.valueOI;
	location_of_data_files=p.Results.location_of_data_files;

	location_of_subbasins=p.Results.location_of_subbasins;
	file_name_prefix=p.Results.file_name_prefix;
	method=p.Results.method;
	rr=p.Results.relief_radius;

	OUT=GRIDobj(DEM);
	OUT=OUT-9999;

	[head_dir,~,~]=fileparts(location_of_data_files);
	if isempty(head_dir)
		file_name_prefix=[pwd filesep file_name_prefix];
	else
		file_name_prefix=[head_dir filesep file_name_prefix];
	end



	switch method
	case 'subdivided'
		%% Build File List
		% Get Basin Numbers
		AllFullFiles=dir([location_of_data_files filesep '*_Data.mat']);
		num_basins=numel(AllFullFiles);
		basin_nums=zeros(num_basins,1);
		for jj=1:num_basins
			FileName=AllFullFiles(jj,1).name;
			basin_nums(jj)=sscanf(FileName,'%*6s %i'); %%%
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

		w1=waitbar(0,'正在构建栅格...');
		for ii=1:num_files
			FileName=[FileList(ii,1).folder filesep FileList(ii,1).name];
			switch valueOI
			case 'ksn'
				load(FileName,'DEMoc','KSNc_stats');
				val=KSNc_stats(:,1);
			case 'gradient'
				load(FileName,'DEMoc','Gc_stats');
				val=Gc_stats(:,1);
			case 'elevation'
				load(FileName,'DEMoc','Zc_stats');
				val=Zc_stats(:,1);
			case 'chir2'
				load(FileName,'DEMoc','Sc','Ac','theta_ref');
				c=chiplot(Sc,DEMoc,Ac,'a0',1,'mn',theta_ref,'plot',false);
				val=c.R2;
			case 'drainage_area'
				load(FileName,'DEMoc','drainage_area');
				val=drainage_area;
			case 'hypsometric_integral'
				load(FileName,'DEMoc','hyps');
				val=abs(trapz((hyps(:,2)-min(hyps(:,2)))/(max(hyps(:,2))-min(hyps(:,2))),hyps(:,1)/100));
			case 'id'
				load(FileName,'DEMoc','RiverMouth');
				val=RiverMouth(:,3);
			case 'theta'
				load(FileName,'DEMoc','Chic');
				val=Chic.mn;
			case 'relief'
				VarList=whos('-file',FileName);
				RLFInd=find(strcmp(cellstr(char(VarList.name)),'rlf'));
				if isempty(RLFInd)
					if isdeployed
						errordlg('确实已为这些流域计算了起伏度。')
					end
					error('确实已经为这些流域计算了高程起伏。')
				end
				load(FileName,'DEMoc','rlf','rlf_stats');
				ix = find(cell2mat(rlf(:,2)) == rr);
				if isempty(ix)
					if isdeployed
						errordlg('输入的高程起伏半径未在高程起伏输出中找到，请检查以确保高程起伏半径正确无误。')
					end
					error('输入的高程起伏半径未在高程起伏输出中找到，请检查以确保高程起伏半径正确无误。')
				end
				val=rlf_stats(ix,1);
			otherwise
				VarList=whos('-file',FileName);
				AGInd=find(strcmp(cellstr(char(VarList.name)),'AGc'));
				ACGInd=find(strcmp(cellstr(char(VarList.name)),'ACGc'));
				if ~isempty(AGInd)
					load(FileName,'DEMoc','AGc');
					AGInd=true;
				end

				if ~isempty(ACGInd)
					load(FileName,'DEMoc','ACGc');
					ACGInd=true;
				end

				if AGInd & any(strcmp(AGc(:,2),valueOI))
					Nix=find(strcmp(AGc(:,2),valueOI));
					load(FileName,'AGc_stats');
					val=AGc_stats(Nix,1);
				elseif ACGInd & any(strcmp(ACGc(:,3),valueOI))
					Nix=find(trcmp(ACGc(:,3),valueOI));
					load(FileName,'ACGc_stats');
					val=ACGc_stats(Nix,1);
				else
					if isdeployed
						errordlg('提供给“valueOI”的名称与其他网格或附加分类网格中的名称不匹配。')
					end
					error('提供给“valueOI”的名称与附加网格或附加分类网格中的任何名称均不匹配。');
				end
			end

			I=~isnan(DEMoc.Z);
			[X,Y]=getcoordinates(DEMoc);
			xmat=repmat(X,numel(Y),1);
			ymat=repmat(Y,1,numel(X));

			xix=xmat(I);
			yix=ymat(I);

			ix=coord2ind(DEM,xix,yix);
			val_list=ones(numel(ix),1).*val;
			OUT.Z(ix)=val_list;
			waitbar(ii/num_files);
		end
		close(w1);

	case 'nested'
		% Build list of indices
		AllFiles=dir([location_of_data_files filesep '*_Data.mat']);
		num_basins=numel(AllFiles);

		ix_cell=cell(num_basins,1);
		basin_list=zeros(num_basins,1);
		FileCell=cell(num_basins,1);
		for jj=1:num_basins
			FileName=AllFiles(jj,1).name;
			FileCell{jj}=FileName;

			load(FileName,'DEMoc');
			[x,y]=getcoordinates(DEMoc);
			xg=repmat(x,numel(y),1);
			yg=repmat(y,1,numel(x));
			xl=xg(~isnan(DEMoc.Z));
			yl=yg(~isnan(DEMoc.Z));
			ix=coord2ind(DEM,xl,yl);

			ix_cell{jj}=ix;

			basin_list(jj,1)=numel(ix);
		end

		% Sort basin size list in descending order
		[~,six]=sort(basin_list,'descend');
		% Apply sorting index to FileCell and ix_cell
		FileCell=FileCell(six);
		ix_cell=ix_cell(six);

		w1=waitbar(0,'正在构建栅格...');
		for ii=1:num_basins
			FileName=FileCell{ii};

			switch valueOI
			case 'ksn'
				load(FileName,'DEMoc','KSNc_stats');
				val=KSNc_stats(:,1);
			case 'gradient'
				load(FileName,'DEMoc','Gc_stats');
				val=Gc_stats(:,1);
			case 'elevation'
				load(FileName,'DEMoc','Zc_stats');
				val=Zc_stats(:,1);
			case 'chir2'
				load(FileName,'DEMoc','Chic');
				val=Chic.R2;
			case 'drainage_area'
				load(FileName,'DEMoc','drainage_area');
				val=drainage_area;
			case 'hypsometric_integral'
				load(FileName,'DEMoc','hyps');
				val=abs(trapz((hyps(:,2)-min(hyps(:,2)))/(max(hyps(:,2))-min(hyps(:,2))),hyps(:,1)/100));
			case 'id'
				load(FileName,'DEMoc','RiverMouth');
				val=RiverMouth(:,3);
			case 'theta'
				load(FileName,'DEMoc','Chic');
				val=Chic.mn;
			case 'relief'
				VarList=whos('-file',FileName);
				RLFInd=find(strcmp(cellstr(char(VarList.name)),'rlf'));
				if isempty(RLFInd)
					if isdeployed
						errordlg('这些流域似乎已计算了高程起伏。')
					end
					error('这些流域似乎已计算了高程起伏。')
				end
				load(FileName,'DEMoc','rlf','rlf_stats');
				ix=find(rlf(:,2)==rr);
				if isempty(ix)
					if isdeployed
						errordlg('输入的高程起伏半径未在高程起伏输出中找到，请检查确保高程起伏半径正确。')
					end
					error('输入的高程起伏半径未在高程起伏输出中找到，请检查确保高程起伏半径正确。')
				end
				val=rlf_stats(ix,1);
			otherwise
				VarList=whos('-file',FileName);
				AGInd=find(strcmp(cellstr(char(VarList.name)),'AGc'));
				ACGInd=find(strcmp(cellstr(char(VarList.name)),'ACGc'));
				if ~isempty(AGInd)
					load(FileName,'DEMoc','AGc');
					AGInd=true;
				end

				if ~isempty(ACGInd)
					load(FileName,'DEMoc','ACGc');
					ACGInd=true;
				end

				if AGInd & any(strcmp(AGc(:,2),valueOI))
					Nix=find(strcmp(AGc(:,2),valueOI));
					load(FileName,'AGc_stats');
					val=AGc_stats(Nix,1);
				elseif ACGInd & any(strcmp(ACGc(:,3),valueOI))
					Nix=find(trcmp(ACGc(:,3),valueOI));
					load(FileName,'ACGc_stats');
					val=ACGc_stats(Nix,1);
				else
					if isdeployed
						errordlg('提供给“valueOI”的名称与附加栅格或附加分类网格中的名称不匹配。')
					end
					error('提供给“valueOI”的名称与附加栅格或附加分类网格中的名称不匹配。');
				end
			end

			I=~isnan(DEMoc.Z);
			[X,Y]=getcoordinates(DEMoc);
			xmat=repmat(X,numel(Y),1);
			ymat=repmat(Y,1,numel(X));

			ix=ix_cell{ii};

			val_list=ones(numel(ix),1).*val;
			OUT.Z(ix)=val_list;
			waitbar(ii/num_files);
		end
		close(w1);
	end

	switch valueOI
	case 'ksn'
		out_file=[file_name_prefix '_ksn.txt'];
		GRIDobj2ascii(OUT,out_file);
	case 'gradient'
		out_file=[file_name_prefix '_gradient.txt'];
		GRIDobj2ascii(OUT,out_file);
	case 'elevation'
		out_file=[file_name_prefix '_elevation.txt'];
		GRIDobj2ascii(OUT,out_file);
	case 'chir2'
		out_file=[file_name_prefix '_chiR2.txt'];
		GRIDobj2ascii(OUT,out_file);
	case 'drainage_area'
		out_file=[file_name_prefix '_drainarea.txt'];
		GRIDobj2ascii(OUT,out_file);
	case 'hypsometric_integral'
		out_file=[file_name_prefix '_hyps_int.txt'];
		GRIDobj2ascii(OUT,out_file);
	case 'id'
		out_file=[file_name_prefix '_id.txt'];
		GRIDobj2ascii(OUT,out_file);
	case 'theta'
		out_file=[file_name_prefix '_theta.txt'];
		GRIDobj2ascii(OUT,out_file);
	case 'relief'
		out_file=[file_name_prefix '_relief.txt'];
		GRIDobj2ascii(OUT,out_file);		
	otherwise
		out_file=[file_name_prefix '_' valueOI '.txt'];
		GRIDobj2ascii(OUT,out_file);
	end

