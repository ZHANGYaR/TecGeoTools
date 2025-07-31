function ProcessRiverBasins(DEM,FD,A,S,river_mouths,basin_dir,varargin)
%
% 使用方法：
%	ProcessRiverBasins(DEM,FD,A,S,river_mouths,basin_dir);
%	ProcessRiverBasins(DEM,FD,A,S,river_mouths,basin_dir,'name',value,...);
%
% 描述：
% 	该函数接收来自MakeStreams脚本的网格对象输出（DEM, FD, A, S），一系列河流口的x,y坐标，
% 	并输出裁剪后的DEM、流域网络、各种地形度量和河流值（ks, ksn, chi）。
%
% 必需输入：
% 		DEM - 区域的数字高程模型（DEM）的GRIDobj，已加载到工作空间中。
% 		FD - 区域的流向（flow direction）的FLOWobj，已加载到工作空间中。
% 		A - 区域的流积累（flow accumulation）的GRID对象，已加载到工作空间中。
% 		S - 区域的河流网络（stream network）的STREAMobj，已加载到工作空间中。
% 		river_mouths - 河流口位置（即汇水点），从这些点上方提取流域。可接受三种形式：
%			1) nx3矩阵，包含河流口的x、y坐标以及标识感兴趣流域的编号（必须与DEM使用相同的投影），
%				从PickBasins输出的矩阵（并保存为Outlets.mat）可以作为RiverMouths输入。
%			2) 单个值，表示一个海拔高度，代码将使用该值自动生成河流口。结合'min_basin_size'使用，以限制提取的流域。
%			3) 指向点Shapefile的完整路径，Shapefile包含一个数字用户输入字段（如默认的'ID'字段，由ArcGIS生成），
%				该字段将作为河流口ID（必须与DEM使用相同的投影）。
%			4) 指向多段线Shapefile的完整路径。代码将假定您希望提取与多段线相交的河流上游的流域。结合'min_basin_size'使用，以限制提取的流域。
% 		basin_dir - 存储流域文件的文件夹位置（如果指定的文件夹不存在，代码将创建它）。
%
% 可选输入：
%		conditioned_DEM [] - 可选提供的平滑DEM，用于此函数（不要将平滑DEM作为主要输入！），
%			该DEM将用于提取海拔。有关制作平滑DEM的选项，请参见'ConditionDEM'函数。如果未提供输入，代码将默认使用mincosthydrocon函数。
%		interp_value [0.1] - 用于mincosthydrocon的插值参数值（介于0和1之间）（如果用户提供了平滑DEM，则不使用）。
% 		threshold_area [1e6] - 用于定义流的最小积累面积（单位：平方米）。
% 		segment_length [1000] - 用于沿ksn平均的平滑距离（单位：米），建议值为1000米。
% 		ref_concavity [0.5] - 用于计算ksn的参考凹度，建议值为0.45。
%		ksn_method [quick] - 切换计算ksn值的方法，选项有'quick'、'trunk'或'trib'，'trib'方法的计算时间比'quick'长3-4倍。
%			在大多数情况下，'quick'方法工作良好，但如果关注支流交汇处的值，'trib'可能更好，因为它为各个渠道段计算ksn值。
%			'trunk'选项则为主要河流计算斜率值（考虑到通过'min_order'指定的流域等级）。
%		min_order [4] - 被视为主要河流的最小流域等级，仅在'ksn_method'设置为'trunk'时使用。
% 		write_arc_files [false] - 设置为true以输出各种网格的ASCII文件和ksn的Shapefile，设置为false则不输出Arc文件。
%		add_grids [] - 可选提供额外网格的单元格数组，这些网格将在选定的河流流域上裁剪。预期输入为nx2单元格数组，
%			其中第一列为GRIDobj，第二列为网格的名称（以便于后续查看输出时记住这些网格的用途，但也用于'Basin2Shape'时的字段名称）。
%			代码将检查任何输入网格是否与输入DEM的尺寸和单元格大小一致，如果不一致，将使用'resample'函数转换输入网格。
%		add_cat_grids [] - 可选提供分类网格（如由'CatPoly2GRIDobj'函数生成的地质图）的单元格数组。预期输入为nx3单元格数组，
%			第一列为GRIDobj，第二列为look_table，第三列为网格名称。假设在预处理这些网格时，使用与主函数相同的DEM GRIDobj。
%		resample_method ['nearest'] - 额外网格的重采样方法（如需要）。可接受的输入为'nearest'、'bilinear'或'bicubic'。
%		gradient_method ['arcslope'] - 用于计算梯度的函数，选项为'arcslope'（默认）或'gradient8'。
%		calc_relief [false] - 计算局部起伏的选项。可以提供'radi'的数组。
%		relief_radii [2500] - 计算局部起伏时使用的半径（单位为地图单位）的1D向量。如果提供多个值，函数将计算所有这些半径的局部起伏。
%		ksn_radius [5000] - 用于平均ksn值的圆形移动区域半径，以制作插值后的ksn网格。如果提供空数组（即[]），则抑制计算（并保存此输出）。
%		min_basin_size [10] - 提取流域的最小面积（单位：km^2），如果'river_mouths'输入为单一值或多段线Shapefile。设置为0以提取满足定义条件的所有流域。
%		precip_AGcol_name [] - 指示'add_grids'提供的降水栅格名称的字符串（即第二列输入）。如果提供有效条目，将使用该降水栅格生成加权流积累栅格。
%
% 注意：
%		- 代码将检查'river_mouths'输入，确认没有重复的ID号（如果有，将丢弃原ID并创建新的ID，并输出包含新ID的河流口位置的文本文件）；
%			同时还会检查提供的河流口是否超出了DEM的边界（如果超出，将删除这些ID）。
%			
%
% 示例：
%		ProcessRiverBasins(DEM,FD,S,RiverMouths);
%		ProcessRiverBasins(DEM,FD,S,RiverMouths,'theta_ref',0.5,'write_arc_files',true);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 由Adam M. Forte编写 - 更新日期：06/18/18 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 解析输入参数
p = inputParser;
p.FunctionName = 'ProcessRiverBasins';
addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
addRequired(p,'A',@(x) isa(x,'GRIDobj'));
addRequired(p,'S',@(x) isa(x,'STREAMobj'));
addRequired(p,'river_mouths',@(x) isnumeric(x) && size(x,2)==3 || isnumeric(x) && isscalar(x) || regexp(x,regexptranslate('wildcard','*.shp')));
addRequired(p,'basin_dir',@(x) ischar(x));

addParameter(p,'ref_concavity',0.5,@(x) isscalar(x) && isnumeric(x));
addParameter(p,'threshold_area',1e6,@(x) isscalar(x) && isnumeric(x));
addParameter(p,'segment_length',1000,@(x) isscalar(x) && isnumeric(x));
addParameter(p,'write_arc_files',false,@(x) isscalar(x));
addParameter(p,'ksn_method','quick',@(x) ischar(validatestring(x,{'quick','trunk','trib'})));
addParameter(p,'min_order',4,@(x) isscalar(x) && isnumeric(x));
addParameter(p,'add_grids',[],@(x) isa(x,'cell') && size(x,2)==2 || isempty(x));
addParameter(p,'add_cat_grids',[],@(x) isa(x,'cell') && size(x,2)==3 || isempty(x));
addParameter(p,'resample_method','nearest',@(x) ischar(validatestring(x,{'nearest','bilinear','bicubic'})));
addParameter(p,'gradient_method','arcslope',@(x) ischar(validatestring(x,{'arcslope','gradient8'})));
addParameter(p,'calc_relief',false,@(x) isscalar(x));
addParameter(p,'relief_radii',[2500],@(x) isnumeric(x) && size(x,2)==1 || size(x,1)==1);
addParameter(p,'conditioned_DEM',[],@(x) isa(x,'GRIDobj') || isempty(x));
addParameter(p,'interp_value',0.1,@(x) isnumeric(x) && x>=0 && x<=1);
addParameter(p,'ksn_radius',5000,@(x) isnumeric(x) && isscalar(x) || isempty(x));
addParameter(p,'min_basin_size',10,@(x) isnumeric(x) && isscalar(x));
addParameter(p,'precip_AGcol_name',[],@(x) ischar(x));
addParameter(p,'out_dir',[],@(x) ischar(x)); % GUI使用的隐藏参数

parse(p,DEM,FD,A,S,river_mouths,basin_dir,varargin{:});
DEM=p.Results.DEM;
FD=p.Results.FD;
A=p.Results.A;
S=p.Results.S;
river_mouths=p.Results.river_mouths;
basin_dir=p.Results.basin_dir;

min_order=p.Results.min_order;
theta_ref=p.Results.ref_concavity;
threshold_area=p.Results.threshold_area;
segment_length=p.Results.segment_length;
write_arc_files=p.Results.write_arc_files;
ksn_method=p.Results.ksn_method;
AG=p.Results.add_grids;
ACG=p.Results.add_cat_grids;
resample_method=p.Results.resample_method;
gradient_method=p.Results.gradient_method;
calc_relief=p.Results.calc_relief;
relief_radii=p.Results.relief_radii;
iv=p.Results.interp_value;
DEMhc=p.Results.conditioned_DEM;
radius=p.Results.ksn_radius;
out_dir=p.Results.out_dir;
mbsz=p.Results.min_basin_size;
prn=p.Results.precip_AGcol_name;

% 设置重做标志
redo_flag=false;

% 处理输出目录
[fp,~,~]=fileparts(basin_dir);
if radius<=0
	radius=[];
end

% 创建流域目录
if isempty(out_dir) && isempty(fp)
	current=pwd;
	bsn_path=fullfile(current,basin_dir);
elseif isempty(out_dir) && ~isempty(fp)
	bsn_path=basin_dir;
else
	current=out_dir;
	bsn_path=fullfile(current,basin_dir);
end
if ~isdir(bsn_path)
	mkdir(bsn_path);
end

% 检查附加网格尺寸并重采样
if ~isempty(AG)
	num_grids=size(AG,1);
	for jj=1:num_grids
		AGoi=AG{jj,1};
		if ~validatealignment(AGoi,DEM);
			disp(['正在将附加网格 ' AG{jj,2} ' 重采样至DEM的分辨率和范围（方法：' resample_method ')']);
			AG{jj,1}=resample(AGoi,DEM,resample_method);
		end
	end
end

% 检查降水网格设置
if ~isempty(prn)
	if isempty(AG)
		warning('未提供附加网格，无法计算ksn-q');
		weight_acc_flag=false;
	else
		AGnames=AG(:,2);
		Pidx=strcmp(AGnames,prn);
		if any(Pidx)
			PRECIP=AG{Pidx,1};
			disp('正在计算加权流积累量用于ksn-q计算');
			WA=flowacc(FD,PRECIP);
			weight_acc_flag=true;
		else
			warning('附加网格中未找到指定降水网格，无法计算ksn-q');
			weight_acc_flag=false;
		end
	end
else
	weight_acc_flag=false;
end

% 检查河段长度参数
if (DEM.cellsize*3)>segment_length
	segment_length=DEM.cellsize*3;
	if isdeployed
		warndlg(['输入的segment_length与DEM分辨率不兼容，已自动调整为' num2str(segment_length)])
	end
	warning(['输入的segment_length与DEM分辨率不兼容，已自动调整为' num2str(segment_length)]);
end

% 处理河流出口输入
if ischar(river_mouths)
	disp('读取Shapefile并捕捉河流出口至河道网络')
	rm_ms=shaperead(river_mouths);
	rm_t=struct2table(rm_ms);
	if strcmp(rm_t.Geometry(1),'Point') % 点要素处理
		fn=rm_t.Properties.VariableNames;
		xi=rm_t.X;
		yi=rm_t.Y;
		riv_nums=rm_t.(fn{4});
		% 检查重复ID
		if numel(riv_nums)~=numel(unique(riv_nums))
			riv_nums=[1:numel(riv_nums)]';
			if isdeployed
				warndlg('检测到重复的河流出口ID，已自动重新编号')
			end
			warning('检测到重复的河流出口ID，已自动重新编号');
			redo_flag=true;
		end	
		% 检查边界范围
		[demx,demy]=getoutline(DEM,true);
		[demix]=inpolygon(xi,yi,demx,demy);	
		xi=xi(demix); yi=yi(demix); riv_nums=riv_nums(demix);
		% 捕捉至河道
		[xn,yn]=snap2stream(S,xi,yi);
		RM=[xn yn riv_nums];
		num_basins=numel(xn);
		if redo_flag
			csvwrite(fullfile(bsn_path,'river_mouths_updated.txt'),RM);
		end
	elseif strcmp(rm_t.Geometry(1),'Line') % 线要素处理
		L=line2GRIDobj(DEM,rm_ms);
		Sup=modify(S,'upstreamto',L);
		outix=streampoi(Sup,'outlets','ix');
		if mbsz>0 % 过滤小流域
			DA=(A.*A.cellsize^2)/(1e6);
			da=DA.Z(outix);
			outix(da<mbsz)=[];
		end
		[rmx rmy]=ind2coord(DEM,outix);
		num_basins=numel(rmx);
		riv_nums=transpose(1:num_basins);
		RM=[rmx rmy riv_nums];
		csvwrite(fullfile(bsn_path,'river_mouths.txt'),RM);
	else
		if isdeployed
			errordlg('输入的Shapefile类型错误，必须为点或线要素')
		end
		error('输入的Shapefile类型错误，必须为点或线要素');
	end
elseif size(river_mouths,2)==3 % 矩阵输入处理
	disp('捕捉河流出口至河道网络')
	xi=river_mouths(:,1);
	yi=river_mouths(:,2);
	riv_nums=river_mouths(:,3);
	% 检查重复ID
	if numel(riv_nums)~=numel(unique(riv_nums))
		riv_nums=[1:numel(riv_nums)]';
		if isdeployed
			warndlg('检测到重复的河流出口ID，已自动重新编号')
		end
		warning('检测到重复的河流出口ID，已自动重新编号');
		redo_flag=true;
	end	
	% 检查边界范围
	[demx,demy]=getoutline(DEM,true);
	[demix]=inpolygon(xi,yi,demx,demy);	
	xi=xi(demix); yi=yi(demix); riv_nums=riv_nums(demix);
	% 捕捉至河道
	[xn,yn]=snap2stream(S,xi,yi);
	RM=[xn yn riv_nums];
	num_basins=numel(xn);
	if redo_flag
		csvwrite(fullfile(bsn_path,'river_mouths_updated.txt'),RM);
	end
elseif isscalar(river_mouths) % 高程阈值处理
	disp('根据输入高程自动生成河流出口')
	sz=getnal(S,DEM);
	ix1=S.IXgrid;
	ix1(sz<river_mouths)=[];
	W=GRIDobj(DEM,'logical');
	W.Z(ix1)=true;
	Stemp=STREAMobj(FD,W);
	outix=streampoi(Stemp,'outlets','ix');
	if mbsz>0 % 过滤小流域
		DA=(A.*A.cellsize^2)/(1e6);
		da=DA.Z(outix);
		outix(da<mbsz)=[];
	end
	[rmx rmy]=ind2coord(DEM,outix);
	num_basins=numel(rmx);
	riv_nums=transpose(1:num_basins);
	RM=[rmx rmy riv_nums];
	csvwrite(fullfile(bsn_path,'river_mouths.txt'),RM);
end

% 处理零值ID
if nnz(RM(:,3))~=numel(RM(:,3))
	if isdeployed
		warndlg('检测到零值ID，已自动重新编号')
	end
	warning('检测到零值ID，已自动重新编号')
	zeroIDX=RM(:,3)==0;
	maxRM=max(RM(:,3));
	numZeros=sum(zeroIDX);
	zeroIX=find(zeroIDX);
	for ii=1:numZeros
		RM(zeroIX(ii),3)=maxRM+ii;
	end
	csvwrite(fullfile(bsn_path,'river_mouths_updated.txt'),RM);
end

% 主处理循环
w1=waitbar(0,['正在处理流域 1/' num2str(num_basins)]);
for ii=1:num_basins
	xx=RM(ii,1);
	yy=RM(ii,2);
	basin_num=RM(ii,3);

	RiverMouth=[xx yy basin_num];

	% 提取子流域
	I=dependencemap(FD,xx,yy);
	DEMoc=crop(DEM,I,nan);
	FDc=crop(FD,I);
	Ac=crop(A,I,nan);

	% 计算流域面积
	dep_map=GRIDobj2mat(I);
	num_pix=sum(sum(dep_map));
	drainage_area=(num_pix*DEMoc.cellsize*DEMoc.cellsize)/(1e6);

	% 计算高程直方图
	[rb,eb]=hypscurve(DEMoc,100);
	hyps=[rb eb];

	% 计算流域质心
	[Cx,Cy]=FindCentroid(DEMoc);
	Centroid=[Cx Cy];

	% 生成子流域河道网络
	Sc=STREAMobj(FDc,'minarea',threshold_area,'unit','mapunits');

	% 检查有效河道网络
	if isempty(Sc.x)
		if isdeployed
			warndlg(['流域 ' num2str(basin_num) ' 的集水面积阈值过大，已自动调整'])
		end
		warning(['流域 ' num2str(basin_num) ' 的集水面积阈值过大，已自动调整']);
		new_thresh=threshold_area;
		while isempty(Sc.x)
			new_thresh=new_thresh/2;
			Sc=STREAMobj(FDc,'minarea',new_thresh,'unit','mapunits');
		end
	end

	% 计算Chi值
	Cc=chitransform(Sc,Ac,'a0',1,'mn',theta_ref);
	ChiOBJc=GRIDobj(DEMoc);
	ChiOBJc.Z(Sc.IXgrid)=Cc;

	% 计算坡度
	switch gradient_method
	case 'gradient8'
		Goc=gradient8(DEMoc);
	case 'arcslope'
		Goc=arcslope(DEMoc);
	end

	% 水文校正DEM
	if isempty(DEMhc)
		zcon=mincosthydrocon(Sc,DEMoc,'interp',iv);
	else
		DEMhcc=crop(DEMhc,I,nan);
		zcon=getnal(Sc,DEMhcc);
	end
	DEMcc=GRIDobj(DEMoc);
	DEMcc.Z(DEMcc.Z==0)=NaN;
	DEMcc.Z(Sc.IXgrid)=zcon;

	% 拟合最佳凹度
	SLc=klargestconncomps(Sc,1);	
	Chic=chiplot(SLc,DEMcc,Ac,'a0',1,'plot',false);

	% 计算ksn值
	switch ksn_method
	case 'quick'
		[MSc]=KSN_Quick(DEMoc,DEMcc,Ac,Sc,Chic.mn,segment_length);
		[MSNc]=KSN_Quick(DEMoc,DEMcc,Ac,Sc,theta_ref,segment_length);
	case 'trunk'
		[MSc]=KSN_Trunk(DEMoc,DEMcc,Ac,Sc,Chic.mn,segment_length,min_order);
		[MSNc]=KSN_Trunk(DEMoc,DEMcc,Ac,Sc,theta_ref,segment_length,min_order);			
	case 'trib'
		if drainage_area>2.5 % 小流域使用快速方法
			[MSc]=KSN_Trib(DEMoc,DEMcc,FDc,Ac,Sc,Chic.mn,segment_length);
			[MSNc]=KSN_Trib(DEMoc,DEMcc,FDc,Ac,Sc,theta_ref,segment_length);
		else
			[MSc]=KSN_Quick(DEMoc,DEMcc,Ac,Sc,Chic.mn,segment_length);
			[MSNc]=KSN_Quick(DEMoc,DEMcc,Ac,Sc,theta_ref,segment_length);
		end
	end

	% 计算加权流积累量（降水）
	if weight_acc_flag
		WAc=crop(WA,I,nan);
		switch ksn_method
		case 'quick'
			[WMSc]=KSN_Quick(DEMoc,DEMcc,WAc,Sc,Chic.mn,segment_length);
			[WMSNc]=KSN_Quick(DEMoc,DEMcc,WAc,Sc,theta_ref,segment_length);
		case 'trunk'
			[WMSc]=KSN_Trunk(DEMoc,DEMcc,WAc,Sc,Chic.mn,segment_length,min_order);
			[WMSNc]=KSN_Trunk(DEMoc,DEMcc,WAc,Sc,theta_ref,segment_length,min_order);			
		case 'trib'
			if drainage_area>2.5
				[WMSc]=KSN_Trib(DEMoc,DEMcc,FDc,WAc,Sc,Chic.mn,segment_length);
				[WMSNc]=KSN_Trib(DEMoc,DEMcc,FDc,WAc,Sc,theta_ref,segment_length);
			else
				[WMSc]=KSN_Quick(DEMoc,DEMcc,WAc,Sc,Chic.mn,segment_length);
				[WMSNc]=KSN_Quick(DEMoc,DEMcc,WAc,Sc,theta_ref,segment_length);
			end
		end

		% 计算流域ksn-q统计量
		min_ksnq=min([WMSNc.ksn],[],'omitnan');
		mean_ksnq=mean([WMSNc.ksn],'omitnan');
		max_ksnq=max([WMSNc.ksn],[],'omitnan');
		std_ksnq=std([WMSNc.ksn],'omitnan');
		se_ksnq=std_ksnq/sqrt(numel(WMSNc)); % 标准误

		KSNQc_stats=[mean_ksnq se_ksnq std_ksnq min_ksnq max_ksnq];
	end

	% 计算流域统计量
	min_ksn=min([MSNc.ksn],[],'omitnan');
	mean_ksn=mean([MSNc.ksn],'omitnan');
	max_ksn=max([MSNc.ksn],[],'omitnan');
	std_ksn=std([MSNc.ksn],'omitnan');
	se_ksn=std_ksn/sqrt(numel(MSNc)); % 标准误

	min_grad=min(Goc.Z(:),[],'omitnan');
	mean_grad=mean(Goc.Z(:),'omitnan');
	max_grad=max(Goc.Z(:),[],'omitnan');
	std_grad=std(Goc.Z(:),'omitnan');
	se_grad=std_grad/sqrt(sum(~isnan(Goc.Z(:)))); % 标准误

	min_z=min(DEMoc.Z(:),[],'omitnan');
	mean_z=mean(DEMoc.Z(:),'omitnan');
	max_z=max(DEMoc.Z(:),[],'omitnan');
	std_z=std(DEMoc.Z(:),'omitnan');
	se_z=std_z/sqrt(sum(~isnan(DEMoc.Z(:)))); % 标准误

	KSNc_stats=[mean_ksn se_ksn std_ksn min_ksn max_ksn];
	Gc_stats=double([mean_grad se_grad std_grad min_grad max_grad]);
	Zc_stats=double([mean_z se_z std_z min_z max_z]);

	% 获取出口高程
	out_ix=coord2ind(DEMoc,xx,yy);
	out_el=double(DEMoc.Z(out_ix));

	% 保存结果文件
	FileName=fullfile(bsn_path,['Basin_' num2str(basin_num) '_Data.mat']);
	save(FileName,'RiverMouth','DEMcc','DEMoc','out_el','drainage_area','hyps','FDc','Ac','Sc','SLc','Chic','Goc','MSc','MSNc','KSNc_stats','Gc_stats','Zc_stats','Centroid','ChiOBJc','ksn_method','gradient_method','theta_ref','-v7.3');

	if weight_acc_flag
		save(FileName,'WAc','KSNQc_stats','-append');
	end

	if strcmp(ksn_method,'trunk')
		save(FileName,'min_order','-append');
	end

	% 生成ksn插值网格
	if ~isempty(radius)
		try 
			[KsnOBJc] = KsnAvg(DEMoc,MSNc,radius);
			save(FileName,'KsnOBJc','radius','-append');
		catch
			if isdeployed
				warndlg(['流域 ' num2str(RiverMouth(:,3)) ' 的KSN插值失败'])
			end
			warning(['流域 ' num2str(RiverMouth(:,3)) ' 的KSN插值失败']);
			save(FileName,'radius','-append');
		end
	else 
		save(FileName,'radius','-append')
    end

		% 如果存在附加网格，将其追加到mat文件中
		if ~isempty(AG)
			num_grids=size(AG,1);
			AGc=cell(size(AG));
			for jj=1:num_grids
				AGcOI=crop(AG{jj,1},I,nan);
				AGc{jj,1}=AGcOI;  % 裁剪后的网格
				AGc{jj,2}=AG{jj,2};  % 保留原始标识
				% 计算网格统计量
				mean_AGc=mean(AGcOI.Z(:),'omitnan');
				min_AGc=min(AGcOI.Z(:),[],'omitnan');
				max_AGc=max(AGcOI.Z(:),[],'omitnan');
				std_AGc=std(AGcOI.Z(:),'omitnan');
				se_AGc=std_AGc/sqrt(sum(~isnan(AGcOI.Z(:))));
				AGc_stats(jj,:)=[mean_AGc se_AGc std_AGc min_AGc max_AGc];
			end
			save(FileName,'AGc','AGc_stats','-append');  % 追加保存结果				
		end

		% 处理累积网格数据
		if ~isempty(ACG)
			num_grids=size(ACG,1);
			ACGc=cell(size(ACG));
			for jj=1:num_grids
				ACGcOI=crop(ACG{jj,1},I,nan);
				ACGc{jj,1}=ACGcOI;  % 裁剪后的累积网格
				ACGc{jj,3}=ACG{jj,3};  % 保留分类标识
				% 生成直方图边界
				edg=ACG{jj,2}.Numbers;
				edg=edg+0.5;
				edg=vertcat(0.5,edg);
				[N,~]=histcounts(ACGcOI.Z(:),edg);
				T=ACG{jj,2};
				T.Counts=N';  % 更新统计表
				ACGc{jj,2}=T;
				ACGc_stats(jj,1)=[mode(ACGcOI.Z(:))];  % 计算众数
			end
			save(FileName,'ACGc','ACGc_stats','-append');	
		end				

		% 计算地形起伏度指标
		if calc_relief
			num_rlf=numel(relief_radii);
			rlf=cell(num_rlf,2);
			rlf_stats=zeros(num_rlf,6);
			for jj=1:num_rlf
				% 计算指定半径的地形起伏度
				radOI=relief_radii(jj);
				rlf{jj,2}=radOI;  % 记录当前半径
				rlfOI=localtopography(DEMoc,radOI);
				rlf{jj,1}=rlfOI;
				% 计算统计指标
				mean_rlf=mean(rlfOI.Z(:),'omitnan');
				min_rlf=min(rlfOI.Z(:),[],'omitnan');
				max_rlf=max(rlfOI.Z(:),[],'omitnan');
				std_rlf=std(rlfOI.Z(:),'omitnan');
				se_rlf=std_rlf/sqrt(sum(~isnan(rlfOI.Z(:))));
				rlf_stats(jj,:)=[mean_rlf se_rlf std_rlf min_rlf max_rlf radOI];
			end
			save(FileName,'rlf','rlf_stats','-append');  % 保存地形起伏度数据
		end

		% 输出ArcGIS兼容文件
		if write_arc_files
			% 将DEM中的NaN替换为-9999
			Didx=isnan(DEMoc.Z);
			DEMoc_temp=DEMoc;
			DEMoc_temp.Z(Didx)=-9999;

			% 生成标准地形文件
			DEMFileName=fullfile(bsn_path,['Basin_' num2str(basin_num) '_DEM.txt']);
			GRIDobj2ascii(DEMoc_temp,DEMFileName);
			% 生成河道陡峭指数文件
			CHIFileName=fullfile(bsn_path,['Basin_' num2str(basin_num) '_CHI.txt']);
			GRIDobj2ascii(ChiOBJc,CHIFileName);
			% 生成KSN形状文件
			KSNFileName=fullfile(bsn_path,['Basin_' num2str(basin_num) '_KSN.shp']);
			shapewrite(MSNc,KSNFileName);

			% 输出地形起伏度文件
			if calc_relief
				for jj=1:num_rlf
					RLFFileName=fullfile(bsn_path,['Basin_' num2str(basin_num) '_RLF_' num2str(rlf{jj,2}) '.txt']);
					GRIDobj2ascii(rlf{jj,1},RLFFileName);
				end
			end

			% 输出附加网格文件
			if ~isempty(AG);
				for jj=1:num_grids
					AGcFileName=fullfile(bsn_path,['Basin_' num2str(basin_num) '_' AGc{jj,2} '.txt']);
					GRIDobj2ascii(AGc{jj,1},AGcFileName);
				end
			end

			% 输出累积网格文件
			if ~isempty(ACG);
				for jj=1:num_grids
					ACGcFileName=fullfile(bsn_path,['Basin_' num2str(basin_num) '_' ACGc{jj,3} '.txt']);
					GRIDobj2ascii(ACGc{jj,1},ACGcFileName);
				end
			end
		end

		% 更新进度条
		waitbar(ii/num_basins,w1,['已完成 ' num2str(ii) ' 个，总计 ' num2str(num_basins) ' 个流域'])
	end

	close(w1)  % 关闭进度条
	
end

function [ksn_ms]=KSN_Quick(DEM,DEMc,A,S,theta_ref,segment_length)
	% 快速计算河道陡峭指数
	g=gradient(S,DEMc);
	G=GRIDobj(DEM);
	G.Z(S.IXgrid)=g;

	Z_RES=DEMc-DEM;  % 计算侵蚀/堆积量

	ksn=G./(A.*(A.cellsize^2)).^(-theta_ref);  % 基本KSN公式

	% 准备空间分布数据
	SD=GRIDobj(DEM);
	SD.Z(S.IXgrid)=S.distance;
	
	% 生成地图结构
	ksn_ms=STREAMobj2mapstruct(S,'seglength',segment_length,'attributes',...
		{'ksn' ksn @mean 'uparea' (A.*(A.cellsize^2)) @mean 'gradient' G @mean 'cut_fill' Z_RES @mean...
		'min_dist' SD @min 'max_dist' SD @max});

	% 计算河段长度并优化字段
	seg_dist=[ksn_ms.max_dist]-[ksn_ms.min_dist];
	distcell=num2cell(seg_dist');
	[ksn_ms(1:end).seg_dist]=distcell{:};
	ksn_ms=rmfield(ksn_ms,{'min_dist','max_dist'});
end

function [ksn_ms]=KSN_Trunk(DEM,DEMc,A,S,theta_ref,segment_length,min_order)
	% 主干河道专用KSN计算

	order_exp=['>=' num2str(min_order)];  % 设置河道分级阈值

    Smax=modify(S,'streamorder',order_exp);  % 提取主干河道
	Smin=modify(S,'rmnodes',Smax);  % 剩余支流

	g=gradient(S,DEMc);
	G=GRIDobj(DEM);
	G.Z(S.IXgrid)=g;

	Z_RES=DEMc-DEM;

	ksn=G./(A.(A.cellsize^2)).^(-theta_ref);

	% 分别处理主干和支流
	SDmax=GRIDobj(DEM);
	SDmin=GRIDobj(DEM);
	SDmax.Z(Smax.IXgrid)=Smax.distance;
	SDmin.Z(Smin.IXgrid)=Smin.distance;

	% 支流处理
	ksn_ms_min=STREAMobj2mapstruct(Smin,'seglength',segment_length,'attributes',...
		{'ksn' ksn @mean 'uparea' (A.*(A.cellsize^2)) @mean 'gradient' G @mean 'cut_fill' Z_RES @mean...
		'min_dist' SDmin @min 'max_dist' SDmin @max});

	% 主干处理
	ksn_ms_max=STREAMobj2mapstruct(Smax,'seglength',segment_length,'attributes',...
		{'ksn' ksn @mean 'uparea' (A.*(A.cellsize^2)) @mean 'gradient' G @mean 'cut_fill' Z_RES @mean...
		'min_dist' SDmax @min 'max_dist' SDmax @max});

	% 合并结果
	ksn_ms=vertcat(ksn_ms_min,ksn_ms_max);
	seg_dist=[ksn_ms.max_dist]-[ksn_ms.min_dist];
	distcell=num2cell(seg_dist');
	[ksn_ms(1:end).seg_dist]=distcell{:};
	ksn_ms=rmfield(ksn_ms,{'min_dist','max_dist'});
end

function [ksn_ms]=KSN_Trib(DEM,DEMc,FD,A,S,theta_ref,segment_length)
	% 支流分段KSN计算方法

	% 定义不相交的河段
	[as]=networksegment_slim(DEM,FD,S);
	seg_bnd_ix=as.ix;
	% 预提取高程参数
	z=getnal(S,DEMc);
	zu=getnal(S,DEM);
	z_res=z-zu;
	g=gradient(S,DEMc);
	c=chitransform(S,A,'a0',1,'mn',theta_ref);  % 计算Chi值
	d=S.distance;
	da=getnal(S,A.*(A.cellsize^2));
	ixgrid=S.IXgrid;
	% 获取有序节点列表
	s_node_list=S.orderednanlist;
	streams_ix=find(isnan(s_node_list));
	streams_ix=vertcat(1,streams_ix);
	% 初始化属性列表
	ksn_nal=zeros(size(d));
	% 主循环处理各河道
	num_streams=numel(streams_ix)-1;
	seg_count=1;
	for ii=1:num_streams
		% 提取当前河道节点
		if ii==1
			snlOI=s_node_list(streams_ix(ii):streams_ix(ii+1)-1);
		else
			snlOI=s_node_list(streams_ix(ii)+1:streams_ix(ii+1)-1);
		end

		% 确定当前河道包含的段
		[~,~,dn]=intersect(snlOI,seg_bnd_ix(:,1));
		[~,~,up]=intersect(snlOI,seg_bnd_ix(:,2));
		seg_ix=intersect(up,dn);

		num_segs=numel(seg_ix);
		dn_up=seg_bnd_ix(seg_ix,:);
		for jj=1:num_segs
			% 定位段节点
			dnix=find(snlOI==dn_up(jj,1));
			upix=find(snlOI==dn_up(jj,2));
			seg_ix_oi=snlOI(upix:dnix);
			% 标准化流程距离
			dOI=d(seg_ix_oi);
			dnOI=dOI-min(dOI);
			num_bins=ceil(max(dnOI)/segment_length);
			bin_edges=[0:segment_length:num_bins*segment_length];
			% 分箱计算KSN
			for kk=1:num_bins
				idx=dnOI>bin_edges(kk) & dnOI<=bin_edges(kk+1);
				bin_ix=seg_ix_oi(idx);
				cOI=c(bin_ix);
				zOI=z(bin_ix);
					if numel(cOI)>2
						[ksn_val,r2]=Chi_Z_Spline(cOI,zOI);
						ksn_nal(bin_ix)=ksn_val;

						% 构建地图要素
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
	% 精简版'networksegment'函数，源自TopoToolbox主库，同时移除长度为零或单节点的段

	%% 识别河道源头、汇合点、分支汇合点和出口
	Vhead = streampoi(S,'channelheads','logical');  ihead=find(Vhead==1);  IXhead=S.IXgrid(ihead);
	Vconf = streampoi(S,'confluences','logical');   iconf=find(Vconf==1);  IXconf=S.IXgrid(iconf);
	Vout = streampoi(S,'outlets','logical');        iout=find(Vout==1);    IXout=S.IXgrid(iout);
	Vbconf = streampoi(S,'bconfluences','logical'); ibconf=find(Vbconf==1);IXbconf=S.IXgrid(ibconf);

	%% 建立汇水区关联
	DB   = drainagebasins(FD,vertcat(IXbconf,IXout));DBhead=DB.Z(IXhead); DBbconf=DB.Z(IXbconf); DBconf=DB.Z(IXconf); DBout=DB.Z(IXout);

	%% 计算流程距离
	D = flowdistance(FD);

	%% 建立河段连接关系
	[~,ind11,ind12]=intersect(DBbconf,DBhead);
	[~,ind21,ind22]=intersect(DBbconf,DBconf);
	[~,ind31,ind32]=intersect(DBout,DBhead);
	[~,ind41,ind42]=intersect(DBout,DBconf);
	IX(:,1) = [ IXbconf(ind11)' IXbconf(ind21)' IXout(ind31)'  IXout(ind41)'  ];   ix(:,1)= [ ibconf(ind11)' ibconf(ind21)' iout(ind31)'  iout(ind41)'  ];
	IX(:,2) = [ IXhead(ind12)'  IXconf(ind22)'  IXhead(ind32)' IXconf(ind42)' ];   ix(:,2)= [ ihead(ind12)'  iconf(ind22)'  ihead(ind32)' iconf(ind42)' ];

	%% 过滤有效河段
	flength=double(abs(D.Z(IX(:,1))-D.Z(IX(:,2))));
	idx=flength>=2*DEM.cellsize;  % 剔除无效短段
	seg.IX=IX(idx,:);
	seg.ix=ix(idx,:);
	seg.flength=flength(idx);
	seg.n=numel(IX(:,1));
end

function [KSN,R2] = Chi_Z_Spline(c,z)
	% 通过样条插值计算KSN和拟合优度

	% 重采样Chi-高程曲线
	[~,minIX]=min(c);
	zb=z(minIX);
	chiF=c-min(c);
	zabsF=z-min(z);
	chiS=linspace(0,max(chiF),numel(chiF)).';
	zS=spline(chiF,zabsF,chiS);

	% 线性回归求KSN
	KSN= chiS\(zS);  % 使用最小二乘法求解

	% 计算R²
	z_pred=chiF.*KSN;
	sstot=sum((zabsF-mean(zabsF)).^2);
	ssres=sum((zabsF-z_pred).^2);
	R2=1-(ssres/sstot);
end

function [KSNGrid] = KsnAvg(DEM,ksn_ms,radius)
	% 空间平均KSN计算

	% 计算像素半径
	radiuspx = ceil(radius/DEM.cellsize);

	% 创建初始栅格
	KSNGrid=GRIDobj(DEM);
	KSNGrid.Z(:,:)=NaN;
	% 填入河道KSN值
	for ii=1:numel(ksn_ms)
		ix=coord2ind(DEM,ksn_ms(ii).X,ksn_ms(ii).Y);
		KSNGrid.Z(ix)=ksn_ms(ii).ksn;
	end

	% 使用圆盘滤波器进行空间平均
	ISNAN=isnan(KSNGrid.Z);
    [~,L] = bwdist(~ISNAN,'e');
    ksng = KSNGrid.Z(L);           
    FLT   = fspecial('disk',radiuspx);
    ksng   = imfilter(ksng,FLT,'symmetric','same','conv');

    % 恢复原始掩膜
    ksng(MASK)=NaN;
    KSNGrid.Z=ksng;
end