function [varargout]=KsnChiBatch(DEM,FD,A,S,product,varargin)
	%
% 用法：
%	KsnChiBatch(DEM,FD,A,S,product);
%	KsnChiBatch(DEM,FD,A,S,product,'name',value...);
%	[outputs]=KsnChiBatch(DEM,FD,A,S,product,'output',true,...);
%
% 描述：
% 	用于生成数字高程模型（DEM）中所有通道的通道陡度（ksn）、chi 图或 chi 网格的函数
% 
% 必需输入：
% 	DEM - DEM 网格对象（假定为未平滑的 DEM）
% 	FD - 流动对象
% 	A - 流积的网格对象
%	S - 流对象
% 	product - 开关，用于确定生成哪些产品
%		'ksn' - 作为 shapefile 的 ksn 图
%		'ksngrid' - 使用移动圆形窗口在网格的所有点上平均插值的 ksn 的 ASCII 文件
%		'chimap' - 在通道网络中计算 chi 的 ASCII 文件
%		'chigrid' - 在网格的所有点上计算 chi 的 ASCII 文件
%		'chi' - chimap 和 chigrid 的结果
%		'all' - ksn、ksngrid、chimap 和 chigrids
%
% 可选输入：
%	conditioned_DEM [] - 可选提供平滑的 DEM 供该函数使用（不要将平滑DEM 提供为主要必需 DEM 输入！），
%		将用于提取高程。请参见 'ConditionDEM' 函数以获取制作平滑 DEM 的选项。如果没有提供输入，代码默认使用 mincosthydrocon 函数。
%	file_name_prefix ['batch'] - 输出的前缀，将附加输出类型，即 'ksn'、'chimap' 等
%	smooth_distance [1000] - 平滑 ksn 测量的距离（地图单位）
% 	ref_concavity [0.50] - 计算 ksn 的参考凹度（正值）
% 	output [false]- 开关，用于将 MATLAB 文件输出到工作区（true），
%		或不输出，只保存指定的文件而不进行工作区输出（false）。输出的数量将取决于 'product' 输入。
%			'ksn' - [KSNG,ksn_ms] 其中 KSNG 是流网络上的 ksn 值的 GRIDobj，ksn_ms 是适合创建 shapefile 的 mapstructure。
%			如果提供了 precip_grid，则输出将为 [KSNG,ksn_ms,KSNQG,ksnq_ms]，后两个为加权 ksn-q 值（例如，Adams et al, 2020）
%			'ksngrid' - [KSNgrid,KSNstdGrid] 其中 KSNgrid 是插值 ksn 值的 GRIDobj，KSNstdGrid 是指定半径内 ksn 的标准偏差的 GRIDobj。
%			如果提供了 precip_grid，则输出将为 [KSNgrid,KSNstdGrid,KSNQgrid,KSNQstdGrid]
%			'chimap' - [ChiMap] 其中 ChiMap 是流网络上的 chi 值的 GRIDobj
%			'chigrid' - [ChiGrid] 其中 ChiGrid 是整个网格上的 chi 值的 GRIDobj
%			'chi' - [ChiMap,ChiGrid]
%			'all' - [KSNG,ksn_ms,KSNGrid,KSNstdGrid,ChiMap,ChiGrid]，如果提供了 precip_grid，则输出将为 
%					[KSNG,ksn_ms,KSNGrid,KSNstdGrid,ChiMap,ChiGrid,KSNQG,ksnq_ms,KSNQGrid,KSNQstdGrid]
%	ksn_method [quick] - 切换计算 ksn 值的方法，选项为 'quick'、'trunk' 或 'trib'，'trib' 方法比 'quick' 方法慢 3-4 倍。
%		在大多数情况下，'quick' 方法效果良好，但如果靠近支流交汇处的值很重要，'trib' 可能更好，
%		因为这会为每个通道段单独计算 ksn 值。'trunk' 选项独立计算大河流的陡度值（被认为是干流的河流由 
%		提供给 'min_order' 的流序值控制）。如果您注意到主干流的通道陡度值异常高，可能是由于值的平均方式。
%	precip_grid [] - 可选的降水 GRIDobj 输入。如果提供了这个参数，将有额外的输出对应 ksn-q 值（请参见 Adams et al, 2020），
%		即降水加权的 ksn 值。请不要同时提供这个可选的降水网格和加权降水的流积网格 "A"，
%		因为这会产生错误的值。降水网格输入应覆盖与提供的 DEM 相同或更大的区域，但不需要与之具有相同的单元大小，
%		因为函数会在必要时对降水网格进行重采样。
%	min_order [4] - 被视为干流的流的最小流序，仅在 'ksn_method' 设置为 'trunk' 时使用
%	outlet_level_method [] - 控制流网络出口水平调整的参数。对出口高程的控制选项包括：
%			'elevation' - 提取仅在给定高程（通过用户提供的 'min_elevation' 参数提供）之上的河流，以确保所有河流的基准高程
%				均匀。如果提供的高程过低（即某些未经修改的河流网络的出口高于此高程），则会显示警告，但代码仍将运行。
%			'max_out_elevation' - 使用所有河流出口的最高高程，仅提取该高程之上的河流，仅适用于仅操作河流线的选项（即不适用于 'ksngrid' 或 'chigrid'）。
%	min_elevation [] - 设置出口水平的最小高程的参数，当 'outlet_level_method' 设置为 'elevation' 时需要提供
%	complete_networks_only [true] - 如果为 true（默认），代码仅填充完整的流网络部分。通常，该选项应该保留为 true 
%			（即如果排水面积不准确，则 chi 将不准确），但这可能在某些 DEM 上过于激进，并且在与 'min_elevation' 
%			一起使用时，计算速度可能较慢，因为这需要重新计算 FLOWobj。
%	interp_value [0.1] - 用于 mincosthydrocon 的插值参数的值（在用户提供平滑 DEM 时不使用）
%	radius [5000] - 用于在设置为 'ksngrid' 或 'all' 时生成插值 ksn 网格的平均 ksn 值的移动区域的半径
%	error_type ['std'] - 控制计算 KsnGrid 时第二个输出是标准偏差（默认）还是标准误差。如果要计算标准误差而不是标准偏差，请提供 'std_error' 作为可选输入。
%
% 注意事项：
%	请注意，生成 ksngrid 和/或 chigrid 可能会耗时，因此请耐心等待...
%
% 示例：
%	KsnChiBatch(DEM,FD,A,S,'ksn');
%	[KSNG,ksn_ms]=KsnChiBatch(DEM,FD,A,S,'ksn','output',true,'theta_ref',0.55);
%	[KSNG,ksn_ms,KSNQG,ksnq_ms]=KsnChiBatch(DEM,FD,A,S,'ksn','output',true,'theta_ref',0.45,'precip_grid',PRECIP);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数编写者：Adam M. Forte - 更新时间：2020 年 10 月 17 日 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	% Parse Inputs
	p = inputParser;         
	p.FunctionName = 'KsnChiBatch';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'FD', @(x) isa(x,'FLOWobj'));
	addRequired(p,'A', @(x) isa(x,'GRIDobj'));
	addRequired(p,'S',@(x) isa(x,'STREAMobj'));
	addRequired(p,'product',@(x) ischar(validatestring(x,{'ksn','ksngrid','chimap','chigrid','chi','all'})));

	addParameter(p,'file_name_prefix','batch',@(x) ischar(x));
	addParameter(p,'smooth_distance',1000,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'ref_concavity',0.50,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'output',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'ksn_method','quick',@(x) ischar(validatestring(x,{'quick','trunk','trib'})));
	addParameter(p,'precip_grid',[],@(x) isa(x,'GRIDobj') || isempty(x));	
	addParameter(p,'min_order',4,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'outlet_level_method',[],@(x) isempty(x) || ischar(validatestring(x,{'elevation','max_out_elevation'})));
	addParameter(p,'min_elevation',[],@(x) isnumeric(x) || isempty(x));
	addParameter(p,'conditioned_DEM',[],@(x) isa(x,'GRIDobj') || isempty(x));
	addParameter(p,'interp_value',0.1,@(x) isnumeric(x) && x>=0 && x<=1);
	addParameter(p,'complete_networks_only',true,@(x) isscalar(x) && islogical(x));
	addParameter(p,'radius',5000,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'error_type','std',@(x) ischar(validatestring(x,{'std','std_error'})));

	parse(p,DEM,FD,A,S,product,varargin{:});
	DEM=p.Results.DEM;
	FD=p.Results.FD;
	A=p.Results.A;
	S=p.Results.S;
	product=p.Results.product;

	file_name_prefix=p.Results.file_name_prefix;
	segment_length=p.Results.smooth_distance;
	theta_ref=p.Results.ref_concavity;
	output=p.Results.output;
	ksn_method=p.Results.ksn_method;
	precip_grid=p.Results.precip_grid;	
	min_order=p.Results.min_order;
	blm=p.Results.outlet_level_method;
	me=p.Results.min_elevation;
	iv=p.Results.interp_value;
	DEMc=p.Results.conditioned_DEM;
	cno=p.Results.complete_networks_only;
	radius=p.Results.radius;
	error_type=p.Results.error_type;

% 检查在设置基准面时是否提供了截止值
	if strcmp(blm,'max_out_elevation') & (strcmp(product,'chigrid') | strcmp(product,'ksngrid') | strcmp(product,'all'))
		if isdeployed
			warndlg('"max_out_elevation" 不是连续网格产品的有效选项，将忽略此输入')
		else
			warning('"max_out_elevation" 不是连续网格产品的有效选项，将忽略此输入')
		end
	elseif strcmp(blm,'elevation') & isempty(me)
		if isdeployed
			warndlg('选择的出口标高调整方法 "elevation" 需要提供参数 "min_elevation" 的输入，将在没有基准面控制的情况下运行')
		else
			warning('选择的出口标高调整方法 "elevation" 需要提供参数 "min_elevation" 的输入，将在没有基准面控制的情况下运行');
		end
		blm=[];
	end


	% Check for precip grid and product, resample if necessary, and make precip weighted flow accumulation GRIDobj
	if ~isempty(precip_grid) & (strcmp(product,'ksn') | strcmp(product,'ksngrid') | strcmp(product,'all'));
		disp('Calculating precpitation weighted flow accumulation grid')
		if ~validatealignment(precip_grid,DEM);
			precip_grid=resample(precip_grid,DEM,'nearest');
		end
		WA=flowacc(FD,precip_grid);
		q_flag=true;
	else
		q_flag=false;
	end

	switch product
	case 'ksn'

		disp('Calculating channel steepness')
		
		if strcmp(blm,'elevation')
			[S]=DTSetOutlet(DEM,FD,A,S,blm,'complete_networks_only',cno,'min_elevation',me);
		elseif strcmp(blm,'max_out_elevation')
			[S]=DTSetOutlet(DEM,FD,A,S,blm,'complete_networks_only',cno);
		elseif isempty(blm) & cno
			[S]=DTSetOutlet(DEM,FD,A,S,'complete_only','complete_networks_only',cno);
		end

		% Hydrologically condition dem
		if isempty(DEMc)
			zc=mincosthydrocon(S,DEM,'interp',iv);
			DEMc=GRIDobj(DEM);
			DEMc.Z(DEMc.Z==0)=NaN;
			DEMc.Z(S.IXgrid)=zc;
		end
		
		switch ksn_method
		case 'quick'
			[ksn_ms]=KSN_Quick(DEM,DEMc,A,S,theta_ref,segment_length);
		case 'trunk'
			[ksn_ms]=KSN_Trunk(DEM,DEMc,A,S,theta_ref,segment_length,min_order);
		case 'trib'
			[ksn_ms]=KSN_Trib(DEM,DEMc,FD,A,S,theta_ref,segment_length);
		end

		disp('Writing ARC files')
		out_file=[file_name_prefix '_ksn.shp'];
		shapewrite(ksn_ms,out_file);

		if q_flag
			disp('Calculating precipitation weighted channel steepness')

			switch ksn_method
			case 'quick'
				[ksnq_ms]=KSN_Quick(DEM,DEMc,WA,S,theta_ref,segment_length);
			case 'trunk'
				[ksnq_ms]=KSN_Trunk(DEM,DEMc,WA,S,theta_ref,segment_length,min_order);
			case 'trib'
				[ksnq_ms]=KSN_Trib(DEM,DEMc,FD,WA,S,theta_ref,segment_length);
			end

			disp('Writing ARC files')
			out_file=[file_name_prefix '_ksnq.shp'];
			shapewrite(ksnq_ms,out_file);			
		end

		switch output
		case true
			if q_flag
				KSNG=GRIDobj(DEM);
				KSNG.Z(:,:)=NaN;
				for ii=1:numel(ksn_ms)
					ix=coord2ind(DEM,ksn_ms(ii).X,ksn_ms(ii).Y);
					KSNG.Z(ix)=ksn_ms(ii).ksn;
				end
				KSNQG=GRIDobj(DEM);
				KSNQG.Z(:,:)=NaN;
				for ii=1:numel(ksnq_ms)
					ix=coord2ind(DEM,ksnq_ms(ii).X,ksnq_ms(ii).Y);
					KSNQG.Z(ix)=ksnq_ms(ii).ksn;
				end
				varargout{1}=KSNG;
				varargout{2}=ksn_ms;
				varargout{3}=KSNQG;
				varargout{4}=ksnq_ms;
			else
				KSNG=GRIDobj(DEM);
				KSNG.Z(:,:)=NaN;
				for ii=1:numel(ksn_ms)
					ix=coord2ind(DEM,ksn_ms(ii).X,ksn_ms(ii).Y);
					KSNG.Z(ix)=ksn_ms(ii).ksn;
				end
				varargout{1}=KSNG;
				varargout{2}=ksn_ms;
			end
		end

	case 'ksngrid'
		disp('Calculating channel steepness')
		if strcmp(blm,'elevation')
			IDX=DEM<me;
			DEM.Z(IDX.Z)=NaN;
		else
			IDX=GRIDobj(DEM,'logical');
		end

		if isempty(DEMc)
			zc=mincosthydrocon(S,DEM,'interp',iv);
			DEMc=GRIDobj(DEM);
			DEMc.Z(DEMc.Z==0)=NaN;
			DEMc.Z(S.IXgrid)=zc;
		else
			DEMc.Z(IDX.Z)=NaN;
		end

		if strcmp(blm,'elevation')
			[S]=DTSetOutlet(DEM,FD,A,S,blm,'complete_networks_only',cno,'min_elevation',me);
		elseif strcmp(blm,'max_out_elevation')
			[S]=DTSetOutlet(DEM,FD,A,S,blm,'complete_networks_only',cno);
		elseif isempty(blm) & cno
			[S]=DTSetOutlet(DEM,FD,A,S,'complete_only','complete_networks_only',cno);
		end
		
		switch ksn_method
		case 'quick'
			[ksn_ms]=KSN_Quick(DEM,DEMc,A,S,theta_ref,segment_length);
		case 'trunk'
			[ksn_ms]=KSN_Trunk(DEM,DEMc,A,S,theta_ref,segment_length,min_order);
		case 'trib'
			[ksn_ms]=KSN_Trib(DEM,DEMc,FD,A,S,theta_ref,segment_length);
		end

		disp('Calculating gridded channel steepness')
		[KSNGrid,KSNstdGrid]=KsnAvg(DEM,ksn_ms,radius,error_type);

		disp('Writing ARC files')
		out_file_ksng=[file_name_prefix '_ksngrid.txt'];
		out_file_ksnstd=[file_name_prefix '_ksnstdgrid.txt'];
		GRIDobj2ascii(KSNGrid,out_file_ksng);
		GRIDobj2ascii(KSNstdGrid,out_file_ksnstd);

		if q_flag
			disp('Calculating precipitation weighted channel steepness')
			switch ksn_method
			case 'quick'
				[ksnq_ms]=KSN_Quick(DEM,DEMc,WA,S,theta_ref,segment_length);
			case 'trunk'
				[ksnq_ms]=KSN_Trunk(DEM,DEMc,WA,S,theta_ref,segment_length,min_order);
			case 'trib'
				[ksnq_ms]=KSN_Trib(DEM,DEMc,FD,WA,S,theta_ref,segment_length);
			end

			disp('Calculating precipitation weighted gridded channel steepness')
			[KSNQGrid,KSNQstdGrid]=KsnAvg(DEM,ksnq_ms,radius,error_type);

			disp('Writing ARC files')
			out_file_ksng=[file_name_prefix '_ksnqgrid.txt'];
			out_file_ksnstd=[file_name_prefix '_ksnqstdgrid.txt'];
			GRIDobj2ascii(KSNQGrid,out_file_ksng);
			GRIDobj2ascii(KSNQstdGrid,out_file_ksnstd);			
		end		

		switch output
		case true
			if q_flag
				varargout{1}=KSNGrid;
				varargout{2}=KSNstdGrid;
				varargout{3}=KSNQGrid;
				varargout{4}=KSNQstdGrid;				
			else
				varargout{1}=KSNGrid;
				varargout{2}=KSNstdGrid;				
			end
		end

	case 'chimap'

		if strcmp(blm,'elevation')
			[S]=DTSetOutlet(DEM,FD,A,S,blm,'complete_networks_only',cno,'min_elevation',me);
		elseif strcmp(blm,'max_out_elevation')
			[S]=DTSetOutlet(DEM,FD,A,S,blm,'complete_networks_only',cno);
		elseif isempty(blm) & cno
			[S]=DTSetOutlet(DEM,FD,A,S,'complete_only','complete_networks_only',cno);
		end

	    disp('Calculating chi map');
		[ChiMap]=MakeChiMap(DEM,FD,A,S,theta_ref);

		disp('Writing ARC files')
		out_file_cm=[file_name_prefix '_chimap.txt'];
		GRIDobj2ascii(ChiMap,out_file_cm);

		switch output
		case true
			varargout{1}=ChiMap;
		end

	case 'chigrid'

	    disp('Calculating chi grid');
		[ChiGrid]=MakeChiGrid(DEM,FD,'theta_ref',theta_ref,'complete_networks_only',cno,'min_elevation',me);

	    disp('Writing ARC files')
		out_file=[file_name_prefix '_chigrid.txt'];
		GRIDobj2ascii(ChiGrid,out_file);

		switch output
		case true
			varargout{1}=ChiGrid;
		end

	case 'chi'

		if strcmp(blm,'elevation')
			[S]=DTSetOutlet(DEM,FD,A,S,blm,'complete_networks_only',cno,'min_elevation',me);
		elseif strcmp(blm,'max_out_elevation')
			[S]=DTSetOutlet(DEM,FD,A,S,blm,'complete_networks_only',cno);
		elseif isempty(blm) & cno
			[S]=DTSetOutlet(DEM,FD,A,S,'complete_only','complete_networks_only',cno);
		end

	    disp('Calculating chi map');
		[ChiMap]=MakeChiMap(DEM,FD,A,S,theta_ref);

	    disp('Calculating chi grid');
		[ChiGrid]=MakeChiGrid(DEM,FD,'theta_ref',theta_ref,'complete_networks_only',cno,'min_elevation',me);

		disp('Writing ARC files')
		out_file_cg=[file_name_prefix '_chigrid.txt'];
		GRIDobj2ascii(ChiGrid,out_file_cg);
		out_file_cm=[file_name_prefix '_chimap.txt'];
		GRIDobj2ascii(ChiMap,out_file_cm);

		switch output
		case true
			varargout{1}=ChiMap;
			varargout{2}=ChiGrid;
		end

	case 'all'

		if strcmp(blm,'elevation')
			IDX=DEM<me;
			DEM.Z(IDX.Z)=NaN;
		else
			IDX=GRIDobj(DEM,'logical');
		end

		disp('Calculating channel steepness')
		if strcmp(blm,'elevation')
			[S]=DTSetOutlet(DEM,FD,A,S,blm,'complete_networks_only',cno,'min_elevation',me);
		elseif strcmp(blm,'max_out_elevation')
			[S]=DTSetOutlet(DEM,FD,A,S,blm,'complete_networks_only',cno);
		elseif isempty(blm) & cno
			[S]=DTSetOutlet(DEM,FD,A,S,'complete_only','complete_networks_only',cno);
		end

		if isempty(DEMc)
			zc=mincosthydrocon(S,DEM,'interp',iv);
			DEMc=GRIDobj(DEM);
			DEMc.Z(DEMc.Z==0)=NaN;
			DEMc.Z(S.IXgrid)=zc;
		else
			DEMc.Z(IDX.Z)=NaN;
		end
		
		switch ksn_method
		case 'quick'
			[ksn_ms]=KSN_Quick(DEM,DEMc,A,S,theta_ref,segment_length);
		case 'trunk'
			[ksn_ms]=KSN_Trunk(DEM,DEMc,A,S,theta_ref,segment_length,min_order);
		case 'trib'
			[ksn_ms]=KSN_Trib(DEM,DEMc,FD,A,S,theta_ref,segment_length);
		end

		disp('Calculating gridded ksn grid')
		[KSNGrid,KSNstdGrid]=KsnAvg(DEM,ksn_ms,radius,error_type);

		if q_flag
			disp('Calculating precipitation weighted channel steepness')

			switch ksn_method
			case 'quick'
				[ksnq_ms]=KSN_Quick(DEM,DEMc,WA,S,theta_ref,segment_length);
			case 'trunk'
				[ksnq_ms]=KSN_Trunk(DEM,DEMc,WA,S,theta_ref,segment_length,min_order);
			case 'trib'
				[ksnq_ms]=KSN_Trib(DEM,DEMc,FD,WA,S,theta_ref,segment_length);
			end

			disp('Calculating precipitation weighted gridded channel steepness')
			[KSNQGrid,KSNQstdGrid]=KsnAvg(DEM,ksnq_ms,radius,error_type);		
		end	

	    disp('Calculating chi map');
		[ChiMap]=MakeChiMap(DEM,FD,A,S,theta_ref);

	    disp('Calculating chi grid');
		[ChiGrid]=MakeChiGrid(DEM,FD,'theta_ref',theta_ref,'complete_networks_only',cno,'min_elevation',me);

		disp('Writing ARC files')
		out_file_ksn=[file_name_prefix '_ksn.shp'];
		shapewrite(ksn_ms,out_file_ksn);
		out_file_ksng=[file_name_prefix '_ksngrid.txt'];
		out_file_ksnstd=[file_name_prefix '_ksnstdgrid.txt'];
		GRIDobj2ascii(KSNGrid,out_file_ksng);
		GRIDobj2ascii(KSNstdGrid,out_file_ksnstd);
		out_file_cg=[file_name_prefix '_chigrid.txt'];
		GRIDobj2ascii(ChiGrid,out_file_cg);
		out_file_cm=[file_name_prefix '_chimap.txt'];
		GRIDobj2ascii(ChiMap,out_file_cm);
		if q_flag
			out_file_ksnqg=[file_name_prefix '_ksnqgrid.txt'];
			out_file_ksnqstd=[file_name_prefix '_ksnqstdgrid.txt'];
			GRIDobj2ascii(KSNQGrid,out_file_ksnqg);
			GRIDobj2ascii(KSNQstdGrid,out_file_ksnqstd);	
		end

		switch output
		case true
			if q_flag
				KSNG=GRIDobj(DEM);
				KSNG.Z(:,:)=NaN;
				for ii=1:numel(ksn_ms)
					ix=coord2ind(DEM,ksn_ms(ii).X,ksn_ms(ii).Y);
					KSNG.Z(ix)=ksn_ms(ii).ksn;
				end
				varargout{1}=KSNG;
				varargout{2}=ksn_ms;
				varargout{3}=KSNGrid;
				varargout{4}=KSNstdGrid;
				varargout{5}=ChiMap;
				varargout{6}=ChiGrid;

				KSNQG=GRIDobj(DEM);
				KSNQG.Z(:,:)=NaN;
				for ii=1:numel(ksnq_ms)
					ix=coord2ind(DEM,ksnq_ms(ii).X,ksnq_ms(ii).Y);
					KSNQG.Z(ix)=ksnq_ms(ii).ksn;
				end
				varargout{7}=KSNQG;
				varargout{8}=ksnq_ms;
				varargout{9}=KSNQGrid;
				varargout{10}=KSNQstdGrid;				
			else
				KSNG=GRIDobj(DEM);
				KSNG.Z(:,:)=NaN;
				for ii=1:numel(ksn_ms)
					ix=coord2ind(DEM,ksn_ms(ii).X,ksn_ms(ii).Y);
					KSNG.Z(ix)=ksn_ms(ii).ksn;
				end
				varargout{1}=KSNG;
				varargout{2}=ksn_ms;
				varargout{3}=KSNGrid;
				varargout{4}=KSNstdGrid;
				varargout{5}=ChiMap;
				varargout{6}=ChiGrid;
			end
		end
	end

% Main Function End
end

function [SC]=DTSetOutlet(DEM,FD,A,S,method,varargin)
	% Clone of Divide Tools 'SetOutlet' function.

	% Parse Inputs
	p = inputParser;         
	p.FunctionName = 'DTSetOutlet';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'FD', @(x) isa(x,'FLOWobj'));
	addRequired(p,'A', @(x) isa(x,'GRIDobj'));
	addRequired(p,'S',@(x) isa(x,'STREAMobj'));
	addRequired(p,'method',@(x) ischar(validatestring(x,{'elevation','max_out_elevation','complete_only'})));

	addParameter(p,'complete_networks_only',true,@(x) islogical(x) & isscalar(x));
	addParameter(p,'min_elevation',[],@(x) isnumeric(x));


	parse(p,DEM,FD,A,S,method,varargin{:});
	DEM=p.Results.DEM;
	FD=p.Results.FD;
	A=p.Results.A;
	S=p.Results.S;
	method=p.Results.method;

	cno=p.Results.complete_networks_only;
	me=p.Results.min_elevation;

	if ~cno & strcmp(method,'complete_networks_only')
		if isdeployed
			errordlg('不能将方法设置为 complete_only 并且将 "仅限完整河流网络" 设置为 false')
		end
		error('不能将方法设置为 complete_only 并且将"仅限完整河流网络"  设置为 false');
	end

	if cno & ~strcmp(method,'complete_networks_only')
		S=removeedgeeffects(S,FD,DEM);
	end

	%% Initiate graphical picker if no values for either min drainage area or min elevation are provided
	if strcmp(method,'elevation') & isempty(me)
		f1=figure(1);
		set(f1,'Units','normalized','Position',[0.1 0.1 0.75 0.75]);
		hold on
		imageschs(DEM,DEM);
		title('Zoom to desired view and press enter')
		pause()
		title('Pick point on DEM to select base level elevation');
		[x,y]=ginput(1);
		hold off
		close(f1)
		el_ix=coord2ind(DEM,x,y);
		me=DEM.Z(el_ix);
		disp(['Selected elevation is : ' num2str(me) ' m']);
	end

	%% Main switch between methods
	switch method
	case 'elevation'
		st_el=getnal(S,DEM);
		idx=st_el>=me;
		IX=S.IXgrid(idx);
		W=GRIDobj(DEM);
		W.Z(IX)=1;
		W.Z=logical(W.Z);
		SC=STREAMobj(FD,W);
		% Check to see if all outlets meet the condition
		coix=streampoi(SC,'outlets','ix');
		coel=DEM.Z(coix);
		max_coel=max(coel);
		if sum(coel>me)~=0 & ~cno
			if isdeployed
				warndlg(['一个或多个河流出口高于提供的高度，最大河流出口高度为 ' num2str(max_coel)])
			else 
				warning(['一个或多个河流出口高于提供的高度，最大河流出口高度为 ' num2str(max_coel)]);
            end
		elseif sum(coel>me)~=0 & cno
			[xo,yo]=getoutline(DEM,true);
			% Control for incosistent output of getoutline
			sz=size(xo);
			if sz(1)==1 & sz(2)>1
				[oxy]=[xo' yo'];
			elseif sz(2)==1 & sz(1)>1
				[oxy]=[xo yo];
			end
			[coxy]=streampoi(SC,'outlets','xy');
			idx2=coel>me & ismember(coxy,oxy,'rows'); % Find streams with outlets greater than min elevation AND along boundary of DEM
			coix(idx2)=[];
			W=GRIDobj(DEM);
			W.Z(coix)=1;
			W.Z=logical(W.Z);
			SC=modify(SC,'upstreamto',W);			
		end
	case 'max_out_elevation'
		coix=streampoi(S,'outlets','ix');
		coel=DEM.Z(coix);
		max_coel=max(coel);
		st_el=getnal(S,DEM);
		idx=st_el>=max_coel;
		IX=S.IXgrid(idx);
		W=GRIDobj(DEM);
		W.Z(IX)=1;
		W.Z=logical(W.Z);
		SC=STREAMobj(FD,W);
	case 'complete_only'
		SC=removeedgeeffects(S,FD,DEM);
	end
end

function [ChiOBJ]=MakeChiMap(DEM,FD,A,S,theta_ref);

	C=chitransform(S,A,'mn',theta_ref,'a0',1);

	% Make Empty GRIDobj
	ChiOBJ=GRIDobj(DEM);
	ChiOBJ.Z(ChiOBJ.Z==0)=NaN;

	% Populate Grid
	ChiOBJ.Z(S.IXgrid)=C;

end

function [ChiOBJ]=MakeChiGrid(DEM,FD,varargin)
	% Clone of DivideTools 'ChiGrid' function with some options disabled.

	% Parse Inputs
	p = inputParser;         
	p.FunctionName = 'ChiGrid';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'FD', @(x) isa(x,'FLOWobj'));

	addParameter(p,'theta_ref',0.5,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'chi_ref_area',1,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'complete_networks_only',true,@(x) isscalar(x) && islogical(x));
	addParameter(p,'min_elevation',[],@(x) isnumeric(x));

	parse(p,DEM,FD,varargin{:});
	DEM=p.Results.DEM;
	FD=p.Results.FD;

	mn=p.Results.theta_ref;
	a0=p.Results.chi_ref_area;
	me=p.Results.min_elevation;
	cno=p.Results.complete_networks_only;

	if isempty(me)
		abl=false;
	else
		abl=true;
	end

	if cno && ~abl
		% Find nodes influenced by edge (From Topotoolbox blog)
		IXE = GRIDobj(DEM,'logical');
		IXE.Z(:,:) = true;
		IXE.Z(2:end-1,2:end-1) = false;
		IXE = influencemap(FD,IXE);
		% Rest is mine
		% Find drainage basins and all outlets
		[DB,oixi]=drainagebasins(FD);
		% Find where these share pixels other than the edge
		db=DB.Z; db=db(3:end-2,3:end-2); % Move 2 pixels in to be conservative
		ixe=IXE.Z; ixe=ixe(3:end-2,3:end-2); % Move 2 pixels in to be conservative
		dbL=db(ixe);
		% Compile list of drainage basins that are influenced by edge pixels
		dbL=unique(dbL);
		% Index list of outlets based on this
		idxi=ismember(DB.Z(oixi),dbL);
		oixi(idxi)=[];
		% Remove drainage basins based on this
		W=dependencemap(FD,oixi);
		% DEM.Z(~mask.Z)=NaN;
		% % Extract info from FLOWobj		
		% W=~isnan(DEM);
		I=W.Z(FD.ix);
		ix=double(FD.ix(I));
		ixc=double(FD.ixc(I));
		% Recalculate indices
		IX        = zeros(FD.size,'uint32');
		IX(W.Z)   = 1:nnz(W.Z);
		ix      = double(IX(ix));
		ixc     = double(IX(ixc));
		I          = ixc == 0;
		ix(I)    = [];
		ixc(I)   = [];
	elseif cno && abl
		% Recalculate flow directions after removing portions below min elevation
		IX=DEM>me;
		DEM.Z(IX.Z==false)=NaN;
		FD=FLOWobj(DEM,'preprocess','carve');
		% Find nodes influenced by edge (From Topotoolbox blog)
		IXE = GRIDobj(DEM,'logical');
		IXE.Z(:,:) = true;
		IXE.Z(2:end-1,2:end-1) = false;
		IXE = influencemap(FD,IXE);
		% Rest is mine
		% Find drainage basins and all outlets
		[DB,oixi]=drainagebasins(FD);
		% Find where these share pixels other than the edge
		db=DB.Z; db=db(3:end-2,3:end-2); % Move 2 pixels in to be conservative
		ixe=IXE.Z; ixe=ixe(3:end-2,3:end-2); % Move 2 pixels in to be conservative
		dbL=db(ixe);
		% Compile list of drainage basins that are influenced by edge pixels
		dbL=unique(dbL);
		% Index list of outlets based on this
		idxi=ismember(DB.Z(oixi),dbL);
		oixi(idxi)=[];
		% Find offending drainage basins and recalculate indicies
		W=dependencemap(FD,oixi);
		I=W.Z(FD.ix);
		ix=double(FD.ix(I));
		ixc=double(FD.ixc(I));
		% Recalculate indices
		IX        = zeros(FD.size,'uint32');
		IX(W.Z)   = 1:nnz(W.Z);
		ix      = double(IX(ix));
		ixc     = double(IX(ixc));
		I          = ixc == 0;
		ix(I)    = [];
		ixc(I)   = [];
	elseif ~cno && ~abl
		% Extract info from FLOWobj	
		W=~isnan(DEM);
		I=W.Z(FD.ix);
		ix=double(FD.ix(I));
		ixc=double(FD.ixc(I));
		% Recalculate indices
		IX        = zeros(FD.size,'uint32');
		IX(W.Z)   = 1:nnz(W.Z);
		ix      = double(IX(ix));
		ixc     = double(IX(ixc));
		I          = ixc == 0;
		ix(I)    = [];
		ixc(I)   = [];
	elseif ~cno && abl
		W=DEM>me;
		I=W.Z(FD.ix);
		ix=double(FD.ix(I));
		ixc=double(FD.ixc(I));
		% Recalculate indices
		IX        = zeros(FD.size,'uint32');
		IX(W.Z)   = 1:nnz(W.Z);
		ix      = double(IX(ix));
		ixc     = double(IX(ixc));
		I          = ixc == 0;
		ix(I)    = [];
		ixc(I)   = [];
	end

	% Generate coordinate list
	IXgrid=find(W.Z);
	[x,y]=ind2coord(DEM,IXgrid);

	% Distance between two nodes throughout grid
	d = nan(size(x));
	dedge = sqrt((x(ixc)-x(ix)).^2 + (y(ixc)-y(ix)).^2);
	d(:) = 0;
	d(ix) = dedge;

	% Cumulative trapezoidal numerical integration of draiange area grid
	DA=flowacc(FD).*(DEM.cellsize^2);
	da=DA.Z(IXgrid);
	c = zeros(size(da));
	da = ((a0) ./da).^mn;
	for r = numel(ix):-1:1;
	    c(ix(r)) = c(ixc(r)) + (da(ixc(r))+(da(ix(r))-da(ixc(r)))/2)*d(ix(r));
	end

	% Generate and populate full chi grid
	ChiOBJ=GRIDobj(DEM);
	ChiOBJ.Z(ChiOBJ.Z==0)=NaN;
	ChiOBJ.Z(IXgrid)=c;
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
    disp(order_exp);
    Smax=modify(S,'streamorder',order_exp);
    if isempty(Smax)
    error('Smax is empty or invalid!');
    end
    disp('Smax is valid');
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

	% Define non-intersecting segments
	w1=waitbar(0,'Finding network segments');
	[as]=networksegment_slim(DEM,FD,S);
	seg_bnd_ix=as.ix;
	% Precompute values or extract values needed for later
	waitbar(1/4,w1,'Calculating hydrologically conditioned stream elevations');
	z=getnal(S,DEMc);
	zu=getnal(S,DEM);
	z_res=z-zu;
	waitbar(2/4,w1,'Calculating chi values');
	g=gradient(S,DEMc);
	c=chitransform(S,A,'a0',1,'mn',theta_ref);
	d=S.distance;
	da=getnal(S,A.*(A.cellsize^2));
	ixgrid=S.IXgrid;
	waitbar(3/4,w1,'Extracting node ordered list');
	% Extract ordered list of stream indices and find breaks between streams
	s_node_list=S.orderednanlist;
	streams_ix=find(isnan(s_node_list));
	streams_ix=vertcat(1,streams_ix);
	waitbar(1,w1,'Pre computations completed');
	close(w1)
	% Generate empty node attribute list for ksn values
	ksn_nal=zeros(size(d));
	% Begin main loop through channels
	num_streams=numel(streams_ix)-1;
	w1=waitbar(0,'Calculating k_{sn} values - 0% Done');
	seg_count=1;
	for ii=1:num_streams
		% Extract node list for stream of interest
		if ii==1
			snlOI=s_node_list(streams_ix(ii):streams_ix(ii+1)-1);
		else
			snlOI=s_node_list(streams_ix(ii)+1:streams_ix(ii+1)-1);
		end

		% Determine which segments are within this stream
		[~,~,dn]=intersect(snlOI,seg_bnd_ix(:,1));
		[~,~,up]=intersect(snlOI,seg_bnd_ix(:,2));
		seg_ix=intersect(up,dn);

		num_segs=numel(seg_ix);
		dn_up=seg_bnd_ix(seg_ix,:);
		for jj=1:num_segs
			% Find positions within node list
			dnix=find(snlOI==dn_up(jj,1));
			upix=find(snlOI==dn_up(jj,2));
			% Extract segment indices of desired segment
			seg_ix_oi=snlOI(upix:dnix);
			% Extract flow distances and normalize
			dOI=d(seg_ix_oi);
			dnOI=dOI-min(dOI);
			num_bins=ceil(max(dnOI)/segment_length);
			bin_edges=[0:segment_length:num_bins*segment_length];
			% Loop through bins
			for kk=1:num_bins
				idx=dnOI>bin_edges(kk) & dnOI<=bin_edges(kk+1);
				bin_ix=seg_ix_oi(idx);
				cOI=c(bin_ix);
				zOI=z(bin_ix);
					if numel(cOI)>2
						[ksn_val,r2]=Chi_Z_Spline(cOI,zOI);
						ksn_nal(bin_ix)=ksn_val;

						% Build mapstructure
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
	perc_of_total=round((ii/num_streams)*1000)/10;
	if rem(perc_of_total,1)==0
		waitbar((ii/num_streams),w1,['Calculating k_{sn} values - ' num2str(perc_of_total) '% Done']);
	end
	
	end
	close(w1);
end

function seg = networksegment_slim(DEM,FD,S)
	% Slimmed down version of 'networksegment' from main TopoToolbox library that also removes zero and single node length segments

	%% Identify channel heads, confluences, b-confluences and outlets
	Vhead = streampoi(S,'channelheads','logical');  ihead=find(Vhead==1);  IXhead=S.IXgrid(ihead);
	Vconf = streampoi(S,'confluences','logical');   iconf=find(Vconf==1);  IXconf=S.IXgrid(iconf);
	Vout = streampoi(S,'outlets','logical');        iout=find(Vout==1);    IXout=S.IXgrid(iout);
	Vbconf = streampoi(S,'bconfluences','logical'); ibconf=find(Vbconf==1);IXbconf=S.IXgrid(ibconf);

	%% Identify basins associated to b-confluences and outlets
	DB   = drainagebasins(FD,vertcat(IXbconf,IXout));DBhead=DB.Z(IXhead); DBbconf=DB.Z(IXbconf); DBconf=DB.Z(IXconf); DBout=DB.Z(IXout);

	%% Compute flowdistance
	D = flowdistance(FD);

	%% Identify river segments
	% links between channel heads and b-confluences
	[~,ind11,ind12]=intersect(DBbconf,DBhead);
	% links between confluences and b-confluences
	[~,ind21,ind22]=intersect(DBbconf,DBconf);
	% links between channel heads and outlets
	[~,ind31,ind32]=intersect(DBout,DBhead);
	% links between channel heads and outlets
	[~,ind41,ind42]=intersect(DBout,DBconf);
	% Connecting links into segments
	IX(:,1) = [ IXbconf(ind11)' IXbconf(ind21)' IXout(ind31)'  IXout(ind41)'  ];   ix(:,1)= [ ibconf(ind11)' ibconf(ind21)' iout(ind31)'  iout(ind41)'  ];
	IX(:,2) = [ IXhead(ind12)'  IXconf(ind22)'  IXhead(ind32)' IXconf(ind42)' ];   ix(:,2)= [ ihead(ind12)'  iconf(ind22)'  ihead(ind32)' iconf(ind42)' ];

	% Compute segment flow length
	flength=double(abs(D.Z(IX(:,1))-D.Z(IX(:,2))));

	% Remove zero and one node length elements
	idx=flength>=2*DEM.cellsize;
	seg.IX=IX(idx,:);
	seg.ix=ix(idx,:);
	seg.flength=flength(idx);

	% Number of segments
	seg.n=numel(IX(:,1));
end

function [KSN,R2] = Chi_Z_Spline(c,z)

	% Resample chi-elevation relationship using cubic spline interpolation
	[~,minIX]=min(c);
	zb=z(minIX);
	chiF=c-min(c);
	zabsF=z-min(z);
	chiS=linspace(0,max(chiF),numel(chiF)).';
	zS=spline(chiF,zabsF,chiS);

	% Calculate ksn via slope
	KSN= chiS\(zS); % mn not needed because a0 is fixed to 1

	% Calculate R^2
	z_pred=chiF.*KSN;
	sstot=sum((zabsF-mean(zabsF)).^2);
	ssres=sum((zabsF-z_pred).^2);
	R2=1-(ssres/sstot);

end

function [KSNGrid,KSNstdGrid] = KsnAvg(DEM,ksn_ms,radius,er_type)

	% Calculate radius
	radiuspx = ceil(radius/DEM.cellsize);
	SE = strel('disk',radiuspx,0);

	% Record mask of current NaNs
	MASK=isnan(DEM.Z);

	disp('Populating channels'); ts=tic;
	% Make grid with values along channels
	KSNGrid=GRIDobj(DEM);
	KSNGrid.Z(:,:)=NaN;
	for ii=1:numel(ksn_ms)
		ix=coord2ind(DEM,ksn_ms(ii).X,ksn_ms(ii).Y);
		KSNGrid.Z(ix)=ksn_ms(ii).ksn;
	end
	tfi=toc(ts); disp(['Completed in ' num2str(tfi) ' secs'])

	disp('Finding average within radius'); ts=tic;
	% Local mean based on radius
	ISNAN=isnan(KSNGrid.Z);
    [~,L] = bwdist(~ISNAN,'e');
    ksng = KSNGrid.Z(L);           
    FLT   = fspecial('disk',radiuspx);
    ksng   = imfilter(ksng,FLT,'symmetric','same','conv');
	tfi=toc(ts); disp(['Completed in ' num2str(tfi) ' secs'])

    disp('Finding standard deviation within radius'); ts=tic;
    nhood   = getnhood(SE);
    ksnstd   = stdfilt(ksng,nhood); 
	tfi=toc(ts); disp(['Completed in ' num2str(tfi) ' secs'])

    switch er_type
    case 'std_error'
    	disp('Finding number of pixels within radius'); ts=tic;
    	II=~MASK; II=single(II);
    	avg_num=imfilter(II,FLT,'symmetric','same','conv');
    	num_nhood_pix=sum(SE.Neighborhood(:));
    	num_pix=avg_num.*num_nhood_pix;
    	num_pix(num_pix==0)=NaN;
		tfi=toc(ts); disp(['Completed in ' num2str(tfi) ' secs'])

		disp('Converting standard deviation to standard error'); ts=tic;
    	ksnstder=ksnstd./sqrt(num_pix);
    	ksnstder(MASK)=NaN;
		tfi=toc(ts); disp(['Completed in ' num2str(tfi) ' secs'])
    end

    % Set original NaN cells back to NaN
    ksng(MASK)=NaN;
    ksnstd(MASK)=NaN;

    % Output
    KSNGrid.Z=ksng;

    switch er_type
    case 'std'
	    KSNstdGrid=GRIDobj(DEM);
	    KSNstdGrid.Z=ksnstd;
	case 'std_error'
	    KSNstdGrid=GRIDobj(DEM);
	    KSNstdGrid.Z=ksnstder;	
    end	
end
