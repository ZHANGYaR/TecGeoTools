function [knl,ksn_master,bnd_list,kn_list,Sc]=KsnProfiler(DEM,FD,A,S,varargin)
%
	% 用法:
	%	[knl,ksn_master,bnd_list,kn_list,Sc]=KsnProfiler(DEM,FD,A,S);
	%	[knl,ksn_master,bnd_list,kn_list,Sc]=KsnProfiler(DEM,FD,A,S,'参数名',参数值,...);
	%
	% 描述:
	% 	该函数用于交互式选择河道源头并定义计算河道陡峭度值的河段。
	% 	该函数设计类似于Profiler_51的操作，但有一些改进。函数将显示带有河流网络的地图，
	%	并期望用户选择感兴趣的河道源头附近位置。然后会提示用户确认所定义的河流是否为所需选择。
	%	最后，将显示所选河流的chi-z和纵剖面图，用户需要在这些图上（通过鼠标点击）定义任何明显的
	%	具有不同河道陡峭指数（或凹度）的河段（参见'pick_method'选项）。完成选择后按回车键。
	%	系统会询问用户是否希望继续选择河流或已完成选择。当完成河流选择后，函数将输出五种不同的结果（见下文），
	%	并生成所选河流的shapefile，包含ksn、凹度、面积和坡度信息。
	%
	% 选择边界:
	%	该函数期望用户在适当的图上（基于'pick_method'参数选择的红色坐标轴图）使用左键点击来定义不同的ksn河段边界。
	%	如果点击其他鼠标按钮或按任何键（除回车键外），代码将识别为选择了一个裂点位置。这与边界不同，
	%	其位置将被记录在'kn_list'输出中，但不会用于定义不同的ksn河段。这是为了提供用户在不希望用作边界的剖面位置上
	%	"标记"位置的能力。边界和这些标记位置都作为配套函数'ClassifyKnicks'的输入（其中边界和这些裂点/标记位置被区分开来）。
	%	
	% 必需输入:
	%	DEM - 数字高程模型，GRIDobj格式，假设为未平滑DEM（例如来自ProcessRiverBasins的DEMoc或MakeStreams的输出）
	%	FD - 流向，FLOWobj格式
	%	S - 河流网络，STREAMobj格式
	%	A - 流量累积，GRIDobj格式
	%
	%%%%%%%%%%%%%%%%%%
	% 可选输入:
	%
	%%% 重新开始选择
	%	restart [] - 提供此参数允许用户重新开始运行，无论是成功完成但想重新开始的运行，还是由于错误或中途退出的失败运行。
	%			代码运行时，会将重启所需的数据保存在名为'*_restart.mat'的mat文件中。如果代码成功完成，此'*_restart.mat'文件将被删除。
	%			代码运行时或失败后希望恢复运行时，请勿删除此文件。如果使用非'interactive'的'input_method'运行代码且代码成功完成
	%			（即您拟合了通过所选输入方法选择的所有河流且未提前停止代码），则使用restart运行将不会执行任何操作。
	%			有效输入为:
	%		'continue' - 将重新开始运行。如果用于已完成或失败的'interactive'运行，将在地图上重新填充已选择的河流，您可以继续选择。
	%			如果用于失败或中止的非交互式输入方法，将从序列中的下一个河流开始。
	%		'skip' - 仅对非交互式运行有意义的输入。这将跳过序列中的下一个河段。这在特定河段导致代码出错时很有用，
	%			这样您可以在不修改河流网络的情况下跳过该河段。
	%
	%%% 主要选项
	% 	input_method ['interactive'] - 控制如何提供感兴趣的河流的参数:
	%		'interactive' - 用户通过在地图上选择河道源头来选取感兴趣的河流，此选项还会随着用户选择更多河流而迭代构建河道陡峭指数图。
	%		'all_streams' - 将使用提供的STREAMobj并遍历所有河道源头。有一个内部参数避免选择太短而无法正确拟合的河流
	%			（主要与'junction method'设置为'check'相关）。默认值约为DEM像元大小的4倍，用户可以通过为可选参数'min_channel_length'提供输入来更改此值，
	%			输入应为地图单位且大于默认值。您可以使用类似'SegmentPicker'的代码来选择STREAMobj的部分。
	%		'stream_length' - 将使用提供的STREAMobj和'min_length_to_extract'的输入来遍历所有长于提供的长度的河流。
	%			有一个内部参数避免选择太短而无法拟合的河流（主要与'junction method'设置为'check'相关）。
	%			默认值约为DEM像元大小的4倍，用户可以通过为可选参数'min_channel_length'提供输入来更改此值，
	%			输入应为地图单位且大于默认值。
	%		'channel_heads' - 将使用提供的河道源头坐标列表来选择和遍历感兴趣的河流。如果使用此选项，
	%			用户必须为可选参数'channel_head_list'提供输入。
	%	pick_method ['chi'] - 选择如何选取河流段。基于选择的图表将以红色轮廓显示。有效输入为:
	%		'chi' - 在chi-z图上选择河段（推荐且默认）
	%		'stream' - 在纵剖面上选择河段
	%		'slope_area' - 在坡度-面积图上选择河段
	%	junction_method ['check'] - 选择如何处理河流交汇点:
	%		'check' - 每次选择后，将检查所选河流的下游部分是否已被拟合，如果已被拟合，
	%			则不会显示或重新拟合该河流的已拟合部分
	%		'ignore' - 每条河流将从其源头到河口独立显示，无论同一河流网络的某些部分是否已被拟合
	%	concavity_method ['ref']- 凹度选项:
	%		'ref' - 使用参考凹度，用户可以通过参考凹度选项指定此值（见下文）
	%		'auto' - 函数为每条选定的河流找到最佳拟合凹度，如果与'junction_method','check'结合使用，
	%			这意味着选取的短河段将自动拟合可能与同一河流下游部分不同的凹度
	%
	%%% 输入方法选项 
	%	min_channel_length [] - 使用'all_streams'输入方法时考虑的最小河道长度，以地图单位提供。
	%	channel_head_list [] - 河道源头x和y坐标的m×2数组，或河道源头点shapefile的名称/位置，
	%			当使用'channel_heads'输入方法时需要，必须与输入DEM等相同的坐标系。
	%			代码将尝试找到最接近您提供的坐标的河道源头，因此提供的用户坐标越接近河道源头，此选择方法越准确。
	%	min_length_to_extract [] - 如果'input_method'设置为'stream_length'，则提取河流的最小长度（地图单位）。
	%
	%%% 重新定义阈值面积选项
	%	redefine_threshold [false] - 逻辑标志，为每条河流启动一个额外步骤，手动定义山坡-河道过渡
	%			（这将覆盖用于生成提供的STREAMobj的阈值面积，并生成具有可变阈值面积的STREAMobj用于河道定义）。
	%			参见附加可选输入'rd_pick_method'。
	%	rd_pick_method ['slope_area'] - 如果'redefine_threshold'设置为true，用于选择新阈值面积的图表。
	%			有效输入为'slopearea'和'chi'。
	%
	%%% 河流网络修改选项
	%	complete_networks_only [false] - 如果为true，代码将在选择河流之前过滤掉不完整的河流网络部分
	%	min_elev [] - 停止提取河道信息的最低高程（如果留空则不执行操作）
	%	max_area [] - 停止提取河道信息的最大集水面积（以平方地图单位，如果留空则不执行操作）
	%
	%%% 平滑选项
	%	conditioned_DEM [] - 提供平滑的DEM用于此函数的选项（不要为主要必需的DEM输入提供平滑过的DEM！），
	%			将用于提取高程。参见'ConditionDEM'函数以获取制作平滑DEM的选项。
	%			如果未提供输入，代码默认使用mincosthydrocon函数。
	%	interp_value [0.1] - mincosthydrocon中插值参数的值（0到1之间）（如果用户提供平滑过的DEM则不使用）。
	%			接近0的值倾向于更多"切割"，而接近1的值倾向于填充。参见'mincosthydrocon'的信息
	%
	%%% 显示选项
	%	display_slope_area [false] - 显示坡度-面积图的逻辑标志。有些人喜欢坡度-面积图（如支持论文的作者之一），
	%			有些人讨厌坡度-面积图（如支持论文的另一位作者），因此您可以选择完全不绘制它们（false - 默认）
	%			或包含它们（true）。如果您选择'slope_area'作为'pick_method'，这将自动设置为true。
	%	plot_type ['vector'] - 期望'vector'或'grid'，默认为'vector'。控制所有河流是绘制为单独的线（'vector'）
	%			还是将河流网络绘制为网格并下采样（'grid'）。'grid'选项在大型数据集上更快，
	%			但可能导致不准确的河道源头选择。'vector'选项更容易查看，但在大型数据集上加载和交互可能非常慢。	
	%
	%%% 常数
	%	ref_concavity [0.50] - 如果'theta_method'设置为'ref'，则使用参考凹度
	%	smooth_distance [1000] - 转换为shapefile时平滑ksn测量的距离（地图单位）
	%	max_ksn [250] - 用于颜色比例的最大ksn，不会影响实际结果，仅用于显示目的
	%	threshold_area [1e6] - 如果'plot_type'设置为'grid'，用于重绘下采样河流网络的阈值面积
	%
	%%% 输出选项
	%	stack_method ['stack'] - 如果'junction_method'设置为'ignore'，此参数将控制函数在生成shapefile时如何处理河流网络的重叠部分。
	%			有效输入为'stack'（默认）和'average'。如果设置为'stack'，输出shapefile将在网络重叠部分具有多个堆叠的多段线。
	%			这与Profiler51的工作方式类似。如果设置为'average'，函数将在节点基础上平均网络的重叠部分。
	%			注意，如果'junction_method'设置为'check'，则忽略此参数。
	%	shape_name ['ksn'] - 要导出的shapefile的名称，不能有空格以成为ArcGIS的有效名称，且不应包含'.shp'
	%	save_figures [false] - 逻辑标志，保存显示ksn拟合的图形（true）或不保存（false - 默认）	
	%
	%%%%%%%%%%	
	% 输出:
	%	knl - 所选河流段的n×12矩阵节点列表，列依次为x坐标、y坐标、集水面积、ksn、负ksn误差、
	%		正ksn误差、参考凹度、最佳拟合凹度、最小阈值面积、坡度、拟合残差和标识号。注意，
	%		如果在'concavity_method','auto'模式下使用代码，则参考凹度和最佳拟合凹度列将相同。
	%	ksn_master - 与knl相同，但作为单元数组，其中单个单元是单个选定的河道
	%	bnd_list - 用于拟合ksn的所选边界的n×4矩阵，列依次为x坐标、y坐标、高程和河流标识号，
	%		 也作为单独的shapefile输出（'_bounds.shp'）。如果x、y和z值显示为NaN，则表示未为此河流选择边界。
	%		
	%	kn_list - 额外识别的裂点（非边界）的n×4矩阵，列依次为x坐标、y坐标、高程和河流标识号，
	%		 也作为单独的shapefile输出（'_knicks.shp'）。如果x、y和z值显示为NaN，则表示未为此河流选择裂点。
	%	Sc - 所选河流的STREAMobj
	%
	% 示例:
	%	[knl,ksn_master,bnd_list,kn_list,Sc]=KSN_Profiler(DEM,FD,A,S);
	%	[knl,ksn_master,bnd_list,kn_list,Sc]=KSN_Profiler(DEM,FD,A,S,'junction_method','ignore','ref_concavity',0.65,'max_ksn',500);
	%	[knl,ksn_master,bnd_list,kn_list,Sc]=KSN_Profiler(DEM,FD,A,S,'input_method','channel_heads','channel_head_list',channel_head_array);
	%	[knl,ksn_master,bnd_list,kn_list,Sc]=KSN_Profiler(DEM,FD,A,S,'input_method','channel_heads','channel_head_list','channel_heads.shp');
	%	重启示例:
	%	[knl,ksn_master,bnd_list,kn_list,Sc]=KSN_Profiler(DEM,FD,A,S,'restart','continue'); % 从上次停止处继续
	%	[knl,ksn_master,bnd_list,kn_list,Sc]=KSN_Profiler(DEM,FD,A,S,'restart','skip'); % 跳过非交互式序列中的下一个河流
	%	
	%
	% 注意:
	%	-如果未为任何选定的河流选择边界/裂点，则不会生成'_bounds.shp'/'_knicks.shp' shapefile。
	%	-保存的'*_profiler.mat'包含除代码正式输出外的其他文件。这些附加变量是能够使用'restart'选项重新开始运行所必需的。
	%	-如果将'save_figures'设置为true，请勿手动关闭图形，否则会导致代码出错。
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% 函数由Adam M. Forte编写 - 更新日期 : 2020年1月9日 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'KsnProfiler';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
	addRequired(p,'A',@(x) isa(x,'GRIDobj'));
	addRequired(p,'S',@(x) isa(x,'STREAMobj'));

	addParameter(p,'shape_name','ksn',@(x) ischar(x));
	addParameter(p,'smooth_distance',1000,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'concavity_method','ref',@(x) ischar(validatestring(x,{'ref','auto'})));
	addParameter(p,'complete_networks_only',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'pick_method','chi',@(x) ischar(validatestring(x,{'chi','stream','slope_area'})));
	addParameter(p,'junction_method','check',@(x) ischar(validatestring(x,{'check','ignore'})));	
	addParameter(p,'ref_concavity',0.50,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'redefine_threshold',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'rd_pick_method','slope_area',@(x) ischar(validatestring(x,{'chi','slope_area'})));
	addParameter(p,'display_slope_area',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'max_ksn',250,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'min_elev',[],@(x) isscalar(x) && isnumeric(x) || isempty(x));
	addParameter(p,'max_area',[],@(x) isscalar(x) && isnumeric(x) || isempty(x));
	addParameter(p,'plot_type','vector',@(x) ischar(validatestring(x,{'vector','grid'})));
	addParameter(p,'threshold_area',1e6,@(x) isnumeric(x));
	addParameter(p,'input_method','interactive',@(x) ischar(validatestring(x,{'interactive','channel_heads','all_streams','stream_length'})));
	addParameter(p,'channel_head_list',[],@(x) isnumeric(x) && size(x,2)==2 || regexp(x,regexptranslate('wildcard','*.shp')) || isempty(x));
	addParameter(p,'min_length_to_extract',[],@(x) isnumeric(x) && isscalar(x) || isempty(x));
	addParameter(p,'min_channel_length',[],@(x) isnumeric(x) && isscalar(x) || isempty(x));
	addParameter(p,'conditioned_DEM',[],@(x) isa(x,'GRIDobj') || isempty(x));
	addParameter(p,'interp_value',0.1,@(x) isnumeric(x) && x>=0 && x<=1);
	addParameter(p,'save_figures',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'restart',[],@(x) ischar(validatestring(x,{'continue','skip'})) || isempty(x));
	addParameter(p,'stack_method','stack',@(x) ischar(validatestring(x,{'average','stack'})));
	addParameter(p,'restart_loc',[],@(x) ischar(x) || isempty(x));

	parse(p,DEM,FD,A,S,varargin{:});
	DEM=p.Results.DEM;
	FD=p.Results.FD;
	S=p.Results.S;
	A=p.Results.A;

	shape_name=p.Results.shape_name;
	smooth_distance=p.Results.smooth_distance;
	theta_method=p.Results.concavity_method;
	cno=p.Results.complete_networks_only;
	pick_method=p.Results.pick_method;
	junction_method=p.Results.junction_method;
	ref_theta=p.Results.ref_concavity;
	display_slope_area=p.Results.display_slope_area;
	plot_type=p.Results.plot_type;
	threshold_area=p.Results.threshold_area;
	input_method=p.Results.input_method;
	chl=p.Results.channel_head_list;
	mlte=p.Results.min_length_to_extract;
	min_channel_length=p.Results.min_channel_length;
	iv=p.Results.interp_value;
	DEMc=p.Results.conditioned_DEM;
	min_elev=p.Results.min_elev;
	max_area=p.Results.max_area;
	save_figures=p.Results.save_figures;
	redefine_thresh=p.Results.redefine_threshold;
	rd_pick_method=p.Results.rd_pick_method;
	mksn=p.Results.max_ksn;
	restart=p.Results.restart;
	stack_method=p.Results.stack_method;
	restart_loc=p.Results.restart_loc; % Hidden parameter for GUI deployed versions

wtb=waitbar(0,'准备输入中...');

	% 设置重启标志
	if ~isempty(restart)
		rf=true;
	else
		rf=false;
	end

	% 检查是否启用了重启并替换参数
	if rf
		if ~isempty(restart_loc)
			restart_file=dir('*_restart.mat');
			out_file=dir('*_profiler.mat');
		else
			restart_file=dir(fullfile(restart_loc,'*_restart.mat'));
			out_file=dir(fullfile(restart_loc,'*_profiler.mat'));			
		end

		if numel(restart_file)>1
			if isdeployed 
				errordlg('发现多个重启文件，请从当前目录或搜索路径中删除非目标重启文件')
			end
			error('发现多个重启文件，请从当前目录或搜索路径中删除非目标重启文件');
		elseif isempty(restart_file) && ~isempty(out_file)		
			if numel(out_file)>1
				if isdeployed
					errordlg('发现多个profiler输出mat文件，请从当前目录或搜索路径中删除非目标文件')
				end
				error('发现多个profiler输出mat文件，请从当前目录或搜索路径中删除非目标文件');
			end
			load(out_file(1,1).name,'input_params');
			r_type='c';
		elseif isempty(restart_file) && isempty(out_file)
			if isdeployed
				errordlg('未找到先前运行文件，如果是首次运行KsnProfiler，请不要使用"restart"参数')
			end
			error('未找到先前运行文件，如果是首次运行KsnProfiler，请不要使用"restart"参数')
		else
			load(restart_file(1,1).name,'input_params');
			r_type='r';
		end

		% Load in parameters from previous run
		shape_name=input_params.shape_name;
		smooth_distance=input_params.smooth_distance;
		theta_method=input_params.concavity_method;
		cno=input_params.complete_networks_only;
		pick_method=input_params.pick_method;
		junction_method=input_params.junction_method;
		ref_theta=input_params.ref_concavity;
		display_slope_area=input_params.display_slope_area;
		plot_type=input_params.plot_type;
		threshold_area=input_params.threshold_area;
		input_method=input_params.input_method;
		chl=input_params.channel_head_list;
		mlte=input_params.min_length_to_extract;
		min_channel_length=input_params.min_channel_length;
		iv=input_params.interp_value;
		DEMc=input_params.conditioned_DEM;
		min_elev=input_params.min_elev;
		max_area=input_params.max_area;
		save_figures=input_params.save_figures;
		redefine_thresh=input_params.redefine_threshold;
		rd_pick_method=input_params.rd_pick_method;
		mksn=input_params.max_ksn;	
	end	

	waitbar(1/4,wtb,'保存参数中...');

	% 保存参数到最终文件和重启文件
	out_mat_name=[shape_name '_profiler.mat'];
	out_restart_name=[shape_name '_restart.mat'];
	if rf
		save(out_mat_name,'input_params','-append');
		if exist(out_restart_name)==2
			save(out_restart_name,'input_params','-append');
		else
			save(out_restart_name,'input_params','-v7.3');
		end
	else
		input_params=p.Results;
		save(out_mat_name,'input_params','-v7.3');
		save(out_restart_name,'input_params','-v7.3');
	end

      % 如果启用标志则移除边缘效应
    if cno
    	S=removeedgeeffects(S,FD,DEM);
    end

	% 查找河道源头
	[ch]=streampoi(S,'channelheads','xy');

	% 创建主KSN色标
	KSN_col=ksncolor(100);

	waitbar(2/4,wtb,'进行参数检查...');

	% 执行参数检查并重新赋值
	if strcmp(input_method,'channel_heads')
		if isempty(chl)
			if isdeployed
				errordlg('选择方法为"channel_heads"，必须为"channel_head_list"参数提供输入')
			end
			error('选择方法为"channel_heads"，必须为"channel_head_list"参数提供输入');
		end

		if ischar(chl)
			% 加载shapefile时的提示
			if logical(regexp(chl,regexptranslate('wildcard','*.shp')))
				ch_ms=shaperead(chl);
				ch_t=struct2table(ch_ms);
				if ~strcmp(ch_t.Geometry(1),'Point')
					if isdeployed
						errordlg('提供的河道源头shapefile不是点类型文件')
					end
					error('提供的河道源头shapefile不是点类型文件');
				end
				xi=ch_t.X;
				yi=ch_t.Y;
				chl=[xi yi];
			end
		end

% 坐标匹配提示
		dists=pdist2(chl,ch);
		[~,s_ch_ix]=min(dists,[],2);
		s_ch=ch(s_ch_ix,:);
		num_ch=size(s_ch,1);

		plot_type='none';
		input_method='preselected';
	elseif strcmp(input_method,'all_streams')
		s_ch=ch;
		num_ch=size(s_ch,1);

		if isempty(min_channel_length)
			min_channel_length=sqrt((4*DEM.cellsize)^2 + (4*DEM.cellsize)^2);
		end

		plot_type='none';
		input_method='preselected';
	elseif strcmp(input_method,'stream_length')
		if isempty(mlte)
			if isdeployed
				errordlg('选择方法为"stream_length"，必须为"mean_length_to_extract"参数提供输入')
			end
			error('选择方法为"stream_length"，必须为"mean_length_to_extract"参数提供输入');
		end

		FlowD=flowdistance(FD);
		ix=coord2ind(DEM,ch(:,1),ch(:,2));
		idx=FlowD.Z(ix)>=mlte;

		s_ch=ch(idx,:);

			if isempty(s_ch)
			if isdeployed
				errordlg('"mean_length_to_extract"参数设置导致无可用河道，请减小该值后重试')
			end
			error('"mean_length_to_extract"参数设置导致无可用河道，请减小该值后重试')
		end

		num_ch=size(s_ch,1);

		if isempty(min_channel_length)
			min_channel_length=sqrt((4*DEM.cellsize)^2 + (4*DEM.cellsize)^2);
		end
		plot_type='none';
		input_method='preselected';
	end


	% 参数冲突提示
	if ~isempty(min_elev) && ~isempty(max_area)
		if isdeployed
			errordlg('不能同时设置"min_elev"和"max_area"参数，请只设置其中一个')
		end
		error('不能同时设置"min_elev"和"max_area"参数，请只设置其中一个')
	end

	% DEM处理提示
if isempty(DEMc)
		zc=mincosthydrocon(S,DEM,'interp',iv);
		DEMc=GRIDobj(DEM);
		DEMc.Z(DEMc.Z==0)=NaN;
		DEMc.Z(S.IXgrid)=zc;
	elseif ~isempty(DEMc) && redefine_thresh
		if isdeployed
			warndlg(['同时提供处理过的DEM和重定义流域面积阈值可能导致意外结果，' ...
				'建议使用极低阈值面积生成的DEM'])
		else
			warning(['同时提供处理过的DEM和重定义流域面积阈值可能导致意外结果，' ...
				'建议使用极低阈值面积生成的DEM'])
		end
 end

	% Make gradient
	G=gradient8(DEMc);
	% Make drainage area
	DA=A.*A.cellsize^2;

	waitbar(3/4,wtb,'生成地形梯度...');

% 根据最小高程或最大流域面积选项修改提供的河流网络
	if ~isempty(min_elev)
		zel=getnal(S,DEMc);
		idx=zel>min_elev;
		new_ix=S.IXgrid(idx);
		W=GRIDobj(DEMc,'logical');
		W.Z(new_ix)=true;
		S=STREAMobj(FD,W);
	elseif ~isempty(max_area)
		DA=A.*(A.cellsize^2);
		zda=getnal(S,DA);
		idx=zda<max_area;
		new_ix=S.IXgrid(idx);
		W=GRIDobj(DEMc,'logical');
		W.Z(new_ix)=true;
		S=STREAMobj(FD,W);		
	end

	% 如果重新定义阈值则生成流距网格
	if redefine_thresh
		FLUS=flowdistance(FD);
	end

	% 生成交互式选择的地图图形
	switch plot_type
	case 'grid'	

		disp('下采样数据以用于显示目的');
		% 重新计算流向	
		DEMr=resample(DEM,DEM.cellsize*4);
		FDr=FLOWobj(DEMr,'preprocess','carve');
		% 真实出水口
		out_T_xy=streampoi(S,'outlets','xy');
		% 下采样后的总河流网络
		Sr_temp=STREAMobj(FDr,'minarea',threshold_area,'unit','mapunits');
		out_D_xy=streampoi(Sr_temp,'outlets','xy');
		out_D_ix=streampoi(Sr_temp,'outlets','ix');
		% 检查出水口列表是否不同
		dists=pdist2(out_T_xy,out_D_xy);
		[~,s_out_ix]=min(dists,[],2);
		out_D_ix=out_D_ix(s_out_ix);
		% 重建下采样网络
		Sr=STREAMobj(FDr,'minarea',threshold_area,'unit','mapunits','outlets',out_D_ix);
		% 转换为栅格
		SG=STREAMobj2GRIDobj(Sr);

		% 初始化地图图形
		f1=figure(1);
		set(f1,'Visible','off');
		hold on
		[RGB]=imageschs(DEMr,SG,'colormap','gray');
		[~,R]=GRIDobj2im(DEMr);
		imshow(flipud(RGB),R);
		axis xy
		colormap(KSN_col);
		caxis([0 mksn])
		c1=colorbar;
		ylabel(c1,'Channel Steepness')
		if ~verLessThan('matlab','9.5')
	        disableDefaultInteractivity(gca);
	    end	
		hold off
		set(f1,'Visible','on','Units','normalized','Position',[0.05 0.1 0.5 0.5],'renderer','painters');

	case 'vector'

		% 初始化矢量地图图形
		f1=figure(1);
		set(f1,'Visible','off');

		[RGB]=imageschs(DEM,DEM,'colormap','gray');
		[~,R]=GRIDobj2im(DEM);

		imshow(flipud(RGB),R);
		axis xy
		hold on
		colormap(KSN_col);
		plot(S,'-w');
		caxis([0 mksn])
		c1=colorbar;
		ylabel(c1,'Channel Steepness')
		if ~verLessThan('matlab','9.5')
	        disableDefaultInteractivity(gca);
	    end			
		hold off
		set(f1,'Visible','on','Units','normalized','Position',[0.05 0.1 0.5 0.5],'renderer','painters');	
	end

	waitbar(1,wtb);
	close(wtb);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% 图形化选择与河道源头列表的主切换 %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	switch input_method
	case 'interactive'
		% 初始化计数器和循环变量
		str1='N';
		str2='Y';

		if rf && strcmp(r_type,'c')
			% 从已完成运行中加载数据
			load(out_mat_name,'count','ksn_master','bnd_master','kn_master','res_master','Sc');
			ii=count+1;
			% 重新生成已绘制的河流
			km=vertcat(ksn_master{:});
			kix=coord2ind(DEM,km(:,1),km(:,2));
			K=GRIDobj(DEM);
			K.Z(kix)=km(:,4);
			figure(1)
			hold on
			plotc(Sc,K);
			hold off
			clear km kix K;
		elseif rf && strcmp(r_type,'r')
			% 从失败或中止的运行中加载数据
			load(out_restart_name,'count','ksn_master','bnd_master','kn_master','res_master','Sc');
			ii=count+1;
			% 重新生成已绘制的河流
			km=vertcat(ksn_master{:});
			kix=coord2ind(DEM,km(:,1),km(:,2));
			K=GRIDobj(DEM);
			K.Z(kix)=km(:,4);
			figure(1)
			hold on
			plotc(Sc,K);
			hold off
			clear km kix K;
		else
			ii=1;
		end

		if strcmp(theta_method,'ref')
			% 自动计算ksn用于对比
			[auto_ksn]=KSN_Quick(DEM,A,S,ref_theta);
		end

		% 开始选择河流
		while strcmpi(str2,'Y')
           			% 选择要拟合的河道
			while strcmpi(str1,'N')	
				str3='R'; % 重置重做标志

	            figure(1)
	            hold on
	            title('缩放或平移至感兴趣区域，然后按回车');
	            hold off
				pause()

	            figure(1)
	            hold on
	            title('在目标河道源头附近选择点');
	            hold off

				[x,y]=ginput(1);
				pOI=[x y];

				% 查找最近的河道源头
				distance=sqrt(sum(bsxfun(@minus, ch, pOI).^2,2));
				chOI=ch(distance==min(distance),:);

				% 构建逻辑栅格
				ix=coord2ind(DEM,chOI(:,1),chOI(:,2));
				IX=GRIDobj(DEM,'logical');
				IX.Z(ix)=true;

				if redefine_thresh

					% 提取目标河流
					Sn=modify(S,'downstreamto',IX);

					figure(f1)
					hold on
					p1=plot(Sn,'-b','LineWidth',2);
					hold off

					qa1=questdlg('这是您想要的河段吗？','河流选择','是','否','是');
					switch qa1
					case '是'

						delete(p1);

						[Sn]=RedefineThreshold(DEM,FD,A,Sn,FLUS,ref_theta,rd_pick_method,smooth_distance,ii,save_figures,shape_name);
						% 更新DEMc
						if any(isnan(getnal(Sn,DEMc)))
							zc=mincosthydrocon(Sn,DEM,'interp',iv);
							DEMc.Z(Sn.IXgrid)=zc;
						end
						% 重新计算自动ksn
						if strcmp(theta_method,'ref')
							[auto_ksn]=KSN_Quick(DEM,A,Sn,ref_theta);
						end

						str1='Y';

						if strcmp(junction_method,'check')
							if ii>1
								[IIXX,~,~,Si]=intersectlocs(Sc,Sn);
								if isempty(IIXX)
									Sn=Sn;
									Sct=union(Sn,Sc,FD);
								else
									Sn=Si;
									Sct=union(Sn,Sc,FD);
								end
							else
								Sct=Sn;
							end
						elseif strcmp(junction_method,'ignore')					
							if ii>1
								Sct=union(Sn,Sc,FD);
							else
								Sct=Sn;
							end
						end

						% 绘制更新后的河流
						figure(f1)
						hold on
						p1=plot(Sn,'-b','LineWidth',2);
						hold off

						Sc=Sct;
					case '否'
						str1='N';
						delete(p1);
					end
				else
					% 提取目标河流
					Sn=modify(S,'downstreamto',IX);

					% 构建已选河流的复合网络
					if strcmp(junction_method,'check')
						if ii>1
							[IIXX,~,~,Si]=intersectlocs(Sc,Sn);
							if isempty(IIXX)
								Sn=Sn;
								Sct=union(Sn,Sc,FD);
							else
								Sn=Si;
								Sct=union(Sn,Sc,FD);
							end
						else
							Sct=Sn;
						end
					elseif strcmp(junction_method,'ignore')					
						if ii>1
							Sct=union(Sn,Sc,FD);
						else
							Sct=Sn;
						end
					end

					figure(f1)
					hold on
					p1=plot(Sn,'-b','LineWidth',2);
					hold off

					qa1=questdlg('这是您想要的河段吗？','河流选择','是','否','是');
					switch qa1
					case '是'
						str1='Y';
						Sc=Sct;
					case '否'
						str1='N';
						delete(p1);
					end					
				end
			end % 结束单河道选择

			%% 提取阈值流域面积
			snchix=streampoi(Sn,'channelheads','ix');
			snda=DA.Z(snchix);

			%% 计算Chi并提取ksn数据
			if strcmp(theta_method,'ref')
				C=ChiCalc(Sn,DEMc,A,1,ref_theta);
				ak=getnal(Sn,auto_ksn);
			elseif strcmp(theta_method,'auto')
				C=ChiCalc(Sn,DEMc,A,1);
				if redefine_thresh
					[auto_ksn]=KSN_Quick(DEM,A,Sn,C.mn);
				else
					[auto_ksn]=KSN_Quick(DEM,A,S,C.mn);
				end
				ak=getnal(Sn,auto_ksn);
			end

			%% 分箱数据
			[DAvg,KsnAvg]=BinAverage(Sn.distance,ak,smooth_distance);
			[~,CAvg]=BinAverage(C.distance,C.chi,smooth_distance);

			%% 开始拟合循环
			while strcmpi(str3,'R')
				%% 初始化边界选择图形
				f2=figure(2);
				set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');
				clf	

				%% 不同选择方法的主切换
				if strcmp(pick_method,'chi')

					if display_slope_area
						[bs,ba,bc,bd,aa,ag,ad,ac]=sa(DEMc,Sn,A,C.chi,smooth_distance);

						ax4=subplot(4,1,4);
						hold on
						scatter(aa,ag,5,ac,'+');
						scatter(ba,bs,20,bc,'filled','MarkerEdgeColor','k');
						xlabel('对数面积');
						ylabel('对数坡度');
						title('坡度-面积图');
						set(ax4,'YScale','log','XScale','log','XDir','reverse');
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax4);
					    end							
						hold off

						ax3=subplot(4,1,3);
						hold on
						pl1=plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
						pl2=plotdz(Sn,DEMc,'dunit','km','Color','k');
						pl3=scatter((C.distance)./1000,C.elev,5,C.chi,'filled');
						xlabel('距河口距离 (km)')
						ylabel('高程 (m)')
						legend([pl1 pl2 pl3],'未平滑DEM','平滑DEM','\chi','location','best');
						title('纵剖面图')
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax3);
					    end							
						hold off

						ax2=subplot(4,1,2);
						hold on
						scatter(CAvg,KsnAvg,20,CAvg,'filled','MarkerEdgeColor','k');
						xlabel('\chi')
						ylabel('k_{sn}');
						title('\chi - k_{sn}');
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax2);
					    end							
						hold off

						ax1=subplot(4,1,1);
						hold on
						plot(C.chi,C.elev,'-k');
						scatter(C.chi,C.elev,10,C.chi,'filled');
						xlabel('\chi')
						ylabel('高程 (m)')
						title(['\chi - Z : 凹度 = ' num2str(C.mn) ' : 选择河段 - 完成后按回车'],'Color','r')
						ax1.XColor='Red';
						ax1.YColor='Red';
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax1);
					    end							
						hold off

						linkaxes([ax1,ax2],'x');
						colormap(ax1,'jet'); colormap(ax2,'jet'); colormap(ax3,'jet'); colormap(ax4,'jet');

					else
						ax3=subplot(3,1,3);
						hold on
						pl1=plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
						pl2=plotdz(Sn,DEMc,'dunit','km','Color','k');
						pl3=scatter((C.distance)./1000,C.elev,5,C.chi,'filled');
						xlabel('距河口距离 (km)')
						ylabel('高程 (m)')
						legend([pl1 pl2 pl3],'未平滑DEM','平滑DEM','\chi','location','best');
						title('纵剖面图')
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax3);
					    end							
						hold off

						ax2=subplot(3,1,2);
						hold on
						scatter(CAvg,KsnAvg,20,CAvg,'filled','MarkerEdgeColor','k');
						xlabel('\chi')
						ylabel('k_{sn}');
						title('\chi - k_{sn}');
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax2);
					    end							
						hold off

						ax1=subplot(3,1,1);
						hold on
						plot(C.chi,C.elev,'-k');
						scatter(C.chi,C.elev,10,C.chi,'filled');
						xlabel('\chi')
						ylabel('高程 (m)')
						title(['\chi - Z : 凹度 = ' num2str(C.mn) ' : 选择河段 - 完成后按回车'],'Color','r')
						ax1.XColor='Red';
						ax1.YColor='Red';
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax1);
					    end							
						hold off

						linkaxes([ax1,ax2],'x');
						colormap(ax1,'jet'); colormap(ax2,'jet'); colormap(ax3,'jet');
					end	

					[cv,~,bttn]=ginput;
					% 检测非边界裂点
					bttn_idx=bttn~=1;
					if any(bttn_idx)
						% 分离非边界点
						cv_kn=cv(bttn_idx);
						cv(bttn_idx)=[];
						% 转换为索引
						rc=C.chi; rx=C.x; ry=C.y;
						kn_ix=zeros(numel(cv_kn),1); 
						for jj=1:numel(cv_kn)
							chidist=sqrt(sum(bsxfun(@minus, rc, cv_kn(jj)).^2,2));
							[~,knbix]=min(chidist);
							knbx=rx(knbix);
							knby=ry(knbix);
							kn_ix(jj)=coord2ind(DEM,knbx,knby);
						end
					else
						kn_ix=NaN;
					end

					if isempty(cv)
						if strcmp(theta_method,'ref')
							Cbf=ChiCalc(Sn,DEMc,A,1);
							ksn_list=[C.x C.y C.area ones(numel(C.x),1)*C.ks ones(numel(C.x),1)*C.ks_neg ones(numel(C.x),1)*C.ks_pos ...
								ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*Cbf.mn ones(numel(C.x),1)*snda];
						elseif strcmp(theta_method,'auto')
							ksn_list=[C.x C.y C.area ones(numel(C.x),1)*C.ks ones(numel(C.x),1)*C.ks_neg ones(numel(C.x),1)*C.ks_pos ...
								ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*snda];
						end

						% 确定ksn值在色阶中的位置并绘图
						ksn_val=C.ks;

						if ksn_val > mksn
							figure(f1)
							hold on
							plot(Sn,'Color',KSN_col(end,:),'LineWidth',2);
							hold off
						else
							edges=linspace(0,mksn,10);
							n=histc(ksn_val,edges);
							figure(f1)
							hold on
							plot(Sn,'Color',KSN_col(logical(n),:),'LineWidth',2);
							hold off	
						end

						[~,lbix]=min(C.chi);
						elbl=C.elev(lbix);

						if display_slope_area
							figure(f2)
							subplot(4,1,1);
							hold on
							plot(C.chi,C.pred+elbl,'-k','LineWidth',2);
							hold off

							subplot(4,1,3);
							hold on
							pl4=plot((C.distance)/1000,C.pred+elbl,'-k','LineWidth',2);
							legend([pl1 pl2 pl3 pl4],'未平滑DEM','平滑DEM','\chi','段拟合','location','best');
							hold off

							subplot(4,1,4);
							hold on
							plot(C.area,ksn_val.*C.area.^(-C.mn),'-k','LineWidth',2);
							hold off

						else
							figure(f2)
							subplot(3,1,1);
							hold on
							plot(C.chi,C.pred+elbl,'-k','LineWidth',2);
							hold off

							subplot(3,1,3);
							hold on
							pl4=plot((C.distance)/1000,C.pred+elbl,'-k','LineWidth',2);
							legend([pl1 pl2 pl3 pl4],'未平滑DEM','平滑DEM','\chi','段拟合','location','best');
							hold off
						end

						res_list=[C.chi C.res];
						bnd_ix=NaN;

					else
						% 排序裂点列表并构建边界列表
						cvs=sortrows(cv);
						bnds=vertcat(0,cvs,C.chi(1));

						num_bnds=numel(bnds);
						rc=C.chi;
						rx=C.x;
						ry=C.y;
						rd=C.distance;
						for jj=1:num_bnds-1
							% 提取边界
							lb=bnds(jj);
							rb=bnds(jj+1);

							% 裁剪河段
							lb_chidist=sqrt(sum(bsxfun(@minus, rc, lb).^2,2));
							rb_chidist=sqrt(sum(bsxfun(@minus, rc, rb).^2,2));

							[~,lbix]=min(lb_chidist);
							[~,rbix]=min(rb_chidist);

							lbx=rx(lbix);
							lby=ry(lbix);

							rbx=rx(rbix);
							rby=ry(rbix);

							lix=coord2ind(DEM,lbx,lby);
							rix=coord2ind(DEM,rbx,rby);

							Seg=modify(Sn,'downstreamto',rix);
							Seg=modify(Seg,'upstreamto',lix);

							% 重建包含下游节点的河段
							WSEG=GRIDobj(DEM,'logical');
							WSEG.Z(Seg.IXgrid)=true;
							WSEG.Z(lix)=true;
							% 如果是河段末端则添加上游节点
							if jj==num_bnds-1
								WSEG.Z(rix)=true;
							end
							Seg=STREAMobj(FD,WSEG);

							% 检查河段长度，小于2个节点时向下游扩展
							lbix_new=lbix+1;
							first_time=true;
							while numel(Seg.IXgrid)<=2
								if first_time
									wrn_mssg=['第' num2str(jj) '段选择的河段太短，已向下游扩展边界'];
									wd=warndlg(wrn_mssg);
									uiwait(wd);
								end
								lix_new=coord2ind(DEM,rx(lbix_new),ry(lbix_new));
								WSEG.Z(lix_new)=true;
								Seg=STREAMobj(FD,WSEG);
								lbix_new=lbix_new+1;
								first_time=false;
							end

							% 构建边界列表
							if jj<num_bnds-1
								bnd_ix(jj,1)=rix;
							end

							% 计算chi获取ksn和最佳凹度
							if strcmp(theta_method,'ref')
								Cseg=ChiCalc(Seg,DEMc,A,1,ref_theta);
							elseif strcmp(theta_method,'auto')
								Cseg=ChiCalc(Seg,DEMc,A,1,C.mn);
							end
							Cbfseg=ChiCalc(Seg,DEMc,A,1);								

							% 确定ksn值在色阶中的位置并绘图
							ksn_val=Cseg.ks;

						if ksn_val > mksn
								figure(f1)
								hold on
								plot(Seg,'Color',KSN_col(end,:),'LineWidth',2);
								hold off
							else
								edges=linspace(0,mksn,10);
								n=histc(ksn_val,edges);
								figure(f1)
								hold on
								plot(Seg,'Color',KSN_col(logical(n),:),'LineWidth',2);
								hold off	
                         end

ksn_nodes{jj,1}=[Cseg.x Cseg.y Cseg.area ones(numel(Cseg.x),1)*Cseg.ks ones(numel(Cseg.x),1)*Cseg.ks_neg ones(numel(Cseg.x),1)*Cseg.ks_pos ...
								ones(numel(Cseg.x),1)*Cseg.mn ones(numel(Cseg.x),1)*Cbfseg.mn ones(numel(Cseg.x),1)*snda];

							% 绘制线性拟合结果	
							rchi=rc(rb_chidist==min(rb_chidist));
							lchi=rc(lb_chidist==min(lb_chidist));
							segChi=linspace(lchi,rchi,numel(Cseg.chi));
							seg0Chi=linspace(0,max(Cseg.chi),numel(Cseg.chi));
							[~,lbix]=min(Cseg.chi);
							elbl=Cseg.elev(lbix);

							if display_slope_area
								figure(f2)
								subplot(4,1,1);
								hold on
								plot(segChi,(seg0Chi.*ksn_val)+elbl,'-k','LineWidth',2);
								hold off

                                seg_st=rd(lb_chidist==min(lb_chidist));
								subplot(4,1,3);
								hold on
								pl4=plot((Cseg.distance+seg_st)/1000,(Cseg.pred)+elbl,'-k','LineWidth',2);
								legend([pl1 pl2 pl3 pl4],'未平滑DEM','平滑DEM','\chi','段拟合','location','best');
								hold off

								subplot(4,1,4);
								hold on
								plot(Cseg.area,ksn_val.*Cseg.area.^(-1*Cseg.mn),'-k','LineWidth',2);
								hold off
							else								
								figure(f2)
								subplot(3,1,1);
								hold on
								plot(segChi,(seg0Chi.*ksn_val)+elbl,'-k','LineWidth',2);
								hold off

								seg_st=rd(lb_chidist==min(lb_chidist));
								subplot(3,1,3);
								hold on
								pl4=plot((Cseg.distance+seg_st)/1000,(Cseg.pred)+elbl,'-k','LineWidth',2);
								legend([pl1 pl2 pl3 pl4],'未平滑DEM','平滑DEM','\chi','段拟合','location','best');
								hold off
							end

							res_list{jj,1}=[Cseg.chi+lchi Cseg.res];	

						end

						ksn_list=vertcat(ksn_nodes{:});
						res_list=vertcat(res_list{:});
					end
	
					%% 绘制结果图
					if display_slope_area
						figure(f2)
						subplot(4,1,2)
						hold on
						plot(res_list(:,1),ksn_list(:,4),'-k','LineWidth',2);
						hold off
					else
						figure(f2)
						subplot(3,1,2)
						hold on
						plot(res_list(:,1),ksn_list(:,4),'-k','LineWidth',2);
						hold off
					end

					f3=figure(3);
					set(f3,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');

					clf
					sbplt2=subplot(2,1,2);
					hold on
					s1=scatter(CAvg,KsnAvg,20,'k','filled');
					p1=plot(res_list(:,1),ksn_list(:,4)-ksn_list(:,5),':k');
					plot(res_list(:,1),ksn_list(:,4)+ksn_list(:,6),':k');
					p2=plot(res_list(:,1),ksn_list(:,4),'-k','LineWidth',2);
					[ksn_vals,~,ksn_ix]=unique(ksn_list(:,4));
					chi_means=accumarray(ksn_ix,res_list(:,1),[],@(x) mean(x,'omitnan'));
					for kk=1:numel(ksn_vals)
						text(chi_means(kk),ksn_vals(kk),['k_{sn} = ' num2str(ksn_vals(kk))],...
							'VerticalAlignment','bottom','HorizontalAlignment','center');
					end
					xlabel('\chi')
					ylabel('k_{sn}')
					legend([s1 p1 p2],{'k_{sn}','k_{sn}误差范围','拟合段k_{sn}'},'location','best');
					title('k_{sn} - \chi 关系图')
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(sbplt2);
				    end						
					hold off

					sbplt1=subplot(2,1,1);
					hold on
					plot([min(res_list(:,1)) max(res_list(:,1))],[0 0],'-k');
					scatter(res_list(:,1),res_list(:,2),10,'k','filled');
					xlabel('\chi')
					ylabel('残差 (m)')
					title('k_{sn}拟合残差分布')
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(sbplt1);
				    end						
					hold off				


					elseif strcmp(pick_method,'stream')	

					if display_slope_area
						[bs,ba,bc,bd,bk,aa,ag,ad,ac]=sa_ksn(DEMc,Sn,A,C.chi,ak,smooth_distance);	

						ax4=subplot(4,1,4);
						hold on
						scatter(aa,ag,5,ad./1000,'+');
						scatter(ba,bs,20,bd./1000,'filled','MarkerEdgeColor','k');
						xlabel('对数面积');
						ylabel('对数坡度');
						title('坡度-面积图');
						set(ax4,'XScale','log','YScale','log','XDir','reverse');
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax4);
					    end							
						hold off

						ax1=subplot(4,1,1);
						hold on
						plot(C.chi,C.elev,'Color',[0.5 0.5 0.5]);
						scatter(C.chi,C.elev,10,(C.distance)./1000,'filled');
						xlabel('\chi')
						ylabel('高程 (m)')
						title(['\chi - Z : 凹度 = ' num2str(C.mn)])
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax1);
					    end							
						hold off

						ax2=subplot(4,1,2);
						hold on
						scatter(DAvg./1000,KsnAvg,20,DAvg./1000,'filled','MarkerEdgeColor','k');
						xlabel('距离 (km)')
						ylabel('k_{sn}');
						title('距离-k_{sn}关系');
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax2);
					    end							
						hold off

						ax3=subplot(4,1,3);
						hold on
						pl1=plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
						pl2=plotdz(Sn,DEMc,'dunit','km','Color','k');
						pl3=scatter((C.distance)./1000,C.elev,5,(C.distance)./1000,'filled');
						xlabel('距河口距离 (km)')
						ylabel('高程 (m)')
						legend([pl1 pl2 pl3],'未平滑DEM','平滑DEM','河道距离','location','best');
						title('纵剖面图 : 选择河段 - 完成后按回车','Color','r')
						ax3.XColor='Red';
						ax3.YColor='Red';
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax3);
					    end							
						hold off

						linkaxes([ax2,ax3],'x');
						colormap(ax1,'jet'); colormap(ax2,'jet'); colormap(ax3,'jet'); colormap(ax4,'jet')
                	else
						ax1=subplot(3,1,1);
						hold on
						plot(C.chi,C.elev,'Color',[0.5 0.5 0.5]);
						scatter(C.chi,C.elev,10,(C.distance)./1000,'filled');
						xlabel('\chi')
						ylabel('高程 (m)')
						title(['\chi - Z : 凹度 = ' num2str(C.mn)])
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax1);
					    end							
						hold off

						ax2=subplot(3,1,2);
						hold on
						scatter(DAvg./1000,KsnAvg,20,DAvg./1000,'filled','MarkerEdgeColor','k');
						xlabel('距离 (km)')
						ylabel('k_{sn}');
						title('距离-k_{sn}关系');
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax2);
					    end							
						hold off

						ax3=subplot(3,1,3);
						hold on
						pl1=plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
						pl2=plotdz(Sn,DEMc,'dunit','km','Color','k');
						pl3=scatter((C.distance)./1000,C.elev,5,(C.distance)./1000,'filled');
						xlabel('距河口距离 (km)')
						ylabel('高程 (m)')
						legend([pl1 pl2 pl3],'未平滑DEM','平滑DEM','河道距离','location','best');
						title('纵剖面图 : 选择河段 - 完成后按回车','Color','r')
						ax3.XColor='Red';
						ax3.YColor='Red';
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax3);
					    end							
						hold off

						linkaxes([ax2,ax3],'x');
						colormap(ax1,'jet'); colormap(ax2,'jet'); colormap(ax3,'jet');	
                    end	

					[d,~,bttn]=ginput;
					d=d*1000; % 转换为米单位;

					% 检测非边界裂点
					bttn_idx=bttn~=1;
					if any(bttn_idx)
						% 分离非边界点
						d_kn=d(bttn_idx);
						d(bttn_idx)=[];
						% 转换为栅格索引
						rd=C.distance; rx=C.x; ry=C.y;
						kn_ix=zeros(numel(d_kn),1); 
						for jj=1:numel(d_kn)
							distdist=sqrt(sum(bsxfun(@minus, rd, d_kn(jj)).^2,2));
							[~,knbix]=min(distdist);
							knbx=rx(knbix);
							knby=ry(knbix);
							kn_ix(jj)=coord2ind(DEM,knbx,knby);
						end
					 else
						kn_ix=NaN;
                     end

					if isempty(d)
						if strcmp(theta_method,'ref')
							Cbf=ChiCalc(Sn,DEMc,A,1);
							ksn_list=[C.x C.y C.area ones(numel(C.x),1)*C.ks ones(numel(C.x),1)*C.ks_neg ones(numel(C.x),1)*C.ks_pos ...
								ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*Cbf.mn ones(numel(C.x),1)*snda];
						elseif strcmp(theta_method,'auto')
							ksn_list=[C.x C.y C.area ones(numel(C.x),1)*C.ks ones(numel(C.x),1)*C.ks_neg ones(numel(C.x),1)*C.ks_pos ...
								ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*snda];
						end

						% 确定ksn值在色阶中的位置并绘图
						ksn_val=C.ks;

						if ksn_val > mksn
							figure(f1)
							hold on
							plot(Sn,'Color',KSN_col(end,:),'LineWidth',2);
							hold off
						else
							edges=linspace(0,mksn,10);
							n=histc(ksn_val,edges);
							figure(f1)
							hold on
							plot(Sn,'Color',KSN_col(logical(n),:),'LineWidth',2);
							hold off	
						end

						[~,lbix]=min(C.chi);
						elbl=C.elev(lbix);
						if display_slope_area
							figure(f2)
							subplot(4,1,1);
							hold on
							plot(C.chi,C.pred+elbl,'-k','LineWidth',2);
							hold off

							subplot(4,1,3);
							hold on
							pl4=plot((C.distance)./1000,C.pred+elbl,'-k','LineWidth',2);
							legend([pl1 pl2 pl3 pl4],'未平滑DEM','平滑DEM','河道距离','段拟合','location','best');
							hold off

							subplot(4,1,4);
							hold on
							plot(C.area,ksn_val.*C.area.^(-C.mn),'-k','LineWidth',2);
							hold off							
						else
							figure(f2)
							subplot(3,1,1);
							hold on
							plot(C.chi,C.pred+elbl,'-k','LineWidth',2);
							hold off

							subplot(3,1,3);
							hold on
							pl4=plot((C.distance)./1000,C.pred+elbl,'-k','LineWidth',2);
							legend([pl1 pl2 pl3 pl4],'未平滑DEM','平滑DEM','河道距离','段拟合','location','best');
							hold off
						end

						res_list=[C.chi C.res];
						bnd_ix=NaN;

					else
						% 排序裂点并构建边界列表
						ds=sortrows(d);
						bnds=vertcat(0,ds,max(Sn.distance));

						num_bnds=numel(bnds);
						rd=C.distance; 
						rx=C.x;
						ry=C.y; 
						rc=C.chi;
						for jj=1:num_bnds-1
							% 提取边界
							lb=bnds(jj);
							rb=bnds(jj+1);

							% 裁剪河段
							lb_dist=sqrt(sum(bsxfun(@minus, rd, lb).^2,2));
							rb_dist=sqrt(sum(bsxfun(@minus, rd, rb).^2,2));

							[~,lbix]=min(lb_dist);
							[~,rbix]=min(rb_dist);

							lbx=rx(lbix);
							lby=ry(lbix);

							rbx=rx(rbix);
							rby=ry(rbix);	

							lix=coord2ind(DEM,lbx,lby);	
							rix=coord2ind(DEM,rbx,rby);

							Seg=modify(Sn,'downstreamto',rix);
							Seg=modify(Seg,'upstreamto',lix);

							% 重建包含下游节点的河段
							WSEG=GRIDobj(DEM,'logical');
							WSEG.Z(Seg.IXgrid)=true;
							WSEG.Z(lix)=true;
							% 河段末端处理
							if jj==num_bnds-1
								WSEG.Z(rix)=true;
							end
							Seg=STREAMobj(FD,WSEG);

							% 河段长度校验与扩展
							lbix_new=lbix+1;
							first_time=true;
							while numel(Seg.IXgrid)<=2
								if first_time
									wrn_mssg=['第' num2str(jj) '段选择的河段太短，已向下游扩展边界'];
									wd=warndlg(wrn_mssg);
									uiwait(wd);
								end
								lix_new=coord2ind(DEM,rx(lbix_new),ry(lbix_new));
								WSEG.Z(lix_new)=true;
								Seg=STREAMobj(FD,WSEG);
								lbix_new=lbix_new+1;
								first_time=false;
                            end

                            % 构建边界索引列表
							if jj<num_bnds-1
								bnd_ix(jj,1)=rix;
							end	
							% 计算chi获取ksn和最佳凹度
							if strcmp(theta_method,'ref')
								Cseg=ChiCalc(Seg,DEMc,A,1,ref_theta);
							elseif strcmp(theta_method,'auto')
								Cseg=ChiCalc(Seg,DEMc,A,1,C.mn);
							end
							Cbfseg=ChiCalc(Seg,DEMc,A,1);

							% 确定ksn值在色阶中的位置并绘图
							ksn_val=Cseg.ks;

							if ksn_val > mksn
								figure(f1)
								hold on
								plot(Seg,'Color',KSN_col(end,:),'LineWidth',2);
								hold off
							else
								edges=linspace(0,mksn,10);
								n=histc(ksn_val,edges);
								figure(f1)
								hold on
								plot(Seg,'Color',KSN_col(logical(n),:),'LineWidth',2);
								hold off	
							end

							ksn_nodes{jj,1}=[Cseg.x Cseg.y Cseg.area ones(numel(Cseg.x),1)*Cseg.ks ones(numel(Cseg.x),1)*Cseg.ks_neg ones(numel(Cseg.x),1)*Cseg.ks_pos ...
								ones(numel(Cseg.x),1)*Cseg.mn ones(numel(Cseg.x),1)*Cbfseg.mn ones(numel(Cseg.x),1)*snda];

							% 绘制线性拟合结果	
							rchi=rc(rb_dist==min(rb_dist));
							lchi=rc(lb_dist==min(lb_dist));
							ld=rd(lb_dist==min(lb_dist));
							segChi=linspace(lchi,rchi,numel(Cseg.chi));
							seg0Chi=linspace(0,max(Cseg.chi),numel(Cseg.chi));
							[~,lbix]=min(Cseg.chi);
							elbl=Cseg.elev(lbix);
							if display_slope_area
								figure(f2)
								subplot(4,1,1);
								hold on
								plot(segChi,(seg0Chi.*ksn_val)+elbl,'-k','LineWidth',2);
								hold off

								seg_st=rd(lb_dist==min(lb_dist));
								subplot(4,1,3);
								hold on
								pl4=plot((Cseg.distance+seg_st)/1000,(Cseg.pred)+elbl,'-k','LineWidth',2);
								legend([pl1 pl2 pl3 pl4],'未平滑DEM','平滑DEM','河道距离','段拟合','location','best');
								hold off

								subplot(4,1,4);
								hold on
								plot(Cseg.area,ksn_val.*Cseg.area.^(-Cseg.mn),'-k','LineWidth',2);
								hold off									
							else
								figure(f2)
								subplot(3,1,1);
								hold on
								plot(segChi,(seg0Chi.*ksn_val)+elbl,'-k','LineWidth',2);
								hold off

								seg_st=rd(lb_dist==min(lb_dist));
								subplot(3,1,3);
								hold on
								pl4=plot((Cseg.distance+seg_st)/1000,(Cseg.pred)+elbl,'-k','LineWidth',2);
								legend([pl1 pl2 pl3 pl4],'未平滑DEM','平滑DEM','河道距离','段拟合','location','best');
								hold off	
							end						

							res_list{jj,1}=[Cseg.distance+ld Cseg.res];
						end

						ksn_list=vertcat(ksn_nodes{:});
						res_list=vertcat(res_list{:});
					end

					%% 绘制拟合结果
					if display_slope_area
						figure(f2)
						subplot(4,1,2)
						hold on
						plot(res_list(:,1)./1000,ksn_list(:,4),'-k','LineWidth',2);
						hold off
					else
						figure(f2)
						subplot(3,1,2)
						hold on
						plot(res_list(:,1)./1000,ksn_list(:,4),'-k','LineWidth',2);
						hold off
					end

					f3=figure(3);
					set(f3,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');

					clf
					sbplt1=subplot(2,1,2);
					hold on
					s1=scatter(DAvg./1000,KsnAvg,20,'k','filled');
					p1=plot(res_list(:,1)./1000,ksn_list(:,4)-ksn_list(:,5),':k');
					plot(res_list(:,1)./1000,ksn_list(:,4)+ksn_list(:,6),':k');
					p2=plot(res_list(:,1)./1000,ksn_list(:,4),'-k','LineWidth',2);
					[ksn_vals,~,ksn_ix]=unique(ksn_list(:,4));
					d_means=accumarray(ksn_ix,res_list(:,1),[],@(x) mean(x,'omitnan'));
					for kk=1:numel(ksn_vals)
						text(d_means(kk)./1000,ksn_vals(kk),['k_{sn} = ' num2str(ksn_vals(kk))],...
							'VerticalAlignment','bottom','HorizontalAlignment','center');
					end
					xlabel('距离 (km)')
					ylabel('k_{sn}')
					title('距离-k_{sn}关系')
					legend([s1 p1 p2],{'k_{sn}','k_{sn}误差范围','拟合段k_{sn}'},'location','best');
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(sbplt1);
				    end						
					hold off

					sbplt2=subplot(2,1,1);
					hold on
					plot([min(res_list(:,1)./1000) max(res_list(:,1)./1000)],[0 0],'-k');
					scatter(res_list(:,1)./1000,res_list(:,2),10,'k','filled');
					xlabel('距离 (km)')
					ylabel('残差 (m)')
					title('k_{sn}拟合残差分布')
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(sbplt2);
				    end						
					hold off

				elseif strcmp(pick_method,'slope_area')

					[bs,ba,bc,bd,bk,aa,ag,ad,ac]=sa_ksn(DEMc,Sn,A,C.chi,ak,smooth_distance);

					ax3=subplot(4,1,3);
					hold on
					pl1=plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
					pl2=plotdz(Sn,DEMc,'dunit','km','Color','k');
					pl3=scatter((C.distance)./1000,C.elev,5,log10(C.area),'filled');
					xlabel('距河口距离 (km)')
					ylabel('高程 (m)')
					legend([pl1 pl2 pl3],'未平滑DEM','平滑DEM','对数流域面积','location','best');
					title('纵剖面图')
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax3);
				    end						
					hold off

ax2=subplot(4,1,2);
hold on
scatter(ba,bk,20,log10(ba),'filled','MarkerEdgeColor','k');
xlabel('对数面积')  % Log Area -> 对数面积
ylabel('k_{sn}');  % Auto k_{sn} -> 自动k_{sn}
title('对数面积 - k_{sn}');  % Log Area - Auto k_{sn} -> 对数面积 - 自动k_{sn}
set(ax2,'XScale','log','XDir','reverse');
if ~verLessThan('matlab','9.5')
    disableDefaultInteractivity(ax2);
end						
hold off

ax1=subplot(4,1,1);
hold on
plot(C.chi,C.elev,'-k');
scatter(C.chi,C.elev,10,log10(C.area),'filled');
xlabel('\chi')
ylabel('高程 (m)')  % Elevation (m) -> 高程 (m)
title('\chi - 高程')  % \chi - Z -> \chi - 高程
if ~verLessThan('matlab','9.5')
    disableDefaultInteractivity(ax1);
end						
hold off

ax4=subplot(4,1,4);
hold on
scatter(aa,ag,5,log10(aa),'+');
scatter(ba,bs,20,log10(ba),'filled','MarkerEdgeColor','k');
xlabel('对数面积');  % Log Area -> 对数面积
ylabel('对数坡度');  % Log Gradient -> 对数坡度
title(['坡度-面积: \theta = ' num2str(C.mn) ' : 选取河段 - 按回车键结束'],'Color','r');  % Slope-Area... -> 坡度-面积... Press Enter When Done -> 按回车键结束
set(ax4,'YScale','log','XScale','log','XDir','reverse');
ax4.XColor='Red';
ax4.YColor='Red';
if ~verLessThan('matlab','9.5')
    disableDefaultInteractivity(ax4);
end						
hold off

linkaxes([ax4,ax2],'x');
colormap(ax1,'jet'); colormap(ax2,'jet'); colormap(ax3,'jet'); colormap(ax4,'jet');

[av,~,bttn]=ginput;
% 确定是否存在非边界节点
bttn_idx=bttn~=1;
if any(bttn_idx)
    % 提取非边界节点
    av_kn=av(bttn_idx);
    av(bttn_idx)=[];
    % 转换为索引
    ra=C.area; rx=C.x; ry=C.y;
    kn_ix=zeros(numel(av_kn),1); 
    for jj=1:numel(av_kn)
        areadist=sqrt(sum(bsxfun(@minus, ra, av_kn(jj)).^2,2));
        [~,knbix]=min(areadist);
        knbx=rx(knbix);
        knby=ry(knbix);
        kn_ix(jj)=coord2ind(DEM,knbx,knby);
    end
else
    kn_ix=NaN;
end					

if isempty(av)
    if strcmp(theta_method,'ref')
        Cbf=ChiCalc(Sn,DEMc,A,1);
        ksn_list=[C.x C.y C.area ones(numel(C.x),1)*C.ks ones(numel(C.x),1)*C.ks_neg ones(numel(C.x),1)*C.ks_pos ...
            ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*Cbf.mn ones(numel(C.x),1)*snda];
    elseif strcmp(theta_method,'auto')
        ksn_list=[C.x C.y C.area ones(numel(C.x),1)*C.ks ones(numel(C.x),1)*C.ks_neg ones(numel(C.x),1)*C.ks_pos ...
            ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*snda];
    end

    % 确定ksn值在色阶中的位置并绘图
    ksn_val=C.ks;

    if ksn_val > mksn
        figure(f1)
        hold on
        plot(Sn,'Color',KSN_col(end,:),'LineWidth',2);
        hold off
    else
        edges=linspace(0,mksn,10);
        n=histc(ksn_val,edges);
        figure(f1)
        hold on
        plot(Sn,'Color',KSN_col(logical(n),:),'LineWidth',2);
        hold off	
    end

    [~,lbix]=min(C.chi);
    elbl=C.elev(lbix);

    figure(f2)
    subplot(4,1,1);
    hold on
    plot(C.chi,C.pred+elbl,'-k','LineWidth',2);
    hold off

    subplot(4,1,3);
    hold on
    pl4=plot((C.distance)/1000,C.pred+elbl,'-k','LineWidth',2);
    legend([pl1 pl2 pl3 pl4],'未平滑DEM','平滑后的DEM','对数集水面积','河段拟合','location','best');  % 更新图例翻译
    hold off

    subplot(4,1,4);
    hold on
    plot(C.area,ksn_val.*C.area.^(-C.mn),'-k','LineWidth',2);
    hold off

    res_list=[C.area C.res];
    bnd_ix=NaN;

else
    % 对节点列表排序并创建边界列表
    avs=sortrows(av,'descend');
    bnds=vertcat(max(C.area,[],'omitnan'),avs,min(C.area,[],'omitnan'));

    num_bnds=numel(bnds);
    rc=C.chi;
    rx=C.x;
    ry=C.y;
    rd=C.distance;
    ra=C.area;
    for jj=1:num_bnds-1
        % 提取边界
        lb=bnds(jj);
        rb=bnds(jj+1);

        % 裁剪河流段
        lb_dadist=sqrt(sum(bsxfun(@minus, ra, lb).^2,2));
        rb_dadist=sqrt(sum(bsxfun(@minus, ra, rb).^2,2));

        [~,lbix]=min(lb_dadist);
        [~,rbix]=min(rb_dadist);

        lbx=rx(lbix);
        lby=ry(lbix);

        rbx=rx(rbix);
        rby=ry(rbix);	

        lix=coord2ind(DEM,lbx,lby);
        rix=coord2ind(DEM,rbx,rby);

        Seg=modify(Sn,'downstreamto',rix);
        Seg=modify(Seg,'upstreamto',lix);

        % 重建包含下游边界节点的河流段
        WSEG=GRIDobj(DEM,'logical');
        WSEG.Z(Seg.IXgrid)=true;
        WSEG.Z(lix)=true;
        % 如果是河段末端则添加上游节点
        if jj==num_bnds-1
            WSEG.Z(rix)=true;
        end
        Seg=STREAMobj(FD,WSEG);

        % 检查河流段的长度，如果节点数小于等于两个，则向下游扩展直到足够长
        lbix_new=lbix+1;
        first_time=true;
        while numel(Seg.IXgrid)<=2
            if first_time
                wrn_mssg=['所选河段' num2str(jj) '太短，已向下游扩展边界'];  % Warning message translation
                wd=warndlg(wrn_mssg);
                uiwait(wd);
            end
            lix_new=coord2ind(DEM,rx(lbix_new),ry(lbix_new));
            WSEG.Z(lix_new)=true;
            Seg=STREAMobj(FD,WSEG);
            lbix_new=lbix_new+1;
            first_time=false;
        end							

        % 构建边界列表
        if jj<num_bnds-1
            bnd_ix(jj,1)=rix;
        end

        % 计算chi来获取ksn和最佳拟合凹度
        if strcmp(theta_method,'ref')
            Cseg=ChiCalc(Seg,DEMc,A,1,ref_theta);
        elseif strcmp(theta_method,'auto')
            Cseg=ChiCalc(Seg,DEMc,A,1,C.mn);
        end
        Cbfseg=ChiCalc(Seg,DEMc,A,1);								

        % 确定ksn值在色阶中的位置并绘图
        ksn_val=Cseg.ks;

        if ksn_val > mksn
            figure(f1)
            hold on
            plot(Seg,'Color',KSN_col(end,:),'LineWidth',2);
            hold off
        else
            edges=linspace(0,mksn,10);
            n=histc(ksn_val,edges);
            figure(f1)
            hold on
            plot(Seg,'Color',KSN_col(logical(n),:),'LineWidth',2);
            hold off	
        end

        ksn_nodes{jj,1}=[Cseg.x Cseg.y Cseg.area ones(numel(Cseg.x),1)*Cseg.ks ones(numel(Cseg.x),1)*Cseg.ks_neg ones(numel(Cseg.x),1)*Cseg.ks_pos ...
            ones(numel(Cseg.x),1)*Cseg.mn ones(numel(Cseg.x),1)*Cbfseg.mn ones(numel(Cseg.x),1)*snda];

        % 绘制线性拟合
        rchi=rc(rb_dadist==min(rb_dadist));
        lchi=rc(lb_dadist==min(lb_dadist));
        segChi=linspace(lchi,rchi,numel(Cseg.chi));
        seg0Chi=linspace(0,max(Cseg.chi),numel(Cseg.chi));
        [~,lbix]=min(Cseg.chi);
        elbl=Cseg.elev(lbix);

        figure(f2)
        subplot(4,1,1);
        hold on
        plot(segChi,(seg0Chi.*ksn_val)+elbl,'-k','LineWidth',2);
        hold off

        seg_st=rd(lb_dadist==min(lb_dadist));
        subplot(4,1,3);
        hold on
        pl4=plot((Cseg.distance+seg_st)/1000,(Cseg.pred)+elbl,'-k','LineWidth',2);
        legend([pl1 pl2 pl3 pl4],'未平滑DEM','平滑DEM','对数集水面积','河段拟合','location','best');  % 更新图例翻译
        hold off

        subplot(4,1,4);
        hold on
        plot(Cseg.area,ksn_val.*Cseg.area.^(-1*Cseg.mn),'-k','LineWidth',2);
        hold off

        res_list{jj,1}=[Cseg.area Cseg.res];
    end

    ksn_list=vertcat(ksn_nodes{:});
    res_list=vertcat(res_list{:});
end

%% 绘制结果图
if display_slope_area
    figure(f2)
    subplot(4,1,2)
    hold on
    plot(res_list(:,1),ksn_list(:,4),'-k','LineWidth',2);
    hold off
else
    figure(f2)
    subplot(3,1,2)
    hold on
    plot(res_list(:,1),ksn_list(:,4),'-k','LineWidth',2);
    hold off
end

f3=figure(3);
set(f3,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');
clf
ax31=subplot(2,1,2);
hold on
s1=scatter(log10(ba),bk,20,'k','filled');
p1=plot(log10(res_list(:,1)),ksn_list(:,4)-ksn_list(:,5),':k');
plot(log10(res_list(:,1)),ksn_list(:,4)+ksn_list(:,6),':k');
p2=plot(log10(res_list(:,1)),ksn_list(:,4),'-k','LineWidth',2);
[ksn_vals,~,ksn_ix]=unique(ksn_list(:,4));
a_means=accumarray(ksn_ix,res_list(:,1),[],@(x) mean(x,'omitnan'));
for kk=1:numel(ksn_vals)
    text(log10(a_means(kk)),ksn_vals(kk),['k_{sn} = ' num2str(ksn_vals(kk))],...
        'VerticalAlignment','bottom','HorizontalAlignment','center');
end
xlabel('对数面积')  % Log Area -> 对数面积
ylabel('k_{sn}')
title('k_{sn} - 对数面积')  % k_{sn} - Log Area -> k_{sn} - 对数面积
legend([s1 p1 p2],{'k_{sn}','k_{sn}不确定度','拟合河段k_{sn}'},'location','best');  % Legend translations
set(ax31,'XScale','log','XDir','reverse');
if ~verLessThan('matlab','9.5')
    disableDefaultInteractivity(ax31);
end						
hold off

ax32=subplot(2,1,1);
hold on
plot([log10(min(res_list(:,1))) log10(max(res_list(:,1)))],[0 0],'-k');
scatter(log10(res_list(:,1)),res_list(:,2),10,'k','filled');
xlabel('对数面积')  % Log Area -> 对数面积
ylabel('残差 (m)')  % Residual (m) -> 残差 (m)
title('k_{sn}拟合残差')  % Residual on k_{sn} fit -> k_{sn}拟合残差
set(ax32,'XScale','log','XDir','Reverse');
if ~verLessThan('matlab','9.5')
    disableDefaultInteractivity(ax32);
end						
hold off

end

qa2=questdlg('请选择操作','河流拟合','停止选取','重新拟合','继续选取','继续选取');  % 对话框翻译

switch qa2
case '继续选取'
    str2 = 'Y';
    str1 = 'N';
    str3 = 'C';
    % 添加坡度、残差和河流编号到节点列表并清除NaN
    kidx=isnan(ksn_list(:,1));
    ksn_list(kidx,:)=[];
    res_list(kidx,:)=[];
    gix=coord2ind(DEM,ksn_list(:,1),ksn_list(:,2));
    kmat=[ksn_list G.Z(gix) res_list(:,2) ones(size(ksn_list,1),1).*ii];
    ksn_master{ii,1}=kmat;
    bnd_master{ii,1}=bnd_ix;
    kn_master{ii,1}=kn_ix;
    res_master{ii,1}=res_list;
    count=ii;
    save(out_restart_name,'ksn_master','bnd_master','kn_master','res_master','Sc','count','-append');
    if save_figures
        f2_name=[shape_name '_stream_fits_' num2str(ii) '.pdf'];
        f3_name=[shape_name '_stream_rsds_' num2str(ii) '.pdf'];
        print(f2,f2_name,'-dpdf','-fillpage');
        print(f3,f3_name,'-dpdf','-fillpage');
    end
    clear ksn_list ksn_nodes res_list bnd_ix kn_ix kidx gix kmat;
    ii=ii+1;
case '停止选取'
    wtb=waitbar(0,'正在清理并生成输出，请勿关闭窗口'); 
    str1 = 'N';
    str2 = 'N';
    str3 = 'C';
    % 添加坡度、残差和河流编号到节点列表并清除NaN
    kidx=isnan(ksn_list(:,1));
    ksn_list(kidx,:)=[];
    res_list(kidx,:)=[];
    gix=coord2ind(DEM,ksn_list(:,1),ksn_list(:,2));
    kmat=[ksn_list G.Z(gix) res_list(:,2) ones(size(ksn_list,1),1).*ii];
    ksn_master{ii,1}=kmat;
    bnd_master{ii,1}=bnd_ix;
    kn_master{ii,1}=kn_ix;
    res_master{ii,1}=res_list;
    count=ii;
    save(out_restart_name,'ksn_master','bnd_master','kn_master','res_master','Sc','count','-append');					
    if save_figures
        f2_name=[shape_name '_stream_fits_' num2str(ii) '.pdf'];
        f3_name=[shape_name '_stream_rsds_' num2str(ii) '.pdf'];
        print(f2,f2_name,'-dpdf','-fillpage');
        print(f3,f3_name,'-dpdf','-fillpage');
    end
    clear ksn_list ksn_nodes res_list bnd_ix kn_ix kidx gix kmat;
case '重新拟合'
    str3 = 'R';
    str2 = 'Y';
    clear ksn_list ksn_nodes res_list bnd_ix kn_ix;
end

close figure 2
close figure 3
end
end

case 'preselected'
% 初始化循环参数
str1='R';
break_flag=false;

if strcmp(theta_method,'ref')
    % 自动计算ksn用于对比
    [auto_ksn]=KSN_Quick(DEM,A,S,ref_theta);
end

if rf && strcmp(r_type,'c') && strcmp(restart,'continue')
    load(out_mat_name,'count','ksn_master','bnd_master','kn_master','res_master','Sc');
    if count>=num_ch
        if isdeployed
            errordlg('运行似乎已完成，无法继续')  % Error dialog translation
        end
        error('运行似乎已完成，无法继续');
    end
    rng=count+1:num_ch;
elseif rf && strcmp(r_type,'c') && strcmp(restart,'skip')
    load(out_mat_name,'count','ksn_master','bnd_master','kn_master','res_master','Sc');
    if count>=num_ch
        if isdeployed
            errordlg('运行似乎已完成，无法继续')
        end
        error('运行似乎已完成，无法继续');
    end
    rng=count+2:num_ch;
elseif rf && strcmp(r_type,'r') && strcmp(restart,'continue')
    load(out_restart_name,'count','ksn_master','bnd_master','kn_master','res_master','Sc');
    if count>=num_ch
        if isdeployed
            errordlg('运行似乎已完成，无法继续')
        end
        error('运行似乎已完成，无法继续');
    end
    rng=count+1:num_ch;
elseif rf && strcmp(r_type,'r') && strcmp(restart,'skip')
    load(out_restart_name,'count','ksn_master','bnd_master','kn_master','res_master','Sc');
    if count>=num_ch
        if isdeployed
            errordlg('运行似乎已完成，无法继续')
        end
        error('运行似乎已完成，无法继续');
    end
    rng=count+2:num_ch;
else
    rng=1:num_ch;
end		

   for ii=rng

		while strcmpi(str1,'R')
           
			chOI=s_ch(ii,:);

			% 构建逻辑栅格
			ix=coord2ind(DEM,chOI(:,1),chOI(:,2));
			IX=GRIDobj(DEM,'logical');
			IX.Z(ix)=true;

			% 提取目标河道
			Sn=modify(S,'downstreamto',IX);

			if redefine_thresh
				[Sn]=RedefineThreshold(DEM,FD,A,Sn,FLUS,ref_theta,rd_pick_method,smooth_distance,ii,save_figures,shape_name);
				% 更新DEMc
				if any(isnan(getnal(Sn,DEMc)));
					zc=mincosthydrocon(Sn,DEM,'interp',iv);
					DEMc.Z(Sn.IXgrid)=zc;
				end

				% 重新计算auto_ksn
				if strcmp(theta_method,'ref')
					[auto_ksn]=KSN_Quick(DEM,A,Sn,ref_theta);
				end
			end

			% 构建选取河道的复合河道网络
			if strcmp(junction_method,'check')
				if ii>1
					[IIXX,~,~,Si]=intersectlocs(Sc,Sn);
					if isempty(IIXX)
						Sn=Sn;
						Sct=union(Sn,Sc,FD);
					else
						Sn=Si;

						if ~isempty(min_channel_length) & max(Sn.distance)<min_channel_length
							disp(['跳过河道头' num2str(ii) '，河道段过短']);
							break;
						end

						Sct=union(Sn,Sc,FD);
					end
				else
					Sct=Sn;
				end
			elseif strcmp(junction_method,'ignore')
				if ii>1
					Sct=union(Sn,Sc,FD);
				else
					Sct=Sn;
				end
			end

			%% 提取阈值集水面积
			snchix=streampoi(Sn,'channelheads','ix');
			snda=DA.Z(snchix);

			%% 计算chi并提取ksn数据
			if strcmp(theta_method,'ref')
				C=ChiCalc(Sn,DEMc,A,1,ref_theta);
				ak=getnal(Sn,auto_ksn);
			elseif strcmp(theta_method,'auto')
				C=ChiCalc(Sn,DEMc,A,1);
				if redefine_thresh
					[auto_ksn]=KSN_Quick(DEM,A,Sn,C.mn);
				else
					[auto_ksn]=KSN_Quick(DEM,A,S,C.mn);
				end
				ak=getnal(Sn,auto_ksn);
			end

			%% 数据分箱
			[DAvg,KsnAvg]=BinAverage(Sn.distance,ak,smooth_distance);
			[~,CAvg]=BinAverage(C.distance,C.chi,smooth_distance);


			%% 初始化图形界面以选取边界
			f2=figure(2);
			set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');
			clf	

			%% 不同选取方法的主逻辑
			if strcmp(pick_method,'chi')

				if display_slope_area
					[bs,ba,bc,bd,aa,ag,ad,ac]=sa(DEMc,Sn,A,C.chi,smooth_distance);

					ax4=subplot(4,1,4);
					hold on
					scatter(aa,ag,5,ac,'+');
					scatter(ba,bs,20,bc,'filled','MarkerEdgeColor','k');
					xlabel('对数面积');
					ylabel('对数坡度');
					title('坡度-面积图');
					set(ax4,'YScale','log','XScale','log','XDir','reverse');
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax4);
				    end							
					hold off

					ax3=subplot(4,1,3);
					hold on
					pl1=plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
					pl2=plotdz(Sn,DEMc,'dunit','km','Color','k');
					pl3=scatter((C.distance)./1000,C.elev,5,C.chi,'filled');
					xlabel('距河口距离 (km)')
					ylabel('高程 (m)')
					legend([pl1 pl2 pl3],'未平滑DEM','平滑DEM','\chi','location','best');
					title('纵剖面')
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax3);
				    end							
					hold off

					ax2=subplot(4,1,2);
					hold on
					scatter(CAvg,KsnAvg,20,CAvg,'filled','MarkerEdgeColor','k');
					xlabel('\chi')
					ylabel('k_{sn}');
					title('\chi - k_{sn}');
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax2);
				    end							
					hold off

					ax1=subplot(4,1,1);
					hold on
					plot(C.chi,C.elev,'-k');
					scatter(C.chi,C.elev,10,C.chi,'filled');
					xlabel('\chi')
					ylabel('高程 (m)')
					title(['\chi - 高程 : \theta = ' num2str(C.mn) ' : 选取河段 - 完成后按回车'],'Color','r')
					ax1.XColor='Red';
					ax1.YColor='Red';
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax1);
				    end							
					hold off

					linkaxes([ax1,ax2],'x');
					colormap(ax1,'jet'); colormap(ax2,'jet'); colormap(ax3,'jet'); colormap(ax4,'jet');

				else
					ax3=subplot(3,1,3);
					hold on
					pl1=plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
					pl2=plotdz(Sn,DEMc,'dunit','km','Color','k');
					pl3=scatter((C.distance)./1000,C.elev,5,C.chi,'filled');
					xlabel('距河口距离 (km)')
					ylabel('高程 (m)')
					legend([pl1 pl2 pl3],'未平滑DEM','平滑后的DEM','\chi','location','best');
					title('纵剖面')
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax3);
				    end							
					hold off

					ax2=subplot(3,1,2);
					hold on
					scatter(CAvg,KsnAvg,20,CAvg,'filled','MarkerEdgeColor','k');
					xlabel('\chi')
					ylabel('k_{sn}');
					title('\chi - k_{sn}');
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax2);
				    end							
					hold off

					ax1=subplot(3,1,1);
					hold on
					plot(C.chi,C.elev,'-k');
					scatter(C.chi,C.elev,10,C.chi,'filled');
					xlabel('\chi')
					ylabel('高程 (m)')
					title(['\chi - 高程 : \theta = ' num2str(C.mn) ' : 选取河段 - 完成后按回车'],'Color','r')
					ax1.XColor='Red';
					ax1.YColor='Red';
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax1);
				    end							
					hold off

					linkaxes([ax1,ax2],'x');
					colormap(ax1,'jet'); colormap(ax2,'jet'); colormap(ax3,'jet');
				end	

				[cv,~,bttn]=ginput;
				% 判断是否存在非边界裂点
				bttn_idx=bttn~=1;
				if any(bttn_idx)
					% 分离非边界点
					cv_kn=cv(bttn_idx);
					cv(bttn_idx)=[];
					% 转换为索引
					rc=C.chi; rx=C.x; ry=C.y;
					kn_ix=zeros(numel(cv_kn),1); 
					for jj=1:numel(cv_kn);
						chidist=sqrt(sum(bsxfun(@minus, rc, cv_kn(jj)).^2,2));
						[~,knbix]=min(chidist);
						knbx=rx(knbix);
						knby=ry(knbix);
						kn_ix(jj)=coord2ind(DEM,knbx,knby);
					end
				else
					kn_ix=NaN;
				end

				if isempty(cv)
					if strcmp(theta_method,'ref')
						Cbf=ChiCalc(Sn,DEMc,A,1);
						ksn_list=[C.x C.y C.area ones(numel(C.x),1)*C.ks ones(numel(C.x),1)*C.ks_neg ones(numel(C.x),1)*C.ks_pos ...
							ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*Cbf.mn ones(numel(C.x),1)*snda];
					elseif strcmp(theta_method,'auto')
						ksn_list=[C.x C.y C.area ones(numel(C.x),1)*C.ks ones(numel(C.x),1)*C.ks_neg ones(numel(C.x),1)*C.ks_pos ...
							ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*snda];
					end

					% 确定ksn值在色标中的位置并绘图
					ksn_val=C.ks;

					[~,lbix]=min(C.chi);
					elbl=C.elev(lbix);
					if display_slope_area
						figure(f2)
						subplot(4,1,1);
						hold on
						plot(C.chi,C.pred+elbl,'-k','LineWidth',2);
						hold off

						subplot(4,1,3);
						hold on
						pl4=plot((C.distance)/1000,C.pred+elbl,'-k','LineWidth',2);
						legend([pl1 pl2 pl3 pl4],'未平滑DEM','平滑后DEM','\chi','拟合河段','location','best');
						hold off

						subplot(4,1,4);
						hold on
						plot(C.area,ksn_val.*C.area.^(-C.mn),'-k','LineWidth',2);
						hold off
					else
						figure(f2)
						subplot(3,1,1);
						hold on
						plot(C.chi,C.pred+elbl,'-k','LineWidth',2);
						hold off

						subplot(3,1,3);
						hold on
						pl4=plot((C.distance)/1000,C.pred+elbl,'-k','LineWidth',2);
						legend([pl1 pl2 pl3 pl4],'未平滑DEM','平滑后的DEM','\chi','拟合河段','location','best');
						hold off
					end

					res_list=[C.chi C.res];
					bnd_ix=NaN;;

				else
					% 对裂点列表排序并构建边界列表
					cvs=sortrows(cv);
					bnds=vertcat(0,cvs,C.chi(1));

					num_bnds=numel(bnds);
					rc=C.chi;
					rx=C.x;
					ry=C.y;
					rd=C.distance;
					for jj=1:num_bnds-1
						% 提取边界
						lb=bnds(jj);
						rb=bnds(jj+1);

						% 裁剪河道段
						lb_chidist=sqrt(sum(bsxfun(@minus, rc, lb).^2,2));
						rb_chidist=sqrt(sum(bsxfun(@minus, rc, rb).^2,2));

						[~,lbix]=min(lb_chidist);
						[~,rbix]=min(rb_chidist);

						lbx=rx(lbix);
						lby=ry(lbix);

						rbx=rx(rbix);
						rby=ry(rbix);

						lix=coord2ind(DEM,lbx,lby);
						rix=coord2ind(DEM,rbx,rby);

						Seg=modify(Sn,'downstreamto',rix);
						Seg=modify(Seg,'upstreamto',lix);

						% 重建包含下游边界节点的河道
						WSEG=GRIDobj(DEM,'logical');
						WSEG.Z(Seg.IXgrid)=true;
						WSEG.Z(lix)=true;
						% 如果是河道末端则添加回上游节点
						if jj==num_bnds-1
							WSEG.Z(rix)=true;
						end
						Seg=STREAMobj(FD,WSEG);

						% 检查河道长度，若节点数不足两个则向下游扩展直至满足条件
						lbix_new=lbix+1;
						first_time=true;
						while numel(Seg.IXgrid)<=2
							if first_time
								wrn_mssg=['第' num2str(jj) '段选取的河段过短，河段边界已向下游扩展'];
								wd=warndlg(wrn_mssg);
								uiwait(wd);
							end
							lix_new=coord2ind(DEM,rx(lbix_new),ry(lbix_new));
							WSEG.Z(lix_new)=true;
							Seg=STREAMobj(FD,WSEG);
							lbix_new=lbix_new+1;
							first_time=false;
						end

						% 构建边界索引列表
						if jj<num_bnds-1
							bnd_ix(jj,1)=rix;
						end

						% 计算chi以获取ksn和最佳拟合凹度 
						if strcmp(theta_method,'ref')
							Cseg=ChiCalc(Seg,DEMc,A,1,ref_theta);
						elseif strcmp(theta_method,'auto')
							Cseg=ChiCalc(Seg,DEMc,A,1,C.mn);
						end
						Cbfseg=ChiCalc(Seg,DEMc,A,1);								

						% 确定ksn值在色标中的位置并绘图
						ksn_val=Cseg.ks;

						ksn_nodes{jj,1}=[Cseg.x Cseg.y Cseg.area ones(numel(Cseg.x),1)*Cseg.ks ones(numel(Cseg.x),1)*Cseg.ks_neg ones(numel(Cseg.x),1)*Cseg.ks_pos ...
							ones(numel(Cseg.x),1)*Cseg.mn ones(numel(Cseg.x),1)*Cbfseg.mn ones(numel(Cseg.x),1)*snda];

						% 绘制线性拟合结果	
						rchi=rc(rb_chidist==min(rb_chidist));
						lchi=rc(lb_chidist==min(lb_chidist));
						segChi=linspace(lchi,rchi,numel(Cseg.chi));
						seg0Chi=linspace(0,max(Cseg.chi),numel(Cseg.chi));
						[~,lbix]=min(Cseg.chi);
						elbl=Cseg.elev(lbix);

						if display_slope_area
							figure(f2)
							subplot(4,1,1);
							hold on
							plot(segChi,(seg0Chi.*ksn_val)+elbl,'-k','LineWidth',2);
							hold off

							seg_st=rd(lb_chidist==min(lb_chidist));
							subplot(4,1,3);
							hold on
							pl4=plot((Cseg.distance+seg_st)/1000,(Cseg.pred)+elbl,'-k','LineWidth',2);
							legend([pl1 pl2 pl3 pl4],'未平滑DEM','平滑DEM','\chi','拟合河段','location','best');
							hold off

							subplot(4,1,4);
							hold on
							plot(Cseg.area,ksn_val.*Cseg.area.^(-Cseg.mn),'-k','LineWidth',2);
							hold off
						else
							figure(f2)
							subplot(3,1,1);
							hold on
							plot(segChi,(seg0Chi.*ksn_val)+elbl,'-k','LineWidth',2);
							hold off

							seg_st=rd(lb_chidist==min(lb_chidist));
							subplot(3,1,3);
							hold on
							pl4=plot((Cseg.distance+seg_st)/1000,(Cseg.pred)+elbl,'-k','LineWidth',2);
							legend([pl1 pl2 pl3 pl4],'未平滑DEM','平滑后DEM','\chi','拟合河段','location','best');
							hold off
						end

						res_list{jj,1}=[Cseg.chi+lchi Cseg.res];	

					end

					ksn_list=vertcat(ksn_nodes{:});
					res_list=vertcat(res_list{:});
				end

				%% 绘制结果图
				if display_slope_area
					figure(f2)
					subplot(4,1,2)
					hold on
					plot(res_list(:,1),ksn_list(:,4),'-k','LineWidth',2);
					hold off
				else
					figure(f2)
					subplot(3,1,2)
					hold on
					plot(res_list(:,1),ksn_list(:,4),'-k','LineWidth',2);
					hold off
				end

				f3=figure(3);
				set(f3,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');

				clf
				sbplt2=subplot(2,1,2);
				hold on
				s1=scatter(CAvg,KsnAvg,20,'k','filled');
				p1=plot(res_list(:,1),ksn_list(:,4)-ksn_list(:,5),':k');
				plot(res_list(:,1),ksn_list(:,4)+ksn_list(:,6),':k');
				p2=plot(res_list(:,1),ksn_list(:,4),'-k','LineWidth',2);
				[ksn_vals,~,ksn_ix]=unique(ksn_list(:,4));
				chi_means=accumarray(ksn_ix,res_list(:,1),[],@(x) mean(x,'omitnan'));
				for kk=1:numel(ksn_vals)
					text(chi_means(kk),ksn_vals(kk),['k_{sn} = ' num2str(ksn_vals(kk))],...
						'VerticalAlignment','bottom','HorizontalAlignment','center','Color','r');
				end
				xlabel('\chi')
				ylabel('k_{sn}')
				legend([s1 p1 p2],{'k_{sn}','k_{sn}误差范围','拟合段k_{sn}'},'location','best');
				title('k_{sn} - \chi 关系图')
				if ~verLessThan('matlab','9.5')
			        disableDefaultInteractivity(sbplt2);
			    end						
				hold off

				sbplt1=subplot(2,1,1);
				hold on
				plot([min(res_list(:,1)) max(res_list(:,1))],[0 0],'-k');
				scatter(res_list(:,1),res_list(:,2),10,'k','filled');
				xlabel('\chi')
				ylabel('残差 (m)')
				title('k_{sn}拟合残差')
				if ~verLessThan('matlab','9.5')
			        disableDefaultInteractivity(sbplt1);
			    end						
				hold off				

			elseif strcmp(pick_method,'stream')	

				if display_slope_area
					[bs,ba,bc,bd,aa,ag,ad,ac]=sa(DEMc,Sn,A,C.chi,smooth_distance);	

					ax4=subplot(4,1,4);
					hold on
					scatter(aa,ag,5,ad./1000,'+');
					scatter(ba,bs,20,bd./1000,'filled','MarkerEdgeColor','k');
					xlabel('对数面积');
					ylabel('对数坡度');
					title('坡度-面积图');
					set(ax4,'XScale','log','YScale','log','XDir','reverse');
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax4);
				    end							
					hold off

					ax1=subplot(4,1,1);
					hold on
					plot(C.chi,C.elev,'Color',[0.5 0.5 0.5]);
					scatter(C.chi,C.elev,10,(C.distance)./1000,'filled');
					xlabel('\chi')
					ylabel('高程 (m)')
					title(['\chi - 高程 : \theta = ' num2str(C.mn)])
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax1);
				    end							
					hold off

					ax2=subplot(4,1,2);
					hold on
					scatter(DAvg./1000,KsnAvg,20,DAvg./1000,'filled','MarkerEdgeColor','k');
					xlabel('距离 (km)')
					ylabel('k_{sn}');
					title('距离 - k_{sn}');
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax2);
				    end							
					hold off

					ax3=subplot(4,1,3);
					hold on
					pl1=plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
					pl2=plotdz(Sn,DEMc,'dunit','km','Color','k');
					pl3=scatter((C.distance)./1000,C.elev,5,(C.distance)./1000,'filled');
					xlabel('距河口距离 (km)')
					ylabel('高程 (m)')
					legend([pl1 pl2 pl3],'未平滑DEM','平滑后DEM','河道距离','location','best');
					title('纵剖面：选取河段 - 完成后按回车','Color','r')
					ax3.XColor='Red';
					ax3.YColor='Red';
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax3);
				    end							
					hold off

					linkaxes([ax2,ax3],'x');
					colormap(ax1,'jet'); colormap(ax2,'jet'); colormap(ax3,'jet'); colormap(ax4,'jet')
				else
					ax1=subplot(3,1,1);
					hold on
					plot(C.chi,C.elev,'Color',[0.5 0.5 0.5]);
					scatter(C.chi,C.elev,10,(C.distance)./1000,'filled');
					xlabel('\chi')
					ylabel('高程 (m)')
					title(['\chi - 高程 : \theta = ' num2str(C.mn)])
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax1);
				    end							
					hold off

					ax2=subplot(3,1,2);
					hold on
					scatter(DAvg./1000,KsnAvg,20,DAvg./1000,'filled','MarkerEdgeColor','k');
					xlabel('距离 (km)')
					ylabel('k_{sn}');
					title('距离 - k_{sn}');
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax2);
				    end							
					hold off

					ax3=subplot(3,1,3);
					hold on
					pl1=plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
					pl2=plotdz(Sn,DEMc,'dunit','km','Color','k');
					pl3=scatter((C.distance)./1000,C.elev,5,(C.distance)./1000,'filled');
					xlabel('距河口距离 (km)')
					ylabel('高程 (m)')
					legend([pl1 pl2 pl3],'未平滑DEM','平滑后DEM','河道距离','location','best');
					title('纵剖面：选取河段 - 完成后按回车','Color','r')
					ax3.XColor='Red';
					ax3.YColor='Red';
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax3);
				    end							
					hold off

					linkaxes([ax2,ax3],'x');
					colormap(ax1,'jet'); colormap(ax2,'jet'); colormap(ax3,'jet');	
				end	

				linkaxes([ax2,ax3],'x');
				colormap(ax1,'jet'); colormap(ax2,'jet'); colormap(ax3,'jet');

				[d,~,bttn]=ginput;
				d=d*1000; % 转换回米制单位

				% 判断是否存在非边界裂点
				bttn_idx=bttn~=1;
				if any(bttn_idx)
					% 分离非边界点
					d_kn=d(bttn_idx);
					d(bttn_idx)=[];
					% 转换为索引
					rd=C.distance; rx=C.x; ry=C.y;
					kn_ix=zeros(numel(d_kn),1); 
					for jj=1:numel(d_kn);
						distdist=sqrt(sum(bsxfun(@minus, rd, d_kn(jj)).^2,2));
						[~,knbix]=min(distdist);
						knbx=rx(knbix);
						knby=ry(knbix);
						kn_ix(jj)=coord2ind(DEM,knbx,knby);
					end
				else
					kn_ix=NaN;
				end

				if isempty(d)
					if strcmp(theta_method,'ref')
						Cbf=ChiCalc(Sn,DEMc,A,1);
						ksn_list=[C.x C.y C.area ones(numel(C.x),1)*C.ks ones(numel(C.x),1)*C.ks_neg ones(numel(C.x),1)*C.ks_pos ...
							ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*Cbf.mn ones(numel(C.x),1)*snda];
					elseif strcmp(theta_method,'auto')
						ksn_list=[C.x C.y C.area ones(numel(C.x),1)*C.ks ones(numel(C.x),1)*C.ks_neg ones(numel(C.x),1)*C.ks_pos ...
							ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*snda];
					end

					% 确定ksn值在色标中的位置并绘图
					ksn_val=C.ks;

					[~,lbix]=min(C.chi);
					elbl=C.elev(lbix);
					if display_slope_area
						figure(f2)
						subplot(4,1,1);
						hold on
						plot(C.chi,C.pred+elbl,'-k','LineWidth',2);
						hold off

						subplot(4,1,3);
						hold on
						pl4=plot((C.distance)./1000,C.pred+elbl,'-k','LineWidth',2);
						legend([pl1 pl2 pl3 pl4],'未平滑DEM','平滑后DEM','河道距离','拟合河段','location','best');
						hold off

						subplot(4,1,4);
						hold on
						plot(C.area,ksn_val.*C.area.^(-C.mn),'-k','LineWidth',2);
						hold off							
					else
						figure(f2)
						subplot(3,1,1);
						hold on
						plot(C.chi,C.pred+elbl,'-k','LineWidth',2);
						hold off

						subplot(3,1,3);
						hold on
						pl4=plot((C.distance)./1000,C.pred+elbl,'-k','LineWidth',2);
						legend([pl1 pl2 pl3 pl4],'未平滑DEM','平滑后DEM','河道距离','拟合河段','location','best');
						hold off
					end

					res_list=[C.chi C.res];
					bnd_ix=NaN;

				else
					% 对裂点列表排序并构建边界列表
					ds=sortrows(d);
					bnds=vertcat(0,ds,max(Sn.distance));

					num_bnds=numel(bnds);
					rd=C.distance; 
					rx=C.x;
					ry=C.y; 
					rc=C.chi;
					for jj=1:num_bnds-1
						% 提取边界
						lb=bnds(jj);
						rb=bnds(jj+1);

						% 裁剪河道段
						lb_dist=sqrt(sum(bsxfun(@minus, rd, lb).^2,2));
						rb_dist=sqrt(sum(bsxfun(@minus, rd, rb).^2,2));

						[~,lbix]=min(lb_dist);
						[~,rbix]=min(rb_dist);

						lbx=rx(lbix);
						lby=ry(lbix);

						rbx=rx(rbix);
						rby=ry(rbix);	

						lix=coord2ind(DEM,lbx,lby);	
						rix=coord2ind(DEM,rbx,rby);

						Seg=modify(Sn,'downstreamto',rix);
						Seg=modify(Seg,'upstreamto',lix);

						% 重建包含下游边界节点的河道
						WSEG=GRIDobj(DEM,'logical');
						WSEG.Z(Seg.IXgrid)=true;
						WSEG.Z(lix)=true;
						% 如果是河道末端则添加回上游节点
						if jj==num_bnds-1
							WSEG.Z(rix)=true;
						end
						Seg=STREAMobj(FD,WSEG);

						% 检查河道长度，若节点数不足两个则向下游扩展直至满足条件
						lbix_new=lbix+1;
						first_time=true;
						while numel(Seg.IXgrid)<=2
							if first_time
								wrn_mssg=['第' num2str(jj) '段选取的河段过短，河段边界已向下游扩展'];
								wd=warndlg(wrn_mssg);
								uiwait(wd);
							end
							lix_new=coord2ind(DEM,rx(lbix_new),ry(lbix_new));
							WSEG.Z(lix_new)=true;
							Seg=STREAMobj(FD,WSEG);
							lbix_new=lbix_new+1;
							first_time=false;
						end

						% 计算chi以获取ksn和最佳拟合凹度 
						if strcmp(theta_method,'ref')
							Cseg=ChiCalc(Seg,DEMc,A,1,ref_theta);
						elseif strcmp(theta_method,'auto')
							Cseg=ChiCalc(Seg,DEMc,A,1,C.mn);
						end
						Cbfseg=ChiCalc(Seg,DEMc,A,1);

						% 确定ksn值在色标中的位置并绘图
						ksn_val=Cseg.ks;

						ksn_nodes{jj,1}=[Cseg.x Cseg.y Cseg.area ones(numel(Cseg.x),1)*Cseg.ks ones(numel(Cseg.x),1)*Cseg.ks_neg ones(numel(Cseg.x),1)*Cseg.ks_pos ...
							ones(numel(Cseg.x),1)*Cseg.mn ones(numel(Cseg.x),1)*Cbfseg.mn ones(numel(Cseg.x),1)*snda];

						% 绘制线性拟合结果	
						rchi=rc(rb_dist==min(rb_dist));
						lchi=rc(lb_dist==min(lb_dist));
						ld=rd(lb_dist==min(lb_dist));
						segChi=linspace(lchi,rchi,numel(Cseg.chi));
						seg0Chi=linspace(0,max(Cseg.chi),numel(Cseg.chi));
						[~,lbix]=min(Cseg.chi);
						elbl=Cseg.elev(lbix);

						if display_slope_area
							figure(f2)
							subplot(4,1,1);
							hold on
							plot(segChi,(seg0Chi.*ksn_val)+elbl,'-k','LineWidth',2);
							hold off

							seg_st=rd(lb_dist==min(lb_dist));
							subplot(4,1,3);
							hold on
							pl4=plot((Cseg.distance+seg_st)/1000,(Cseg.pred)+elbl,'-k','LineWidth',2);
							legend([pl1 pl2 pl3 pl4],'未平滑DEM','平滑后DEM','河道距离','拟合河段','location','best');
							hold off

							subplot(4,1,4);
							hold on
							plot(Cseg.area,ksn_val.*Cseg.area.^(-Cseg.mn),'-k','LineWidth',2);
							hold off
						else
							figure(f2)
							subplot(3,1,1);
							hold on
							plot(segChi,(seg0Chi.*ksn_val)+elbl,'-k','LineWidth',2);
							hold off

							seg_st=rd(lb_dist==min(lb_dist));
							subplot(3,1,3);
							hold on
							pl4=plot((Cseg.distance+seg_st)/1000,(Cseg.pred)+elbl,'-k','LineWidth',2);
							legend([pl1 pl2 pl3 pl4],'未平滑DEM','平滑后DEM','河道距离','拟合河段','location','best');
							hold off
						end

						res_list{jj,1}=[Cseg.distance+ld Cseg.res];
					end

					ksn_list=vertcat(ksn_nodes{:});
					res_list=vertcat(res_list{:});
				end

				%% 绘制拟合结果
				if display_slope_area
					figure(f2)
					subplot(4,1,2)
					hold on
					plot(res_list(:,1)./1000,ksn_list(:,4),'-k','LineWidth',2);
					hold off
				else
					figure(f2)
					subplot(3,1,2)
					hold on
					plot(res_list(:,1)./1000,ksn_list(:,4),'-k','LineWidth',2);
					hold off
				end

				f3=figure(3);
				set(f3,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');

				clf
				sbplt2=subplot(2,1,2);
				hold on
				s1=scatter(DAvg./1000,KsnAvg,20,'k','filled');
				p1=plot(res_list(:,1)./1000,ksn_list(:,4)-ksn_list(:,5),':k');
				plot(res_list(:,1)./1000,ksn_list(:,4)+ksn_list(:,6),':k');
				p2=plot(res_list(:,1)./1000,ksn_list(:,4),'-k','LineWidth',2);
				[ksn_vals,~,ksn_ix]=unique(ksn_list(:,4));
				d_means=accumarray(ksn_ix,res_list(:,1),[],@(x) mean(x,'omitnan'));
				for kk=1:numel(ksn_vals)
					text(d_means(kk)./1000,ksn_vals(kk),['k_{sn} = ' num2str(ksn_vals(kk))],...
						'VerticalAlignment','bottom','HorizontalAlignment','center');
				end
				xlabel('距离 (km)')
				ylabel('k_{sn}')
				title('k_{sn} - 距离关系图')
				legend([s1 p1 p2],{'k_{sn}','k_{sn}误差范围','拟合段k_{sn}'},'location','best');
				if ~verLessThan('matlab','9.5')
			        disableDefaultInteractivity(sbplt2);
			    end						
				hold off

				sbplt1=subplot(2,1,1);
				hold on
				plot([min(res_list(:,1)./1000) max(res_list(:,1)./1000)],[0 0],'-k');
				scatter(res_list(:,1)./1000,res_list(:,2),10,'k','filled');
				xlabel('距离 (km)')
				ylabel('残差 (m)')
				title('k_{sn}拟合残差')
				if ~verLessThan('matlab','9.5')
			        disableDefaultInteractivity(sbplt1);
			    end						
				hold off

				elseif strcmp(pick_method,'slope_area')

				[bs,ba,bc,bd,bk,aa,ag,ad,ac]=sa_ksn(DEMc,Sn,A,C.chi,ak,smooth_distance);

				ax3=subplot(4,1,3);
				hold on
				pl1=plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
				pl2=plotdz(Sn,DEMc,'dunit','km','Color','k');
				pl3=scatter((C.distance)./1000,C.elev,5,log10(C.area),'filled');
				xlabel('距河口距离 (km)')
				ylabel('高程 (m)')
				legend([pl1 pl2 pl3],'未平滑DEM','平滑后DEM','对数集水面积','location','best');
				title('纵剖面')
				if ~verLessThan('matlab','9.5')
			        disableDefaultInteractivity(ax3);
			    end						
				hold off

				ax2=subplot(4,1,2);
				hold on
				scatter(ba,bk,20,log10(ba),'filled','MarkerEdgeColor','k');
				xlabel('对数面积')
				ylabel('k_{sn}');
				title('对数面积 - k_{sn}');
				set(ax2,'XScale','log','XDir','reverse');
				if ~verLessThan('matlab','9.5')
			        disableDefaultInteractivity(ax2);
			    end						
				hold off

				ax1=subplot(4,1,1);
				hold on
				plot(C.chi,C.elev,'-k');
				scatter(C.chi,C.elev,10,log10(C.area),'filled');
				xlabel('\chi')
				ylabel('高程 (m)')
				title('\chi - 高程')
				if ~verLessThan('matlab','9.5')
			        disableDefaultInteractivity(ax1);
			    end						
				hold off

				ax4=subplot(4,1,4);
				hold on
				scatter(aa,ag,5,log10(aa),'+');
				scatter(ba,bs,20,log10(ba),'filled','MarkerEdgeColor','k');
				xlabel('对数面积');
				ylabel('对数坡度');
				title(['坡度-面积图: \theta = ' num2str(C.mn) ' : 选取河段 - 完成后按回车'],'Color','r');
				set(ax4,'YScale','log','XScale','log','XDir','reverse');
				ax4.XColor='Red';
				ax4.YColor='Red';
				if ~verLessThan('matlab','9.5')
			        disableDefaultInteractivity(ax4);
			    end						
				hold off

				linkaxes([ax4,ax2],'x');
				colormap(ax1,'jet'); colormap(ax2,'jet'); colormap(ax3,'jet'); colormap(ax4,'jet');

				[av,~,bttn]=ginput;
				% 判断是否存在非边界裂点
				bttn_idx=bttn~=1;
				if any(bttn_idx)
					% 分离非边界点
					av_kn=av(bttn_idx);
					av(bttn_idx)=[];
					% 转换为索引
					ra=C.area; rx=C.x; ry=C.y;
					kn_ix=zeros(numel(av_kn),1); 
					for jj=1:numel(av_kn);
						areadist=sqrt(sum(bsxfun(@minus, ra, av_kn(jj)).^2,2));
						[~,knbix]=min(areadist);
						knbx=rx(knbix);
						knby=ry(knbix);
						kn_ix(jj)=coord2ind(DEM,knbx,knby);
					end
				else
					kn_ix=NaN;
				end	

				if isempty(av)
					if strcmp(theta_method,'ref')
						Cbf=ChiCalc(Sn,DEMc,A,1);
						ksn_list=[C.x C.y C.area ones(numel(C.x),1)*C.ks ones(numel(C.x),1)*C.ks_neg ones(numel(C.x),1)*C.ks_pos ...
							ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*Cbf.mn ones(numel(C.x),1)*snda];
					elseif strcmp(theta_method,'auto')
						ksn_list=[C.x C.y C.area ones(numel(C.x),1)*C.ks ones(numel(C.x),1)*C.ks_neg ones(numel(C.x),1)*C.ks_pos ...
							ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*snda];
					end

					% 确定ksn值在色标中的位置并绘图
					ksn_val=C.ks;

					[~,lbix]=min(C.chi);
					elbl=C.elev(lbix);

					figure(f2)
					subplot(4,1,1);
					hold on
					plot(C.chi,C.pred+elbl,'-k','LineWidth',2);
					hold off

					subplot(4,1,3);
					hold on
					pl4=plot((C.distance)/1000,C.pred+elbl,'-k','LineWidth',2);
					legend([pl1 pl2 pl3 pl4],'未平滑DEM','平滑后DEM','对数集水面积','拟合河段','location','best');
					hold off

					subplot(4,1,4);
					hold on
					plot(C.area,ksn_val.*C.area.^(-C.mn),'-k','LineWidth',2);
					hold off

					res_list=[C.area C.res];
					bnd_ix=NaN;;

				else
					% 对裂点列表排序并构建边界列表
					avs=sortrows(av,'descend');
					bnds=vertcat(max(C.area,[],'omitnan'),avs,min(C.area,[],'omitnan'));

					num_bnds=numel(bnds);
					rc=C.chi;
					rx=C.x;
					ry=C.y;
					rd=C.distance;
					ra=C.area;
					for jj=1:num_bnds-1
						% 提取边界
						lb=bnds(jj);
						rb=bnds(jj+1);

						% 裁剪河道段
						lb_dadist=sqrt(sum(bsxfun(@minus, ra, lb).^2,2));
						rb_dadist=sqrt(sum(bsxfun(@minus, ra, rb).^2,2));

						[~,lbix]=min(lb_dadist);
						[~,rbix]=min(rb_dadist);

						lbx=rx(lbix);
						lby=ry(lbix);

						rbx=rx(rbix);
						rby=ry(rbix);	

						lix=coord2ind(DEM,lbx,lby);
						rix=coord2ind(DEM,rbx,rby);

						Seg=modify(Sn,'downstreamto',rix);
						Seg=modify(Seg,'upstreamto',lix);

						% 重建包含下游边界节点的河道
						WSEG=GRIDobj(DEM,'logical');
						WSEG.Z(Seg.IXgrid)=true;
						WSEG.Z(lix)=true;
						% 如果是河道末端则添加回上游节点
						if jj==num_bnds-1
							WSEG.Z(rix)=true;
						end
						Seg=STREAMobj(FD,WSEG);

						% 检查河道长度，若节点数不足两个则向下游扩展直至满足条件
						lbix_new=lbix+1;
						first_time=true;
						while numel(Seg.IXgrid)<=2
							if first_time
								wrn_mssg=['第' num2str(jj) '段选取的河段过短，河段边界已向下游扩展'];
								wd=warndlg(wrn_mssg);
								uiwait(wd);
							end
							lix_new=coord2ind(DEM,rx(lbix_new),ry(lbix_new));
							WSEG.Z(lix_new)=true;
							Seg=STREAMobj(FD,WSEG);
							lbix_new=lbix_new+1;
							first_time=false;
						end

						% 构建边界索引列表
						if jj<num_bnds-1
							bnd_ix(jj,1)=rix;
						end

						% 计算chi以获取ksn和最佳拟合凹度 
						if strcmp(theta_method,'ref')
							Cseg=ChiCalc(Seg,DEMc,A,1,ref_theta);
						elseif strcmp(theta_method,'auto')
							Cseg=ChiCalc(Seg,DEMc,A,1,C.mn);
						end
						Cbfseg=ChiCalc(Seg,DEMc,A,1);								

						% 确定ksn值在色标中的位置并绘图
						ksn_val=Cseg.ks;

						ksn_nodes{jj,1}=[Cseg.x Cseg.y Cseg.area ones(numel(Cseg.x),1)*Cseg.ks ones(numel(Cseg.x),1)*Cseg.ks_neg ones(numel(Cseg.x),1)*Cseg.ks_pos ...
							ones(numel(Cseg.x),1)*Cseg.mn ones(numel(Cseg.x),1)*Cbfseg.mn ones(numel(Cseg.x),1)*snda];

						% 绘制线性拟合结果	
						rchi=rc(rb_dadist==min(rb_dadist));
						lchi=rc(lb_dadist==min(lb_dadist));
						segChi=linspace(lchi,rchi,numel(Cseg.chi));
						seg0Chi=linspace(0,max(Cseg.chi),numel(Cseg.chi));
						[~,lbix]=min(Cseg.chi);
						elbl=Cseg.elev(lbix);

						figure(f2)
						subplot(4,1,1);
						hold on
						plot(segChi,(seg0Chi.*ksn_val)+elbl,'-k','LineWidth',2);
						hold off

						seg_st=rd(lb_dadist==min(lb_dadist));
						subplot(4,1,3);
						hold on
						pl4=plot((Cseg.distance+seg_st)/1000,(Cseg.pred)+elbl,'-k','LineWidth',2);
						legend([pl1 pl2 pl3 pl4],'未平滑DEM','平滑后的DEM','对数集水面积','拟合河段','location','best');
						hold off

						subplot(4,1,4);
						hold on
						plot(Cseg.area,ksn_val.*Cseg.area.^(-1*Cseg.mn),'-k','LineWidth',2);
						hold off


						res_list{jj,1}=[Cseg.area Cseg.res];


					end

					ksn_list=vertcat(ksn_nodes{:});
					res_list=vertcat(res_list{:});
				end


				%% 绘制结果图
				if display_slope_area
					figure(f2)
					subplot(4,1,2)
					hold on
					plot(res_list(:,1),ksn_list(:,4),'-k','LineWidth',2);
					hold off
				else
					figure(f2)
					subplot(3,1,2)
					hold on
					plot(res_list(:,1),ksn_list(:,4),'-k','LineWidth',2);
					hold off
				end

				f3=figure(3);
				set(f3,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');
				clf
				ax31=subplot(2,1,2);
				hold on
				s1=scatter(log10(ba),bk,20,'k','filled');
				p1=plot(log10(res_list(:,1)),ksn_list(:,4)-ksn_list(:,5),':k');
				plot(log10(res_list(:,1)),ksn_list(:,4)+ksn_list(:,6),':k');
				p2=plot(log10(res_list(:,1)),ksn_list(:,4),'-k','LineWidth',2);
				[ksn_vals,~,ksn_ix]=unique(ksn_list(:,4));
				a_means=accumarray(ksn_ix,res_list(:,1),[],@(x) mean(x,'omitnan'));
				for kk=1:numel(ksn_vals)
					text(log10(a_means(kk)),ksn_vals(kk),['k_{sn} = ' num2str(ksn_vals(kk))],...
						'VerticalAlignment','bottom','HorizontalAlignment','center');
				end
				xlabel('对数面积')
				ylabel('k_{sn}')
				title('k_{sn} - 对数面积')
				legend([s1 p1 p2],{'k_{sn}','k_{sn}误差范围','拟合段k_{sn}'},'location','best');
				set(ax31,'XScale','log','XDir','reverse');
				if ~verLessThan('matlab','9.5')
			        disableDefaultInteractivity(ax31);
			    end						
				hold off

				ax32=subplot(2,1,1);
				hold on
				plot([log10(min(res_list(:,1))) log10(max(res_list(:,1)))],[0 0],'-k');
				scatter(log10(res_list(:,1)),res_list(:,2),10,'k','filled');
				xlabel('对数面积')
				ylabel('残差 (m)')
				title('k_{sn}拟合残差')
				set(ax32,'XScale','log','XDir','Reverse');
				if ~verLessThan('matlab','9.5')
			        disableDefaultInteractivity(ax32);
			    end						
				hold off
			% 结束选取方法分支	
			end

			if ii<num_ch

				ignore_str=['忽略剩余' num2str(num_ch-ii) '条河道'];
				qa2=questdlg('请选择操作','河道拟合',ignore_str,'重新拟合','继续选取','继续选取');

				switch qa2
				case '继续选取'
					str2 = 'Y';
					str1 = [];
					% 添加坡度、残差及河道编号到节点列表并清除空值
					kidx=isnan(ksn_list(:,1));
					ksn_list(kidx,:)=[];
					res_list(kidx,:)=[];
					gix=coord2ind(DEM,ksn_list(:,1),ksn_list(:,2));
					kmat=[ksn_list G.Z(gix) res_list(:,2) ones(size(ksn_list,1),1).*ii];
					ksn_master{ii,1}=kmat;
					bnd_master{ii,1}=bnd_ix;
					kn_master{ii,1}=kn_ix;
					res_master{ii,1}=res_list;		
					Sc=Sct;
					count=ii;
					save(out_restart_name,'ksn_master','bnd_master','kn_master','res_master','Sc','count','-append');						
					if save_figures
						f2_name=[shape_name '_stream_fits_' num2str(ii) '.pdf'];
						f3_name=[shape_name '_stream_rsds_' num2str(ii) '.pdf'];
						print(f2,f2_name,'-dpdf','-fillpage');
						print(f3,f3_name,'-dpdf','-fillpage');
					end
					clear ksn_list ksn_nodes res_list bnd_ix kn_ix kidx gix kmat;
				case '重新拟合'
					str2 = 'R';
					str1 = 'R';
					clear ksn_list ksn_nodes res_list bnd_ix kn_ix;
				case ignore_str
					wtb=waitbar(0,'正在生成最终输出，请勿关闭窗口');
					str1=[];
					str2=[];
					% 添加坡度、残差及河道编号到节点列表并清除空值
					kidx=isnan(ksn_list(:,1));
					ksn_list(kidx,:)=[];
					res_list(kidx,:)=[];
					gix=coord2ind(DEM,ksn_list(:,1),ksn_list(:,2));
					kmat=[ksn_list G.Z(gix) res_list(:,2) ones(size(ksn_list,1),1).*ii];
					ksn_master{ii,1}=kmat;
					bnd_master{ii,1}=bnd_ix;
					kn_master{ii,1}=kn_ix;
					res_master{ii,1}=res_list;		
					Sc=Sct;
					count=ii;
					save(out_restart_name,'ksn_master','bnd_master','kn_master','res_master','Sc','count','-append');
					if save_figures
						f2_name=[shape_name '_stream_fits_' num2str(ii) '.pdf'];
						f3_name=[shape_name '_stream_rsds_' num2str(ii) '.pdf'];
						print(f2,f2_name,'-dpdf','-fillpage');
						print(f3,f3_name,'-dpdf','-fillpage');
					end
					clear ksn_list ksn_nodes res_list bnd_ix kn_ix kidx gix kmat;	
					break_flag=true;
				end

				close figure 2
				close figure 3
			else
				qa2=questdlg('请选择操作','河道拟合','重新拟合','完成处理','完成处理');

				switch qa2
				case '完成处理'
					wtb=waitbar(0,'正在生成最终输出，请勿关闭窗口');
					str2 = 'Y';
					str1 = [];
					% 添加坡度、残差及河道编号到节点列表并清除空值
					kidx=isnan(ksn_list(:,1));
					ksn_list(kidx,:)=[];
					res_list(kidx,:)=[];
					gix=coord2ind(DEM,ksn_list(:,1),ksn_list(:,2));
					kmat=[ksn_list G.Z(gix) res_list(:,2) ones(size(ksn_list,1),1).*ii];
					ksn_master{ii,1}=kmat;
					bnd_master{ii,1}=bnd_ix;
					kn_master{ii,1}=kn_ix;
					res_master{ii,1}=res_list;		
					Sc=Sct;
					count=ii;
					save(out_restart_name,'ksn_master','bnd_master','res_master','Sc','count','-append');
					if save_figures
						f2_name=[shape_name '_stream_fits_' num2str(ii) '.pdf'];
						f3_name=[shape_name '_stream_rsds_' num2str(ii) '.pdf'];
						print(f2,f2_name,'-dpdf','-fillpage');
						print(f3,f3_name,'-dpdf','-fillpage');
					end
					clear ksn_list ksn_nodes res_list bnd_ix kn_ix kidx gix kmat;
				case '重新拟合'
					str2 = 'R';
					str1 = 'R';
					clear ksn_list ksn_nodes res_list bnd_ix kn_ix;
				end

				close figure 2
				close figure 3					
			end
			end
			% 重置循环变量
str1='R';
			
if break_flag
	break
end

end

% 结束交互模式与预设模式的分支判断
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 数据整合、可视化与导出 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch input_method
case 'interactive'
	close figure 1
end

% 为边界点添加河道编号并转换格式
bnd_list=cell(numel(bnd_master),1);
for jj=1:numel(bnd_master)
	bnd=bnd_master{jj,1};
	if ~isnan(bnd)
		% 提取高程值
		bz=DEM.Z(bnd);
		% 获取坐标
		[bx,by]=ind2coord(DEM,bnd);
		% 保留线性索引
		bix=bnd;
		% 添加河道编号
		brm=ones(size(bnd)).*jj;
		% 整合存储
		bnd_list{jj,1}=[bx by bz brm bix];
	elseif isempty(bnd)
		bnd_list{jj,1}=[];
	else
		bnd_list{jj,1}=[0 0 0 jj bnd];
	end
end

% 合并所有边界点数据
bnd_list=vertcat(bnd_list{:});

% 为裂点添加河道编号并转换格式
kn_list=cell(numel(kn_master),1);
for jj=1:numel(kn_master)
	kn=kn_master{jj,1};
	if ~isnan(kn)
		knz=DEM.Z(kn);
		[knx,kny]=ind2coord(DEM,kn);
		knix=kn;
		krm=ones(size(kn)).*jj;
		kn_list{jj,1}=[knx kny knz krm knix];
	elseif isempty(kn)
		kn_list{jj,1}=[];
	else
		kn_list{jj,1}=[0 0 0 jj kn];
	end
end
kn_list=vertcat(kn_list{:});

waitbar(1/5,wtb);

% 构建KSN空间数据结构
if strcmp(junction_method,'check')
	% 整合所有观测点数据
	knl=vertcat(ksn_master{:});
	ix=coord2ind(DEM,knl(:,1),knl(:,2));

	% 创建空栅格
	ksnR=GRIDobj(DEM);
	ksnRn=GRIDobj(DEM);
	ksnRp=GRIDobj(DEM);
	thetaR=GRIDobj(DEM);
	segthetaR=GRIDobj(DEM);
	threshaR=GRIDobj(DEM);
	resR=GRIDobj(DEM);
	rivnumR=GRIDobj(DEM);

	% 填充栅格数据
	ksnR.Z(ix)=knl(:,4);
	ksnRn.Z(ix)=knl(:,5);
	ksnRp.Z(ix)=knl(:,6);
	thetaR.Z(ix)=knl(:,7);
	segthetaR.Z(ix)=knl(:,8);
	threshaR.Z(ix)=knl(:,9);
	resR.Z(ix)=knl(:,11);
	rivnumR.Z(ix)=knl(:,12);

	% 创建河道矢量结构
	KSN=STREAMobj2mapstruct(Sc,'seglength',smooth_distance,'attributes',...
		{'ksn' ksnR @mean 'ksn_neg' ksnRn @mean 'ksn_pos' ksnRp @mean 'uparea' (A.*(A.cellsize^2)) @mean...
		'gradient' G @mean 'theta' thetaR @mean 'seg_theta' segthetaR @mean 'thrsh_ar' threshaR @mean...
		'resid' resR @mean 'riv_num' rivnumR @median});
elseif strcmp(junction_method,'ignore') & strcmp(stack_method,'stack')
	% 分河道存储模式
	num_streams=numel(ksn_master);
	KSN=cell(num_streams,1);
	for ii=1:num_streams
		knl=ksn_master{ii};
		ix=coord2ind(DEM,knl(:,1),knl(:,2));
		WW=GRIDobj(DEM,'logical');
		WW.Z(ix)=true;
		ScT=STREAMobj(FD,WW);

		% 创建空栅格
		ksnR=GRIDobj(DEM);
		ksnRn=GRIDobj(DEM);
		ksnRp=GRIDobj(DEM);
		thetaR=GRIDobj(DEM);
		segthetaR=GRIDobj(DEM);
		threshaR=GRIDobj(DEM);
		resR=GRIDobj(DEM);
		rivnumR=GRIDobj(DEM);

		% 填充数据
		ksnR.Z(ix)=knl(:,4);
		ksnRn.Z(ix)=knl(:,5);
		ksnRp.Z(ix)=knl(:,6);
		thetaR.Z(ix)=knl(:,7);
		segthetaR.Z(ix)=knl(:,8);
		threshaR.Z(ix)=knl(:,9);
		resR.Z(ix)=knl(:,11);
		rivnumR.Z(ix)=knl(:,12);

		KSN{ii}=STREAMobj2mapstruct(ScT,'seglength',smooth_distance,'attributes',...
			{'ksn' ksnR @mean 'ksn_neg' ksnRn @mean 'ksn_pos' ksnRp @mean 'uparea' (A.*(A.cellsize^2)) @mean...
			'gradient' G @mean 'theta' thetaR @mean 'seg_theta' segthetaR @mean 'thrsh_ar' threshaR @mean...
			'resid' resR @mean 'riv_num' rivnumR @median});	
	end
	KSN=vertcat(KSN{:});
elseif strcmp(junction_method,'ignore') & strcmp(stack_method,'average')
	% 平均模式整合数据
	knl=vertcat(ksn_master{:});
	ix=coord2ind(DEM,knl(:,1),knl(:,2));

	% 创建空栅格
	ksnR=GRIDobj(DEM);
	ksnRn=GRIDobj(DEM);
	ksnRp=GRIDobj(DEM);
	thetaR=GRIDobj(DEM);
	segthetaR=GRIDobj(DEM);
	threshaR=GRIDobj(DEM);
	resR=GRIDobj(DEM);
	rivnumR=GRIDobj(DEM);

	% 累加计算平均值
	ksnr=accumarray(ix,knl(:,4),[],@mean);
	ksnrn=accumarray(ix,knl(:,5),[],@mean);
	ksnrp=accumarray(ix,knl(:,6),[],@mean);
	thetar=accumarray(ix,knl(:,7),[],@mean);
	segthetar=accumarray(ix,knl(:,8),[],@mean);
	threshar=accumarray(ix,knl(:,9),[],@mean);
	resr=accumarray(ix,knl(:,11),[],@mean);
	rivnumr=accumarray(ix,knl(:,12),[],@mode);

	% 填充栅格数据
	ksnR.Z(ix)=ksnr;
	ksnRn.Z(ix)=ksnrn;
	ksnRp.Z(ix)=ksnrp;
	thetaR.Z(ix)=thetar;
	segthetaR.Z(ix)=segthetar;
	threshaR.Z(ix)=threshar;
	resR.Z(ix)=resr;
	rivnumR.Z(ix)=rivnumr;

	% 创建矢量结构
	KSN=STREAMobj2mapstruct(Sc,'seglength',smooth_distance,'attributes',...
		{'ksn' ksnR @mean 'ksn_neg' ksnRn @mean 'ksn_pos' ksnRp @mean 'uparea' (A.*(A.cellsize^2)) @mean...
		'gradient' G @mean 'theta' thetaR @mean 'seg_theta' segthetaR @mean 'thrsh_ar' threshaR @mean...
		'resid' resR @mean 'riv_num' rivnumR @mode});
end	

waitbar(2/5,wtb);	

% 创建裂点矢量数据并整理边界输出
idx=~isnan(bnd_list(:,5));
bnd_strc=bnd_list(idx,:);
bnd_list(~idx,1)=NaN; bnd_list(~idx,2)=NaN; bnd_list(~idx,3)=NaN;
bnd_list=bnd_list(:,[1:4]);

if ~isempty(bnd_strc)
	KNK=struct;
	for jj=1:numel(bnd_strc(:,1))
		KNK(jj,1).Geometry='Point';
		KNK(jj,1).X=double(bnd_strc(jj,1));
		KNK(jj,1).Y=double(bnd_strc(jj,2));
		KNK(jj,1).Elev=double(bnd_strc(jj,3));
		KNK(jj,1).StrNum=double(bnd_strc(jj,4));
	end
	out_bound_name=[shape_name '_bounds.shp'];
	shapewrite(KNK,out_bound_name);
end

% 处理裂点数据
idx=~isnan(kn_list(:,5));
kn_strc=kn_list(idx,:);
kn_list(~idx,1)=NaN; kn_list(~idx,2)=NaN; kn_list(~idx,3)=NaN;
kn_list=kn_list(:,[1:4]);

if ~isempty(kn_strc)
	XKNK=struct;
	for jj=1:numel(kn_strc(:,1))
		XKNK(jj,1).Geometry='Point';
		XKNK(jj,1).X=double(kn_strc(jj,1));
		XKNK(jj,1).Y=double(kn_strc(jj,2));
		XKNK(jj,1).Elev=double(kn_strc(jj,3));
		XKNK(jj,1).StrNum=double(kn_strc(jj,4));
	end
	out_knick_name=[shape_name '_knicks.shp'];
	shapewrite(XKNK,out_knick_name);
end

waitbar(3/5,wtb);

% 导出河道矢量文件
out_shape_name=[shape_name '.shp'];
shapewrite(KSN,out_shape_name);

waitbar(4/5,wtb);

% 保存MAT文件
save(out_mat_name,'knl','ksn_master','bnd_list','kn_list','Sc','bnd_master','kn_master','res_master','count','-append');
% 删除临时重启文件
delete(out_restart_name);

waitbar(5/5,wtb);
close(wtb);

% 函数结束
end

function [OUT]=ChiCalc(S,DEM,A,a0,varargin)
	% 改进的Chi剖面计算方法：通过样条插值均衡chi间距，避免高集水区数据点聚集对ksn拟合的影响

	% 获取河道节点总数
	nrc = numel(S.x);
	M   = sparse(double(S.ix),double(S.ixc),true,nrc,nrc);
	% 确定出口点
	outlet = sum(M,2) == 0 & sum(M,1)'~=0;
	if nnz(outlet)>1
	    error('河道网络不能存在多个出口');
	end

	% 节点高程值
	zx   = double(DEM.Z(S.IXgrid));
	% 出口高程
	if nnz(outlet)==0
		zb=min(zx);
	else
		zb   = double(DEM.Z(S.IXgrid(outlet)));
	end
	% 计算Chi积分参数
	a    = double(a0./(A.Z(S.IXgrid)*(A.cellsize.^2)));
	x    = S.distance;
	Lib = true(size(x));

	% 确定最优凹度值
	if isempty(varargin)
	    mn0  = 0.5; % 初始值
	    mn   = fminsearch(@mnfit,mn0);
	else
		mn=varargin{1};
	end

	% 计算Chi值
	chi = netcumtrapz(x,a.^mn,S.ix,S.ixc);

	% 样条插值均衡chi间距
	chiF=chi(Lib);
	zabsF=zx(Lib)-zb;

	% 关闭插值过程中的警告
	warning off

	chiS=linspace(0,max(chiF),numel(chiF)).';
	try
		zS=spline(chiF,zabsF,chiS);
	catch
		zS=zabsF;
		chiS=chiF;
	end

	% 构建输出结构体
	OUT=struct;
	try
		% 线性拟合
		ft=fittype('a*x');
		fobj=fit(chiS,zS,ft,'StartPoint',chiS\zS);
		BETA=coeffvalues(fobj);
		BETA_UNC=confint(fobj);
		OUT.ks   = BETA*a0^mn;
		OUT.ks_neg = (BETA*a0^mn)-(min(BETA_UNC)*a0^mn);
		OUT.ks_pos = (max(BETA_UNC)*a0^mn)-(BETA*a0^mn);
		OUT.mn   = mn;
	catch
		% 拟合失败时的处理
		BETA = chiS\(zS);
		OUT.ks   = BETA*a0^mn;
		OUT.ks_neg = 0;
		OUT.ks_pos = 0;
		OUT.mn   = mn;
	end

	warning on

	% 将结果映射回空间对象
	[OUT.x,...
	 OUT.y,...
	 OUT.chi,...
	 OUT.elev,...
	 OUT.elevbl,...
	 OUT.distance,...
	 OUT.pred,...
	 OUT.area] = STREAMobj2XY(S,chi,DEM,zx-zb,S.distance,BETA*chi,A.*(A.cellsize^2));
	 OUT.res = OUT.elevbl - OUT.pred;

	% 凹度优化的嵌套函数
	function sqres = mnfit(mn)
		CHI = netcumtrapz(x(Lib),a(Lib).^mn,S.ix,S.ixc);
		CHI = CHI ./ max(CHI);
		z   = zx(Lib)-zb;
		z   = z./max(z);
        sqres = sum((CHI - z).^2);
	end
end

function z = netcumtrapz(x,y,ix,ixc)
	% 沿河道网络进行积分计算
	z = zeros(size(x));
	for lp = numel(ix):-1:1;
	    z(ix(lp)) = z(ixc(lp)) + (y(ixc(lp))+(y(ix(lp))-y(ixc(lp)))/2) *(abs(x(ixc(lp))-x(ix(lp))));
	end
end

function [Xavg,Yavg]=BinAverage(X,Y,bin_size);
	% 数据分箱平均计算
	ix=~isnan(X);
	X=X(ix); Y=Y(ix);

	minX=min(X);
	maxX=max(X);

	b=[minX:bin_size:maxX+bin_size];

	try
		[idx]=discretize(X,b);
	catch
		[~,idx]=histc(X,b);
	end

	Xavg=accumarray(idx(:),X,[],@mean);
	Yavg=accumarray(idx(:),Y,[],@mean);
end

function [ksn]=KSN_Quick(DEM,A,S,theta_ref)
	% 快速计算ksn值
	zc=mincosthydrocon(S,DEM,'interp',0.1);
	g=gradient(S,zc);
	G=GRIDobj(DEM);
	G.Z(G.Z==0)=NaN;
	G.Z(S.IXgrid)=g;
	ksn=G./(A.*(A.cellsize^2)).^(-theta_ref);
end

function [bs,ba,bc,bd,a,g,d,C]=sa(DEM,S,A,C,bin_size)
	% 坡度-面积关系分析函数
	minX=min(S.distance);
	maxX=max(S.distance);
	b=[minX:bin_size:maxX+bin_size];

	numbins=round(max([numel(b) numel(S.IXgrid)/10]));

	an=getnal(S,A.*A.cellsize^2);
	z=getnal(S,DEM);
	gn=gradient(S,z,'unit','tangent','method','robust','drop',20);
	gn=smooth(gn,3);

	[~,~,a,g,d]=STREAMobj2XY(S,an,gn,S.distance);
	a(isnan(a))=[];
	g(isnan(g))=[];
	d(isnan(d))=[];
	C(isnan(C))=[];

	mina=min(a);
	maxa=max(a);

    edges = logspace(log10(mina-0.1),log10(maxa+1),numbins+1);
    try
    	[ix]=discretize(a,edges);
    catch
	    [~,ix] = histc(a,edges);
	end

	ba=accumarray(ix,a,[numbins 1],@median,nan);
	bs=accumarray(ix,g,[numbins 1],@(x) mean(x(~isnan(x))),nan);
	bd=accumarray(ix,d,[numbins 1],@mean,nan);
	bc=accumarray(ix,C,[numbins 1],@mean,nan);

	idx=bs>=0 & ba>=0 & bc>=0 & bd>=0;
	bs=bs(idx);
	ba=ba(idx);
	bc=bc(idx);
	bd=bd(idx);

	idx=a>=0 & g>=0 & d>=0 & C>=0;
	a=a(idx);
	g=g(idx);
	d=d(idx);
	C=C(idx);
end

function [bs,ba,bc,bd,bk,a,g,d,C]=sa_ksn(DEM,S,A,C,ak,bin_size);
	% 含ksn的坡度-面积分析函数
	minX=min(S.distance);
	maxX=max(S.distance);
	b=[minX:bin_size:maxX+bin_size];

	numbins=round(max([numel(b) numel(S.IXgrid)/10]));

	an=getnal(S,A.*A.cellsize^2);
	z=getnal(S,DEM);
	gn=gradient(S,z,'unit','tangent','method','robust','drop',20);
	gn=smooth(gn,3);

	[~,~,a,g,d,k]=STREAMobj2XY(S,an,gn,S.distance,ak);
	a(isnan(a))=[];
	g(isnan(g))=[];
	d(isnan(d))=[];
	C(isnan(C))=[];
	k(isnan(k))=[];

	mina=min(a);
	maxa=max(a);

    edges = logspace(log10(mina-0.1),log10(maxa+1),numbins+1);
    try
    	[ix]=discretize(a,edges);
    catch
	    [~,ix] = histc(a,edges);
	end

	ba=accumarray(ix,a,[numbins 1],@median,nan);
	bs=accumarray(ix,g,[numbins 1],@(x) mean(x(~isnan(x))),nan);
	bd=accumarray(ix,d,[numbins 1],@mean,nan);
	bc=accumarray(ix,C,[numbins 1],@mean,nan);
	bk=accumarray(ix,k,[numbins 1],@mean,nan);

	idx=bs>=0 & ba>=0 & bc>=0 & bd>=0 & bk>=0;
	bs=bs(idx);
	ba=ba(idx);
	bc=bc(idx);
	bd=bd(idx);
	bk=bk(idx);

	idx=a>=0 & g>=0 & d>=0 & C>=0 & k>=0;
	a=a(idx);
	g=g(idx);
	d=d(idx);
	C=C(idx);
	k=k(idx);
end

function [Sn]=RedefineThreshold(DEM,FD,A,S,FLUS,ref_theta,pick_method,bin_size,count,figure_flag,shape_name)
	% 动态调整集水面积阈值

	chix=streampoi(S,'channelheads','ix');
	DA=A.*(DEM.cellsize^2);

	UP=dependencemap(FD,chix);
	FLDSt=DEM.*UP;
	[~,ix]=max(FLDSt);
	IX=influencemap(FD,ix);
	ST=STREAMobj(FD,IX);
	z=mincosthydrocon(ST,DEM,'interp',0.1);

	C=chiplot(ST,z,A,'a0',1,'mn',ref_theta,'plot',false);
	[bs,ba,bc,bd,aa,ag,ad,ac]=sa(DEM,ST,A,C.chi,bin_size);

	str11='R';

	while strcmp(str11,'R');
		f4=figure(4);
		set(f4,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');
		clf
		colormap(jet);

		switch pick_method
		case 'chi'
			% Chi剖面选取模式
			ax2=subplot(2,1,2);
			hold on 
			scatter(aa,ag,5,ac,'+');
			scatter(ba,bs,20,bc,'filled','MarkerEdgeColor','k');
			xlabel('对数集水面积');
			ylabel('对数坡度');
			caxis([0 max(C.chi)]);
			set(ax2,'YScale','log','XScale','log','XDir','reverse');
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(ax2);
		    end	
			hold off

			ax1=subplot(2,1,1);
			hold on
			plot(C.chi,C.elev,'-k');
			scatter(C.chi,C.elev,10,C.chi,'filled');
			xlabel('\chi');
			ylabel('高程 (m)');
			title('选择坡面-河道转换点');
			caxis([0 max(C.chi)]);
			ax1.XColor='Red';
			ax1.YColor='Red';
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(ax1);
		    end				
			hold off

			% 获取用户选择的阈值
			[c,~]=ginput(1);
			if ~isempty(c)
				[~,cix]=min(abs(C.chi-c),[],'omitnan');
				a=C.area(cix);
			else
				a=1e6;
			end

			% 可视化筛选结果
			chi_idx=C.area<a;
			sl_idx=ba<a;
			aa_idx=aa<a;

			subplot(2,1,2)
			hold on
			scatter(aa(aa_idx),ag(aa_idx),10,'k','+');
			scatter(ba(sl_idx),bs(sl_idx),30,'k','filled');
			hold off

			subplot(2,1,1)
			hold on
			scatter(C.chi(chi_idx),C.elev(chi_idx),20,'k','filled');
			title('黑色点将被排除在河道定义之外');
			hold off

		case 'slope_area'
			% 坡度-面积图选取模式
			ax2=subplot(2,1,2);
			hold on
			plot(C.chi,C.elev,'-k');
			scatter(C.chi,C.elev,10,C.chi,'filled');
			xlabel('\chi');
			ylabel('高程 (m)');
			caxis([0 max(C.chi)]);
			xlim([0 max(C.chi)+0.5]);
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(ax2);
		    end				
			hold off

			ax1=subplot(2,1,1);
			hold on 
			scatter(aa,ag,5,ac,'+');
			scatter(ba,bs,20,bc,'filled','MarkerEdgeColor','k');
			xlabel('对数集水面积');
			ylabel('对数坡度');
			title('选择坡面-河道转换点');
			caxis([0 max(C.chi)]);
			set(ax1,'YScale','log','XScale','log','XDir','reverse');
			ax1.XColor='Red';
			ax1.YColor='Red';
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(ax1);
		    end				
			hold off

			% 获取阈值并可视化
			[a,~]=ginput(1);
			if isempty(a)
				a=1e6;
			end

			chi_idx=C.area<a;
			sl_idx=ba<a;
			aa_idx=aa<a;

			subplot(2,1,2)
			hold on
			scatter(C.chi(chi_idx),C.elev(chi_idx),20,'k','filled');
			hold off

			subplot(2,1,1)
			hold on
			scatter(aa(aa_idx),ag(aa_idx),10,'k','+');
			scatter(ba(sl_idx),bs(sl_idx),30,'k','filled');
			title('黑色点将被排除在河道定义之外');
			hold off
		end

		% 确认阈值选择
		qa3=questdlg('是否接受新阈值？','设置阈值','重新选择','接受','接受');
		switch qa3
		case '接受'
			str11='C';
		case '重新选择'
			str11='R';
		end
	end

	% 保存阈值设置图
	if figure_flag
		f4_name=[shape_name '_stream_thresh_' num2str(count) '.pdf'];
		print(f4,f4_name,'-dpdf','-fillpage');
	end

	close(f4);

	% 根据新阈值重建河道网络
	da=getnal(ST,A.*A.cellsize^2);
	nix=ST.IXgrid(da>=a);
	IX=GRIDobj(DEM,'logical');
	IX.Z(nix)=true;
	Sn=STREAMobj(FD,IX);
end