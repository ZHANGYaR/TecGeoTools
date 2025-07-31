function [DEMc]=ConditionDEM(DEM,FD,S,method,varargin)
	% 用法：
%	[DEMc]=ConditionDEM(DEM,FD,S,method);
%	[DEMc]=ConditionDEM(DEM,FD,S,method,'name',value,...);
%
% 描述：
%	封装了TopoToolbox中多种河流剖面平滑方法。除'quantc_grid'和'mingrad'外，所有方法仅修改河流网络处的高程。
%	各方法详细说明请参考原始函数。本函数生成对比原始DEM与处理后DEM纵剖面的图形，便于快速评估结果。
%	不同方法在计算复杂度与耗时上差异较大，建议根据需求选择。'mincost'方法可作为初步尝试，
%	注意除'mincost'和'mingrad'外，多数方法需优化工具箱，并行处理可加速运算。若无优化工具箱，
%	可考虑使用编译版本并将输出转换为GRIDobj格式。
%
% 必需输入：
%	DEM  - 数字高程模型（GRIDobj格式），建议为原始数据（如ProcessRiverBasins输出的DEMoc）
%	FD   - 水流方向对象（FLOWobj格式）
%	S    - 河流网络对象（STREAMobj格式）
%	method - DEM修正方法，可选：
%		'mincost'     - 使用mincosthydrocon函数，可选参数'mc_method'和'fillp'
%		'mingrad'     - 使用imposemin函数，可选参数'ming'（注意过大梯度可能导致河床下切异常）
%		'quantc'      - 使用quantcarve（STREAMobj版），可选'tau'、'ming'、'split'，需优化工具箱
%		'quantc_grid' - 使用quantcarve（GRIDobj版），可选'tau'，需优化工具箱，计算量大慎用
%		'smooth'      - 使用smooth函数，可选'sm_method'、'split'、'stiffness'、'stiff_tribs'、'positive'
%		'crs'         - 使用crs函数，可选'stiffness'、'tau'、'ming'、'stiff_tribs'、'knicks'、'split'
%		'crslin'      - 使用crslin函数，可选'stiffness'、'stiff_tribs'、'ming'等系列参数
%
% 可选输入：
%	mc_method [interp] - 'mincost'方法选择，'minmax'或'interp'
%	fillp [0.1]       - 雕刻与填充的平衡参数（0-1），用于'mincost'
%	ming [0]          - 最小下坡梯度（m/m），用于'mingrad'、'quantc'、'crs'、'crslin'
%	tau [0.5]         - 雕刻分位数（0-1），用于'quantc'、'quantc_grid'、'crs'
%	split [true]      - 是否并行处理支流，用于'quantc_grid'、'smooth'、'crs'
%	sm_method ['regularization'] - 'smooth'方法选择，'regularization'或'movmean'
%	stiffness [10]    - 刚度惩罚系数（正值），用于'smooth'、'crs'、'crslin'
%	stiff_tribs [true] - 支流汇合处放松刚度，用于'smooth'、'crs'、'crslin'
%	knicks []         - 陡坎点坐标矩阵（nx2），用于'crs'、'crslin'
%	imposemin [false] - 是否预处理DEM保证最小梯度，用于'crslin'
%	attachtomin [false] - 防止低于剖面最小值，用于'crslin'
%	attachheads [false] - 固定源头高程，用于'crslin'
%	discardflats [false] - 剔除平坦段，用于'crslin'
%	maxcurvature []   - 最大凸曲率限制，用于'crslin'
%	precisecoords []  - 必经点坐标矩阵（nx3），用于'crslin'
%
% 示例：
%		[DEMc]=ConditionDEM(DEM,FD,S,'mincost');
%		[DEMc]=ConditionDEM(DEM,FD,S,'quantc','tau',0.6);
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% 函数作者：Adam M. Forte - 最后更新：2018/06/18        %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	% 解析输入参数
	p = inputParser;
	p.FunctionName = 'ConditionDEM';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
	addRequired(p,'S',@(x) isa(x,'STREAMobj'));
	addRequired(p,'method',@(x) ischar(validatestring(x,{'mincost','mingrad','quantc','quantc_grid','smooth','crs','crslin'})))

	addParameter(p,'mc_method','interp',@(x) ischar(validatestring(x,{'minmax','interp'})));
	addParameter(p,'fillp',0.1,@(x) isscalar(x) && isnumeric(x) && x>=0 && x<=1);
	addParameter(p,'ming',0,@(x) isscalar(x) && isnumeric(x) && x>=0);
	addParameter(p,'tau',0.5,@(x) isscalar(x) && isnumeric(x) && x>=0 && x<=1);
	addParameter(p,'split',true,@(x) islogical(x))
	addParameter(p,'sm_method','regularization',@(x) ischar(validatestring(x,{'regularization','movmean'})));
	addParameter(p,'stiffness',10,@(x) isscalar(x) && isnumeric(x) && x>=0);
	addParameter(p,'stiff_tribs',true,@(x) islogical(x));
	addParameter(p,'positive',true,@(x) islogical(x));
	addParameter(p,'knicks',[],@(x) isnumeric(x) && size(x,2)==2);
	addParameter(p,'imposemin',false,@(x) islogical(x));
	addParameter(p,'attachtomin',false,@(x) islogical(x));
	addParameter(p,'attachheads',false,@(x) islogical(x));
	addParameter(p,'discardflats',false,@(x) islogical(x));
	addParameter(p,'maxcurvature',[],@(x) isnumeric(x) && isscalar(x) || isempty(x));
	addParameter(p,'precisecoords',[],@(x) isnumeric(x) && size(x,2)==3);

	parse(p,DEM,FD,S,method,varargin{:});
	DEM=p.Results.DEM;
	FD=p.Results.FD;
	S=p.Results.S;
	method=p.Results.method;

	mc_method=p.Results.mc_method;
	fillp=p.Results.fillp;
	ming=p.Results.ming;
	tau=p.Results.tau;
	split=p.Results.split;
	sm_method=p.Results.sm_method;
	K=p.Results.stiffness;
	st=p.Results.stiff_tribs;
	po=p.Results.positive;
	knicks=p.Results.knicks;
	% CRSLIN专用参数
	p1=p.Results.imposemin;
	p2=p.Results.attachtomin;
	p3=p.Results.attachheads;
	p4=p.Results.discardflats;
	p5=p.Results.maxcurvature;
	p6=p.Results.precisecoords;


	switch method
	case 'mincost'
		% 最小成本法处理河网高程
		zc=mincosthydrocon(S,DEM,mc_method,fillp);
		DEMc=GRIDobj(DEM);
		DEMc.Z(DEMc.Z==0)=NaN;
		DEMc.Z(S.IXgrid)=zc;
	case 'mingrad'
		% 强制最小梯度处理
		DEMc=imposemin(FD,DEM,ming);
	case 'quantc'
		% 分位数雕刻法（河网版）
		if split
			sp=2;  % 启用并行支流处理
		elseif ~split
			sp=1;  % 禁用并行
		end
		[zc]=quantcarve(S,DEM,tau,'split',sp,'mingradient',ming);
		DEMc=GRIDobj(DEM);
		DEMc.Z(DEMc.Z==0)=NaN;
		DEMc.Z(S.IXgrid)=zc;	
	case 'quantc_grid'
		% 分位数雕刻法（全栅格版）
		DEMc=quantcarve(FD,DEM,tau);
	case 'smooth'
		% 平滑处理法
		zc=smooth(S,DEM,'method',sm_method,'split',split,'K',K,'nstribs',st,'positive',po);
		DEMc=GRIDobj(DEM);
		DEMc.Z(DEMc.Z==0)=NaN;
		DEMc.Z(S.IXgrid)=zc;
	case 'crs'
		% 曲率限制平滑法
		if split
			sp=2;  % 启用并行
		elseif ~split
			sp=0;  % 禁用并行
		end	

		% 处理陡坎点坐标转换
		if isempty(knicks)
			knicksix=[];
		else
			knicksix=coord2ind(DEM,knicks(:,1),knicks(:,2));
		end

		[zc]=crs(S,DEM,'K',K,'tau',tau,'mingradient',ming,'split',sp,'nonstifftribs',st,'knickpoints',knicksix);
		DEMc=GRIDobj(DEM);
		DEMc.Z(DEMc.Z==0)=NaN;
		DEMc.Z(S.IXgrid)=zc;		
	case 'crslin'
		% 线性曲率限制法
		[zc,~,~]=crslin(S,DEM,'K',K,'mingradient',ming,'nonstifftribs',st,'knickpoints',knicks,'imposemin',p1,'attachtomin',p2,'attachheads',p3,'maxcurvature',p4,'precisecoords',p5);
		DEMc=GRIDobj(DEM);
		DEMc.Z(DEMc.Z==0)=NaN;
		DEMc.Z(S.IXgrid)=zc;	
	end

	% 绘制对比图形
	f1=figure(1);
	set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');
	clf

	% 提取主河道
	SL=trunk(klargestconncomps(S,1));

	% 纵剖面对比
	subplot(2,2,1:2)
	hold on
	plotdz(SL,DEM,'color','k');      % 原始DEM剖面
	plotdz(SL,DEMc,'color','r');     % 处理后的DEM剖面
	legend('原始DEM','处理后DEM','最佳位置');
	hold off

	% 高程差剖面
	subplot(2,2,3)
	hold on
	plotdz(SL,DEM-DEMc,'color','k');
	legend('高程差异剖面','最佳位置');
	hold off

	% 高程差空间分布
	subplot(2,2,4)
	hold on
	imagesc(DEM-DEMc);
	colorbar;
	title('高程差异空间分布');
	hold off