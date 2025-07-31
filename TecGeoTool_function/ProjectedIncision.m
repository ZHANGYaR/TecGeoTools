function [S,zpOUT,inOUT]=ProjectedIncision(DEM,A,S,Sc,OUT,basin_num,varargin)
	% 用法:
	%	[S,zpOUT,inOUT]=ProjectionIncision(DEM,A,S,Sc,basin_num,OUT);
	%	[S,zpOUT,inOUT]=ProjectionIncision(DEM,A,S,Sc,OUT,basin_num,'name',value);
	%
	% 功能描述:
	% 	根据SegmentProjector的河流剖面投影结果，生成整个网络的投影下切图。本函数利用拟合河段的陡峭指数和全网的χ值，
	%	计算原始河流网络的投影高程及隐含下切量。假设在SegmentProjector中拟合的河段代表被抬升切割的前低起伏地形，
	%	本代码首先计算若无下切（仅有地表抬升）时的河流高程（存储在'zpOUT'节点属性数组中）。随后通过当前河网高程与投影高程的差值
	%	计算下切量（存储在'inOUT'节点属性数组中）。若投影了多条河流（存储在OUT元胞数组中），代码将分别计算每条河流的投影高程和下切量，
	%	最终输出全网的均值、标准差、最小值和最大值。
	%
	%	注意：若Sc输入包含多个出口点，投影高程和下切量仅计算满足以下条件的河段：
	%	1) S网络中位于Sc出口上游的部分 
	%	2) 与用于投影的河道相连的部分
	%	例如，若Sc定义了两个子网络A和B，且分别投影了A中的2条、B中的3条河流，则网络A部分的下切量仅基于A中的投影结果计算。
	%
	%	注：负的下切值表示现代河网高程高于投影高程。本方法假设地表抬升期间流域面积/网络拓扑未发生改变。
	% 
	% 输入参数:
	%	DEM - 数字高程模型，GRIDobj对象（建议使用ProcessRiverBasins输出的未处理DEMoc）
	%	A - 汇流累积量GRIDobj
	%	S - 与DEM和A对应的完整STREAMobj
	%	Sc - 传递给SegmentProjector的子集STREAMobj
	%	OUT - SegmentProjector的输出元胞数组
	%	basin_num - 流域编号，用于输出文件命名标识
	%
	% 可选参数:
	%	display_figure [true] - 逻辑标志，是否显示计算下切量均值、最大值、最小值及标准差的图表
	%	exlcude_streams [] - 排除指定编号的投影河流（例如[2 15]表示排除OUT中第2和第15条投影河流）
	%	conditioned_DEM [] - 可选提供经过水文校正的DEM（主DEM参数请勿使用校正后的DEM！）
	%	interp_value [0.1] - mincosthydrocon函数的插值参数（未提供conditioned_DEM时生效）
	%	out_dir [] - 指定输出目录（默认当前目录）
	%
	% 输出参数:
	%	S - 经裁剪的STREAMobj，仅保留Sc出口上游的河网
	%	zpOUT - n×4数组，投影高程的均值、标准差、最小值和最大值
	%	inOUT - n×4数组，下切量的均值、标准差、最小值和最大值
	%
	% 输出文件:
	%	生成两个shapefile：
	%	'*_Pnts_Used.shp' - 点状文件，记录用于计算的投影河道源头位置
	%	'*_Proj_Incision.shp' - 线状文件，包含zpOUT和inOUT的属性字段
	%
	% 示例:
	%	[S,zpOUT,inOUT]=ProjectedIncision(DEM,A,S,Sc,OUT);
	%	[S,zpOUT,inOUT]=ProjectedIncision(DEM,A,S,Sc,OUT,'exclude_streams',[2 15 40]);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% 函数作者: Adam M. Forte - 最后更新: 2019/04/02        %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	p = inputParser;
	p.FunctionName = 'ProjectedIncision';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'A',@(x) isa(x,'GRIDobj'));
	addRequired(p,'S',@(x) isa(x,'STREAMobj'));
	addRequired(p,'Sc',@(x) isa(x,'STREAMobj'));
	addRequired(p,'OUT',@(x) iscell(x));
	addRequired(p,'basin_num',@(x) isnumeric(x));	

	addParameter(p,'display_figure',true,@(x) isscalar(x) && islogical(x));
	addParameter(p,'exclude_streams',[],@(x) isnumeric(x) || isempty(x));
	addParameter(p,'conditioned_DEM',[],@(x) isa(x,'GRIDobj') || isempty(x));
	addParameter(p,'interp_value',0.1,@(x) isnumeric(x) && x>=0 && x<=1);
	addParameter(p,'out_dir',[],@(x) isdir(x));

	parse(p,DEM,A,S,Sc,OUT,basin_num,varargin{:});
	DEM=p.Results.DEM;
	A=p.Results.A;
	S=p.Results.S;
	Sc=p.Results.Sc;
	OUT=p.Results.OUT;
	basin_num=p.Results.basin_num;

	display_figure=p.Results.display_figure;
	exclude_streams=p.Results.exclude_streams;
	DEMc=p.Results.conditioned_DEM;
	iv=p.Results.interp_value;
	out_dir=p.Results.out_dir;

	if isempty(out_dir)
		out_dir=pwd;
	end

	% 如果没有提供经过条件处理的DEM，则进行处理
	if isempty(DEMc)
		zc=mincosthydrocon(S,DEM,'interp',iv);
		DEMc=GRIDobj(DEM);
		DEMc.Z(DEMc.Z==0)=NaN;
		DEMc.Z(S.IXgrid)=zc;
	end

	% 根据exclude_streams参数过滤投影河流
	ref_num=[1:size(OUT,2)];
	if ~isempty(exclude_streams)
		idx=logical(ones(1,size(OUT,2)));
		idx(exclude_streams)=false;
		OUT=OUT(:,idx);
		ref_num=ref_num(idx);
	end

	% 获取选定河流的出口点
	outix=streampoi(Sc,'outlets','ix');
	% 裁剪主河流网络
	S=modify(S,'upstreamto',outix);
	% 创建连接分量标识网格
	[L,nc]=conncomps(S);
	LG=GRIDobj(DEM);
	LG.Z(S.IXgrid)=L;

	% 初始化存储数组
	num_proj=size(OUT,2);
	num_nodes=numel(S.x);
	zpM=zeros(num_nodes,num_proj);
	inM=zeros(num_nodes,num_proj);
	ch=struct;

	% 主计算循环：逐个处理投影河流
	for ii=1:num_proj
		% 提取投影数据
		chx=OUT{1,ii}(:,1);
		chy=OUT{1,ii}(:,2);
		c=OUT{2,ii}(:,5);
		mn=OUT{2,ii}(:,6); mn=mn(1);
		zp=OUT{2,ii}(:,8);

		% 确定所属子网络
		chix=coord2ind(DEM,chx,chy);
		loi=LG.Z(chix);

		% 清理无效数据
		idx=~isnan(c);
		c=c(idx); zp=zp(idx);

		% 计算投影出口高程
		outzp=min(zp);

		% 计算等效陡峭指数
		ksn=c\(zp-outzp);

		% 计算全网投影高程
		cnal=chitransform(S,A,'a0',1,'mn',mn);
		zpnal=cnal.*ksn;
		zpnal=zpnal+outzp;

		% 计算下切量
		inc=zpnal-getnal(S,DEMc);

		% 屏蔽非相关网络数据
		lidx=L~=loi;
		zpnal(lidx)=NaN;
		inc(lidx)=NaN;

		% 存储结果
		zpM(:,ii)=zpnal;
		inM(:,ii)=inc;

		% 构建河道起点结构
		ch(ii,1).Geometry='Point';
		ch(ii,1).X=double(chx);
		ch(ii,1).Y=double(chy);
		ch(ii,1).rivID=double(ref_num(ii));
	end

	% 计算统计量
	mean_in=mean(inM,2,'omitnan');
	std_in=std(inM,0,2,'omitnan');
	min_in=min(inM,[],2,'omitnan');
	max_in=max(inM,[],2,'omitnan');

	mean_zp=mean(zpM,2,'omitnan');
	std_zp=std(zpM,0,2,'omitnan');
	min_zp=min(zpM,[],2,'omitnan');
	max_zp=max(zpM,[],2,'omitnan');

	zpOUT=[mean_zp std_zp min_zp max_zp];
	inOUT=[mean_in std_in min_in max_in];

	% 可视化模块
	if display_figure
		f1=figure(1);
		set(f1,'unit','normalized','position',[0.1 0.1 0.8 0.8]);

		sbplt1=subplot(2,2,1);
		hold on
		[RGB]=imageschs(DEM,DEM,'colormap','gray');
		[~,R]=GRIDobj2im(DEM);
		imshow(flipud(RGB),R);
		axis xy
		plotc(S,mean_in);
		colorbar;
		title('平均下切量');
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(sbplt1);
        end 
		hold off

		sbplt2=subplot(2,2,2);
		hold on
		imshow(flipud(RGB),R);
		axis xy
		plotc(S,std_in);
		colorbar;
		title('下切量标准差');
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(sbplt2);
        end 		
		hold off

		sbplt3=subplot(2,2,3);
		hold on 
		imshow(flipud(RGB),R);
		axis xy
		plotc(S,min_in);
		colorbar;
		title('最小下切量');
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(sbplt3);
        end 		
		hold off

		sbplt4=subplot(2,2,4);
		hold on 
		imshow(flipud(RGB),R);
		axis xy
		plotc(S,max_in);
		colorbar;
		title('最大下切量');
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(sbplt4);
        end 		
		hold off
	end	

	% 输出shapefile
	ms=STREAMobj2mapstruct(S,'seglength',DEM.cellsize*3,'attributes',{'mean_inc' mean_in @mean 'std_inc' std_in @mean...
		'min_inc' min_in @mean 'max_inc' max_in @mean 'mean_zp' mean_zp @mean 'std_zp' std_zp @mean 'min_zp' min_zp @mean...
		'max_zp' max_zp @mean});

	pnts_name=fullfile(out_dir,['ProjectedIncision_' num2str(basin_num) '_Pnts_Used.shp']);
	strm_name=fullfile(out_dir,['ProjectedIncision_' num2str(basin_num) '_Map.shp']);

	shapewrite(ch,pnts_name);
	shapewrite(ms,strm_name);

end