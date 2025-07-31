function [Sn,thresh_list,xd_list]=FindThreshold(DEM,FD,A,S,num_streams,varargin)
	%
	% 用法:
	%	[Sn,thresh_list,xd_list]=FindThreshold(DEM,FD,A,S,num_streams);
	%	[Sn,thresh_list,xd_list]=FindThreshold(DEM,FD,A,S,num_streams,'name',value,...);
	% 
	% 功能描述:
	% 	该函数用于交互式选择适合的河道提取阈值面积。可通过以下两种方式操作：
	%	1. 当num_streams为数值时，遍历指定数量的河道（从分水岭开始）进行阈值选择，最终使用用户选择阈值的平均值生成新河道网络
	%	2. 当num_streams设为'all'时，遍历所有河道单独选择阈值，可根据'remake_network'参数决定使用个体阈值还是平均阈值
	%	3. 当num_streams设为'auto'时，使用ischange函数在标准化Chi-高程曲线上自动检测冲积-基岩过渡带（需要MATLAB 2017a+）
	%	支持在Chi-高程图或坡度-面积图上进行阈值选择，输出包含所选阈值列表及河源距分水岭距离列表
	%
	% 必需输入参数:
	%	DEM - 数字高程模型GRID对象（建议使用原始DEM）
	%	FD - 流向对象FLOWobj
	%	A - 汇流累积量GRID对象
	%	S - 现有河道网络STREAMobj（原始阈值不限）
	%	num_streams - 选择模式：数值型指定河道数量/'all'处理所有河道/'auto'自动检测模式
	%
	% 可选输入参数:
	%	ref_concavity [0.50] - Chi图参考凹度参数，自动检测时用于计算Chi值
	%	pick_method ['slope_area'] - 选择界面类型：'chi'使用Chi-高程图，'slope_area'使用坡度-面积图
	%	remake_network [false] - 是否使用平均阈值重建网络（当num_streams为'all'或'auto'时）
	%	max_threshold [] - 自动检测时的最大允许阈值面积（平方米）
	%
	% 输出参数:
	%	Sn - 新生成的STREAMobj河道网络
	%	thresh_list - 所有选择的阈值面积列表
	%	xd_list - 所有河源距分水岭的距离列表
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% 函数作者: Adam M. Forte - 最近更新: 2019/02/09        %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% 解析输入参数
	p = inputParser;
	p.FunctionName = 'FitThreshold';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
	addRequired(p,'A',@(x) isa(x,'GRIDobj'));
	addRequired(p,'S',@(x) isa(x,'STREAMobj'));
	addRequired(p,'num_streams',@(x) isnumeric(x) & isscalar(x) || ischar(validatestring(x,{'all','auto'})));

	addParameter(p,'pick_method','slope_area',@(x) ischar(validatestring(x,{'chi','slope_area'})));
	addParameter(p,'ref_concavity',0.50,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'remake_network',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'max_threshold',[],@(x) isscalar(x) && isnumeric(x) || isempty(x));

	parse(p,DEM,FD,A,S,num_streams,varargin{:});
	DEM=p.Results.DEM;
	FD=p.Results.FD;
	A=p.Results.A;
	S=p.Results.S;
	num_streams=p.Results.num_streams;

	pick_method=p.Results.pick_method;
	ref_theta=p.Results.ref_concavity;
	remake_network=p.Results.remake_network;
	mt=p.Results.max_threshold;

	% 获取河源头及流向距离
	chix=streampoi(S,'channelheads','ix');
	FLUS=flowdistance(FD);
	FLDS=flowdistance(FD,'downstream');
	DA=A.*(DEM.cellsize^2);

	% 设置空间分箱尺寸
	bin_size=15*DEM.cellsize;

	% 判断操作模式
	if isnumeric(num_streams)
		op=1;
	elseif strcmp(num_streams,'all')
		op=2;
	elseif strcmp(num_streams,'auto')
		op=3;
	end

	% 版本兼容性检查
	if verLessThan('matlab','9.3') & op==3;
		warning('自动检测功能需要MATLAB 2017a或更高版本，已切换为全河道模式');
		op=2;
	end

	switch op
	case 1
		% 数值模式：选择指定数量河道计算平均阈值

		% 验证河道数量
		num_ch=numel(chix);
		if num_streams>num_ch;
			num_streams=num_ch;
			if isdeployed
				warndlg('拟拟合河道数量超过现有河源数，将处理全部河道')
			else
				warning('拟拟合河道数量超过现有河源数，将处理全部河道')
			end
		end

		% 按河道长度排序（从最长开始处理）
		fl=FLUS.Z(chix);
		[fl,six]=sort(fl,'descend');
		chix=chix(six);

		for ii=1:num_streams
			chOI=chix(ii);

			% 提取当前河道
			UP=dependencemap(FD,chOI);
			FLDSt=FLUS.*UP;
			[~,ix]=max(FLDSt);
			IX=influencemap(FD,ix);
			St=STREAMobj(FD,IX);
			z=mincosthydrocon(St,DEM,'interp',0.1);

			% 计算Chi和坡度-面积参数
			C=chiplot(St,z,A,'a0',1,'mn',ref_theta,'plot',false);
			[bs,ba,bc,bd,aa,ag,ac]=sa(DEM,St,A,C.chi,bin_size);

			% 创建交互式图形界面
			f1=figure(1);
			set(f1,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');
			clf
			colormap(jet);

			switch pick_method
			case 'chi'
				% Chi-高程图模式
				ax2=subplot(2,1,2);
				hold on 
				scatter(aa,ag,5,ac,'+');
				scatter(ba,bs,20,bc,'filled','MarkerEdgeColor','k');
				xlabel('对数集水区面积');
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
				ylabel('高程（米）');
				title(['选择坡面-河道过渡点：剩余 ' num2str(num_streams-ii) ' 条河道待选择']);
				caxis([0 max(C.chi)]);
				ax1.XColor='Red';
				ax1.YColor='Red';
		        if ~verLessThan('matlab','9.5')
		            disableDefaultInteractivity(ax1);
		        end 				
				hold off

				% 获取用户选择的Chi值
				[c,~]=ginput(1);
				[~,cix]=min(abs(C.chi-c),[],'omitnan');
				a=C.area(cix);

				% 计算距分水岭距离
				cx=C.x(cix);
				cy=C.y(cix);
				ccix=coord2ind(DEM,cx,cy);
				xd=FLDS.Z(ccix);

			case 'slope_area'
				% 坡度-面积图模式
				ax2=subplot(2,1,2);
				hold on
				plot(C.chi,C.elev,'-k');
				scatter(C.chi,C.elev,10,C.chi,'filled');
				xlabel('\chi');
				ylabel('高程（米）');
				caxis([0 max(C.chi)]);
		        if ~verLessThan('matlab','9.5')
		            disableDefaultInteractivity(ax2);
		        end 
				hold off

				ax1=subplot(2,1,1);
				hold on 
				scatter(aa,ag,5,ac,'+');
				scatter(ba,bs,20,bc,'filled','MarkerEdgeColor','k');
				xlabel('对数集水区面积');
				ylabel('对数坡度');
				title(['选择坡面-河道过渡点：剩余 ' num2str(num_streams-ii) ' 条河道待选择']);
				caxis([0 max(C.chi)]);
				set(ax1,'YScale','log','XScale','log','XDir','reverse');
				ax1.XColor='Red';
				ax1.YColor='Red';
		        if ~verLessThan('matlab','9.5')
		            disableDefaultInteractivity(ax1);
		        end 
				hold off

				% 获取用户选择的面积阈值
				[a,~]=ginput(1);

				% 计算距分水岭距离
				[~,cix]=min(abs(C.area-a),[],'omitnan');
				cx=C.x(cix);
				cy=C.y(cix);
				ccix=coord2ind(DEM,cx,cy);
				xd=FLDS.Z(ccix);
			end
			close(f1);
			thresh_list(ii,1)=a;
			xd_list(ii,1)=xd;
		end

		% 计算平均阈值并生成新河道网络
		mean_thresh=mean(thresh_list);
		mean_xd=mean(xd_list);

		% 显示统计结果
		f1=figure(1);
		set(f1,'Units','normalized','Position',[0.5 0.1 0.4 0.8],'renderer','painters');
		sbplt1=subplot(2,1,1);
		hold on
		edges=logspace(log10(min(thresh_list)),log10(max(thresh_list)),10);
		histogram(thresh_list,edges);
		plot([mean_thresh,mean_thresh],[0,max(histcounts(thresh_list,edges))],'-k','LineWidth',2);
		xlabel('选择的阈值面积（平方米）');
		a_str=sprintf('%0.1e',mean_thresh);
		title(['平均阈区面积为 ' a_str]);
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(sbplt1);
        end 
		hold off

		sbplt2=subplot(2,1,2);
		hold on
		histogram(xd_list,10);
		plot([mean_xd,mean_xd],[0,max(histcounts(xd_list,10))],'-k','LineWidth',2);
		xlabel('河源头至分水岭的平均距离（米）');
		title(['平均距离 = ' num2str(round(mean_xd,1))]);
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(sbplt2);
        end 		
		hold off

		Sn=STREAMobj(FD,'minarea',mean_thresh,'unit','mapunits');

	case 2
		% 全河道模式：处理所有河道

		% 排序河道并初始化变量
		fl=FLUS.Z(chix);
		[fl,six]=sort(fl,'descend');
		chix=chix(six);
		xd_list=zeros(numel(chix),1);
		thresh_list=zeros(numel(chix),1);
		ix_list=cell(numel(chix),1);

		for ii=1:numel(chix)
			chOI=chix(ii);

			% 提取当前河道
			UP=dependencemap(FD,chOI);
			FLDSt=FLUS.*UP;
			[~,ix]=max(FLDSt);
			IX=influencemap(FD,ix);
			St=STREAMobj(FD,IX);
			z=mincosthydrocon(St,DEM,'interp',0.1);

			% 计算地形参数
			C=chiplot(St,z,A,'a0',1,'mn',ref_theta,'plot',false);
			[bs,ba,bc,bd,aa,ag,ac]=sa(DEM,St,A,C.chi,bin_size);

			% 创建交互界面
			f1=figure(1);
			set(f1,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');
			clf
			colormap(jet);

			switch pick_method
			case 'chi'
				% Chi-高程图交互
				ax2=subplot(2,1,2);
				hold on 
				scatter(aa,ag,5,ac,'+');
				scatter(ba,bs,20,bc,'filled','MarkerEdgeColor','k');
				xlabel('对数集水区面积');
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
				ylabel('高程（米）');
				title(['选择坡面-河道过渡点：剩余 ' num2str(numel(chix)-ii) ' 条河道待选择']);
				caxis([0 max(C.chi)]);
				ax1.XColor='Red';
				ax1.YColor='Red';
		        if ~verLessThan('matlab','9.5')
		            disableDefaultInteractivity(ax1);
		        end 				
				hold off

				% 获取用户选择
				[c,~]=ginput(1);
				[~,cix]=min(abs(C.chi-c),[],'omitnan');
				a=C.area(cix);

				% 记录河道节点
				allx=C.x(C.area>=a);
				ally=C.y(C.area>=a);
				ix_list{ii,1}=coord2ind(DEM,allx,ally);

			case 'slope_area'
				% 坡度-面积图交互
				ax2=subplot(2,1,2);
				hold on
				plot(C.chi,C.elev,'-k');
				scatter(C.chi,C.elev,10,C.chi,'filled');
				xlabel('\chi');
				ylabel('高程（米）');
				caxis([0 max(C.chi)]);
		        if ~verLessThan('matlab','9.5')
		            disableDefaultInteractivity(ax2);
		        end 
				hold off

				ax1=subplot(2,1,1);
				hold on 
				scatter(aa,ag,5,ac,'+');
				scatter(ba,bs,20,bc,'filled','MarkerEdgeColor','k');
				xlabel('对数集水区面积');
				ylabel('对数坡度');
				title(['选择坡面-河道过渡点：剩余 ' num2str(numel(chix)-ii) ' 条河道待选择']);
				caxis([0 max(C.chi)]);
				set(ax1,'YScale','log','XScale','log','XDir','reverse');
				ax1.XColor='Red';
				ax1.YColor='Red';
		        if ~verLessThan('matlab','9.5')
		            disableDefaultInteractivity(ax1);
		        end 
				hold off

				% 获取用户选择
				[a,~]=ginput(1);

				[~,cix]=min(abs(C.area-a),[],'omitnan');
				cx=C.x(cix);
				cy=C.y(cix);
				ccix=coord2ind(DEM,cx,cy);
				xd=FLDS.Z(ccix);
                
				% 记录河道节点
				allx=C.x(C.area>=a);
				ally=C.y(C.area>=a);
				ix_list{ii,1}=coord2ind(DEM,allx,ally);
			end
			close(f1);
			thresh_list(ii,1)=a;
			xd_list(ii,1)=xd;
		end

		% 显示统计直方图
		f1=figure(1);
		set(f1,'Units','normalized','Position',[0.5 0.1 0.4 0.8],'renderer','painters');
		subplot(2,1,1);
		hold on
		edges=logspace(log10(min(thresh_list)),log10(max(thresh_list)),10);
		histogram(thresh_list,edges);
		xlabel('选择的阈值面积（平方米）');
		hold off

		subplot(2,1,2);
		hold on
		histogram(xd_list,10);
		xlabel('河源头至分水岭的距离（米）');
		hold off

		% 生成最终河道网络
		if remake_network
			Sn=STREAMobj(FD,'minarea',mean(thresh_list),'unit','mapunits');
		else
			W=GRIDobj(DEM,'logical');
			W.Z(unique(vertcat(ix_list{:})))=true;
			Sn=STREAMobj(FD,W);
		end

	case 3
		% 自动检测模式
		xd_list=zeros(numel(chix),1);
		thresh_list=zeros(numel(chix),1);
		ix_list=cell(numel(chix),1);

		w1=waitbar(0,'基于Chi-高程拐点自动寻找河道阈值...');
		for ii=1:numel(chix)
			% 提取当前河道
			chOI=chix(ii);
			UP=dependencemap(FD,chOI);
			FLDSt=FLUS.*UP;
			[~,ix]=max(FLDSt);
			IX=influencemap(FD,ix);
			St=STREAMobj(FD,IX);
			z=mincosthydrocon(St,DEM,'interp',0.1);

			% 计算Chi参数
			c=chitransform(St,A,'mn',ref_theta,'a0',1);
			da=getnal(St,DA);

			% 标准化处理
			[cs,six]=sort(c);
			zs=z(six);
			das=da(six);
			zn=(zs-min(zs))/(max(zs)-min(zs));
			cn=cs/max(cs);

			% 拐点检测
			TF=ischange(zn,'linear','SamplePoints',cn,'Threshold',0.05);
			chi_ix=find(TF,1,'last');

			% 阈值确定
			if isempty(chi_ix)
				ta=DA.Z(chOI);
			else 
				ta=das(chi_ix);
				if ~isempty(mt) && ta>mt
					ta=mt;
				end
			end

			% 记录结果
			thresh_list(ii,1)=ta;
			[~,loc_ix]=min(abs(da-ta));
			ix_list{ii,1}=St.IXgrid(da>=ta);			
			waitbar(ii/numel(chix),w1,['正在处理第 ' num2str(ii) ' 个河流源头...']);
		end
		close(w1);

		% 显示统计结果
		f1=figure(1);
		set(f1,'Units','normalized','Position',[0.5 0.1 0.4 0.8],'renderer','painters');
		subplot(2,1,1);
		histogram(thresh_list,logspace(log10(min(thresh_list)),log10(max(thresh_list)),10));
		set(gca,'XScale','log');
		xlabel('自动检测的阈值面积（平方米）');

		subplot(2,1,2);
		histogram(xd_list,50);
		xlabel('河源头至分水岭的距离（米）');

		% 生成最终河道网络
		if remake_network
			Sn=STREAMobj(FD,'minarea',mean(thresh_list),'unit','mapunits');
		else
			W=GRIDobj(DEM,'logical');
			W.Z(unique(vertcat(ix_list{:})))=true;
			Sn=STREAMobj(FD,W);
		end
	end
end

function [bs,ba,bc,bd,a,g,C]=sa(DEM,S,A,C,bin_size)
	% 坡度-面积分析辅助函数
	% 输入参数：
	% DEM - 高程模型
	% S - 河道对象
	% A - 汇流面积
	% C - Chi值
	% bin_size - 空间分箱尺寸

	minX=min(S.distance);
	maxX=max(S.distance);
	b=[minX:bin_size:maxX+bin_size];
	numbins=round(max([numel(b) numel(S.IXgrid)/10]));

	% 计算地形参数
	an=getnal(S,A.*A.cellsize^2);
	z=getnal(S,DEM);
	gn=gradient(S,z,'unit','tangent','method','robust','drop',20);
	gn=smooth(gn,3);

	% 数据预处理
	[~,~,a,g,d]=STREAMobj2XY(S,an,gn,S.distance);
	a(isnan(a))=[];
	g(isnan(g))=[];
	d(isnan(d))=[];
	C(isnan(C))=[];

	% 分箱统计
	edges = logspace(log10(min(a)),log10(max(a)),numbins+1);
	[~,ix] = histc(a,edges);
	ba=accumarray(ix,a,[numbins 1],@median,nan);
	bs=accumarray(ix,g,[numbins 1],@(x) mean(x(~isnan(x))),nan);
	bd=accumarray(ix,d,[numbins 1],@mean,nan);
	bc=accumarray(ix,C,[numbins 1],@mean,nan);

	% 数据过滤
	valid=bs>=0 & ba>=0 & bc>=0 & bd>=0;
	bs=bs(valid);
	ba=ba(valid);
	bc=bc(valid);
	bd=bd(valid);
end