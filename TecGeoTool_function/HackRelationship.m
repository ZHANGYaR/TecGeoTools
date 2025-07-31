function [C,h,drain_area,outs,globC,globh]=HackRelationship(DEM,FD,A,S,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 该函数正在积极开发中，请注意！%%%
%% 函数文档在手册中尚未完全完善 %%
%%%%%% 可能会有更改和故障 %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 用法：
%   [C, h, drain_area, outs, globC, globh] = HackRelationship(DEM, FD, A, S);
%   [C, h, drain_area, outs, globC, globh] = HackRelationship(DEM, FD, A, S, 'name', value, ...);
%
% 描述：
%   该函数用于计算景观中各流域的Hack关系参数。假设输入的STREAMobj在流域出口处终止，
%   可计算Hack系数（C）和指数（h）。
%
% 必需输入：
%   DEM - 数字高程模型（GRIDobj格式）
%   FD  - 水流方向（FLOWobj格式）
%   A   - 流量累积（GRIDobj格式）
%   S   - 河流网络（STREAMobj格式）
%
% 可选输入：
%   method ['trunks']       - 数据拟合方法：'trunks'(主干流)/'streams'(全河网)/'grids'(全栅格)
%   relation ['original']   - Hack关系式：'original'(L=CA^h) 或 'inverse'(A=CL^h)
%   include_hillslope [false] - 是否包含河道源头至分水岭的坡面部分
%   measure_from ['divide'] - 距离测量起点：'divide'(分水岭) 或 'channelheads'(河道源头)
%   draw_fig [false]        - 是否绘制参数随流域面积变化图
%
% 输出：
%   C          - 各流域Hack系数（n×1数组）
%   h          - 各流域Hack指数（n×1数组）
%   drain_area - 各流域总面积（地图单位平方）
%   outs       - 出口位置索引（同streampoi输出的河道出口）
%   globC      - 全局Hack系数（全数据集拟合）
%   globh      - 全局Hack指数（全数据集拟合）
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 作者：Adam M. Forte - 最后更新：2019/05/02       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% 解析输入参数
	p=inputParser;
	p.FunctionName='HackRelationship';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
	addRequired(p,'A', @(x) isa(x,'GRIDobj'));
	addRequired(p,'S',@(x) isa(x,'STREAMobj'));

	addParameter(p,'method','trunks',@(x) ischar(validatestring(x,{'trunks','streams','grids'})));
	addParameter(p,'relation','original',@(x) ischar(validatestring(x,{'original','inverse'})));
	addParameter(p,'include_hillslope',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'measure_from','divide',@(x) ischar(validatestring(x,{'divide','channelheads'})));
	addParameter(p,'draw_fig',false,@(x) isscalar(x) && islogical(x));

	parse(p,DEM,FD,A,S,varargin{:});
	DEM=p.Results.DEM;
	FD=p.Results.FD;
	A=p.Results.A;
	S=p.Results.S;

	method=p.Results.method;
	relation=p.Results.relation;
	include_hillslope=p.Results.include_hillslope;
	measure_from=p.Results.measure_from;
	draw_fig=p.Results.draw_fig;

	% 预计算所需参数
	outs=streampoi(S,'outlets','ix'); % 获取出口索引
	FLDS=flowdistance(FD,'downstream'); % 计算下游水流距离
	DA=A.*A.cellsize^2; % 换算实际流域面积
	DB=drainagebasins(FD,outs); % 划分流域盆地
	num_basins=numel(outs); % 流域数量

	% 根据方法提取节点属性
	switch method
	case 'trunks'
		% 处理主干流数据
		S=trunk(S);
		% 若包含山坡部分则重新计算流网
		if include_hillslope
			FLUS=flowdistance(FD); % 上游水流距离
			chix=streampoi(S,'channelheads','ix'); % 获取河道源头
			ix=zeros(numel(chix),1);
			% 寻找各源头的最远点
			for ii=1:numel(chix)
				chOI=chix(ii);
				UP=dependencemap(FD,chOI); % 依赖区域
				FLUSt=FLUS.*UP;
				[~,ix(ii,1)]=max(FLUSt);
			end
			IX=influencemap(FD,ix); % 影响区域
			S=STREAMobj(FD,IX); % 重建流网
		end
		% 获取节点属性
		danal=getnal(S,DA); % 流域面积
		dnal=getnal(S,FLDS); % 水流距离
		dbnal=getnal(S,DB); % 流域编号
	case 'streams'
		% 处理全河网数据（逻辑同上）
		if include_hillslope
			FLUS=flowdistance(FD);
			chix=streampoi(S,'channelheads','ix');
			ix=zeros(numel(chix),1);
			for ii=1:numel(chix)
				chOI=chix(ii);
				UP=dependencemap(FD,chOI);
				FLUSt=FLUS.*UP;
				[~,ix(ii,1)]=max(FLUSt);
			end
			IX=influencemap(FD,ix);
			S=STREAMobj(FD,IX);
		end
		danal=getnal(S,DA);
		dnal=getnal(S,FLDS);
		dbnal=getnal(S,DB);
	end

	% 分流域计算Hack参数
	for ii=1:num_basins
		switch method
		case 'grids' % 全栅格模式
			IDX=DB==ii; % 当前流域掩膜
			da=DA.Z(IDX.Z); % 流域面积数据
			l=double(FLDS.Z(IDX.Z)); % 水流距离数据
		case {'streams','trunks'} % 河网模式
			IDX=dbnal==ii; % 当前流域节点
			da=danal(IDX);
			l=double(dnal(IDX));
		end

		% 拟合Hack关系
		switch relation
		case 'original' % 标准形式 L=CA^h
			f=fit(da,l,'power1'); % 幂律拟合
		case 'inverse'  % 逆向形式 A=CL^h
			nzidx=l>0; % 过滤零值
			f=fit(l(nzidx),da(nzidx),'power1');
		end

		% 提取拟合参数
		cf=coeffvalues(f);
		C(ii,1)=cf(1); % 系数C
		h(ii,1)=cf(2); % 指数h
		drain_area(ii,1)=max(da); % 最大流域面积
	end

	% 全局拟合
	switch relation
	case 'original'
		switch method
		case 'grids'
			globfit=fit(DA.Z(:),double(FLDS.Z(:)),'power1');
		case {'streams','trunks'}
			globfit=fit(danal,double(dnal),'power1');
		end
	case 'inverse'
		switch method
		case 'grids'
			flds=double(FLDS.Z(:));
			da=DA.Z(:);
			nzidx=flds>0;
			globfit=fit(flds(nzidx),da(nzidx),'power1');
		case {'streams','trunks'}
			nzidx=dnal>0;
			globfit=fit(double(dnal(nzidx)),danal(nzidx),'power1');
		end	
	end

	% 提取全局参数
	globcf=coeffvalues(globfit);
	globC=globcf(1);
	globh=globcf(2);

	% 可视化结果
	if draw_fig
		f1=figure(1);
		clf

		% Hack系数分布
		subplot(2,1,1);
		hold on 
		scatter(drain_area,C,20,'k','filled'); % 单流域结果
		set(gca,'XScale','log'); % 对数坐标
		plot(xlim,[globC globC],'--b'); % 全局拟合
		plot(xlim,[mean(C) mean(C)],'--r'); % 平均值
		xlabel('流域面积');
		ylabel('Hack系数 C');
		legend('单流域拟合','全局拟合','流域均值','Location','best');
		hold off

		% Hack指数分布
		subplot(2,1,2);
		hold on 
		scatter(drain_area,h,20,'k','filled');
		set(gca,'XScale','log');
		plot(xlim,[globh globh],'--b');
		plot(xlim,[mean(h) mean(h)],'--r');
		xlabel('流域面积');
		ylabel('Hack指数 h');
		hold off
	end
end