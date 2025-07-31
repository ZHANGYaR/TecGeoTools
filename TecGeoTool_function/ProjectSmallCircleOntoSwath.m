function [ds,db,dab]=ProjectSmallCircleOntoSwath(SW,x,y,cx,cy)
	% 用法：
	%	[ds,db,dab]=ProjectSmallCircleOntoSwath(SW,x,y,cx,cy);
	%
	% 描述：
	%	该函数将输入数据沿共享中心点的小圆族投影到扫掠路径上，
	%	实现基于小圆投影的空间数据映射。
	%
	% 必需输入：
	%	SW - SWATH对象，包含扫掠路径信息
	%	x - 待投影数据的x坐标（nx1数组）
	%	y - 待投影数据的y坐标（nx1数组）
	%	cx - 小圆族中心点的x坐标
	%	cy - 小圆族中心点的y坐标
	%	
	% 输出：
	%	ds - nx1数组，投影点沿扫掠路径的距离。包含NaN表示该点投影位置不在扫掠线上
	%	db - nx1数组，投影点到扫掠基线的垂直距离（地图单位），带符号
	%	dab - nx1数组，沿小圆投影的角距离（弧度），带符号
	%
	% 示例：
	%	[ds,db,dab]=ProjectSmallCircleOntoSwath(SW,x,y,cx,cy);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% 函数作者：Adam M. Forte - 最后更新：2019/04/02 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% 从SWATHobj中提取参数
	xypoints=SW.xy;      % 扫掠路径坐标点
	swdist=SW.distx;     % 沿扫掠线的累积距离
	proj=SW.georef.mstruct;  % 地理投影参数

	% 在投影坐标系中计算小圆
	[dlat,dlon]=projinv(proj,x,y);        % 将输入点反投影为地理坐标
	[clat,clon]=projinv(proj,cx,cy);      % 将中心点反投影为地理坐标
	% 生成小圆坐标（500米间隔）
	[sclat,sclon]=scircle2(clat,clon,dlat,dlon,[],[],500);
	[scx,scy]=projfwd(proj,sclat,sclon);  % 将小圆坐标转回投影坐标系

	% 沿扫掠线寻找最近交点
	d=zeros(size(xypoints,1),numel(x));   % 初始化距离矩阵
	w1=waitbar(0,'正在沿小圆查找交点...');  % 创建进度条
	for ii=1:size(xypoints,1)
		% 计算扫掠路径点与小圆点的最小距离
		d(ii,:)=min(hypot(xypoints(ii,1)-scx,xypoints(ii,2)-scy));
		waitbar(ii/size(xypoints,1));     % 更新进度
	end
	close(w1);  % 关闭进度条
	[~,ix]=min(d,[],1);  % 获取最小距离索引
	ds=swdist(ix);       % 转换索引为沿扫掠线距离

	% 计算扫掠中心线与投影点的弧长差
	xsw=xypoints(ix,1);  % 获取最近交点x坐标
	ysw=xypoints(ix,2);  % 获取最近交点y坐标
	% 转换为极坐标
	[tsw,rsw]=cart2pol(xsw-cx,ysw-cy);  % 交点相对中心点的极坐标
	[tp,rp]=cart2pol(x-cx,y-cy);        % 原始点相对中心点的极坐标
	dab=tsw-tp;    % 计算角度差（弧度）
	db=dab.*rp;    % 转换为线性距离（带符号）

	% 边界点过滤：处理超出扫掠范围的点
	if any(ix==1) || any(ix==size(xypoints,1))
		% 获取扫掠线端点地理坐标
		[lat1,lon1]=projinv(proj,xypoints(end,1),xypoints(end,2));
		[lat0,lon0]=projinv(proj,xypoints(1,1),xypoints(1,2));
		% 生成边界小圆
		[eslat,eslon]=scircle2(clat,clon,vertcat(lat1,lat0),vertcat(lon1,lon0));
		[esx,esy]=projfwd(proj,eslat,eslon);
		% 构建复杂多边形区域
		esx=vertcat(esx(:,1),NaN,flipud(esx(:,2)));
		esy=vertcat(esy(:,1),NaN,flipud(esy(:,2)));
		% 判断点是否在有效区域内
		[inp,onp]=inpolygon(x,y,esx,esy);
		in = inp | onp;
		% 过滤外部点
		ds(~in)=NaN;
		db(~in)=NaN;
	end
end