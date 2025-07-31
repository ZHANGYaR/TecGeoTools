function [SW,SwathMat,xypoints,bends]=MakeTopoSwath(DEM,points,width,varargin)
%
% 使用方法：
% [SW,SwathMat,xypoints,bends]=MakeTopoSwath(DEM,points,width);
% [SW,SwathMat,xypoints,bends]=MakeTopoSwath(DEM,points,width,'属性名',值);
%
% 功能描述：
% 封装TopToolbox的SWATHobj功能，用于生成地形剖面分析
%
% 必需输入参数：
% DEM - 用于生成地形剖面的DEM网格对象（GRIDobj类）
% points - n×2矩阵，包含剖面路径点的x,y坐标，至少需要两个点（起点和终点）。
% 首行为起点坐标，后续行表示剖面的转折点。除首尾外的其他点均视为路径转折点。
% 坐标必须与DEM使用相同坐标系且位于DEM范围内（不能是DEM边界坐标）。
% 若输入空数组，将调用SWATHobj的交互界面显示DEM影像供用户手动选点。
% width - 剖面宽度（地图单位）
%
% 可选输入参数：
% sample [] - 沿剖面线的重采样距离（地图单位）。若未指定，默认使用DEM像元大小，
% 此时不进行重采样。
% smooth [0] - 平滑距离（滤波器宽度，地图单位），默认0表示不进行平滑
% vex [10] - 图形显示的垂直夸张系数
% plot_figure [false] - 逻辑标志，控制是否绘制结果图形
% plot_as_points [false] - 逻辑标志，切换为散点图形式显示剖面
% plot_as_heatmap [false] - 逻辑标志，切换为热力图形式显示剖面
% save_figure [false] - 逻辑标志，控制是否将图形保存为PDF（设为true时自动启用plot_figure）
%
% 输出参数：
% SW - TopoToolbox剖面对象（结构体），包含各类剖面信息。可通过plot(SW)绘制路径，
% plotdz(SW)绘制高程剖面图
% SwathMat - n×4矩阵，列依次为：沿剖面距离、最小高程、平均高程、最大高程
% xypoints - n×2矩阵，记录剖面中心线上各采样点的x,y坐标
% bends - 剖面转折点处的累计距离，若无转折点则返回0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 编写：Adam M. Forte - 最后更新日期：2018年6月18日 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 解析输入参数
p = inputParser;
p.FunctionName = 'MakeTopoSwath';
addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
addRequired(p,'points',@(x) isempty(x) || isnumeric(x) & size(x,1)>=2 && size(x,2)==2);
addRequired(p,'width',@(x) isscalar(x) && isnumeric(x));

addParameter(p,'sample',[],@(x) isscalar(x) && isnumeric(x) || isempty(x));
addParameter(p,'smooth',0,@(x) isscalar(x) && isnumeric(x));
addParameter(p,'vex',10,@(x) isscalar(x) && isnumeric(x));
addParameter(p,'plot_figure',false,@(x) isscalar(x) && islogical(x));
addParameter(p,'plot_as_points',false,@(x) isscalar(x) && islogical(x));
addParameter(p,'plot_as_heatmap',false,@(x) isscalar(x) && islogical(x));
addParameter(p,'save_figure',false,@(x) isscalar(x) && islogical(x));
addParameter(p,'make_shape',true,@(x) islogical(x)); % 隐藏参数，用于解决编译版本生成shape文件的问题
addParameter(p,'out_dir',[],@(x) isdir(x));

parse(p,DEM,points,width,varargin{:});
DEM=p.Results.DEM;
points=p.Results.points;
wdth=p.Results.width;

sample=p.Results.sample;
smth=p.Results.smooth;
vex=p.Results.vex;
plot_figure=p.Results.plot_figure;
plot_as_points=p.Results.plot_as_points;
plot_as_heatmap=p.Results.plot_as_heatmap;
save_figure=p.Results.save_figure;
make_shape=p.Results.make_shape;
out_dir=p.Results.out_dir;

% 设置默认采样距离
if isempty(sample)
	sample=DEM.cellsize;
end

% 设置输出目录
if isempty(out_dir)
	out_dir=pwd;
end

% 检查绘图模式冲突
if plot_as_points && plot_as_heatmap
	if isdeployed
		errordlg('请仅设置"plot_as_points"和"plot_as_heatmap"中的一个为真')
	end
	error('请仅设置"plot_as_points"和"plot_as_heatmap"中的一个为真');
end

% 自动开启绘图标志
if save_figure
	plot_figure=true;
end

% 计算剖面转折点距离
num_points=size(points,1);

if num_points>2
	kk=1;
	while kk<num_points-1
		% 计算相邻点间距离
		bx=points(kk,1);
		by=points(kk,2);
		ex=points(kk+1,1);
		ey=points(kk+1,2);
		xx=ex-bx;
		yy=ey-by;
		dist_to_bend(kk)=sqrt((xx^2)+(yy^2));
		kk=kk+1;
	end
	bends=cumsum(dist_to_bend); % 累计转折点距离
else
	bends=0; % 无转折点时设为0
end

% 创建剖面对象
% 处理TopoToolbox不同版本的差异
if isempty(points)
	SWt=SWATHobj(DEM);
	points=SWt.xy0;
	fig=gcf;
	%close(fig); % 关闭交互选点窗口
end

% 尝试不同版本的构造函数
try
	SW=SWATHobj(DEM,points,'width',wdth,'dx',sample,'smooth',smth); % 旧版本语法
catch
	SW=SWATHobj(DEM,points(:,1),points(:,2),'width',wdth,'dx',sample,'smooth',smth); % 新版本语法
end

% 从剖面对象中提取有用值
try
	elevs=cell2mat(SW.Z); % 旧版本数据存储方式
catch 
	elevs=SW.Z; % 新版本直接获取数据
end

% 计算高程统计量
mean_elevs=mean(elevs,'omitnan');
min_elevs=min(elevs,[],'omitnan');
max_elevs=max(elevs,[],'omitnan');

% 获取坐标点信息
try
	xypoints=cell2mat(SW.xy); % 旧版本数据提取
	swdist=cell2mat(SW.distx);
catch
	xypoints=SW.xy; % 新版本直接获取
	swdist=SW.distx;
end

% 构建输出矩阵
SwathMat=[swdist min_elevs.' mean_elevs.' max_elevs.'];

% 绘图处理
if plot_figure
	f1=figure(1);
	clf 
	% 设置图形属性
	set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.4],'renderer','painters');
	hold on

	% 散点图模式
	if plot_as_points
		for ii=1:size(elevs,1)
			scatter(swdist,elevs(ii,:),1,'k','.');
		end

	% 热力图模式
	elseif plot_as_heatmap
		el_range=linspace(min(min_elevs)-1,max(max_elevs)+1,101);
		el_range_p=linspace(min(min_elevs)-1,max(max_elevs)+1,100);
		C=zeros(100,numel(swdist));

		cmap = jet(256);
		cmap(1,:) = 1;
		colormap(cmap);

		% 构建热力矩阵
		for ii=1:numel(swdist)
			[N,~]=histcounts(SW.Z(:,ii),el_range);
			N=N';
			mi=min_elevs(ii);
			ma=max_elevs(ii);
			idx=el_range_p>ma | el_range_p<mi;
			N(idx)=-1;
			C(:,ii)=N;
		end

		imagesc(swdist,el_range_p,C);
		plot(swdist,min_elevs,'-k');
		plot(swdist,max_elevs,'-k');		
	% 默认填充模式
	else
		xx=vertcat(swdist,flipud(swdist));
		yy=horzcat(min_elevs,fliplr(max_elevs));
		patch(xx,yy,[0.8 0.8 0.8]); % 填充高程范围

		plot(swdist,min_elevs,'-k'); % 绘制最小高程线
		plot(swdist,max_elevs,'-k'); % 绘制最大高程线
		plot(swdist,mean_elevs,'-k','LineWidth',2); % 绘制平均高程线
	end

	% 设置坐标比例
	daspect([vex 1 1])

	% 标记转折点位置
	yl=ylim;
	for jj=1:numel(bends)
		plot([bends(jj),bends(jj)],yl,'-k');
	end

	% 设置坐标标签
	xlabel(['沿剖面距离（米） : 垂直夸张度 = ' num2str(vex)]);
	ylabel('高程（米）');
	xlim([0 max(swdist)]);
	% 禁用新版MATLAB的交互缩放
	if ~verLessThan('matlab','9.5')
        disableDefaultInteractivity(gca);
    end 		
	hold off
end

% 图形保存处理
if save_figure
	orient(f1,'Landscape') % 设置为横向
	if plot_as_points
		print(f1,'-dtiff',fullfile(out_dir,'Swath.tif'));
	else
		print(f1,'-dpdf','-bestfit',fullfile(out_dir,'Swath.pdf'));
	end
end

% 生成剖面边界Shape文件
if make_shape
	% 构建shape结构体
	ms=struct;
	ms(1,1).Geometry='Line';
	ms(1,1).X=SW.xy0(:,1);
	ms(1,1).Y=SW.xy0(:,2);
	ms(1,1).Type='中心线';

	% 生成剖面多边形（需MATLAB R2018a及以上）
	if ~verLessThan('matlab','9.4')
		verts=SwathPolygon(SW,wdth);
		ms(2,1).Geometry='Line';
		ms(2,1).X=verts(:,1);
		ms(2,1).Y=verts(:,2);
		ms(2,1).Type='地形宽度';
	end

	shapewrite(ms,fullfile(out_dir,'SwathBounds.shp'));
end

function [verts]=SwathPolygon(SW,w)

	% 提取中心线坐标
cx=SW.xy(:,1);
cy=SW.xy(:,2);

% 计算坐标差分
dx=diff(cx);
dy=diff(cy);

% 计算半宽
w=w/2;

% 计算法线方向
[sw_angle,~]=cart2pol(dx,dy);
sw_angle=vertcat(sw_angle,sw_angle(end));
[px,py]=pol2cart(sw_angle+(pi/2),w);  % 右侧偏移
[mx,my]=pol2cart(sw_angle-(pi/2),w); % 左侧偏移

% 构建多边形边线
swx=[cx+px cx+mx];
swy=[cy+py cy+my];

% 创建并清理多边形
warning off
P=polyshape(vertcat(swx(:,1),flipud(swx(:,2))),vertcat(swy(:,1),flipud(swy(:,2))));
warning on
P=rmholes(P); % 移除孔洞
verts=P.Vertices;
verts=vertcat(verts,verts(1,:)); % 闭合多边形

end

end
