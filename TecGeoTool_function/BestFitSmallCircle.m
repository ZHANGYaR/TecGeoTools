function [clat,clon,crad,cpx,cpy]=BestFitSmallCircle(x,y,proj,varargin)
% 用法：
%   [clat, clon, crad, cpx, cpy] = BestFitSmallCircle(x, y, proj)
%
% 描述：
%   该函数根据一系列x, y坐标找到最佳拟合的小圆。
%
% 必需的输入：
%   x - 拟合点的x坐标，nx1数组
%   y - 拟合点的y坐标，nx1数组
%   proj - 输入x, y坐标的投影（例如，存储在DEM.georef.mstruct中的数据，
%          或旧版TopoToolbox中的DEM.georef）。如果为'proj'参数提供一个空数组，
%          则假定x和y坐标分别是经度和纬度。在这种情况下，输出的'cpx'和'cpy'也将是
%          经度和纬度。
% 
% 可选输入：
%   trim_circle [false] - 逻辑标志，是否修剪圆到提供的x, y点的近似范围，
%                         如果增加num_points超过默认的100，性能会有所提高。
%   num_points [100] - 在生成的小圆上创建的点数。
%
% 输出：
%   clat - 小圆中心的纬度
%   clon - 小圆中心的经度
%   crad - 小圆的半径（以度为单位）
%   cpx - 圆周的x坐标（投影坐标系）
%   cpy - 圆周的y坐标（投影坐标系）
%
% 示例：
%   [clat, clon, crad, cpx, cpy] = BestFitSmallCircle(x, y, DEM.georef.mstruct);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数作者：Adam M. Forte - 最后更新日期：2021年5月25日 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    p = inputParser;
    p.FunctionName = 'BestFitSmallCircle';
    addRequired(p,'x',@(x) isnumeric(x));
    addRequired(p,'y',@(x) isnumeric(x));
    addRequired(p,'proj',@(x) isstruct(x));

    addParameter(p,'trim_circle',false,@(x) islogical(x) && isscalar(x));
    addParameter(p,'num_points',100,@(x) isscalar(x) && isnumeric(x));

    parse(p,x,y,proj,varargin{:});
    x=p.Results.x;
    y=p.Results.y;
    proj=p.Results.proj;

    trim_circle=p.Results.trim_circle;
    num_points=p.Results.num_points;

	% 移除所有NaN值
	idx=~isnan(x) | ~isnan(y);
	x=x(idx);
	y=y(idx);

	% 将输入坐标转换为经纬度
	if ~isempty(proj)
		[lat,lon]=projinv(proj,x,y);
	else
		lon=x; lat=y;
	end

	% 寻找合理的初始点
	mlat=mean(lat);
	mlon=mean(lon);
	d=hypot(mlat-lat,mlon-lon);
	md=max(d);

	% 执行最小化计算
	model=@bfc;
	x0=[mlat mlon md];
	est=fminsearch(model,x0);

	% 提取计算结果
	clat=est(1);
	clon=est(2);
	crad=est(3);

	% 生成投影坐标系下的圆周点
	[circ_lat,circ_lon]=scircle1(clat,clon,crad,[],[],[],num_points);
	if ~isempty(proj)
		[cpx,cpy]=projfwd(proj,circ_lat,circ_lon);
	else
		cpx=circ_lon; cpy=circ_lat;
	end

	% 根据输入点范围修剪圆周
	if trim_circle
		x0=x(1); y0=y(1);
		x1=x(end); y1=y(end);

		d0=hypot((x0-cpx),(y0-cpy));
		d1=hypot((x1-cpx),(y1-cpy));

		[~,ix0]=min(d0);
		[~,ix1]=min(d1);

		if ix0<ix1
			cpx=cpx(ix0:ix1);
			cpy=cpy(ix0:ix1);
		elseif ix0>ix1
			cpx=cpx(ix1:ix0);
			cpy=cpy(ix1:ix0);
		end
	end

	% 最小化目标函数
	function [ssq]=bfc(params)
		p1=params(1);
		p2=params(2);
		p3=params(3);

		% 生成当前参数的圆周
		[lat1,lon1]=scircle1(p1,p2,p3);

		% 计算每个输入点到圆周的最小距离
		for ii=1:numel(lat)
			mind(ii,1)=min(hypot(lat1-lat(ii),lon1-lon(ii)));
		end

		% 计算残差平方和
		ssq=sum(mind.^2);
	end
end