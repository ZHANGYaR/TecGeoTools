function [ds,db]=ProjectOntoSwath(SW,x,y,data_width,varargin)
%
% 用法：
%   [ds,db]=ProjectOntoSwath(SW,x,y,data_width);
%   [ds,db]=ProjectOntoSwath(SW,x,y,data_width,'参数名',参数值);   
%
% 描述：
%   该函数将空间点投影到SWATHobj扫掠对象上，计算沿扫掠线的距离（ds）和
%   到中心线的垂直距离（db），主要用于'MakeCombinedSwath'功能。
%
% 必需输入：
%   SW - SWATHobj扫掠对象
%   x - 待投影点的x坐标数组
%   y - 待投影点的y坐标数组
%   data_width - 从扫掠中心线算起的最大采样宽度（地图单位）
%
% 可选参数：
%   signed [默认false] - 控制距离符号的标志：
%       false：db返回绝对值
%       true：db带符号，正值表示扫掠线右侧（沿节点方向），负值表示左侧
%   include_concave_bend_regions [默认true] - 控制是否包含内凹弯曲三角区域的标志
%
% 输出：
%   ds - 各点沿扫掠线的投影距离
%   db - 各点到中心线的垂直距离（地图单位）
% 
% 注意：
%   无法投影的点将返回NaN
%   扫掠线的急转弯处可能出现异常结果（尤其在计算符号距离时）
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数作者：Adam M. Forte - 最后更新：2020/09/19 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% 解析输入参数
	p = inputParser;
	p.FunctionName = 'ClassifyKnicks';
	addRequired(p,'SW',@(x) isa(x,'SWATHobj'));
	addRequired(p,'x',@(x) isnumeric(x));
	addRequired(p,'y',@(x) isnumeric(x));
	addRequired(p,'data_width',@(x) isnumeric(x) && isscalar(x));

	addParameter(p,'signed',false,@(x) islogical(x) && isscalar(x));
	addParameter(p,'include_concave_bend_regions',true,@(x) islogical(x) && isscalar(x));

	parse(p,SW,x,y,data_width,varargin{:});
	SW=p.Results.SW;
	x=p.Results.x;
	y=p.Results.y;
	data_width=p.Results.data_width;

	signed=p.Results.signed;  % 符号距离标志
	in_cn_bnds=p.Results.include_concave_bend_regions;  % 包含凹区标志

	% 从SWATHobj中提取扫掠参数
	xypoints=SW.xy;     % 扫掠路径坐标点
	swdist=SW.distx;    % 沿扫掠线累积距离
	xy0=SW.xy0;        % 扫掠拐点坐标

	% 寻找拐点索引
	try
		% 优先使用统计工具箱的高效方法
		for kk=1:numel(SW.xy0(:,1))
			[~,bend_ix(kk,1)]=min(pdist2(SW.xy,SW.xy0(kk,:)));
		end
	catch 
		% 备用方法：手动计算欧氏距离
		for kk=1:numel(SW.xy0(:,1))
			d=zeros(numel(SW.xy(:,1)),1);
			for ll=1:numel(SW.xy(:,1))
				d(ll)=hypot(SW.xy(ll,1)-SW.xy0(kk,1),SW.xy(ll,2)-SW.xy0(kk,2));
			end
			[~,bend_ix(kk,1)]=min(d);
		end
	end

	% 处理扫掠分段
	swxB=xy0(:,1); swyB=xy0(:,2);  % 拐点坐标
	num_segs=numel(swxB)-1;        % 计算扫掠段数量

	% 计算段长度和累积距离
	if num_segs>1
		kk=1;
		while kk<=num_segs
			% 计算段端点坐标差
			x0=swxB(kk,1); y0=swyB(kk,1);
			x1=swxB(kk+1,1); y1=swyB(kk+1,1);
			xx=x1-x0; yy=y1-y0;
			dist_to_bend(kk,1)=sqrt((xx^2)+(yy^2));  % 单段长度
			kk=kk+1;
		end
		dist_to_bend=vertcat(0,dist_to_bend);  % 添加起始点
		bends=cumsum(dist_to_bend);            % 累积距离
	else
		% 单段处理
		bends=[0;max(swdist)];
		dist_to_bend=bends;
	end

	% 获取扫掠线全部点坐标
	swx=xypoints(:,1); swy=xypoints(:,2);

	% 初始化输出矩阵
	dist_in_swath=zeros(numel(x),num_segs);  % 沿扫掠距离
	dist_from_base=zeros(numel(x),num_segs); % 基线距离

	% 主处理循环：逐段处理
	for ii=1:num_segs
		% 获取当前段端点
		swx0=swxB(ii); swx1=swxB(ii+1);
		swy0=swyB(ii); swy1=swyB(ii+1);

		% 计算段内各点距离
		x_dist=swx-swx0; y_dist=swy-swy0;
		seg_dist=sqrt((x_dist.^2)+(y_dist.^2));
		% 处理段间重叠区域
		if ii>1 && ii<num_segs
			seg_dist([1:bend_ix(ii) bend_ix(ii+1)+1:end])=NaN;
		elseif ii==1
			seg_dist(bend_ix(ii+1)+1:end)=NaN;
		else
			seg_dist(1:bend_ix(ii))=NaN;
		end

		% 计算到结束点的距离
		xn_dist=x-swx1; yn_dist=y-swy1;
		n_dist=sqrt((xn_dist.^2)+(yn_dist.^2));

		% 计算段方向角（弧度）
		sw_angle=-1*atan((swy1-swy0)/(swx1-swx0));

		% 坐标旋转变换
		[n_swx,n_swy]=RotCoord(swx,swy,sw_angle,swx0,swy0);  % 扫掠线旋转
		[n_x,n_y]=RotCoord(x,y,sw_angle,swx0,swy0);         % 数据点旋转

		% 遍历所有数据点
		num_points=numel(n_x);
		pd_in_seg=zeros(num_points,1);       % 段内投影距离
		DistFromBaseLine=zeros(num_points,1); % 基线距离
		for jj=1:num_points
			xoi=n_x(jj);  % 当前点旋转后x坐标

			% 寻找最近投影点
			if abs(xoi)<=max(abs(n_swx))
				distances=n_swx-xoi;
				[~,I]=min(abs(distances)); 
				pd_in_seg(jj)=seg_dist(I);  % 记录段内距离
				% 计算基线距离（带符号或绝对值）
				BaseLine=n_swy(I);
				if signed
					DistFromBaseLine(jj)=BaseLine-n_y(jj);  % 带符号距离
				else
					DistFromBaseLine(jj)=abs(BaseLine-n_y(jj));  % 绝对距离
				end
			else 
				pd_in_seg(jj)=NaN;
				DistFromBaseLine(jj)=NaN;
			end

			% 处理弯曲重叠区域
			if in_cn_bnds
				if isnan(pd_in_seg(jj)) && n_dist(jj)<=data_width && ii~=num_segs
					pd_in_seg(jj)=seg_dist(bend_ix(ii+1));
					if signed
						% 确定符号方向
						[n_swx1,n_swy1]=RotCoord(swx1,swy1,sw_angle,swx0,swy0);
						if n_swy1>n_y(jj)
							DistFromBaseLine(jj)=n_dist(jj);  % 右侧正值
						else 
							DistFromBaseLine(jj)=-1*n_dist(jj);  % 左侧负值
						end
					else
						DistFromBaseLine(jj)=n_dist(jj);  % 记录绝对值
					end
				end
			end 
		end

		% 过滤无效点
		idx=isnan(pd_in_seg) | pd_in_seg<=0;	
		pd_in_seg(idx)=NaN;
		DistFromBaseLine(idx)=NaN;

		% 存储结果
		dist_in_swath(:,ii)=pd_in_seg(:)+bends(ii);  % 累积总距离
		dist_from_base(:,ii)=DistFromBaseLine(:);     % 当前段基线距离
	end

	% 确定最近段结果
	[~,c]=min(abs(dist_from_base),[],2,'omitnan');
	r=[1:numel(c)]; r=r(:);
	ix=sub2ind(size(dist_from_base),r,c);
	ds=dist_in_swath(ix);  % 沿扫掠线最终距离
	db=dist_from_base(ix);  % 最终基线距离

	% 过滤超出扫掠范围的点
	idx=single(ds)>=max(swdist);
	db(idx)=NaN;
	ds(idx)=NaN;
end

function [n_x,n_y]=RotCoord(x,y,theta,x0,y0)
	% 坐标旋转变换函数
	% 输入：
	%   x,y - 原始坐标
	%   theta - 旋转角度（弧度）
	%   x0,y0 - 旋转中心
	% 输出：
	%   n_x,n_y - 旋转后的新坐标
	n_x=(x-x0).*cos(theta)-(y-y0).*sin(theta);
	n_y=(x-x0).*sin(theta)+(y-y0).*cos(theta);
end