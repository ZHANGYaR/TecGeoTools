function [ds,db,mag,unc,nc0,ec0]=ProjectGPSOntoSwath(SW,x,y,data_width,nc,ec,nu,eu)
	% 用法：
%   [ds,db,mag,unc,nc0,ec0]=ProjectGPSOntoSwath(SW,x,y,data_width,nc,ec,nu,eu);
%
% 描述：
%   该函数是ProjectOntoSwath的定制版本，专门处理GPS速度数据，
%   可计算扫掠方向上的矢量模量（北向和东向分量）及其不确定性。
%
% 必需输入：
%   SW - 扫掠对象，用于投影GPS数据
%   x - GPS站点的x坐标（nx1数组）
%   y - GPS站点的y坐标（nx1数组）
%   data_width - 采样宽度（地图单位），从扫掠基线开始测量
%   nc - 速度北向分量
%   ec - 速度东向分量
%   nu - 北向分量的不确定性
%   eu - 东向分量的不确定性
%
% 输出：
%   ds - 站点沿扫掠线的投影距离
%   db - 站点到扫掠基线的垂直距离
%   mag - 沿扫掠线的矢量模量，符号表示方向：
%         正值表示与扫掠方向一致（扫掠距离增加方向）
%         负值表示与扫掠方向相反
%   unc - 沿扫掠线的不确定性（误差椭圆分量）
%   nc0 - 投影矢量的北向分量（适用于MATLAB箭图）
%   ec0 - 投影矢量的东向分量（适用于MATLAB箭图）
%
% 示例：
%   [ds,db,mag,unc,nc0,ec0]=ProjectGPSOntoSwath(SW,x,y,10000,nc,ec,nu,eu);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数作者：Adam M. Forte - 最后更新：2019/04/02 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% 从SWATHobj中提取参数
	xypoints=SW.xy;
	swdist=SW.distx;
	xy0=SW.xy0;

	% 寻找弯曲点索引
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

	% 提取弯曲点并计算分段信息
	swxB=xy0(:,1); swyB=xy0(:,2);
	num_segs=numel(swxB)-1;  % 计算扫掠段数量

	% 计算各扫掠段长度
	if num_segs>1
		kk=1;
		while kk<=num_segs
			% 计算段端点坐标差
			x0=swxB(kk,1); y0=swyB(kk,1);
			x1=swxB(kk+1,1); y1=swyB(kk+1,1);
			xx=x1-x0; yy=y1-y0;
			% 计算段长度
			dist_to_bend(kk,1)=sqrt((xx^2)+(yy^2));
			kk=kk+1;
		end
		dist_to_bend=vertcat(0,dist_to_bend);  % 添加起始点
		bends=cumsum(dist_to_bend);  % 累积距离
	else
		% 单段处理
		bends=[0;max(swdist)];
		dist_to_bend=bends;
	end

	% 获取扫掠线所有点坐标
	swx=xypoints(:,1); swy=xypoints(:,2);

	% 初始化输出矩阵
	dist_in_swath=zeros(numel(x),num_segs);
	dist_from_base=zeros(numel(x),num_segs);
	mag_on_swath=zeros(numel(x),num_segs);
	unc_on_swath=zeros(numel(x),num_segs);
	nc0_on_swath=zeros(numel(x),num_segs);
	ec0_on_swath=zeros(numel(x),num_segs);

	% 主处理循环：逐个扫掠段处理
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

		% 计算段方向角
		sw_angle=-1*atan((swy1-swy0)/(swx1-swx0));

		% 坐标旋转变换
		[n_swx,n_swy]=RotCoord(swx,swy,sw_angle,swx0,swy0);
		[n_x,n_y]=RotCoord(x,y,sw_angle,swx0,swy0);

		% 遍历所有GPS点
		num_points=numel(n_x);
		pd_in_seg=zeros(num_points,1);
		DistFromBaseLine=zeros(num_points,1);
		for jj=1:num_points
			xoi=n_x(jj);
			
			% 寻找最近投影点
			if abs(xoi)<=max(abs(n_swx))
				distances=n_swx-xoi;
				[~,I]=min(abs(distances)); 
				pd_in_seg(jj)=seg_dist(I);
				% 计算基线距离
				BaseLine=n_swy(I);
				DistFromBaseLine(jj)=abs(BaseLine-n_y(jj));
			else 
				pd_in_seg(jj)=NaN;
				DistFromBaseLine(jj)=NaN;
			end

			% 处理弯曲重叠区域
			if isnan(pd_in_seg(jj)) && n_dist(jj)<=data_width && ii~=num_segs
				pd_in_seg(jj)=seg_dist(bend_ix(ii+1));
				DistFromBaseLine(jj)=n_dist(jj);
			end 
		end

		% 过滤无效点
		idx=isnan(pd_in_seg) | pd_in_seg<=0;	
		pd_in_seg(idx)=NaN;
		DistFromBaseLine(idx)=NaN;

		% 计算GPS矢量参数
		[mag_seg,unc_seg,ec0_seg,nc0_seg]=AngVecSw(swx0,swy0,swx1,swy1,nc,ec,nu,eu);
		% 应用过滤
		mag_seg(idx)=NaN; unc_seg(idx)=NaN;
		ec0_seg(idx)=NaN; nc0_seg(idx)=NaN;
		% 存储结果
		mag_on_swath(:,ii)=mag_seg;
		unc_on_swath(:,ii)=unc_seg;
		ec0_on_swath(:,ii)=ec0_seg;
		nc0_on_swath(:,ii)=nc0_seg;		

		dist_in_swath(:,ii)=pd_in_seg(:)+bends(ii);  % 累积距离
		dist_from_base(:,ii)=DistFromBaseLine(:);
	end

	% 寻找最近段的结果
	[db,c]=min(dist_from_base,[],2,'omitnan');
	r=[1:numel(c)]; r=r(:);
	ix=sub2ind(size(dist_from_base),r,c);
	ds=dist_in_swath(ix);
	mag=mag_on_swath(ix);
	unc=unc_on_swath(ix);
	ec0=ec0_on_swath(ix);
	nc0=nc0_on_swath(ix);

	% 过滤超出扫掠范围的点
	idx=single(ds)>=max(swdist);
	db(idx)=NaN; ds(idx)=NaN;
	mag(idx)=NaN; unc(idx)=NaN;
	ec0(idx)=NaN; nc0(idx)=NaN;
end

function [n_x,n_y]=RotCoord(x,y,theta,x0,y0)
	% 坐标旋转变换
	n_x=(x-x0).*cos(theta)-(y-y0).*sin(theta);
	n_y=(x-x0).*sin(theta)+(y-y0).*cos(theta);
end

function [proj_mag,unc,ec0,nc0]=AngVecSw(swx0,swy0,swx1,swy1,nc,ec,nu,eu)
	% 将扫掠段定向到北象限
	if swx0>swx1 && swy0>=swy1
		OX=swx1; OY=swy1; HX=swx0; HY=swy0;
	elseif swx0<swx1 && swy0<=swy1
		OX=swx0; OY=swy0; HX=swx1; HY=swy1;
	elseif swx0>=swx1 && swy0<swy1
		OX=swx0; OY=swy0; HX=swx1; HY=swy1;
	elseif swx0<=swx1 && swy0>swy1
		OX=swx1; OY=swy1; HX=swx0; HY=swy0;
	end

	% 计算方向角
	[swT,~]=cart2pol(HX-OX,HY-OY);  % 扫掠段方向角
	[pT,~]=cart2pol(ec,nc);          % GPS矢量方向角

	% 计算投影模量
	theta=swT-pT;
	mag=hypot(ec,nc);
	proj_mag=mag.*cos(theta);  % 矢量投影

	% 分解投影分量
	[ec0,nc0]=pol2cart(swT,proj_mag);

	% 计算误差椭圆半径
	unc=ellipserad(eu,nu,swT);
end

function [r]=ellipserad(a,b,theta)
	% 误差椭圆半径计算
	% 输入：
	%   a - 长轴（东向不确定度）
	%   b - 短轴（北向不确定度） 
	%   theta - 扫掠方向角
	r=(a.*b)./sqrt((b.*cos(theta)).^2 + (a.*sin(theta)).^2);
end