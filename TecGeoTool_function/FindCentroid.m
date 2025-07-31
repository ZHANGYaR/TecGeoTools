function [Cx,Cy]=FindCentroid(DEM)
	%
	% 用法:
	%	[Cx,Cy]=FindCentroid(DEM);
	%
	% 描述:
	% 	该函数用于计算流域的质心坐标
	%
	% 必需输入: 
	%	DEM - GRIDobj格式的数字高程模型
	%
	% 输出:
	%	[Cx,Cy] - 质心坐标（X,Y），坐标系与输入DEM一致
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% 函数作者：Adam M. Forte - 最后更新：2018/06/18      %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% 创建非NaN区域的逻辑矩阵
	I=isnan(DEM);
	I=~I;
	
	% 转换逻辑矩阵为双精度矩阵
	Imat=GRIDobj2mat(I);
	Imat=double(Imat);

	% 获取DEM的坐标网格
	[X,Y]=getcoordinates(DEM);

	% 计算加权坐标
	Ix=bsxfun(@times,Imat,X); % X方向加权
	Iy=bsxfun(@times,Imat,Y); % Y方向加权

	% 提取非零值
	Ixs=nonzeros(Ix); % 有效X坐标集合
	Iys=nonzeros(Iy); % 有效Y坐标集合

	% 计算质心坐标均值
	Cxi=mean(Ixs); % X方向质心
	Cyi=mean(Iys); % Y方向质心

	% 寻找最近网格坐标
	x_dist=abs(X-Cxi); % X方向距离差
	y_dist=abs(Y-Cyi); % Y方向距离差

	[~,Ixx]=min(x_dist); % 最近X索引
	[~,Iyy]=min(y_dist); % 最近Y索引

	% 获取最终质心坐标
	Cx=X(Ixx); % 质心X坐标
	Cy=Y(Iyy); % 质心Y坐标

% 函数结束
end