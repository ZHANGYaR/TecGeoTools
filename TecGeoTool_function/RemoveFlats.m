function [DEMn,MASK]=RemoveFlats(dem,strength)
	%
	% 用法:
	%	[DEMn,MASK]=RemoveFlats(dem,strength);
	%
	% 描述:
	% 	本函数采用半自动化流程处理数字高程模型（DEM），识别并移除平坦区域。
	% 	处理效果可能不如GIS软件手动编辑精确，但执行速度显著更快。
	%
	% 必需输入参数:
	% 	dem - 支持以下格式：
	%           - DEM文件完整路径（ASCII文本或GeoTIFF格式）
	%           - 工作空间中的GRIDobj对象名称
	%	strength - 平坦区域识别强度参数（整型1-4）
	%           强度与邻域大小对应关系：
	%           1:3x3  | 2:5x5 | 3:7x7 | 4:9x9
	%           若结果欠佳可增大强度值，若包含过多非平坦区域可减小强度值
	%
	% 输出:
	%	DEMn - 平坦区域被置为NaN的新DEM
	%	MASK - 逻辑型GRIDobj，标记识别出的平坦区域（true表示平坦区）
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% 作者：Adam M. Forte - 最后更新：2018/06/18 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% 解析输入参数
	p = inputParser;
	p.FunctionName = 'RemoveFlats';
	addRequired(p,'dem',@(x) isa(x,'GRIDobj') | ischar(x));
	addRequired(p,'strength',@(x) isnumeric(x) && x<=4 && mod(x,1)==0);

	parse(p,dem,strength);
	dem=p.Results.dem;
	strength=p.Results.strength;	

	disp('正在准备网格数据...')

	% 定义邻域卷积核
	nhood1=ones(3,3);  % 基础邻域
	% 根据强度选择扩展邻域
	if strength==1
		nhood2=ones(3,3);
	elseif strength==2
		nhood2=ones(5,5);
	elseif strength==3
		nhood2=ones(7,7);
	elseif strength==4
		nhood2=ones(9,9);
	else
		if isdeployed
			errordlg('强度参数无效！请输入1-4之间的整数')
		end
		error('强度参数无效！请输入1-4之间的整数')
	end

	% 处理输入数据类型
	if isa(dem,'GRIDobj')
		DEM=dem;
	elseif ischar(dem)
		disp('正在加载DEM文件')
		DEM=GRIDobj(dem);
	else
		if isdeployed
			errordlg('DEM输入类型错误：应为GRIDobj对象或文件路径')
		end
		error('DEM输入类型错误：应为GRIDobj对象或文件路径')
	end

	% 填充地形凹陷
	[DEMs]=fillsinks(DEM);

	% 识别平坦区域
	FL=GRIDobj(DEMs);
	FL=erode(DEMs,nhood1)==DEMs;  % 通过腐蚀操作识别平坦区
	CFL=dilate(FL,nhood2);       % 扩展平坦区域形成连续区域

	% 标记连续平坦区
	L=bwlabel(CFL.Z);  % 连通区域标记
	[xdim,ydim]=getcoordinates(DEM);
	L=GRIDobj(xdim,ydim,L);

	% 生成地形坡度图用于可视化
	G=gradient8(DEM);  % 计算八方向坡度

	% 用户交互选择洼地区域
	f1=figure(1);
	hold on 
	title('DEM坡度图：请框选洼地区域后按回车键确认')
	imageschs(DEM,G,'colormap','jet','caxis',[0 1]);  % 可视化坡度
	if ~verLessThan('matlab','9.5')
        disableDefaultInteractivity(gca);  % 禁用新版MATLAB的默认交互
    end 	
	hold off
	[x,y]=ginput;  % 获取用户选择的坐标点
	close(f1)

	% 转换坐标点为索引
	[ix]=coord2ind(DEM,x,y);
	l=L.Z(ix);
	l=unique(l);  % 获取唯一标记值

	% 生成逻辑掩膜
	MASK=GRIDobj(DEM);
	for ii=1:numel(l)
		MASK.Z(L.Z==l(ii))=1;  % 标记选定区域
	end
	MASK.Z=logical(MASK.Z);

	% 生成处理后的DEM
	DEMn=DEM;
	DEMn.Z(MASK.Z)=nan;  % 掩膜区置为NaN

	% 可视化验证结果
	f1=figure(1);
	hold on
	imageschs(DEMn,MASK);  % 显示处理后的DEM与掩膜
	if ~verLessThan('matlab','9.5')
        disableDefaultInteractivity(gca);
    end 	
	hold off
end