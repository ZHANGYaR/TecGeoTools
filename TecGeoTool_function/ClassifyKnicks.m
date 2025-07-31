function ClassifyKnicks(DEM,FD,A,Sc,ksn_master,bnd_list,kn_list,varargin)
	%
% 用法：
%	ClassifyKnicks(DEM, FD, A, Sc, ksn_master, bnd_list, kn_list);
%	ClassifyKnicks(DEM, FD, A, Sc, ksn_master, bnd_list, kn_list, 'shape_name', '名称');
%
% 描述：
% 	本函数用于迭代处理通过'KsnProfiler'生成的边界点与裂点数据集，
% 	展示单个河段的长剖面及χ-高程图，支持对KsnProfiler输出的每个边界/裂点进行分类标注。
% 	用户需通过弹窗输入数字或字符进行分类（单次运行需保持输入类型一致），
% 	建议字符分类使用无空格短字符串（如'裂点'或'边界'），便于Shapefile属性表存储。
%
% 必需输入：
%	DEM - 数字高程模型（GRIDobj格式）
%	FD - 流向对象（FLOWobj格式）
%	A - 流量累积对象（GRIDobj格式）
%	Sc - 河网对象（STREAMobj格式）
%	ksn_master - 河道陡度分析主数据（单元格数组）
%	bnd_list - 边界点矩阵（KsnProfiler输出）
%	kn_list - 裂点矩阵（KsnProfiler输出）
%
% 可选输入：
%	shape_name ['ksn'] - 输出矢量文件名（不含.shp后缀，需符合ArcGIS命名规范）
%
% 输出：
%	生成包含分类信息的裂点矢量文件（Shapefile格式）
%
% 注意：
%	即使bnd_list或kn_list包含NaN值，也需完整提供KsnProfiler的输出结果
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数作者：Adam M. Forte - 更新日期：2020年1月10日 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% 解析输入参数
	p = inputParser;
	p.FunctionName = 'ClassifyKnicks';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
	addRequired(p,'A',@(x) isa(x,'GRIDobj'));
	addRequired(p,'Sc',@(x) isa(x,'STREAMobj'));
	addRequired(p,'ksn_master',@(x) iscell(x));
	addRequired(p,'bnd_list',@(x) ismatrix(x));
	addRequired(p,'kn_list',@(x) ismatrix(x));

	addParameter(p,'shape_name','ksn',@(x) ischar(x));

	parse(p,DEM,FD,A,Sc,ksn_master,bnd_list,kn_list,varargin{:});
	DEM=p.Results.DEM;
	FD=p.Results.FD;
	Sc=p.Results.Sc;
	A=p.Results.A;
	ksn_master=p.Results.ksn_master;
	bnd_list=p.Results.bnd_list;
	kn_list=p.Results.kn_list;

	shape_name=p.Results.shape_name;

	% 计算出口点上游流向距离
	OUTIX=GRIDobj(DEM,'logical');
	outix=streampoi(Sc,'outlets','ix');
	OUTIX.Z(outix)=true;
	FLDS=flowdistance(FD,OUTIX); 

	% 过滤无效边界点
	idx=isnan(bnd_list(:,1));
	bnd_list=bnd_list(~idx,:);
	bnd_x=bnd_list(:,1);
	bnd_y=bnd_list(:,2);
	bnd_ix=coord2ind(DEM,bnd_x,bnd_y);
	bnd_el=DEM.Z(bnd_ix);

	% 过滤无效裂点
	idx=isnan(kn_list(:,1));
	kn_list=kn_list(~idx,:);
	kn_x=kn_list(:,1);
	kn_y=kn_list(:,2);
	kn_ix=coord2ind(DEM,kn_x,kn_y);
	kn_el=DEM.Z(kn_ix);	
 
	% 识别包含裂点的河道
	str_list=unique(vertcat(bnd_list(:,4),kn_list(:,4)));
	num_streams=numel(str_list);

	% 初始化分类标记数组
	b_is_classified=false(numel(bnd_ix),1);
	k_is_classified=false(numel(kn_ix),1);

	for ii=1:num_streams
		% 重建完整河道网络
		ix=coord2ind(DEM,ksn_master{str_list(ii),1}(:,1),ksn_master{str_list(ii),1}(:,2));
		ref_theta=ksn_master{str_list(ii),1}(1,7);
		IX=influencemap(FD,ix);
		ST=STREAMobj(FD,IX);
		ST=intersect(ST,Sc); % 与原始河网校正
		z=mincosthydrocon(ST,DEM,'interp',0.1); % 高程插值

		% 计算χ参数
		c=chitransform(ST,A,'a0',1,'mn',ref_theta);
		C=GRIDobj(DEM);
		C.Z(ST.IXgrid)=c;

		% 定位当前河道的特征点
		bidx=ismember(bnd_ix,ST.IXgrid);
		kidx=ismember(kn_ix,ST.IXgrid);
		bidx2= bidx & ~b_is_classified;
		kidx2= kidx & ~k_is_classified;

		% 更新分类标记
		b_is_classified(bidx2)=true;
		k_is_classified(kidx2)=true;

		% 提取待分类点
		bnd_ixOI=bnd_ix(bidx2);
		bnd_xOI=bnd_x(bidx2);
		bnd_yOI=bnd_y(bidx2);
		bnd_class=zeros(numel(bnd_ixOI),1);

		kn_ixOI=kn_ix(kidx2);
		kn_xOI=kn_x(kidx2);
		kn_yOI=kn_y(kidx2);
		kn_class=ones(numel(kn_ixOI),1);

		% 合并边界点与裂点
		cm_ixOI=vertcat(bnd_ixOI,kn_ixOI);
		cm_xOI=vertcat(bnd_xOI,kn_xOI);
		cm_yOI=vertcat(bnd_yOI,kn_yOI);
		cm_class=vertcat(bnd_class,kn_class);

		if ~isempty(cm_ixOI)
			% 创建可视化窗口
			f1=figure(1);
			set(f1,'Units','normalized','Position',[0.5 0.1 0.5 0.5],'renderer','painters');
			clf
			
			% 绘制长剖面图
			sbplt1=subplot(2,1,1);
			hold on
			plotdz(ST,DEM,'Color',[0.5 0.5 0.5]); % 原始高程
			plotdz(ST,z,'Color','k'); % 插值后高程
			cnt=0;
			
			% 添加图例元素
			if ~isempty(bnd_ixOI)
				cnt=cnt+1;
				s1(cnt)=scatter(FLDS.Z(bnd_ixOI),DEM.Z(bnd_ixOI),30,'k','filled');
				leg{cnt,1}='河道边界';
			end
			if ~isempty(kn_ixOI)
				cnt=cnt+1;
				s1(cnt)=scatter(FLDS.Z(kn_ixOI),DEM.Z(kn_ixOI),30,'k','s','LineWidth',2);
				leg{cnt,1}='裂点位置';
			end
			if ~verLessThan('matlab','9.5')
				disableDefaultInteractivity(sbplt1);
			end	
			hold off

			% 绘制χ-高程图
			sbplt2=subplot(2,1,2);
			hold on
			cvec=getnal(ST,C);
			[cvec,six]=sort(cvec);
			evec=z(six);
			plot(cvec,evec,'-k','LineWidth',2);
			if ~isempty(bnd_ixOI)
				scatter(C.Z(bnd_ixOI),DEM.Z(bnd_ixOI),30,'k','filled');
			end
			if ~isempty(kn_ixOI)
				scatter(C.Z(kn_ixOI),DEM.Z(kn_ixOI),30,'k','s','LineWidth',2);
			end
			if ~verLessThan('matlab','9.5')
				disableDefaultInteractivity(sbplt2);
			end				
			hold off

			% 交互式分类处理
			for jj=1:numel(cm_ixOI)
				subplot(2,1,1)
				hold on
				if cm_class(jj)==0
					p1=scatter(FLDS.Z(cm_ixOI(jj)),DEM.Z(cm_ixOI(jj)),50,'r','filled');
				else
					p1=scatter(FLDS.Z(cm_ixOI(jj)),DEM.Z(cm_ixOI(jj)),50,'r','filled','s');
				end
				% 动态更新图例
				if numel(unique(leg))==1
					legend(s1(1),leg{1},'location','best');
				else
					legend(s1,leg,'location','best');
				end
				hold off

				subplot(2,1,2)
				hold on
				if cm_class(jj)==0
					p2=scatter(C.Z(cm_ixOI(jj)),DEM.Z(cm_ixOI(jj)),50,'r','filled');
				else
					p2=scatter(C.Z(cm_ixOI(jj)),DEM.Z(cm_ixOI(jj)),50,'r','filled','s');
				end
				xlabel('\chi');
				ylabel('高程（米）');
				hold off

				% 获取用户输入
				c=inputdlg('请输入选定特征的分类：','特征分类');
				cn=str2num(c{1});
				char_flag=isempty(cn);

				% 存储分类结果
				new_bnds(jj,:)=[cm_xOI(jj) cm_yOI(jj) FLDS.Z(cm_ixOI(jj)) DEM.Z(cm_ixOI(jj)) str_list(ii)];
				if char_flag
					bnd_cat{jj,1}=c;
				else
					bnd_cat(jj,1)=cn;
				end
				delete(p1);
				delete(p2);
			end

			new_bnds_c{ii,1}=new_bnds;
			bnd_cats{ii,1}=bnd_cat;
		else
			new_bnds_c{ii,1}=NaN;
			bnd_cats{ii,1}=NaN;
		end
	end

	% 整合输出结果
	new_bnds=vertcat(new_bnds_c{:});
	[new_bnds,bix,~]=unique(new_bnds,'rows');
	bnd_cats=vertcat(bnd_cats{:});
	bnd_cats=bnd_cats(bix);
	if exist('char_flag','var') && char_flag
		bnd_cats=vertcat(bnd_cats{:});
	end

	% 构建Shapefile结构
	KNK=struct;
	for jj=1:numel(new_bnds(:,1))
		KNK(jj,1).Geometry='Point';
		KNK(jj,1).X=double(new_bnds(jj,1));
		KNK(jj,1).Y=double(new_bnds(jj,2));
		KNK(jj,1).Dist=double(new_bnds(jj,3));
		KNK(jj,1).Elev=double(new_bnds(jj,4));
		KNK(jj,1).StrNum=double(new_bnds(jj,5));
		if exist('char_flag','var') && char_flag
			KNK(jj,1).Category=bnd_cats{jj,1};
		else 
			KNK(jj,1).Category=double(bnd_cats(jj,1));
		end
	end
	out_knick_name=[shape_name '_knicks_classified.shp'];
	shapewrite(KNK,out_knick_name);

	close(f1)
end