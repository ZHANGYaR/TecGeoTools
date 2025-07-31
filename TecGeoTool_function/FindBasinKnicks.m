function [KnickTable]=FindBasinKnicks(Basin_Data_File,plot_result,varargin)
	%
% 用法:
%	[KnickTable]=FindBasinKnicks(Basin_Data_File,plot_result);
%	[KnickTable]=FindBasinKnicks(Basin_Data_File,plot_result,'name',value,...);
%
% 描述:
% 	该函数用于在ProcessRiverBasins生成的结果文件中手动选择河道裂点。用户通过在Chi-高程图
% 	上点击鼠标选取裂点，完成选择后按回车确认。已选裂点会显示为红色标记。自动检测裂点可参考
% 	TopoToolbox的'knickpointfinder'。启用分类功能时（'classify_knicks' = true），需为每个
% 	高亮裂点输入分类标识（建议统一使用数字或简短字符串，避免混合类型导致错误）。
%
% 必需输入:
% 	Basin_Data_File - ProcessRiverBasins输出文件的完整路径
% 	plot_result     - 是否可视化结果（true显示图形，false不显示）
%
% 可选输入:
%	classify_knicks [false] - 是否对裂点进行分类标注
% 	ref_concavity [0.5]    - Chi计算使用的参考凹度值
%	save_mat [true]        - 是否保存.mat结果文件（文件名格式：Knicks_NUM.mat）
%	shape_name []          - 输出shapefile的文件名（不含扩展名）
%
% 输出:
%	KnickTable - 裂点信息表，包含坐标、高程、距离和Chi值。启用分类时包含第六列分类标识
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 作者：Adam M. Forte - 最后更新：2018/07/02       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% 解析输入参数
	p = inputParser;
	p.FunctionName = 'FindBasinKnicks';
	addRequired(p,'Basin_Data_File',@(x) ischar(x));
	addRequired(p,'plot_result',@(x) islogical(x));

	addParameter(p,'classify_knicks',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'ref_concavity',0.5,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'shape_name',[],@(x) ischar(x));
	addParameter(p,'save_mat',true,@(x) isscalar(x) && islogical(x));
	addParameter(p,'out_dir',[],@(x) isdir(x));

	parse(p,Basin_Data_File,plot_result,varargin{:});
	Basin_Data_File=p.Results.Basin_Data_File;
	plot_result=p.Results.plot_result;

	classify_knicks=p.Results.classify_knicks;
	theta_ref=p.Results.ref_concavity;
	shape_name=p.Results.shape_name;
	save_mat=p.Results.save_mat;
	out_dir=p.Results.out_dir;

	if isempty(out_dir)
		out_dir=pwd; % 默认输出到当前目录
	end

	% 加载流域数据文件
	load(Basin_Data_File);
    
	% 简化河网密度
	if drainage_area>20
		S=modify(Sc,'streamorder','>1'); % 保留二级以上河道
		if isempty(S.x)
			S=Sc; % 若无二级河道则保留原始
		end
	else
		S=Sc; % 小流域保留全部河道
	end

	% 获取河道源头坐标
	ChXY=streampoi(S,'channelheads','xy');
	ChIX=coord2ind(DEMoc,ChXY(:,1),ChXY(:,2));
	NumHeads=numel(ChIX); % 源头数量

	% 创建逻辑种子栅格
	SEED=GRIDobj(DEMoc,'logical');

	% 初始化图形窗口
	f1=figure(1);
	set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');
	f2=figure(2);
	set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');

	% 遍历各河道进行裂点选择
	Channels=cell(NumHeads,1);
	Chi=cell(NumHeads,1);
	Knicks=cell(NumHeads,1);

	% 显示操作提示
	uiwait(msgbox('在 χ-高程图上使用鼠标点击选择裂点，完成选择后按回车确认'));

	for ii=1:NumHeads
		ChOI=ChIX(ii); % 当前河道源头索引
		ChR=SEED;
		ChR.Z(ChOI)=true; % 标记当前河道起点

		% 提取下游河道网络
		SS=modify(S,'downstreamto',ChR);
		Channels{ii}=SS;

		% 计算Chi参数
		C=chiplot(SS,DEMcc,Ac,'a0',1,'mn',theta_ref,'plot',false);
		Chi{ii}=C;

		% 提取坐标参数
		chi=C.chi;
		elev=C.elev;
		x=C.x;
		y=C.y;
		d=C.distance;

		% 绘制当前河道网络
		figure(f2);
		clf
		hold on
		imageschs(DEMoc,DEMoc,'colormap',[0.5 0.5 0.5]); % 灰度地形图
		plot(S,'-w'); % 全河网白色显示
		plot(SS,'-r'); % 当前河道红色高亮
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(gca); % 禁用新版MATLAB默认交互
        end 
		hold off

		% 绘制Chi-高程图
		figure(f1);
		clf
		hold on
		scatter(C.chi,C.elev,10,'k','filled'); % 黑色散点

        % 显示已选裂点（红色标记）
        KC=vertcat(Knicks{1:ii-1});        
		if ii>1 && ~isempty(KC)
			xx=KC(:,1); yy=KC(:,2); cc=KC(:,5); ee=KC(:,3);
			currentRiv=[x y];
			pastKnicks=[xx yy];
			cidx=ismember(pastKnicks,currentRiv,'rows'); % 查找重复裂点
			cc=cc(cidx); ee=ee(cidx);
			scatter(cc,ee,30,'r','filled'); % 红色标记已选点
		end

		xlabel('χ值');
		ylabel('高程（米）');
		title(['剩余需处理河道数：' num2str(NumHeads-ii)]);
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(gca);
        end 
		hold off

		% 获取用户输入
		[c,e]=ginput; % 图形坐标采集

		if numel(c)>=1
			% 匹配最近河道点
			for jj=1:numel(c)
				coi=c(jj);
				[~,idx]=min(abs(chi-coi)); % 寻找最近χ值索引
				knp(jj,:)=[x(idx) y(idx) elev(idx) d(idx) chi(idx)]; % 存储裂点信息
			end
			Knicks{ii}=knp;
		end

		% 裂点分类处理
		if classify_knicks && numel(c)>=1
			figure(f1);
			for jj=1:numel(c)
				hold on
				s1=scatter(knp(jj,5),knp(jj,3),40,'g'); % 临时绿色标记
				cl=inputdlg('请输入裂点分类标识','分类输入'); % 弹出分类对话框
				hold off
				delete(s1); % 删除临时标记

				% 分类数据类型处理
				cl_num=str2num(cl{1});
				if isempty(cl_num)
					knpClass{jj,1}=cl{1}; % 字符型分类
				else 
					knpClass{jj,1}=double(cl_num); % 数值型分类
				end
			end	
			knpClasses{ii}=knpClass;
		end
	end

	% 整合裂点数据
	KnickPoints=vertcat(Knicks{:});
	[KnickPoints,idx,~]=unique(KnickPoints,'rows'); % 去重

	% 创建输出表格
	KnickTable=array2table(KnickPoints,'VariableNames',...
		{'x坐标','y坐标','高程','距离','χ值'});

	% 添加分类列
	if classify_knicks && exist('knpClasses','var')
		knpClasses=vertcat(knpClasses{:});
		knpClasses=knpClasses(idx);
		KnickTable. classification=knpClasses;
	end

	close(f1);
	close(f2);

	% 结果可视化
	if plot_result
		f1=figure(1);
		set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.8],'renderer','painters');
       
		% 空间分布子图
		sbplt1=subplot(1,2,1);
		hold on
		[RGB]=imageschs(DEMoc,DEMoc,'colormap','gray');
		[~,R]=GRIDobj2im(DEMoc);
		imshow(flipud(RGB),R);
		axis xy
		colormap(jet);
        caxis([0 max(KnickPoints(:,3))]); % 高程色标
		plot(S,'-w');
		scatter(KnickPoints(:,1),KnickPoints(:,2),20,KnickPoints(:,3),'filled');
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(sbplt1);
        end 
		
		% 剖面展示子图
		sbplt2=subplot(1,2,2);
		hold on
        plotdz(S,DEMoc,'Color',[0.5 0.5 0.5]); % 原始DEM剖面
		plotdz(S,DEMcc,'color','k'); % 填洼后剖面
        caxis([0 max(KnickPoints(:,3))]);
		scatter(KnickPoints(:,4),KnickPoints(:,3),20,KnickPoints(:,3),'filled');
		c1=colorbar;
		ylabel(c1,'裂点高程（米）');
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(sbplt2);
        end 		
	end

	% 结果输出
	if save_mat
		save(fullfile(out_dir,['Knicks_' num2str(RiverMouth(:,3)) '.mat']),...
			'KnickTable','-v7.3');
	end

	% 生成Shapefile
	if ~isempty(shape_name)
		MS=struct;
		for ii=1:size(KnickPoints,1)
			MS(ii,1).ID=ii;
			MS(ii,1).Geometry='Point';
			MS(ii,1).X=KnickPoints(ii,1);
			MS(ii,1).Y=KnickPoints(ii,2);
			MS(ii,1).Elev=KnickPoints(ii,3);
			MS(ii,1).Distance=KnickPoints(ii,4);
			MS(ii,1).Chi=KnickPoints(ii,5);
			if classify_knicks
				MS(ii,1).Class=knpClasses{ii};
			end
		end
		shp_out=fullfile(out_dir,[shape_name '.shp']);
		shapewrite(MS,shp_out); % 写入shapefile
	end
end