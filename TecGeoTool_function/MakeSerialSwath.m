function [SWcell,points]=MakeSerialSwath(DEM,points,divisions,sw_length,varargin)
% 用法：
% [SWcell,points]=MakeSerialSwath(DEM,points,divisions,sw_length);
% [SWcell,points]=MakeSerialSwath(DEM,points,divisions,sw_length,'name',value);
%
% 描述：
%    创建一系列垂直于给定线的地形剖面
%
% 必要输入：
%   DEM - 数字高程模型GRID对象，用于生成地形剖面
%   points - n×2矩阵，包含定义主测线路径的x,y坐标点。至少需要两个端点坐标。
%       坐标需与DEM坐标系一致且位于DEM范围内（不能位于边缘）。空数组将触发交互式绘图选点。
%   divisions - 剖面划分参数，解释方式取决于'div_type'参数：
%       'number'模式：生成的剖面总数
%       'width'模式：单个剖面的宽度（地图单位）
%   sw_length - 单个剖面的垂直延伸长度（地图单位）
%
% 可选输入：
%   div_type ['number'] - 划分类型：'number'（数量）/'width'（宽度）
%   alignment ['center'] - 剖面定位方式：
%       'center'：剖面中心线与主测线重合
%       'right'：剖面位于测线行进方向右侧
%       'left'：剖面位于测线行进方向左侧
%       'between'：在两条测线之间生成剖面（需配合points2参数）
%   points2 [] - 第二条边界线坐标点（仅用于'between'模式）
%   add_grids [] - 附加数据网格单元数组（包含GRID对象及其标签）
%   sample [DEM像元大小] - 沿剖面的重采样距离（地图单位）
%   smooth [0] - 平滑滤波窗口宽度（地图单位）
%   plot_map [true] - 是否显示剖面位置示意图
%   plot_individual [false] - 是否单独显示每个剖面图
%
% 输出：
%   SWcell - 包含SWATHobj对象的单元数组，每行对应一个剖面
%   points - 使用的测线坐标点（交互式选点时返回用户绘制的点）
%
% 示例：
%   [SWcell,points]=MakeSerialSwath(DEM,points,100,1000);
%   [SWcell,points]=MakeSerialSwath(DEM,points,1000,1000,'div_type','width');
%   [SWcell,points]=MakeSerialSwath(DEM,points,100,1000,'alignment','center');	
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 作者：Adam M. Forte - 最后更新：2019年4月2日 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	% 解析输入参数
	p = inputParser;
	p.FunctionName = 'MakeSerialSwath';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'points',@(x) isempty(x) || isnumeric(x) && size(x,1)>=2 && size(x,2)==2);
	addRequired(p,'divisions',@(x) isscalar(x) && isnumeric(x));
	addRequired(p,'sw_length',@(x) isscalar(x) && isnumeric(x));

	addParameter(p,'points2',[],@(x) isnumeric(x) && size(x,1)>=2 && size(x,2)==2 || isempty(x));
	addParameter(p,'div_type','number',@(x) ischar(validatestring(x,{'number','width'})));
	addParameter(p,'alignment','center',@(x) ischar(validatestring(x,{'center','right','left','between'})));
	addParameter(p,'add_grids',[],@(x) isa(x,'cell') && size(x,2)==2 || isempty(x));
	addParameter(p,'sample',[],@(x) isscalar(x) && isnumeric(x) || isempty(x));
	addParameter(p,'smooth',0,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'plot_map',true,@(x) isscalar(x) && islogical(x));
	addParameter(p,'plot_individual',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'make_shape',true,@(x) islogical(x)); % 隐藏参数，用于编译版本处理shapefile生成问题
	addParameter(p,'out_dir',[],@(x) isdir(x));

	parse(p,DEM,points,divisions,sw_length,varargin{:});
	DEM=p.Results.DEM;
	points=p.Results.points;
	divisions=p.Results.divisions;
	sw_length=p.Results.sw_length;

	points2=p.Results.points2;
	div_type=p.Results.div_type;
	alignment=p.Results.alignment;
	AG=p.Results.add_grids;
	sample=p.Results.sample;
	smth=p.Results.smooth;
	plot_map=p.Results.plot_map;
	plot_individual=p.Results.plot_individual;
	make_shape=p.Results.make_shape;
	out_dir=p.Results.out_dir;

	if isempty(sample)
		sample=DEM.cellsize;
	end

	if isempty(out_dir)
		out_dir=pwd;
	end	

	% 检查between模式参数完整性
	if ~isempty(points) && strcmp(alignment,'between') && isempty(points2)
		if isdeployed
			errordlg('当选择"between"模式时必须提供第二条边界线"points2"')
		end
		error('当选择"between"模式时必须提供第二条边界线"points2"')
	end

	% 交互式获取测线坐标点
	if isempty(points) && ~strcmp(alignment,'between')
		f1=figure(1);
		set(f1,'Units','normalized','Position',[0.05 0.1 0.9 0.9]);
		clf
        imagesc(DEM)
    	if ~verLessThan('matlab','9.5')
        disableDefaultInteractivity(gca);
	    end 
        title('绘制连续剖面基准线（双击结束绘制）')
        [points] = getline;
        close(f1);
    elseif isempty(points) && strcmp(alignment,'between')
		f1=figure(1);
		set(f1,'Units','normalized','Position',[0.05 0.1 0.9 0.9]);
		clf
		hold on
        imagesc(DEM)
		if ~verLessThan('matlab','9.5')
	        disableDefaultInteractivity(gca);
	    end         
        title('绘制第一条边界线（双击结束）')
        [points1] = getline;
        plot(points1(:,1),points1(:,2),'-r');
        title('绘制第二条边界线（双击结束）')
        [points2] = getline;
        hold off
        close(f1);

        points={points1,points2};
    elseif ~isempty(points) && ~isempty(points2) && strcmp(alignment,'between')
    	points1=points;
    	points={points1,points2};
    end

	% 获取DEM范围
	[demx,demy]=getoutline(DEM,true);	

    % 处理不同对齐模式
    if strcmp(alignment,'between')
    	num_bounds=2;

    	% 检查边界线是否在DEM范围内
    	in_bnd1=inpolygon(points1(:,1),points1(:,2),demx,demy);
    	in_bnd2=inpolygon(points2(:,1),points2(:,2),demx,demy);

    	if nnz(in_bnd1)~=numel(in_bnd1)
    		if isdeployed
    			errordlg('第一条边界线超出DEM范围')
    		end
    		error('第一条边界线超出DEM范围');
    	elseif nnz(in_bnd2)~=numel(in_bnd2)
    		if isdeployed
    			errordlg('第二条边界线超出DEM范围')
    		end
    		error('第二条边界线超出DEM范围');
    	end

    else
    	num_bounds=1;
    end

    switch num_bounds
    case 1 % 单线模式（左/右/中对齐）

		% 计算测线总长度
		xd=diff(points(:,1)); yd=diff(points(:,2));
		dst=sqrt((xd.^2)+(yd.^2));
		dst=cumsum(dst);
		dst=vertcat(0,dst);
		tot_dst=dst(end);

		% 确定划分方式
		if strcmp(div_type,'number')
			number=divisions;
			trim_flag=false;
			if tot_dst/number < DEM.cellsize*2
				if isdeployed
					warndlg('剖面宽度接近单个像元尺寸，请检查是否误将宽度参数作为数量输入')
				end
				warning('剖面宽度接近单个像元尺寸，请检查是否误将宽度参数作为数量输入')
			end
		elseif strcmp(div_type,'width')
			rmndr=rem(tot_dst,divisions);
			number=(tot_dst-rmndr)/divisions;
			trim_flag=(rmndr~=0);
			if trim_flag
				tot_dst=tot_dst-rmndr;
			end
		end

		% 计算剖面边缘和中心点
		edg=linspace(0,tot_dst,number+1);
		wdth=unique(round(diff(edg)));
		mp=edg(2:end)-(wdth/2);

		% 坐标插值
		[nx,ny,nd]=interpline(points(:,1),points(:,2),dst,1);

		if trim_flag
			nx=nx(nd<=tot_dst);
			ny=ny(nd<=tot_dst);
			nd=nd(nd<=tot_dst);
		end

		% 定位剖面端点
		edg_x=zeros(numel(edg),1);
		edg_y=zeros(numel(edg),1);
		for ii=1:numel(edg)
			[~,ix]=min(abs(edg(ii)-nd));
			edg_x(ii)=nx(ix);
			edg_y(ii)=ny(ix);
		end

		mp_x=zeros(numel(mp),1);
		mp_y=zeros(numel(mp),1);
		for ii=1:numel(mp)
			[~,ix]=min(abs(mp(ii)-nd));
			mp_x(ii)=nx(ix);
			mp_y(ii)=ny(ix);
		end

		% 计算剖面起点终点
		dx=diff(edg_x); dy=diff(edg_y);
		[ang_alng,~]=cart2pol(dx,dy);

		switch alignment
		case 'center'
			[px,py]=pol2cart(ang_alng+(pi/2),sw_length/2);
			[mx,my]=pol2cart(ang_alng-(pi/2),sw_length/2);

			start_x=mp_x+px; start_y=mp_y+py;
			stop_x=mp_x+mx; stop_y=mp_y+my;
		case 'left'
			[px,py]=pol2cart(ang_alng+(pi/2),sw_length);
			start_x=mp_x+px; start_y=mp_y+py;
			stop_x=mp_x; stop_y=mp_y;
		case 'right'
			[mx,my]=pol2cart(ang_alng-(pi/2),sw_length);
			start_x=mp_x; start_y=mp_y;
			stop_x=mp_x+mx; stop_y=mp_y+my;
		end

		% 检查剖面端点是否在DEM内
		start_in=inpolygon(start_x,start_y,demx,demy);
		stop_in=inpolygon(stop_x,stop_y,demx,demy);

		if nnz(start_in)~=numel(start_in)
			if isdeployed
				errordlg('部分剖面起点超出DEM范围，请调整长度或对齐方式')
			end
			error('部分剖面起点超出DEM范围，请调整长度或对齐方式');
		elseif nnz(stop_in)~=numel(stop_in)
			if isdeployed
				errordlg('部分剖面终点超出DEM范围，请调整长度或对齐方式')
			end
			error('部分剖面终点超出DEM范围，请调整长度或对齐方式');
		end

	case 2 % 双线间模式

		% 计算双线长度
		xd1=diff(points1(:,1)); yd1=diff(points1(:,2));
		dst1=sqrt((xd1.^2)+(yd1.^2));
		dst1=cumsum(dst1);
		dst1=vertcat(0,dst1);
		tot_dst1=dst1(end);

		xd2=diff(points2(:,1)); yd2=diff(points2(:,2));
		dst2=sqrt((xd2.^2)+(yd2.^2));
		dst2=cumsum(dst2);
		dst2=vertcat(0,dst2);
		tot_dst2=dst2(end);	

		% 确定划分方式
		[min_dst,min_ix]=min([tot_dst1 tot_dst2]);
		[max_dst]=max([tot_dst1 tot_dst2]);

		if strcmp(div_type,'number')
			number=divisions;
			if min_dst/number < DEM.cellsize*2
				if isdeployed
					warndlg('剖面宽度接近单个像元尺寸，请检查参数设置')
				end
				warning('剖面宽度接近单个像元尺寸，请检查参数设置')
			end

			% 生成划分点
			edg_min=linspace(0,min_dst,number+1);
			wdth=unique(round(diff(edg_min)));
			mp_min=edg_min(2:end)-(wdth/2);

			edg_max=linspace(0,max_dst,number+1);
			wdth_max=unique(round(diff(edg_max)));
			mp_max=edg_max(2:end)-(wdth_max/2);	

			% 坐标插值
			[nx1,ny1,nd1]=interpline(points1(:,1),points1(:,2),dst1,1);
			[nx2,ny2,nd2]=interpline(points2(:,1),points2(:,2),dst2,1);

			% 确定起点终点
			if min_ix==1
				start_x=arrayfun(@(x) nx1(findclosest(nd1,x)), mp_min);
				stop_x=arrayfun(@(x) nx2(findclosest(nd2,x)), mp_max);
				start_y=arrayfun(@(x) ny1(findclosest(nd1,x)), mp_min);
				stop_y=arrayfun(@(x) ny2(findclosest(nd2,x)), mp_max);
			else
				start_x=arrayfun(@(x) nx1(findclosest(nd1,x)), mp_max);
				stop_x=arrayfun(@(x) nx2(findclosest(nd2,x)), mp_min);
				start_y=arrayfun(@(x) ny1(findclosest(nd1,x)), mp_max);
				stop_y=arrayfun(@(x) ny2(findclosest(nd2,x)), mp_min);
			end

		elseif strcmp(div_type,'width')
			rmndr=rem(min_dst,divisions);
			number=(min_dst-rmndr)/divisions;
			trim_flag=(rmndr~=0);
			if trim_flag
				min_dst=min_dst-rmndr;
			end

			% 生成划分点
			edg_min=linspace(0,min_dst,number+1);
			wdth=unique(round(diff(edg_min)));
			mp_min=edg_min(2:end)-(wdth/2);

			edg_max=linspace(0,max_dst,number+1);
			wdth_max=unique(round(diff(edg_max)));
			mp_max=edg_max(2:end)-(wdth_max/2);	

			% 坐标插值
			[nx1,ny1,nd1]=interpline(points1(:,1),points1(:,2),dst1,1);
			[nx2,ny2,nd2]=interpline(points2(:,1),points2(:,2),dst2,1);

			% 截断超长测线
			if min_ix==1 && trim_flag
				nx1=nx1(nd1<=min_dst);
				ny1=ny1(nd1<=min_dst);
				nd1=nd1(nd1<=min_dst);
			elseif trim_flag
				nx2=nx2(nd2<=min_dst);
				ny2=ny2(nd2<=min_dst);
				nd2=nd2(nd2<=min_dst);
			end

			% 确定起点终点
			if min_ix==1
				start_x=arrayfun(@(x) nx1(findclosest(nd1,x)), mp_min);
				stop_x=arrayfun(@(x) nx2(findclosest(nd2,x)), mp_max);
				start_y=arrayfun(@(x) ny1(findclosest(nd1,x)), mp_min);
				stop_y=arrayfun(@(x) ny2(findclosest(nd2,x)), mp_max);
			else
				start_x=arrayfun(@(x) nx1(findclosest(nd1,x)), mp_max);
				stop_x=arrayfun(@(x) nx2(findclosest(nd2,x)), mp_min);
				start_y=arrayfun(@(x) ny1(findclosest(nd1,x)), mp_max);
				stop_y=arrayfun(@(x) ny2(findclosest(nd2,x)), mp_min);
			end
		end
	end

	% 生成剖面对象
	if ~isempty(AG)
		num_grids=size(AG,1);
	else
		num_grids=0;
	end

	SWcell=cell(number,1+num_grids);
	PLOTcell=cell(number,1+num_grids);
	if ~verLessThan('matlab','9.4')
		VERTcell=cell(number,1);
	end

	w1=waitbar(0,'正在生成剖面...');
	for ii=1:number
		SWcell{ii,1}=SWATHobj(DEM,vertcat(start_x(ii),stop_x(ii)),vertcat(start_y(ii),stop_y(ii)),'width',wdth,'dx',sample,'smooth',smth);
		if ~isempty(AG)
			for jj=1:num_grids
				AGoi=AG{jj,1};
				SWcell{ii,jj+1}=SWATHobj(AGoi,vertcat(start_x(ii),stop_x(ii)),vertcat(start_y(ii),stop_y(ii)),'width',wdth);
			end
		end

		PLOTcell{ii,1}=[SWcell{ii,1}.distx min(SWcell{ii,1}.Z,[],'omitnan').' mean(SWcell{ii,1}.Z,'omitnan').' max(SWcell{ii,1}.Z,[],'omitnan').'];
		if ~isempty(AG)
			for jj=1:num_grids
				PLOTcell{ii,jj+1}=[SWcell{ii,jj+1}.distx min(SWcell{ii,jj+1}.Z,[],'omitnan').' mean(SWcell{ii,jj+1}.Z,'omitnan').' max(SWcell{ii,jj+1}.Z,[],'omitnan').'];
			end
		end

		if ~verLessThan('matlab','9.4')
			VERTcell{ii}=SwathPolygon(SWcell{ii,1},wdth);
		end
		waitbar(ii/number,w1);
	end
	close(w1);

	% 可视化输出
	if plot_map
		f1=figure(1);
		clf
		set(f1,'Units','normalized','Position',[0.05 0.1 0.6 0.9]);
		cmap=jet(number);

		switch num_bounds
		case 1
			sbplt1=subplot(3,1,1);
			hold on
			imagesc(DEM);
			for ii=1:number
				if ~verLessThan('matlab','9.4')
					plot(VERTcell{ii}(:,1),VERTcell{ii}(:,2),'Color',cmap(ii,:),'LineWidth',2);
				else
					plot(SWcell{ii,1}.xy0(:,1),SWcell{ii,1}.xy0(:,2),'Color',cmap(ii,:),'LineWidth',2);
				end
			end
			plot(points(:,1),points(:,2),'-k','LineWidth',2);
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(sbplt1);
		    end 
			hold off

			sbplt2=subplot(3,1,2);
			hold on 
			for ii=1:number
				plot(PLOTcell{ii,1}(:,1),PLOTcell{ii,1}(:,3),'LineWidth',2,'Color',cmap(ii,:));
				plot(PLOTcell{ii,1}(:,1),PLOTcell{ii,1}(:,2),'LineWidth',0.5,'Color',cmap(ii,:));
				plot(PLOTcell{ii,1}(:,1),PLOTcell{ii,1}(:,4),'LineWidth',0.5,'Color',cmap(ii,:));
			end
			title('所有高程剖面')
			xlabel('剖面沿线距离（米）');
			ylabel('高程（米）');
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(sbplt2);
		    end 			
			hold off

			% 控制沿x方向采样点数的轻微变化
			for ii=1:number
				num_dist(ii,1)=numel(PLOTcell{ii,1}(:,3));
			end
			mnd=min(num_dist);

			% 组装数组
			for ii=1:number
				mean_comp(:,ii)=PLOTcell{ii,1}(1:mnd,3);
				max_comp(:,ii)=PLOTcell{ii,1}(1:mnd,4);
				min_comp(:,ii)=PLOTcell{ii,1}(1:mnd,2);
			end

			sbplt3=subplot(3,1,3);
			hold on
			plt(1)=plot(PLOTcell{1,1}(1:mnd,1),mean(mean_comp,2),'LineWidth',2,'Color','k');
			plt(2)=plot(PLOTcell{1,1}(1:mnd,1),max(mean_comp,[],2),'LineWidth',0.5,'Color','k');
			plot(PLOTcell{1,1}(1:mnd,1),min(mean_comp,[],2),'LineWidth',0.5,'Color','k');
			plt(3)=plot(PLOTcell{1,1}(1:mnd,1),mean(mean_comp,2)+std(mean_comp,0,2),'--','LineWidth',0.5,'Color','k');
			plot(PLOTcell{1,1}(1:mnd,1),mean(mean_comp,2)-std(mean_comp,0,2),'--','LineWidth',0.5,'Color','k');
			plt(4)=plot(PLOTcell{1,1}(1:mnd,1),max(max_comp,[],2),':','LineWidth',0.5,'Color','k');
			plot(PLOTcell{1,1}(1:mnd,1),min(min_comp,[],2),':','LineWidth',0.5,'Color','k');
			title('沿剖面的平均地形')
			legend(plt,{'平均值','均值的极值','均值标准差','全体极值'},'location','best');
			xlabel('沿剖面距离（米）');
			ylabel('沿剖面高程（米）');
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(sbplt3);
		    end 			
			hold off

		case 2

			sbplt1=subplot(2,1,1);
			hold on
			imagesc(DEM);
			for ii=1:number
				if ~verLessThan('matlab','9.4');
					plot(VERTcell{ii}(:,1),VERTcell{ii}(:,2),'Color',cmap(ii,:),'LineWidth',2);
				else
					plot(SWcell{ii,1}.xy0(:,1),SWcell{ii,1}.xy0(:,2),'Color',cmap(ii,:),'LineWidth',2);
				end
			end
			plot(points1(:,1),points1(:,2),'-k','LineWidth',2);
			plot(points2(:,1),points2(:,2),'-k','LineWidth',2);
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(sbplt1);
		    end 			
			hold off

			sbplt2=subplot(2,1,2);
			hold on 
			for ii=1:number
				plot(PLOTcell{ii,1}(:,1),PLOTcell{ii,1}(:,3),'LineWidth',2,'Color',cmap(ii,:));
				plot(PLOTcell{ii,1}(:,1),PLOTcell{ii,1}(:,2),'LineWidth',0.5,'Color',cmap(ii,:));
				plot(PLOTcell{ii,1}(:,1),PLOTcell{ii,1}(:,4),'LineWidth',0.5,'Color',cmap(ii,:));
			end
			title('全部高程剖面')
			xlabel('沿剖面距离（米）');
			ylabel('沿剖面高程（米）');
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(sbplt2);
		    end 			
			hold off

		end

	end

	if plot_individual
		disp('正在生成独立图表...')
		for ii=1:number
			f=figure(ii+1);
			clf
			set(f,'Units','normalized','Position',[0.05 0.1 0.6 0.6]);

			subplot(2+num_grids,1,1)
			hold on
			imagesc(DEM);
			if ~verLessThan('matlab','9.4')
				plot(VERTcell{ii}(:,1),VERTcell{ii}(:,2),'Color','k','LineWidth',2);
			else
				plot(SWcell{ii,1}.xy0(:,1),SWcell{ii,1}.xy0(:,2),'Color','k','LineWidth',2);
			end

			switch num_bounds
			case 1
				plot(points(:,1),points(:,2),'-k','LineWidth',2);
			case 2
				plot(points1(:,1),points1(:,2),'-k','LineWidth',2);
				plot(points2(:,1),points2(:,2),'-k','LineWidth',2);
			end
			title(['第' num2str(ii) '号剖面']);
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(gca);
		    end 			
			hold off

			subplot(2+num_grids,1,2)
			hold on 
			plot(PLOTcell{ii,1}(:,1),PLOTcell{ii,1}(:,3),'LineWidth',2,'Color','k');
			plot(PLOTcell{ii,1}(:,1),PLOTcell{ii,1}(:,2),'LineWidth',0.5,'Color','k');
			plot(PLOTcell{ii,1}(:,1),PLOTcell{ii,1}(:,4),'LineWidth',0.5,'Color','k');
			xlabel('沿剖面距离（米）');
			ylabel('沿剖面高程（米）');
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(gca);
		    end 			
			hold off

			if ~isempty(AG)
				for jj=1:num_grids
					subplot(2+num_grids,1,2+jj)
					hold on 
					plot(PLOTcell{ii,1+jj}(:,1),PLOTcell{ii,1+jj}(:,3),'LineWidth',2,'Color','k');
					plot(PLOTcell{ii,1+jj}(:,1),PLOTcell{ii,1+jj}(:,2),'LineWidth',0.5,'Color','k');
					plot(PLOTcell{ii,1+jj}(:,1),PLOTcell{ii,1+jj}(:,4),'LineWidth',0.5,'Color','k');
					xlabel('沿剖面距离（米）');
					ylabel(AG{jj,2});
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(gca);
				    end 
					hold off
				end
			end
		end	
	end


	if make_shape
		% 创建输出形状文件
		ms=struct;
		switch num_bounds
		case 1
			ms(1,1).Geometry='Line';
			ms(1,1).X=points(:,1);
			ms(1,1).Y=points(:,2);
			ms(1,1).Type='基准线';

			for ii=1:number
				if ~verLessThan('matlab','9.4')	
					ms(ii+1,1).Geometry='Line';
					ms(ii+1,1).X=VERTcell{ii}(:,1);
					ms(ii+1,1).Y=VERTcell{ii}(:,2);
					ms(ii+1,1).Type=['剖面' num2str(ii) '边界'];
				else
					ms(ii+1,1).Geometry='Line';
					ms(ii+1,1).X=SWcell{ii,1}.xy0(:,1);
					ms(ii+1,1).Y=SWcell{ii,1}.xy0(:,2);
					ms(ii+1,1).Type=['剖面' num2str(ii) '中线'];
				end
			end
		case 2
			ms(1,1).Geometry='Line';
			ms(1,1).X=points1(:,1);
			ms(1,1).Y=points1(:,2);
			ms(1,1).Type='边界线1';

			ms(2,1).Geometry='Line';
			ms(2,1).X=points2(:,1);
			ms(2,1).Y=points2(:,2);
			ms(2,1).Type='边界线2';
			for ii=1:number
				if ~verLessThan('matlab','9.4')	
					ms(ii+2,1).Geometry='Line';
					ms(ii+2,1).X=VERTcell{ii}(:,1);
					ms(ii+2,1).Y=VERTcell{ii}(:,2);
					ms(ii+2,1).Type=['剖面' num2str(ii) '边界'];
				else
					ms(ii+2,1).Geometry='Line';
					ms(ii+2,1).X=SWcell{ii,1}.xy0(:,1);
					ms(ii+2,1).Y=SWcell{ii,1}.xy0(:,2);
					ms(ii+2,1).Type=['剖面' num2str(ii) '中线'];
				end
			end
		end
		shapewrite(ms,fullfile(out_dir,'SerialSwathBounds.shp'));
	end

end

function [verts]=SwathPolygon(SW,w)

	cx=SW.xy(:,1);
	cy=SW.xy(:,2);

	dx=diff(cx);
	dy=diff(cy);

	w=w/2;

	[sw_angle,~]=cart2pol(dx,dy);
	sw_angle=vertcat(sw_angle,sw_angle(end));
	[px,py]=pol2cart(sw_angle+(pi/2),w);
	[mx,my]=pol2cart(sw_angle-(pi/2),w);

	swx=[cx+px cx+mx];
	swy=[cy+py cy+my];

	warning off
	P=polyshape(vertcat(swx(:,1),flipud(swx(:,2))),vertcat(swy(:,1),flipud(swy(:,2))));
	warning on
	P=rmholes(P);
	verts=P.Vertices;
	verts=vertcat(verts,verts(1,:));

end