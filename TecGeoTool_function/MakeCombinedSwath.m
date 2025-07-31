function [SW,SwathMat,xypoints,outData]=MakeCombinedSwath(DEM,points,width,data_type,data,data_width,varargin)
	% 用法：
%   [SW, SwathMat, xypoints, outData] = MakeCombinedSwath(DEM, points, width, data_type, data, data_width);
%   [SW, SwathMat, xypoints, outData] = MakeCombinedSwath(DEM, points, width, data_type, data, data_width, 'name', value, ...);
%
% 描述：
%   该函数用于在剖面图上绘制附加数据。
%
% 必要输入：
%   DEM - DEM网格对象，用于生成拓扑剖面
%   points - n x 2 矩阵，包含用于生成剖面的x,y坐标，至少需要两个点（起点和终点）。
%            第一行包含起点坐标，后续行是弯道点。坐标必须与DEM处于同一坐标系中，并且必须位于DEM范围内。
%            如果提供空数组，将显示DEM图像以选择点。
%   width - 剖面宽度（以地图单位为单位）
%   data_type - 提供的附加数据类型，支持的输入有：
%       'points3' - 通用点数据集，期望一个 n x 3 矩阵，包含x, y, z值
%       'points4' - 通用点数据集，期望一个 n x 4 矩阵，包含x, y, z和额外的值。点将根据额外值着色
%       'points5' - 通用点数据集，期望一个 n x 5 矩阵，包含x, y, z和两个额外值。点根据第一个额外值（第四列）着色，按第二个额外值（第五列）缩放
%       'eqs' - 地震数据，期望一个 n x 4 矩阵，包含x, y, 深度和震级。点将根据震级缩放，并根据距离剖面线的距离着色。深度应为正值。
%       'gps' - GPS速度向量，期望一个 n x 6 矩阵，包含x, y, 北向分量，东向分量，北向不确定性，东向不确定性。详见 'ProjectGPSOntoSwath' 。
%       'STREAMobj' - 将流域剖面（作为点）投影到剖面线上。期望一个基于提供的DEM生成的STREAMobj。
%       'ksn_chandata' - 将通过ksn值绘制剖面，期望旧的Profiler51代码生成的chandata文件。
%       'ksn_batch' - 将通过ksn值绘制剖面，期望 'KsnChiBatch' 函数的输出map结构（即运行 'KsnChiBatch' 时设置产品为 'ksn'，输出设置为true，第二个输出即为data）。
%       'ksn_profiler' - 将通过ksn值绘制剖面，期望 'KsnProfiler' 函数的'knl'输出。
%       'basin_stats' - 通过从 'ProcessRiverBasins' 计算的平均流域值绘制剖面，期望来自 'CompileBasinStats' 的输出，并需要在可选输入 'basin_value' 中指定要根据该值着色的点。
%       'basin_knicks' - 通过 'FindBasinKnicks' 选择的河流断点绘制剖面。提供包含河流断点数据的文件夹或路径。
%   data - 输入数据，根据选择的 data_type 其格式不同
%   data_width - 提供数据的剖面宽度，以地图单位为单位。大于 data_width/2 的值将不绘制。
%
% 可选输入：
%   small_circ_center [] - 提供一个 1 x 2 数组，包含用于投影数据的小圆中心坐标，使用 'ProjectSmallCircleOntoSwath' 进行投影。
%   dist_type ['mapdist'] - 控制 'data_width' 如何解释的选项。可选 'mapdist' 或 'angle'，默认值为 'mapdist'。仅当提供了 'small_circ_center' 时，'angle' 选项有效。
%   sample [] - 沿着地形剖面的重采样距离（地图单位）。如果未提供输入，代码将使用DEM的单元格大小（即不进行重采样）。
%   smooth [0] - 平滑距离，滤波器的宽度（地图单位），默认为 0，表示不进行平滑。
%   vex [10] - 垂直夸张度，用于拓扑剖面图。由于Matlab在控制物理轴尺寸上的限制，垂直夸张度的控制对某些图表（如'ksn_batch'，'ksn_profiler'）可能无效。
%   basin_value [] - 对于 'basin_stats' 选项，指定提供的数据表中的值名称，根据该值着色点。
%   basin_scale [] - 对于 'basin_stats' 选项，指定提供的数据表中的值名称，根据该值缩放点。
%   plot_map [true] - 逻辑标志，指示是否绘制地图，显示拓扑剖面位置及其上附加数据的位置。
%   cmap ['parula'] - 使用的颜色映射，支持颜色映射名称（如 'jet'）或n x 3的颜色映射数组。
%   save_figure [false] - 逻辑标志，指示是否将剖面图保存为PDF。
%
% 输出：
%   SW - TopoToolbox Swath对象，包含各种信息，作为结构体。可以使用plot(SW)绘制路径和剖面框，使用plotdz(SW)绘制剖面图。
%   SwathMat - n x 4 矩阵，包含沿剖面的距离、最小海拔、高程均值和最大海拔。
%   xypoints - n x 2 矩阵，包含每个剖面采样点的x,y坐标。
%   outData - 用于绘制附加数据的输出。根据data_type的不同，其形式不同：
%       'points3' - 距离、高程、距离基线的距离、x坐标、y坐标
%       'points4' - 距离、高程、值、距离基线的距离、x坐标、y坐标
%       'eqs' - 距离、深度、震级、距离基线的距离、x坐标、y坐标
%       'STREAMobj' - 距离、高程、距离基线的距离、x坐标、y坐标
%       'ksn_chandata' - 距离、高程、ksn、距离基线的距离、x坐标、y坐标
%       'ksn_batch' - 距离、ksn、距离基线的距离、x坐标、y坐标
%       'ksn_profiler' - 距离、ksn、距离基线的距离、x坐标、y坐标
%       'basin_stats' - 距离、平均流域海拔、'basin_value'、'basin_scale'（如果提供）、距离基线的距离、x坐标、y坐标
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数编写者：Adam M. Forte - 更新日期：06/18/18 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% 解析输入参数
	p = inputParser;
	p.FunctionName = 'MakeCombinedSwath';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'points',@(x) isempty(x) || isnumeric(x) & size(x,1)>=2 && size(x,2)==2);
	addRequired(p,'width',@(x) isscalar(x) && isnumeric(x));
	addRequired(p,'data_type',@(x) ischar(validatestring(x,{'points3','points4','points5','eqs','gps','STREAMobj','ksn_chandata','ksn_batch','ksn_profiler','basin_stats','basin_knicks'})));
	addRequired(p,'data');
	addRequired(p,'data_width',@(x) isnumeric(x) && isscalar(x));

	addParameter(p,'file_name_prefix','Combined',@(x) ischar(x));
	addParameter(p,'small_circ_center',[],@(x) isnumeric(x) && numel(x)==2 || isempty(x));
	addParameter(p,'dist_type','mapdist',@(x) ischar(validatestring(x,{'mapdist','angle'})));	
	addParameter(p,'sample',[],@(x) isscalar(x) && isnumeric(x) || isempty(x));
	addParameter(p,'smooth',0,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'vex',10,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'basin_value',[],@(x) ischar(x));
	addParameter(p,'basin_scale',[],@(x) ischar(x) || isempty(x));
	addParameter(p,'plot_map',true,@(x) isscalar(x) && islogical(x));
	addParameter(p,'cmap','parula',@(x) ischar(x) || isnumeric(x) & size(x,2)==3);
	addParameter(p,'save_figure',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'out_dir',[],@(x) isdir(x));
	addParameter(p,'sw_num',[],@(x) isscalar(x) && isnumeric(x)); % 循环运行时控制参数

	parse(p,DEM,points,width,data_type,data,data_width,varargin{:});
	DEM=p.Results.DEM;
	points=p.Results.points;
	wdth=p.Results.width;
	data_type=p.Results.data_type;
	data=p.Results.data;
	data_width=p.Results.data_width;

	fnp=p.Results.file_name_prefix;
	small_circ_center=p.Results.small_circ_center;
	dist_type=p.Results.dist_type;
	sample=p.Results.sample;
	smth=p.Results.smooth;
	vex=p.Results.vex;
	bv=p.Results.basin_value;
	bs=p.Results.basin_scale;
	plot_map=p.Results.plot_map;
	cmap=p.Results.cmap;
	save_figure=p.Results.save_figure;
	out_dir=p.Results.out_dir;
	sw_num=p.Results.sw_num;

	if isempty(sample)
		sample=DEM.cellsize;
	end

	if isempty(out_dir)
		out_dir=pwd;
	end

	if isempty(sw_num)
		fo=1;
		fe=2;
	else
		fo=(sw_num*2)-1;
		fe=sw_num*2;
	end

	if isempty(small_circ_center)
		proj_flag=1;
	else
		proj_flag=2;
		cx=small_circ_center(1);
		cy=small_circ_center(2);
	end

	% 生成地形剖面及相关数据集
	if isempty(points)
		SWt=SWATHobj(DEM);
		points=SWt.xy0;
		fig=gcf;
		close(fig);
	end
	
	[SW,SwathMat,xypoints,bends]=MakeTopoSwath(DEM,points,wdth,'sample',sample,'smooth',smth,'make_shape',false);
	swdist=SwathMat(:,1);
	min_elevs=SwathMat(:,2);
	mean_elevs=SwathMat(:,3);
	max_elevs=SwathMat(:,4);

	% 获取DEM范围
	[demx,demy]=getoutline(DEM,true);	

	% 设置颜色映射
	colormap(cmap);

	% 根据数据类型执行不同处理流程

	switch data_type
		case 'points3'
			x_coord=data(:,1);
			y_coord=data(:,2);
			z=data(:,3);

			% 移除超出DEM范围的数据点
			[demix]=inpolygon(x_coord,y_coord,demx,demy);	
			x_coord=x_coord(demix); y_coord=y_coord(demix); z=z(demix);

			switch proj_flag
			case 1
				[ds,db]=ProjectOntoSwath(SW,x_coord,y_coord,data_width);
				outData=[ds z db x_coord y_coord];
				idx=outData(:,3)<=(data_width/2) & ~isnan(ds);
			case 2
				[ds,db,dab]=ProjectSmallCircleOntoSwath(SW,x_coord,y_coord,cx,cy);
				switch dist_type
				case 'mapdist'
					outData=[ds z db x_coord y_coord];
					idx=abs(outData(:,3))<=(data_width/2) & ~isnan(ds);
				case 'angle'
					outData=[ds z dab x_coord y_coord];
					idx=abs(outData(:,3))<=(data_width/2) & ~isnan(ds);
				end
			end

			% 绘制地形剖面图
			f1=figure(fo);
			clf 
			set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.4],'renderer','painters');

			hold on
			plot(swdist,min_elevs,'-k');
			plot(swdist,max_elevs,'-k');
			plot(swdist,mean_elevs,'-k','LineWidth',2);

			daspect([vex 1 1])

			yl=ylim;
			for jj=1:numel(bends)
				plot([bends(jj),bends(jj)],yl,'-k');
			end

			scatter(outData(idx,1),outData(idx,2),30,'k','filled');

			xlabel('沿剖面线距离 (m)');
			ylabel('高程 (m)');
			xlim([0 max(swdist)]);
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(gca);
		    end	
			hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
		case 'points4'
			x_coord=data(:,1);
			y_coord=data(:,2);
			z=data(:,3);
			col=data(:,4);

			% 移除超出DEM范围的数据点
			[demix]=inpolygon(x_coord,y_coord,demx,demy);	
			x_coord=x_coord(demix); y_coord=y_coord(demix); z=z(demix); col=col(demix);

			switch proj_flag
			case 1
				[ds,db]=ProjectOntoSwath(SW,x_coord,y_coord,data_width);
				outData=[ds z col db x_coord y_coord];
				idx=outData(:,4)<=(data_width/2) & ~isnan(ds);
			case 2
				[ds,db,dab]=ProjectSmallCircleOntoSwath(SW,x_coord,y_coord,cx,cy);
				switch dist_type
				case 'mapdist'
					outData=[ds z col db x_coord y_coord];
					idx=abs(outData(:,4))<=(data_width/2) & ~isnan(ds);
				case 'angle'
					outData=[ds z col dab x_coord y_coord];
					idx=abs(outData(:,4))<=(data_width/2) & ~isnan(ds);
				end
			end			

			% 绘制地形剖面图
			f1=figure(fo);
			clf 
			set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.4],'renderer','painters');

			hold on
			plot(swdist,min_elevs,'-k');
			plot(swdist,max_elevs,'-k');
			plot(swdist,mean_elevs,'-k','LineWidth',2);

			daspect([vex 1 1])

			yl=ylim;
			for jj=1:numel(bends)
				plot([bends(jj),bends(jj)],yl,'-k');
			end

			scatter(outData(idx,1),outData(idx,2),30,outData(idx,3),'filled');

			c1=colorbar;
			xlabel(c1,'用户自定义值')

			xlabel('沿剖面线距离 (m)');
			ylabel('高程 (m)');
			xlim([0 max(swdist)]);
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(gca);
		    end				
			hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
		case 'points5'
			x_coord=data(:,1);
			y_coord=data(:,2);
			z=data(:,3);
			col=data(:,4)
			scle=data(:,5);

			% 移除超出DEM范围的数据点
			[demix]=inpolygon(x_coord,y_coord,demx,demy);	
			x_coord=x_coord(demix); y_coord=y_coord(demix); z=z(demix); col=col(demix); scle=scle(demix);

			switch proj_flag
			case 1
				[ds,db]=ProjectOntoSwath(SW,x_coord,y_coord,data_width);
				outData=[ds z col scle db x_coord y_coord];
				idx=outData(:,5)<=(data_width/2) & ~isnan(ds);
			case 2
				[ds,db,dab]=ProjectSmallCircleOntoSwath(SW,x_coord,y_coord,cx,cy);
				switch dist_type
				case 'mapdist'
					outData=[ds z col scle db x_coord y_coord];
					idx=abs(outData(:,5))<=(data_width/2) & ~isnan(ds);
				case 'angle'
					outData=[ds z col scle dab x_coord y_coord];
					idx=abs(outData(:,5))<=(data_width/2) & ~isnan(ds);
				end
			end

			% 绘制地形剖面图
			f1=figure(fo);
			clf 
			set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.4],'renderer','painters');

			hold on
			plot(swdist,min_elevs,'-k');
			plot(swdist,max_elevs,'-k');
			plot(swdist,mean_elevs,'-k','LineWidth',2);

			daspect([vex 1 1])

			% 缩放尺寸向量
			sz_val=outData(idx,4);
			sz=(sz_val/max(sz_val))*100;
			% 创建尺寸图例
			sz_sizes=linspace(min(sz),max(sz),5);
			sz_val_scale=(sz_sizes/100)*max(sz_val);
			for ii=1:5
				sz_leg(ii)=plot(0,0,'ko','MarkerSize',sqrt(sz_sizes(ii)),'MarkerFaceColor','k','LineStyle','none');
				set(sz_leg(ii),'visible','off');
				leg_ent{ii}=num2str(sz_val_scale(ii));
			end	

			yl=ylim;
			for jj=1:numel(bends)
				plot([bends(jj),bends(jj)],yl,'-k');
			end

			scatter(outData(idx,1),outData(idx,2),sz,outData(idx,3),'filled');
			c1=colorbar('southoutside');
			xlabel(c1,'用户自定义值1')
			legend(sz_leg,leg_ent);
			xlabel('沿剖面线距离 (m)');
			ylabel('高程 (m)');
			xlim([0 max(swdist)]);
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(gca);
		    end				
			hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
		case 'eqs'
			x_coord=data(:,1);
			y_coord=data(:,2);
			depth=data(:,3);
			magnitude=data(:,4);

			% 移除超出DEM范围的数据点
			[demix]=inpolygon(x_coord,y_coord,demx,demy);	
			x_coord=x_coord(demix); y_coord=y_coord(demix); depth=depth(demix); magnitude=magnitude(demix);

			switch proj_flag
			case 1
				[ds,db]=ProjectOntoSwath(SW,x_coord,y_coord,data_width);
				outData=[ds depth magnitude db x_coord y_coord];
				idx=outData(:,4)<=(data_width/2) & ~isnan(ds);
			case 2
				[ds,db,dab]=ProjectSmallCircleOntoSwath(SW,x_coord,y_coord,cx,cy);
				switch dist_type
				case 'mapdist'
					outData=[ds depth magnitude db x_coord y_coord];
					idx=abs(outData(:,4))<=(data_width/2) & ~isnan(ds);
				case 'angle'
					outData=[ds depth magnitude dab x_coord y_coord];
					idx=abs(outData(:,4))<=(data_width/2) & ~isnan(ds);
				end
			end				

			% 绘制地形剖面图
			f1=figure(fo);
			clf 
			set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.4],'renderer','painters');

			ax1=subplot(2,1,1);
			hold on
			plot(swdist,min_elevs,'-k');
			plot(swdist,max_elevs,'-k');
			plot(swdist,mean_elevs,'-k','LineWidth',2);

			yl=ylim;
			for jj=1:numel(bends)
				plot([bends(jj),bends(jj)],yl,'-k');
			end
			ylabel('高程 (m)');
			xlim([0 max(swdist)]);
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(ax1);
		    end				
			hold off

			ax2=subplot(2,1,2);
			hold on
			% 缩放尺寸向量
			sz_val=outData(idx,3);
			sz=(sz_val/max(sz_val))*100;
			% 创建尺寸图例
			sz_sizes=linspace(min(sz),max(sz),5);
			sz_val_scale=(sz_sizes/100)*max(sz_val);
			for ii=1:5
				sz_leg(ii)=plot(0,0,'ko','MarkerSize',sqrt(sz_sizes(ii)),'MarkerFaceColor','k','LineStyle','none');
				set(sz_leg(ii),'visible','off');
				leg_ent{ii}=num2str(sz_val_scale(ii));
			end		

			scatter(outData(idx,1),outData(idx,2),sz,outData(idx,4),'filled');
			xlabel('沿剖面线距离 (m)');
			ylabel('深度 (km')
			legend(sz_leg,leg_ent);
			xlim([0 max(swdist)]);
			c1=colorbar(ax2,'southoutside');
			xlabel(c1,'距剖面线距离')
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(ax2);
		    end				
			hold off

			set(ax2,'YDir','reverse');
			linkaxes([ax1,ax2],'x')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
		case 'gps'
			x_coord=data(:,1);
			y_coord=data(:,2);
			nc=data(:,3);
			ec=data(:,4);
			nu=data(:,5);
			eu=data(:,6);

			% 移除超出DEM范围的数据点
			[demix]=inpolygon(x_coord,y_coord,demx,demy);	
			x_coord=x_coord(demix); y_coord=y_coord(demix); nc=nc(demix); ec=ec(demix); nu=nu(demix); eu=eu(demix);

			switch proj_flag
			case 1
				[ds,db,mag,unc,nc0,ec0]=ProjectGPSOntoSwath(SW,x_coord,y_coord,data_width,nc,ec,nu,eu);
				outData=[ds mag unc nc0 ec0 db x_coord y_coord];
				idx=outData(:,6)<=(data_width/2) & ~isnan(ds);
			case 2
				[~,~,mag,unc,nc0,ec0]=ProjectGPSOntoSwath(SW,x_coord,y_coord,data_width,nc,ec,nu,eu);
				[ds,db,dab]=ProjectSmallCircleOntoSwath(SW,x_coord,y_coord,cx,cy);
				switch dist_type
				case 'mapdist'
					outData=[ds mag unc nc0 ec0 db x_coord y_coord];
					idx=abs(outData(:,6))<=(data_width/2) & ~isnan(ds);
				case 'angle'
					outData=[ds mag unc nc0 ec0 dab x_coord y_coord];
					idx=abs(outData(:,6))<=(data_width/2) & ~isnan(ds);
				end
			end

			% 绘制地形剖面图
			f1=figure(fo);
			clf 
			set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.4],'renderer','painters');

			ax1=subplot(2,1,1);
			hold on
			plot(swdist,min_elevs,'-k');
			plot(swdist,max_elevs,'-k');
			plot(swdist,mean_elevs,'-k','LineWidth',2);

			yl=ylim;
			for jj=1:numel(bends)
				plot([bends(jj),bends(jj)],yl,'-k');
			end
			ylabel('高程 (m)');
			xlim([0 max(swdist)]);
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(ax1);
		    end				
			hold off

			ax2=subplot(2,1,2);
			hold on	
			e1=errorbar(outData(idx,1),outData(idx,2),outData(idx,3),'.');
			e1.CapSize=0;
			e1.Color='k';
			scatter(outData(idx,1),outData(idx,2),30,outData(idx,6),'filled');
			xlabel('沿剖面线距离 (m)');
			ylabel('剖面线方向速度 (mm/yr)');
			xlim([0 max(swdist)]);
			c1=colorbar(ax2,'southoutside');
			xlabel(c1,'距剖面线距离')
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(ax2);
		    end				
			hold off

			linkaxes([ax1,ax2],'x')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		case 'STREAMobj'
			x_coord=data.x;
			y_coord=data.y;
			z=getnal(data,DEM);

			% 移除超出DEM范围的数据点
			[demix]=inpolygon(x_coord,y_coord,demx,demy);	
			x_coord=x_coord(demix); y_coord=y_coord(demix); z=z(demix);

			switch proj_flag
			case 1
				[ds,db]=ProjectOntoSwath(SW,x_coord,y_coord,data_width);
				outData=[ds z db x_coord y_coord];
				idx=outData(:,3)<=(data_width/2) & ~isnan(ds);
			case 2
				[ds,db,dab]=ProjectSmallCircleOntoSwath(SW,x_coord,y_coord,cx,cy);
				switch dist_type
				case 'mapdist'
					outData=[ds z db x_coord y_coord];
					idx=abs(outData(:,3))<=(data_width/2) & ~isnan(ds);
				case 'angle'
					outData=[ds z dab x_coord y_coord];
					idx=abs(outData(:,3))<=(data_width/2) & ~isnan(ds);
				end
			end

			% 绘制地形剖面图
			f1=figure(fo);
			clf 
			set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.4],'renderer','painters');

			hold on
			plot(swdist,min_elevs,'-k');
			plot(swdist,max_elevs,'-k');
			plot(swdist,mean_elevs,'-k','LineWidth',2);

			daspect([vex 1 1])

			yl=ylim;
			for jj=1:numel(bends)
				plot([bends(jj),bends(jj)],yl,'-k');
			end

			scatter(outData(idx,1),outData(idx,2),5,'k','filled');

			xlabel('沿剖面线距离 (m)');
			ylabel('高程 (m)');
			xlim([0 max(swdist)]);
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(gca);
		    end				
			hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		case 'ksn_chandata'
			x_coord=data(:,9);
			y_coord=data(:,10);
			ksn=data(:,8);
			elev=data(:,4);

			% 移除超出DEM范围的数据点
			[demix]=inpolygon(x_coord,y_coord,demx,demy);	
			x_coord=x_coord(demix); y_coord=y_coord(demix); ksn=ksn(demix); elev=elev(demix);

			switch proj_flag
			case 1
				[ds,db]=ProjectOntoSwath(SW,x_coord,y_coord,data_width);
				outData=[ds elev ksn db x_coord y_coord];
				idx=outData(:,4)<=(data_width/2) & ~isnan(ds);
			case 2
				[ds,db,dab]=ProjectSmallCircleOntoSwath(SW,x_coord,y_coord,cx,cy);
				switch dist_type
				case 'mapdist'
					outData=[ds elev ksn db x_coord y_coord];
					idx=abs(outData(:,4))<=(data_width/2) & ~isnan(ds);
				case 'angle'
					outData=[ds elev ksn dab x_coord y_coord];
					idx=abs(outData(:,4))<=(data_width/2) & ~isnan(ds);
				end
			end			

			% 绘制地形剖面图
			f1=figure(fo);
			clf 
			set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.4],'renderer','painters');

			ax1=subplot(2,1,1);
			hold on

			plot(swdist,min_elevs,'-k');
			plot(swdist,max_elevs,'-k');
			plot(swdist,mean_elevs,'-k','LineWidth',2);

			yl=ylim;
			for jj=1:numel(bends)
				plot([bends(jj),bends(jj)],yl,'-k');
			end

			xlabel('沿剖面线距离 (m)');
			ylabel('高程 (m)');
			xlim([0 max(swdist)]);
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(ax1);
		    end				
			hold off

			ax2=subplot(2,1,2);
			hold on
			scatter(outData(idx,1),outData(idx,3),30,'k','filled');
			xlabel('沿剖面线距离 (m)');
			ylabel('KSN');
			xlim([0 max(swdist)]);
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(ax2);
		    end	
			hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
case 'basin_stats'
			x_coord=data.center_x;
			y_coord=data.center_y;
			z=data.mean_el;
			col=data.(bv);
			if ~isempty(bs)
				scl=data.(bs);
				if ~isnumeric(scl)
					if isdeployed
						errordlg('用于缩放点的值必须是数值型') % 错误提示中文化
					end
					error('用于缩放点的值必须是数值型')
				end
			end

			if ~isnumeric(col);
				if isdeployed
					errordlg('用于着色的值必须是数值型') % 错误提示中文化
				end
				error('用于着色的值必须是数值型')
			end

			% 移除超出DEM范围的数据点
			[demix]=inpolygon(x_coord,y_coord,demx,demy);	
			x_coord=x_coord(demix); y_coord=y_coord(demix); z=z(demix); col=col(demix);
			if ~isempty(bs)
				scl=scl(demix);
			end

			switch proj_flag
			case 1
				% 数据投影转换
				[ds,db]=ProjectOntoSwath(SW,x_coord,y_coord,data_width);

				% 组装输出数据
				if isempty(bs)
					outData=[ds z col db x_coord y_coord];
				else 
					outData=[ds z col scl db x_coord y_coord];	
				end			

				% 根据数据宽度过滤
				if isempty(bs)
					idx=outData(:,4)<=(data_width/2) & ~isnan(ds);
				else
					idx=outData(:,5)<=(data_width/2) & ~isnan(ds);
				end
			case 2
				[ds,db,dab]=ProjectSmallCircleOntoSwath(SW,x_coord,y_coord,cx,cy);

				switch dist_type
				case 'mapdist'
					% 组装输出数据
					if isempty(bs)
						outData=[ds z col db x_coord y_coord];
					else 
						outData=[ds z col scl db x_coord y_coord];	
					end			

					% 根据数据宽度过滤
					if isempty(bs)
						idx=abs(outData(:,4))<=(data_width/2) & ~isnan(ds);
					else
						idx=abs(outData(:,5))<=(data_width/2) & ~isnan(ds);
					end
				case 'angle'
					% 组装输出数据
					if isempty(bs)
						outData=[ds z col dab x_coord y_coord];
					else 
						outData=[ds z col scl dab x_coord y_coord];	
					end			

					% 根据数据宽度过滤
					if isempty(bs)
						idx=abs(outData(:,4))<=(data_width/2) & ~isnan(ds);
					else
						idx=abs(outData(:,5))<=(data_width/2) & ~isnan(ds);
					end
				end
			end			

			% 绘制地形剖面图
			f1=figure(fo);
			clf 
			set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.4],'renderer','painters');

			hold on
			plot(swdist,min_elevs,'-k');
			plot(swdist,max_elevs,'-k');
			plot(swdist,mean_elevs,'-k','LineWidth',2);

			daspect([vex 1 1])

			yl=ylim;
			for jj=1:numel(bends)
				plot([bends(jj),bends(jj)],yl,'-k');
			end

			if isempty(bs)
				scatter(outData(idx,1),outData(idx,2),30,outData(idx,3),'filled');
			else
				% 缩放尺寸向量
				sz_val=outData(idx,4);
				sz=(sz_val/max(sz_val))*100;
				% 创建尺寸图例
				sz_sizes=linspace(min(sz),max(sz),5);
				sz_val_scale=(sz_sizes/100)*max(sz_val);
				for ii=1:5
					sz_leg(ii)=plot(0,0,'ko','MarkerSize',sqrt(sz_sizes(ii)),'MarkerFaceColor','k','LineStyle','none');
					set(sz_leg(ii),'visible','off');
					leg_ent{ii}=num2str(sz_val_scale(ii));
				end
				scatter(outData(idx,1),outData(idx,2),sz,outData(idx,3),'filled');
				legend(sz_leg,leg_ent);
				bs_n=strrep(bs,'_',' ');
				title(['点按 ' bs_n ' 缩放']); % 标题中文化
			end

			c1=colorbar;
			bv_n=strrep(bv,'_',' ');
			ylabel(c1,bv_n);

			xlabel('沿剖面线距离 (m)'); % 坐标轴标签中文化
			ylabel('高程 (m)');
			xlim([0 max(swdist)]);
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(gca);
		    end				
			hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		case 'ksn_batch'
			% 统计KSN地图结构中的分段数量	
			numSegs=numel(data);
			% 遍历河道分段提取坐标和ksn值
			streamData=zeros(numSegs,3);
			for kk=1:numSegs
				xx=mean(data.ksn_ms(kk,1).X);
				yy=mean(data.ksn_ms(kk,1).Y);
				ksn=mean(data.ksn_ms(kk,1).ksn);
				streamData(kk,:)=[xx yy ksn];
			end	

			x_coord=streamData(:,1);
			y_coord=streamData(:,2);
			ksn=streamData(:,3);

			% 移除超出DEM范围的数据点
			[demix]=inpolygon(x_coord,y_coord,demx,demy);	
			x_coord=x_coord(demix); y_coord=y_coord(demix); ksn=ksn(demix);

			switch proj_flag
			case 1
				[ds,db]=ProjectOntoSwath(SW,x_coord,y_coord,data_width);
				outData=[ds ksn db x_coord y_coord];
				idx=outData(:,3)<=(data_width/2) & ~isnan(ds);
			case 2
				[ds,db,dab]=ProjectSmallCircleOntoSwath(SW,x_coord,y_coord,cx,cy);
				switch dist_type
				case 'mapdist'
					outData=[ds ksn db x_coord y_coord];
					idx=abs(outData(:,3))<=(data_width/2) & ~isnan(ds);
				case 'angle'
					outData=[ds ksn dab x_coord y_coord];
					idx=abs(outData(:,3))<=(data_width/2) & ~isnan(ds);
				end
			end			

			% 绘制地形剖面图
			f1=figure(fo);
			clf 
			set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.4],'renderer','painters');

			subplot(2,1,1);
			hold on
			plot(swdist,min_elevs,'-k');
			plot(swdist,max_elevs,'-k');
			plot(swdist,mean_elevs,'-k','LineWidth',2);

			yl=ylim;
			for jj=1:numel(bends)
				plot([bends(jj),bends(jj)],yl,'-k');
			end

			xlabel('沿剖面线距离 (m)'); % 坐标轴标签中文化
			ylabel('高程 (m)');
			xlim([0 max(swdist)]);
			hold off

			subplot(2,1,2);
			hold on
			scatter(outData(idx,1),outData(idx,2),20,'k','filled');
			xlabel('沿剖面线距离 (m)');
			ylabel('KSN');
			xlim([0 max(swdist)]);
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(gca);
		    end				
			hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%			
		case 'ksn_profiler'
			x_coord=data.knl(:,1);
			y_coord=data.knl(:,2);
			ksn=data.knl(:,4);

			% 移除超出DEM范围的数据点
			[demix]=inpolygon(x_coord,y_coord,demx,demy);	
			x_coord=x_coord(demix); y_coord=y_coord(demix); ksn=ksn(demix);			

			switch proj_flag
			case 1
				[ds,db]=ProjectOntoSwath(SW,x_coord,y_coord,data_width);
				outData=[ds ksn db x_coord y_coord];
				idx=outData(:,3)<=(data_width/2) & ~isnan(ds);
			case 2
				[ds,db,dab]=ProjectSmallCircleOntoSwath(SW,x_coord,y_coord,cx,cy);
				switch dist_type
				case 'mapdist'
					outData=[ds ksn db x_coord y_coord];
					idx=abs(outData(:,3))<=(data_width/2) & ~isnan(ds);
				case 'angle'
					outData=[ds ksn dab x_coord y_coord];
					idx=abs(outData(:,3))<=(data_width/2) & ~isnan(ds);
				end
			end				

			% 绘制地形剖面图
			f1=figure(fo);
			clf 
			set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.4],'renderer','painters');

			ax1=subplot(2,1,1);
			hold on
			plot(swdist,min_elevs,'-k');
			plot(swdist,max_elevs,'-k');
			plot(swdist,mean_elevs,'-k','LineWidth',2);

			yl=ylim;
			for jj=1:numel(bends)
				plot([bends(jj),bends(jj)],yl,'-k');
			end

			xlabel('沿剖面线距离 (m)');
			ylabel('高程 (m)');
			xlim([0 max(swdist)]);
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(ax1);
		    end				
			hold off

			ax2=subplot(2,1,2);
			hold on
			scatter(outData(idx,1),outData(idx,2),20,'k','filled');
			xlabel('沿剖面线距离 (m)');
			ylabel('KSN');
			xlim([0 max(swdist)]);
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(ax2);
		    end				
			hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		case 'basin_knicks'
			% 验证输入目录有效性
			if ~isdir(data)
				if isdeployed
					errordlg('使用"basin_knicks"时必须提供包含"Knicks_*.mat"文件的有效目录') % 错误提示中文化
				end
				error('使用"basin_knicks"时必须提供包含"Knicks_*.mat"文件的有效目录')
			end

			current=pwd;
			cd(data);

			% 加载断点数据文件
			fileList=dir('Knicks_*.mat');
			if isempty(fileList)
				if isdeployed
					errordlg('指定目录中未找到Knicks_*.mat文件') % 错误提示中文化
				end
				error('指定目录中未找到Knicks_*.mat文件')
			end
			knps=cell(numel(fileList),1);
			for jj=1:numel(fileList)
				load(fileList(jj,1).name);
				knps{jj}=[KnickTable.x_coord KnickTable.y_coord KnickTable.elevation];
			end

			knps=vertcat(knps{:});

			x_coord=knps(:,1);
			y_coord=knps(:,2);
			z=knps(:,3);

			% 移除超出DEM范围的数据点
			[demix]=inpolygon(x_coord,y_coord,demx,demy);	
			x_coord=x_coord(demix); y_coord=y_coord(demix); z=z(demix);

			switch proj_flag
			case 1
				[ds,db]=ProjectOntoSwath(SW,x_coord,y_coord,data_width);
				outData=[ds z db x_coord y_coord];
				idx=outData(:,3)<=(data_width/2) & ~isnan(ds);
			case 2
				[ds,db,dab]=ProjectSmallCircleOntoSwath(SW,x_coord,y_coord,cx,cy);
				switch dist_type
				case 'mapdist'
					outData=[ds z db x_coord y_coord];
					idx=abs(outData(:,3))<=(data_width/2) & ~isnan(ds);
				case 'angle'
					outData=[ds z dab x_coord y_coord];
					idx=abs(outData(:,3))<=(data_width/2) & ~isnan(ds);
				end
			end				

			% 绘制地形剖面图
			f1=figure(fo);
			clf 
			set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.4],'renderer','painters');

			hold on
			plot(swdist,min_elevs,'-k');
			plot(swdist,max_elevs,'-k');
			plot(swdist,mean_elevs,'-k','LineWidth',2);

			daspect([vex 1 1])

			yl=ylim;
			for jj=1:numel(bends)
				plot([bends(jj),bends(jj)],yl,'-k');
			end

			scatter(outData(idx,1),outData(idx,2),20,'k','filled');

			xlabel('沿剖面线距离 (m)');
			ylabel('高程 (m)');
			xlim([0 max(swdist)]);
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(gca);
		    end	
			hold off
	end

	% 生成边界形状文件
	ms=struct;
	ms(1,1).Geometry='Line';
	ms(1,1).X=SW.xy0(:,1);
	ms(1,1).Y=SW.xy0(:,2);
	ms(1,1).Type='中心线';

	if ~verLessThan('matlab','9.4')
		Tverts=SwathPolygon(SW,wdth);
		ms(2,1).Geometry='Line';
		ms(2,1).X=Tverts(:,1);
		ms(2,1).Y=Tverts(:,2);
		ms(2,1).Type='地形宽度';

		Dverts=SwathPolygon(SW,data_width);
		ms(3,1).Geometry='Line';
		ms(3,1).X=Dverts(:,1);
		ms(3,1).Y=Dverts(:,2);
		ms(3,1).Type='数据宽度';
	end

	shapewrite(ms,fullfile(out_dir,'SwathBounds.shp'));

	% 绘制空间位置图
	if plot_map
		f2=figure(fe);
		set(f2,'Units','normalized','Position',[0.05 0.1 0.6 0.6]);
		hold on
		imageschs(DEM,DEM,'colormap','gray');
		plot(SW.xy0(:,1),SW.xy0(:,2),'-g','LineWidth',0.5);
		if ~verLessThan('matlab','9.4')
			plot(Tverts(:,1),Tverts(:,2),'-g','LineWidth',0.5);
			plot(Dverts(:,1),Dverts(:,2),'-r','LineWidth',0.5);	
		end	
		scatter(x_coord(idx),y_coord(idx),20,'r','filled');
		% 避免显示过多非相关点
		if ~any([strcmp(data_type,'STREAMobj') strcmp(data_type,'ksn_profiler') strcmp(data_type,'ksn_batch')])
			scatter(x_coord(~idx),y_coord(~idx),20,'w','filled');
		end
		if ~verLessThan('matlab','9.5')
	        disableDefaultInteractivity(gca);
	    end	
		hold off
	end

	% 保存输出图形
	if save_figure && isempty(sw_num)
		orient(f1,'Landscape')
		print(f1,'-dpdf','-bestfit',fullfile(out_dir,[fnp '_Swath.pdf']));
	elseif save_figure && ~isempty(sw_num)
		orient(f1,'Landscape')
		print(f1,'-dpdf','-bestfit',fullfile(out_dir,[fnp '_Swath_' num2str(sw_num) '.pdf']));
	end

end

function [verts]=SwathPolygon(SW,w)
	% 生成剖面多边形顶点
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






