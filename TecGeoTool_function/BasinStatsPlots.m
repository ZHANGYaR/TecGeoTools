function BasinStatsPlots(basin_table,plots,varargin)
	%
% 使用方法：
%   BasinStatsPlots(basin_table,plots);
%   BasinStatsPlots(basin_table,plots,'name',value,...);
%
% 描述：
%   本函数用于处理来自 'ProcessRiverBasins' 和 'SubDivideBigBasins' 的输出，并生成各种汇总流域值的绘图。
%
% 必需输入：
%   basin_table - 来自 'CompileBasinStats' 的表格输出
%   plots - 要生成的绘图类型，有效输入包括：
%       'grd_ksn' - 绘制平均流域坡度与平均流域河道陡度的图（例如，参见 Forte 等，2016，地球与行星科学快报，讨论了这些图的用途）
%       'grd_rlf' - 类似于 'grd_ksn'，但使用局部高程替代 ksn，要求在运行 'ProcessRiverBasins' 时计算了高程，假定高程半径为 2500（可以通过 'rlf_radius' 可选参数设置替代半径）
%       'rlf_ksn' - 绘制平均流域高程与平均流域河道陡度的图
%       'compare_filtered' - 比较平均值与筛选后的平均值的图，如果在运行 'CompileBasinStats' 时使用了筛选
%       'category_mean_hist' - 如果在运行 'CompileBasinStats' 时计算了 'means_by_category'，可以绘制按类别分布的均值直方图，使用此选项。需要输入 'cat_mean1'
%       'category_mean_compare' - 如果为多个值计算了 'means_by_category'（例如，坡度和 ksn），可以通过此图比较按类别的均值。需要输入 'cat_mean1'（绘制在 x 轴上的值）和 'cat_mean2'（绘制在 y 轴上的值）
%       'stacked_hypsometry' - 绘制流域的高程分布图
%       'compare_mean_and_dist' - 绘制选定流域或所有流域中的某统计量的直方图，以便与均值进行比较。需要输入 'statistic_of_interest' 和 'basin_num'
%       'scatterplot_matrix' - 散点图和直方图矩阵，类似于 R 中的 'lattice' 图。提供一个计算了筛选均值的表格，并将 'use_filtered' 设置为 false，可能会生成一个大矩阵
%       'xy' - 通用绘图，需要输入可选的 'xval' 和 'yval'
%
%%%%%%%%%%%%%%%%%%
% 可选输入：
%
%%% 一般参数
%   uncertainty ['se'] - 绘图中使用的误差类型，有效选项为 'se'（标准误差）、'std'（标准差）或 'none'。如果提供 'none'，则表示不希望绘制误差条。此选项的行为取决于您如何运行 'ProcessRiverBasins'，例如，如果在运行 'ProcessRiverBasins' 时只计算了标准差，而在此提供 'se'，代码将忽略您的选择，使用标准差值。
%   use_filtered [false] - 逻辑标志，指示是否使用筛选后的值进行 'grd_ksn'、'grd_rlf'、'rlf_ksn' 或 'scatterplot_matrix' 绘图。仅在运行 'CompileBasinStats' 时计算了筛选值时有效。
%   color_by [] - 用于着色点的值，仅适用于 'grd_ksn'、'grd_rlf'、'rlf_ksn' 和 'xy'，可以是提供表格中的列名或与提供表格长度相同的 m x 1 数字值数组
%   cmap [] - 如果提供了 'color_by'，则使用的颜色图，可以是标准颜色图的名称，也可以是 m x 3 数组的 RGB 值，用作颜色图。
%   define_region [] - 一组坐标，用于定义要绘制数据的矩形区域，期望输入四个元素的矩阵（行或列），定义最小 x、最大 x、最小 y 和最大 y 坐标，或者定义为 true 来绘制所有流域中心图，供您通过绘制矩形选择区域。适用于所有图形。
%   rlf_radius [2500] - 绘制与高程相关的值时使用的高程半径
%   save_figure [false] - 逻辑标志，指示是否保存所有绘制的图形为 PDF 文件
%
%%% xy 绘图
%   xval [] - xy 绘图中绘制在 x 轴上的值，可以是表格中列的名称，或者是与提供表格长度相同的 m x 1 数字值数组
%   yval [] - xy 绘图中绘制在 y 轴上的值，可以是表格中列的名称，或者是与提供表格长度相同的 m x 1 数字值数组
%
%%% compare_mean_and_dist 绘图
%   statistic_of_interest ['ksn'] - 用于绘制直方图并与均值进行比较的统计量。有效输入为 'ksn'、'gradient'、'elevation'、'relief'（如果提供了高程，代码将查找在指定的 'rlf_radius' 半径计算的高程），或者提供给 'ProcessRiverBasins' 的其他网格名称，例如，如果提供了降水网格，并提供了名称 'precip' 且表格中存在名为 'mean_precip' 的列，那么 'precip' 就是此参数的有效输入。
%   basin_num [] - 要用于 'compare_mean_and_dist' 的流域编号（如表格中的 ID 列所示），如果为空，则 'compare_mean_and_dist' 将使用所有流域。
%
%%% category_mean_hist 或 category_mean_compare
%   cat_mean1 [] - 用于绘制的类别，请参见 'category_mean_hist' 或 'category_mean_compare'，有效输入为 'ksn'、'rlf'、'gradient'，或 'ProcessRiverMeans' 提供的其他网格名称。
%   cat_mean2 [] - 用于绘制的类别，仅适用于 'category_mean_compare'，有效输入为 'ksn'、'rlf'、'gradient'，或 'ProcessRiverMeans' 提供的其他网格名称。
%
%%% 拟合坡度-Ksn 关系
%   fit_grd_ksn [false] - 逻辑标志，启动拟合坡度 - ksn 关系。仅在绘图类型设置为 'grd_ksn' 时产生结果。该关系使用与侵蚀率和河道陡度以及侵蚀率和平均坡度之间的幂律关系进行拟合。参见 Forte 等，2016，地球与行星科学快报以获取更多讨论。拟合优化坡度扩散率（D）、河流侵蚀性（K）和阈值坡度（Sc）的值。最佳拟合值将在控制台中打印。
%   start_diffusivity [0.01] - 用于优化坡度扩散率参数的初始值。
%   start_erodibility [1e-7] - 用于优化河流侵蚀性参数的初始值。
%   start_threshold_gradient [0.8] - 用于优化阈值坡度参数的初始值。
%   n_val [2] - 坡度参数的 n 值，这是拟合中的固定参数。
%
%%% 拟合高程-Ksn 关系
%   fit_rlf_ksn [false] - 逻辑标志，启动高程 - ksn 关系的简单线性拟合（预期是线性关系）。仅在绘图类型设置为 'rlf_ksn' 时产生结果。
%
%%% 拟合筛选数据
%   fit_filtered [false] - 逻辑标志，启动对筛选数据的简单线性拟合。仅在绘图类型设置为 'compare_filtered' 时产生结果。
%
%%%%%%%%%%%%%%%%%%
% 示例：
%
%   % 绘制平均流域坡度与 2500 m^2 高程的图，使用默认的高程半径
%   BasinStatsPlots(T,'grd_rlf');
%
%   % 绘制平均流域坡度与平均流域河道陡度的图，并按平均高程着色
%   BasinStatsPlots(T,'grd_ksn','color_by','mean_el');
%
%   % 绘制平均流域坡度与平均流域河道陡度的图，并按地质模式着色，颜色图已缩放为提供不同颜色的单元
%   cmap=colorcube(numel(unique(T.mode_geology)));
%   BasinStatsPlots(T,'grd_ksn','color_by','mode_geology','cmap',cmap);
%
%   % 绘制平均流域河道陡度与流域排水面积的图
%   BasinStatsPlots(T,'xy','xval','drainage_area','yval','mean_ksn');
%
%   % 按单独类别（如地质）绘制均值 2500 m^2 高程的直方图
%   BasinStatsPlots(T,'category_mean_hist','cat_mean1','rlf');
%
%   % 绘制平均流域坡度与平均流域河道陡度的图，并按类别绘制
%   BasinStatsPlots(T,'category_mean_compare','cat_mean1','ksn','cat_mean2','gradient');
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数编写者：Adam M. Forte - 更新于：06/18/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	% 解析输入
	p = inputParser;
	p.FunctionName = 'BasinStatsPlots';
	addRequired(p,'basin_table',@(x) isa(x,'table'));
	addRequired(p,'plots',@(x) ischar(validatestring(x,{'grd_ksn','grd_rlf','rlf_ksn',...
		'compare_filtered','category_mean_hist','category_mean_compare','xy','stacked_hypsometry',...
		'compare_mean_and_dist','scatterplot_matrix'})));

	addParameter(p,'uncertainty','se',@(x) ischar(validatestring(x,{'se','std','none'})));
	addParameter(p,'use_filtered',false,@(x) islogical(x) && isscalar(x));
	addParameter(p,'color_by',[],@(x) ischar(x) || isnumeric(x) && size(x,2)==1 || isempty(x));
	addParameter(p,'cmap',[],@(x) ischar(x) || isnumeric(x) && size(x,2)==3);
	addParameter(p,'xval',[],@(x) ischar(x) || isnumeric(x) && size(x,2)==1);
	addParameter(p,'yval',[],@(x) ischar(x) || isnumeric(x) && size(x,2)==1);
	addParameter(p,'define_region',[],@(x) isnumeric(x) && numel(x)==4 || islogical(x));
	addParameter(p,'statistic_of_interest','ksn',@(x) ischar(x));
	addParameter(p,'basin_num',[],@(x) isnumeric(x) && isscalar(x));
	addParameter(p,'rlf_radius',2500,@(x) isnumeric(x) && isscalar(x));
	addParameter(p,'cat_mean1',[],@(x) ischar(x));
	addParameter(p,'cat_mean2',[],@(x) ischar(x));
	addParameter(p,'fit_grd_ksn',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'start_diffusivity',0.01,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'start_erodibility',1e-7,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'start_threshold_gradient',0.8,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'n_val',2,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'fit_rlf_ksn',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'fit_filtered',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'save_figure',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'out_dir',[],@(x) isfolder(x));

	parse(p,basin_table,plots,varargin{:});
	T=p.Results.basin_table;
	plts=p.Results.plots;

	uncertainty=p.Results.uncertainty;
	use_filtered=p.Results.use_filtered;
	color_by=p.Results.color_by;
	cmap=p.Results.cmap;
	xval=p.Results.xval;
	yval=p.Results.yval;
	regionOI=p.Results.define_region;
	basin_num=p.Results.basin_num;
	stOI=p.Results.statistic_of_interest;
	rr=p.Results.rlf_radius;
	cm1=p.Results.cat_mean1;
	cm2=p.Results.cat_mean2;
	fit_grd_ksn=p.Results.fit_grd_ksn;
	fit_rlf_ksn=p.Results.fit_rlf_ksn;
	fit_filtered=p.Results.fit_filtered;
	D_in=p.Results.start_diffusivity;
	K_in=p.Results.start_erodibility;
	s_in=p.Results.start_threshold_gradient;
	n_val=p.Results.n_val;
	save_figure=p.Results.save_figure;
	out_dir=p.Results.out_dir;

	if isempty(out_dir)
		out_dir=pwd;
	end
 
	if isempty(cmap)
		cmap=parula(50);
	end

	% 处理变量输入
	if ~isempty(color_by) && isnumeric(color_by)
		cval=color_by;
		color_by='color_by';
		T.color_by=cval;
	end

	if ~isempty(xval) && isnumeric(xval)
		xv=xval;
		xval='xval';
		T.xval=xv;
	end

	if ~isempty(yval) && isnumeric(yval)
		yv=yval;
		yval='yval';
		T.yval=yv;
	end	

	% 处理区域选择
	if ~isempty(regionOI) && ~islogical(regionOI)
		rIDX=T.center_x>=regionOI(1) & T.center_x<=regionOI(2) & T.center_y>=regionOI(3) & T.center_y<=regionOI(4);
		T=T(rIDX,:);
		if isempty(T)
			if isdeployed
				errordlg('提供的区域坐标已排除表中所有条目，请确认坐标正确')
			end
			error('提供的区域坐标已排除表中所有条目，请确认坐标正确');
		end
	elseif ~isempty(regionOI) && islogical(regionOI) && regionOI

		f1=figure(1);
		clf
		hold on
		scatter(T.center_x,T.center_y,20,'k','filled');
		title('在需要选择的数据周围绘制矩形');
		if ~verLessThan('matlab','9.5')
			disableDefaultInteractivity(gca);
		end		
		hold off

		rgn=getrect;

		close(f1);

		regionOI=zeros(4,1);
		regionOI(1)=rgn(1);
		regionOI(2)=rgn(1)+rgn(3);
		regionOI(3)=rgn(2);
		regionOI(4)=rgn(2)+rgn(4);

		rIDX=T.center_x>=regionOI(1) & T.center_x<=regionOI(2) & T.center_y>=regionOI(3) & T.center_y<=regionOI(4);
		T=T(rIDX,:);
	end



	% 生成图形
	switch plts
	case 'grd_ksn'

		if use_filtered
			g=T.mean_gradient_f;
			k=T.mean_ksn_f;
		else
			g=T.mean_gradient;
			k=T.mean_ksn;
		end

		if use_filtered
			if ismember('std_ksn_f',T.Properties.VariableNames) && strcmp(uncertainty,'std')
				sk=T.std_ksn_f;
				sg=T.std_gradient_f;
			elseif ismember('se_ksn_f',T.Properties.VariableNames) && strcmp(uncertainty,'se')
				sk=T.se_ksn_f;
				sg=T.se_gradient_f;
			elseif ismember('std_ksn_f',T.Properties.VariableNames) && ~ismember('se_ksn_f',T.Properties.VariableNames)
				sk=T.std_ksn_f;
				sg=T.std_gradient_f;
			elseif ismember('se_ksn_f',T.Properties.VariableNames) && ~ismember('std_ksn_f',T.Properties.VariableNames)
				sk=T.se_ksn_f;
				sg=T.se_gradient_f;
			end
		else
			if ismember('std_ksn',T.Properties.VariableNames) && strcmp(uncertainty,'std')
				sk=T.std_ksn;
				sg=T.std_gradient;
			elseif ismember('se_ksn',T.Properties.VariableNames) && strcmp(uncertainty,'se')
				sk=T.se_ksn;
				sg=T.se_gradient;
			elseif ismember('std_ksn',T.Properties.VariableNames) && ~ismember('se_ksn',T.Properties.VariableNames)
				sk=T.std_ksn;
				sg=T.std_gradient;
			elseif ismember('se_ksn',T.Properties.VariableNames) && ~ismember('std_ksn',T.Properties.VariableNames)
				sk=T.se_ksn;
				sg=T.se_gradient;
			end
		end


		if fit_grd_ksn
			[pf]=KGF(k,g,D_in,K_in,s_in,n_val);
			[mg,mk,~,~]=KGG(pf.D,pf.K,pf.Sc,pf.n,k);
		end			


		f=figure(1);
		clf
		set(f,'Units','normalized','Position',[0.05 0.1 0.5 0.5],'renderer','painters');
		hold on 

		if ~strcmp(uncertainty,'none')
			errorbar(k,g,sg,sg,sk,sk,'.k','CapSize',0);
        end
       



		if ~isempty(color_by) && isnumeric(T.(color_by))
			colormap(ksncolor(20));

			scatter(k,g,30,T.(color_by),'filled');
			cb=colorbar;
			% 去除下划线
			color_by_label=strrep(color_by,'_',' ');
			ylabel(cb,color_by_label);
		elseif ~isempty(color_by) && isa(T.(color_by),'cell')
			colormap(ksncolor(20));
            
			scatter(k,g,30,categorical(T.(color_by)),'filled');
			cb=colorbar('Ticks',[1:numel(unique(T.(color_by)))],'YTickLabel',unique(T.(color_by)));
			% 去除下划线
			color_by_label=strrep(color_by,'_',' ');
			ylabel(cb,color_by_label);
		else
			scatter(k,g,30,'k','filled');
		end

		if fit_grd_ksn
			p1=plot(mk,mg,'-r','LineWidth',2);
			disp(['河流侵蚀性 (K) = ' num2str(pf.K)]);
			disp(['扩散率 (D) = ' num2str(pf.D)]);
			disp(['阈值坡度 (Sc) = ' num2str(pf.Sc)]);
			leg_text=['河流侵蚀性 (K) = ' num2str(pf.K) newline '扩散率 (D) = ' num2str(pf.D) newline '阈值坡度 (Sc) = ' num2str(pf.Sc)];
			legend(p1,leg_text,'location','best');
		end

		xlabel('流域平均 k_{sn}');
		ylabel('流域平均坡度');

		if ~verLessThan('matlab','9.5')
			disableDefaultInteractivity(gca);
		end		
		hold off
	case 'grd_rlf'
		% 验证高程条目
		if use_filtered
			m_rlfN=['mean_rlf' num2str(rr) '_f'];
			se_rlfN=['se_rlf' num2str(rr) '_f'];
			std_rlfN=['std_rlf' num2str(rr) '_f'];
			if ismember(m_rlfN,T.Properties.VariableNames)
				r=T.(m_rlfN);
			else
				if isdeployed
					errordlg('无法识别的高程半径，请确认在运行“ProcessRiverBasins”时计算了该半径的高程')
				end
				error('无法识别的高程半径，请确认在运行“ProcessRiverBasins”时计算了该半径的高程')
			end

			g=T.mean_gradient_f;
		else
			m_rlfN=['mean_rlf' num2str(rr)];
			se_rlfN=['se_rlf' num2str(rr)];
			std_rlfN=['std_rlf' num2str(rr)];
			if ismember(m_rlfN,T.Properties.VariableNames)
				r=T.(m_rlfN);
			else
				if isdeployed
					errordlg('无法识别的高程半径，请确认在运行“ProcessRiverBasins”时计算了该半径的高程')
				end				
				error('无法识别的高程半径，请确认在运行“ProcessRiverBasins”时计算了该半径的高程')
			end

			g=T.mean_gradient;
		end

		if use_filtered
			if ismember('std_gradient_f',T.Properties.VariableNames) && strcmp(uncertainty,'std')
				sr=T.(std_rlfN);
				sg=T.std_gradient_f;
			elseif ismember('se_gradient_f',T.Properties.VariableNames) && strcmp(uncertainty,'se')
				sr=T.(se_rlfN);
				sg=T.se_gradient_f;
			elseif ismember('std_gradient_f',T.Properties.VariableNames) && ~ismember('se_gradient_f',T.Properties.VariableNames)
				sr=T.(std_rlfN);
				sg=T.std_gradient_f;
			elseif ismember('se_gradient_f',T.Properties.VariableNames) && ~ismember('std_gradient_f',T.Properties.VariableNames)
				sr=T.(se_rlfN);
				sg=T.se_gradient_f;
			end
		else
			if ismember('std_gradient',T.Properties.VariableNames) && strcmp(uncertainty,'std')
				sr=T.(std_rlfN);
				sg=T.std_gradient;
			elseif ismember('se_gradient',T.Properties.VariableNames) && strcmp(uncertainty,'se')
				sr=T.(se_rlfN);
				sg=T.se_gradient;
			elseif ismember('std_gradient',T.Properties.VariableNames) && ~ismember('se_gradient',T.Properties.VariableNames)
				sr=T.(std_rlfN);
				sg=T.std_gradient;
			elseif ismember('se_gradient',T.Properties.VariableNames) && ~ismember('std_gradient',T.Properties.VariableNames)
				sr=T.(se_rlfN);
				sg=T.se_gradient;
			end
		end

		f=figure(1);
		set(f,'Units','normalized','Position',[0.05 0.1 0.5 0.5],'renderer','painters');
		clf
		hold on 

		if ~strcmp(uncertainty,'none')
			errorbar(r,g,sg,sg,sr,sr,'.k','CapSize',0);
		end

		if ~isempty(color_by) && isnumeric(T.(color_by))
			colormap(cmap);
			scatter(r,g,30,T.(color_by),'filled');
			cb=colorbar;
			% 去除下划线
			color_by_label=strrep(color_by,'_',' ');
			ylabel(cb,color_by_label);
		elseif ~isempty(color_by) && isa(T.(color_by),'cell')
			colormap(cmap);
			scatter(r,g,30,categorical(T.(color_by)),'filled');
			cb=colorbar('Ticks',[1:numel(unique(T.(color_by)))],'YTickLabel',unique(T.(color_by)));
			% 去除下划线
			color_by_label=strrep(color_by,'_',' ');
			ylabel(cb,color_by_label);
		else
			scatter(r,g,30,'k','filled');
		end

		xlabel(['流域平均 ' num2str(rr) ' m^2 高程']);
		ylabel('流域平均坡度');

		if ~verLessThan('matlab','9.5')
			disableDefaultInteractivity(gca);
		end	

		hold off
		
	case 'rlf_ksn'
			% 验证高程条目
		if use_filtered
			m_rlfN=['mean_rlf' num2str(rr) '_f'];
			se_rlfN=['se_rlf' num2str(rr) '_f'];
			std_rlfN=['std_rlf' num2str(rr) '_f'];
			if ismember(m_rlfN,T.Properties.VariableNames)
				r=T.(m_rlfN);
			else
				if isdeployed
					errordlg('无法识别的高程半径，请确认在运行“ProcessRiverBasins”时计算了该半径的高程')
				end				
				error('无法识别的高程半径，请确认在运行“ProcessRiverBasins”时计算了该半径的高程')
			end

			k=T.mean_ksn_f;
		else
			m_rlfN=['mean_rlf' num2str(rr)];
			se_rlfN=['se_rlf' num2str(rr)];
			std_rlfN=['std_rlf' num2str(rr)];
			if ismember(m_rlfN,T.Properties.VariableNames)
				r=T.(m_rlfN);
			else
				if isdeployed
					errordlg('无法识别的高程半径，请确认在运行“ProcessRiverBasins”时计算了该半径的高程')
				end				
				error('无法识别的高程半径，请确认在运行“ProcessRiverBasins”时计算了该半径的高程')
			end

			k=T.mean_ksn;
		end

		if use_filtered
			if ismember('std_ksn_f',T.Properties.VariableNames) && strcmp(uncertainty,'std')
				sk=T.std_ksn_f;
				sr=T.(std_rlfN);
			elseif ismember('se_ksn_f',T.Properties.VariableNames) && strcmp(uncertainty,'se')
				sk=T.se_ksn_f;
				sr=T.(se_rlfN);
			elseif ismember('std_ksn_f',T.Properties.VariableNames) && ~ismember('se_ksn_f',T.Properties.VariableNames)
				sk=T.std_ksn_f;
				sr=T.(std_rlfN);
			elseif ismember('se_ksn_f',T.Properties.VariableNames) && ~ismember('std_ksn_f',T.Properties.VariableNames)
				sk=T.se_ksn_f;
				sr=T.(se_rlfN);
			end		
		else
			if ismember('std_ksn',T.Properties.VariableNames) && strcmp(uncertainty,'std')
				sk=T.std_ksn;
				sr=T.(std_rlfN);
			elseif ismember('se_ksn',T.Properties.VariableNames) && strcmp(uncertainty,'se')
				sk=T.se_ksn;
				sr=T.(se_rlfN);
			elseif ismember('std_ksn',T.Properties.VariableNames) && ~ismember('se_ksn',T.Properties.VariableNames)
				sk=T.std_ksn;
				sr=T.(std_rlfN);
			elseif ismember('se_ksn',T.Properties.VariableNames) && ~ismember('std_ksn',T.Properties.VariableNames)
				sk=T.se_ksn;
				sr=T.(se_rlfN);
			end
		end

		if fit_rlf_ksn
			ft=fittype('a*x');
			[fobj,gof]=fit(k,r,ft,'StartPoint',k\r);
			krs=coeffvalues(fobj);
			krs_unc=confint(fobj);
		end	

		f=figure(1);
		set(f,'Units','normalized','Position',[0.05 0.1 0.5 0.5],'renderer','painters');
		clf
		hold on 

		if ~strcmp(uncertainty,'none')
			errorbar(k,r,sr,sr,sk,sk,'.k','CapSize',0);
		end

		if ~isempty(color_by) && isnumeric(T.(color_by))
			colormap(cmap);
			scatter(k,r,30,T.(color_by),'filled');
			cb=colorbar;
			% 去除下划线
			color_by_label=strrep(color_by,'_',' ');
			ylabel(cb,color_by_label);
		elseif ~isempty(color_by) && isa(T.(color_by),'cell')
			colormap(cmap);
			scatter(k,r,30,categorical(T.(color_by)),'filled');
			cb=colorbar('Ticks',[1:numel(unique(T.(color_by)))],'YTickLabel',unique(T.(color_by)));
			% 去除下划线
			color_by_label=strrep(color_by,'_',' ');
			ylabel(cb,color_by_label);
		else
			scatter(k,r,30,'k','filled');
		end

		if fit_rlf_ksn
			ksn_vec=linspace(0,max(k)+50,100);
			rlf_vec=krs*ksn_vec;
			rlf_vec_pos=(krs_unc(2))*ksn_vec;
			rlf_vec_neg=(krs_unc(1))*ksn_vec;	
			p1=plot(ksn_vec,rlf_vec,'-r','LineWidth',2);
			plot(ksn_vec,rlf_vec_pos,'--r');
			plot(ksn_vec,rlf_vec_neg,'--r');
			disp(['关系斜率 = ' num2str(krs)]);
			disp(['斜率不确定性区间 = ' num2str(krs_unc(1)) ' : ' num2str(krs_unc(2))]);
			disp(['R2 = ' num2str(gof.rsquare)]);
			leg_text=['关系斜率 = ' num2str(krs) newline '斜率不确定性区间 = ' num2str(krs_unc(1)) ' : ' num2str(krs_unc(2)) newline 'R平方 = ' num2str(gof.rsquare)];
			legend(p1,leg_text,'location','best');
		end		

		ylabel(['流域平均 ' num2str(rr) ' m^2 地形起伏']);
		xlabel('流域平均 k_{sn}');

		if ~verLessThan('matlab','9.5')
			disableDefaultInteractivity(gca);
		end	

		hold off

	case 'stacked_hypsometry'
		hypsCell=T.hypsometry;
		HI=T.hyp_integral;

		numHyps=numel(hypsCell);
		normEl=zeros(100,numHyps);
		normF=zeros(100,numHyps);

		El=normEl;
		F=normF;
		for ii=1:numHyps
			normEl(:,ii)=(hypsCell{ii}(:,2)-min(hypsCell{ii}(:,2)))/(max(hypsCell{ii}(:,2))-min(hypsCell{ii}(:,2)));
			normF(:,ii)=(hypsCell{ii}(:,1))/100;

			El(:,ii)=hypsCell{ii}(:,2);
			F(:,ii)=hypsCell{ii}(:,1);
		end

		histIm=zeros(100,100);
		jj=100:-1:1;
		for ii=1:100
			[N,~]=histcounts(normF(ii,:),linspace(0,1,101));
			histIm(jj(ii),:)=log10(N);
		end

		if ~isempty(color_by) && isnumeric(T.(color_by))
			col_val=T.(color_by);
			f(1)=figure(1);
			set(f(1),'Units','normalized','Position',[0.05 0.5 0.4 0.4],'renderer','painters');
			clf
			hold on
			colormap(ksncolor);
			cm=colormap;
			num_col=size(cm,1);
			[cix,ed]=discretize(col_val,num_col);
			for ii=1:num_col
				idx=cix==ii;
				plot(normF(:,idx),normEl(:,idx),'Color',cm(ii,:));
			end
			caxis([min(ed) max(ed)]);
			cb=colorbar;
			% 去除下划线
			color_by_label=strrep(color_by,'_',' ');
			ylabel(cb,color_by_label);
			axis equal
			xlabel('归一化面积');
			ylabel('归一化高程');
			xlim([0 1]);
			ylim([0 1]);

			if ~verLessThan('matlab','9.5')
				disableDefaultInteractivity(gca);
			end	
			hold off			
		else
			f(1)=figure(1);
			set(f(1),'Units','normalized','Position',[0.05 0.5 0.4 0.4],'renderer','painters');
			clf
			hold on
			plot(normF,normEl,'-k');
			axis equal
			xlabel('归一化面积');
			ylabel('归一化高程');
			xlim([0 1]);
			ylim([0 1]);
			if ~verLessThan('matlab','9.5')
				disableDefaultInteractivity(gca);
			end	
			hold off
		end


		if ~isempty(color_by) && isnumeric(T.(color_by))
			col_val=T.(color_by);
			f(2)=figure(2);
			set(f(2),'Units','normalized','Position',[0.05 0.05 0.4 0.4],'renderer','painters');
			clf
			hold on
			colormap(cmap);
			cm=colormap;
			num_col=size(cm,1);
			[cix,ed]=discretize(col_val,num_col);
			for ii=1:num_col
				idx=cix==ii;
				plot(F(:,idx),El(:,idx),'Color',cm(ii,:));
			end
			caxis([min(ed) max(ed)]);			
			cb=colorbar;
			% 去除下划线
			color_by_label=strrep(color_by,'_',' ');
			ylabel(cb,color_by_label);
			axis square
			xlabel('面积百分比');
			ylabel('高程');

			if ~verLessThan('matlab','9.5')
				disableDefaultInteractivity(gca);
			end	
			hold off
		else
			f(2)=figure(2);
			set(f(2),'Units','normalized','Position',[0.05 0.05 0.4 0.4],'renderer','painters');
			clf
			hold on
			plot(F,El,'-k');
			axis square
			xlabel('面积百分比');
			ylabel('高程');

			if ~verLessThan('matlab','9.5')
				disableDefaultInteractivity(gca);
			end	
			hold off
		end

		num_bins=20;
		[hix,hed]=discretize(HI,linspace(0,1,num_bins+1));

		f(3)=figure(3);
		set(f(3),'Units','normalized','Position',[0.5 0.5 0.4 0.4],'renderer','painters');
		clf 

		for ii=1:num_bins
			sbplt(ii)=subplot(4,5,ii);
			hold on
			idx=hix==ii;
			perc(ii,1)=round((nnz(idx)/numel(idx))*100,1);
			nF=normF(:,idx);
			nE=normEl(:,idx);
			mF{ii,1}=mean(nF,2);
			mE{ii,1}=mean(nE,2);
			plot(nF,nE,'LineWidth',0.5,'Color',[0.4 0.4 0.4]);
			plot(mF{ii,1},mE{ii,1},'-r','LineWidth',2);
			if ii<=num_bins/2
				text(0.75,0.9,[num2str(perc(ii,1)) '%']);
			else
				text(0.1,0.1,[num2str(perc(ii,1)) '%']);
			end
			axis equal
			title(['HI ' num2str(hed(ii)) ' 至 ' num2str(hed(ii+1))])
			xlabel('归一化面积');
			ylabel('归一化高程');
			xlim([0 1]);
			ylim([0 1]);

			if ~verLessThan('matlab','9.5')
				disableDefaultInteractivity(sbplt(ii));
			end				
			hold off
		end				

		f(4)=figure(4);
		set(f(4),'Units','normalized','Position',[0.5 0.05 0.4 0.4]);
		clf 
		X=linspace(0,1,100);
		Y=linspace(0,1,100);
		hold on
		colormap(jet);
		imagesc(X,Y,histIm);
		axis equal
		xlabel('归一化面积');
		ylabel('归一化高程');
		xlim([0 1]);
		ylim([0 1]);

		if ~verLessThan('matlab','9.5')
			disableDefaultInteractivity(gca);
		end	
		hold off

		f(5)=figure(5);
		set(f(5),'Units','normalized','Position',[0.25 0.25 0.4 0.4],'renderer','painters');
		clf 
		hold on
		colormap(jet);
		idx=perc>0;
		pf=perc(idx);
		jc=jet(100);
		for ii=1:num_bins
			ix=round((perc(ii,1)/max(pf))*100);
			if ix==0
				cl=[1 1 1];
			else
				cl=jc(ix,:);
			end
			plot(mF{ii,1},mE{ii,1},'Color',cl,'LineWidth',2);
		end
		caxis([min(pf) max(pf)]);
		cb=colorbar;
		ylabel(cb,'流域百分比')
		axis equal
		xlabel('归一化面积');
		ylabel('归一化高程');
		xlim([0 1]);
		ylim([0 1]);
		if ~verLessThan('matlab','9.5')
			disableDefaultInteractivity(gca);
		end			
		hold off

	case 'scatterplot_matrix'
		% 查找包含均值的变量
		if use_filtered
			VN=T.Properties.VariableNames;
			ix=regexp(VN,regexptranslate('wildcard','mean_*_f'));
			ix=cellfun(@any,ix);
			VNoi=VN(ix);
		else
			VN=T.Properties.VariableNames;
			ix=regexp(VN,regexptranslate('wildcard','mean_*'));
			ix=cellfun(@any,ix);
			VNoi=VN(ix);
		end

		% 添加方位角计算（如果存在）
		aix=regexp(VN,regexptranslate('wildcard','dist_along*'));
		aix=cellfun(@any,aix);
		if ~isempty(aix)
			VNoi=horzcat(VNoi,VN(aix));
		end

		num_sc=numel(VNoi);

		f=figure(1);
		set(f,'Units','normalized','Position',[0.05 0.5 0.9 0.9],'renderer','painters');		
		clf
		pos=1;

		for ii=1:num_sc
			yval=T.(VNoi{ii});
			yname=VNoi{ii};
			yname=strrep(yname,'_',' ');
			for jj=1:num_sc
				subplot(num_sc,num_sc,pos)
				hold on
				if ii==jj
					histogram(yval,25,'FaceColor','k');
					if jj==num_sc && ii==num_sc
						xlabel(yname);
					end
					axis square
				else
					xval=T.(VNoi{jj});
					xname=VNoi{jj};
					xname=strrep(xname,'_',' ');
					scatter(xval,yval,5,'k','filled');
					warning off
					f=fit(xval,yval,'poly2');
					warning on
					xx=linspace(min(xval),max(xval),50);
					yy=f(xx);
					plot(xx,yy,'-r');
					if ii==num_sc
						xlabel(xname);
					end

					if jj==1
						ylabel(yname);
					end

					axis square
				end

				if ~verLessThan('matlab','9.5')
					disableDefaultInteractivity(gca);
				end	
				hold off
				pos=pos+1;
			end
		end

	case 'compare_filtered'

		VN=T.Properties.VariableNames;
		ix=regexp(VN,regexptranslate('wildcard','mean_*_f'));
		ix=cellfun(@any,ix);
		VNoi=VN(ix);
		if isempty(VNoi)
			if isdeployed
				errordlg('提供的表格中未找到筛选后的值')
			end
			error('提供的表格中未找到筛选后的值');
		end

		if fit_filtered
			ft=fittype('a*x');
		end		

		num_filt=numel(VNoi);
		for ii=1:num_filt
			fN=VNoi{ii};
			N=strrep(fN,'_f','');

			% 解析参数名称
			if strcmp(N,'mean_ksn')
				t='平均河道陡度';
			elseif strcmp(N,'mean_gradient')
				t='平均坡度';
			elseif regexp(N,regexptranslate('wildcard','mean_rlf*'))
				t=['平均' strrep(N,'mean_rlf','') 'm²地形起伏'];
			else
				t=['平均' strrep(N,'mean_','')];
			end

			max_val=max([max(T.(fN)) max(T.(N))]);
			max_vec=[0 max_val];

			slp=round(T.(N)./T.(fN),1);

			idx1=slp==1;
			idx2=slp<1;
			idx3=slp>1;

			if fit_filtered
				idx=isnan(T.(fN)) | isnan(T.(N));
				x=double(T.(fN)(~idx));
				y=double(T.(N)(~idx));
				[fobj,gof]=fit(x,y,ft,'StartPoint',x\y);
				cf=coeffvalues(fobj);
				cf_unc=confint(fobj);
			end	

			f(ii)=figure(ii);
			set(gcf,'Units','normalized','Position',[0.05 0.1 0.5 0.5],'renderer','painters');
			clf
			hold on
			plot(max_vec,max_vec,'-k');
			sp(1)=scatter(T.(fN)(idx1),T.(N)(idx1),30,'k','filled');
			sp(2)=scatter(T.(fN)(idx2),T.(N)(idx2),30,'r','filled');
			sp(3)=scatter(T.(fN)(idx3),T.(N)(idx3),30,'b','filled');

			if fit_filtered
				fvec=linspace(0,max(T.(fN)),100);
				vec=cf*fvec;
				vec_pos=(cf_unc(2))*fvec;
				vec_neg=(cf_unc(1))*fvec;	
				pl=plot(fvec,vec,'-g','LineWidth',2);
				plot(fvec,vec_pos,'--g');
				plot(fvec,vec_neg,'--g');
				disp(['拟合结果 ' t ':']);
				disp(['	斜率 = ' num2str(cf)]);
				disp(['	斜率不确定区间 = ' num2str(cf_unc(1)) ' : ' num2str(cf_unc(2))]);
				disp(['	R2 = ' num2str(gof.rsquare)]);
				leg_text={'筛选后 = 未筛选','筛选后 > 未筛选','筛选后 < 未筛选',['最佳拟合:' newline '	斜率 = ' num2str(cf) newline '	不确定区间 = ' num2str(cf_unc(1)) ' : ' num2str(cf_unc(2)) newline '	R平方 = ' num2str(gof.rsquare)]};
				legend([sp pl],leg_text,'location','northwest')
			else
				legend(sp,{'筛选后 = 未筛选','筛选后 > 未筛选','筛选后 < 未筛选'},'location','northwest')			
			end

			xlabel('筛选后均值');
			ylabel('原始均值');
			title(t);
			if ~verLessThan('matlab','9.5')
				disableDefaultInteractivity(gca);
			end	
			hold off
		end

	case 'category_mean_hist'

		if isempty(cm1)
			if isdeployed
				errordlg('对于绘图选项 "category_mean_hist"，必须为 "cat_mean1" 提供输入')
			end
			error('对于绘图选项 "category_mean_hist"，必须为 "cat_mean1" 提供输入');
		elseif strcmp(cm1,'ksn')
			VN=T.Properties.VariableNames;
			ix=regexp(VN,regexptranslate('wildcard','mksn*'));
			ix=cellfun(@any,ix);
			VNoi=VN(ix);
			Cat_Names=cellfun(@(x) strrep(x,'mksn_',''),VNoi,'UniformOutput',false);
			Main_Title='平均k_{sn}（';
		elseif strcmp(cm1,'gradient')
			VN=T.Properties.VariableNames;
			ix=regexp(VN,regexptranslate('wildcard','mgrad*'));
			ix=cellfun(@any,ix);
			VNoi=VN(ix);
			Cat_Names=cellfun(@(x) strrep(x,'mgrad_',''),VNoi,'UniformOutput',false);
			Main_Title='平均坡度（';
		elseif strcmp(cm1,'rlf')
			VN=T.Properties.VariableNames;
			srch_strng=['mr' num2str(rr) '*'];
			ix=regexp(VN,regexptranslate('wildcard',srch_strng));
			ix=cellfun(@any,ix);
			VNoi=VN(ix);
			if isempty(VNoi)
				if isdeployed
					errordlg('参数 "cat_mean1" 无法识别，请检查地形起伏半径')
				end
				error('参数 "cat_mean1" 无法识别，请检查地形起伏半径');
			end
			srch_strng=['mr' num2str(rr) '_'];
			Cat_Names=cellfun(@(x) strrep(x,srch_strng,''),VNoi,'UniformOutput',false);	
			Main_Title=['平均' num2str(rr) 'm²地形起伏（'];	
		else
			VN=T.Properties.VariableNames;
			srch_strng=['m' cm1 '*'];
			ix=regexp(VN,regexptranslate('wildcard',srch_strng));
			ix=cellfun(@any,ix);
			VNoi=VN(ix);
			if isempty(VNoi)
				if isdeployed
					errordlg('参数 "cat_mean1" 无法识别')
				end
				error('参数 "cat_mean1" 无法识别');
			end	
			srch_strng=['m' cm1 '_'];			
			Cat_Names=cellfun(@(x) strrep(x,srch_strng,''),VNoi,'UniformOutput',false);	
			Main_Title=['平均' cm1 '（'];
		end	

		for ii=1:numel(VNoi)
			vals=T.(VNoi{ii});
			vals(isnan(vals))=[];

			if ~isempty(vals)
				val_list{ii,1}=vals;
			end
		end
		val_list=vertcat(val_list{:});
		[~,edges]=discretize(val_list,100);

		for ii=1:numel(VNoi)
			vals=T.(VNoi{ii});
			vals(isnan(vals))=[];

			if ~isempty(vals)
				f(ii)=figure(ii);
				set(gcf,'Units','normalized','Position',[0.05 0.1 0.5 0.5],'renderer','painters');
				clf
				hold on
				histogram(vals,edges);
				title([Main_Title Cat_Names{ii} '）']);
				if ~verLessThan('matlab','9.5')
					disableDefaultInteractivity(gca);
				end	
				hold off
			end
		end

	case 'category_mean_compare'
    disp(cm1);
		if isempty(cm1)
			if isdeployed
				errordlg('对于绘图选项 "category_mean_compare"，必须为 "cat_mean1" 提供输入')
			end
			error('对于绘图选项 "category_mean_compare"，必须为 "cat_mean1" 提供输入');
		elseif strcmp(cm1,'ksn')
			VN=T.Properties.VariableNames;
			ix=regexp(VN,regexptranslate('wildcard','mksn*'));
			ix=cellfun(@any,ix);
			VNoi1=VN(ix);
			Cat_Names=cellfun(@(x) strrep(x,'mksn_',''),VNoi1,'UniformOutput',false);
			axis1='平均k_{sn}（';
		elseif strcmp(cm1,'gradient')
			VN=T.Properties.VariableNames;
			ix=regexp(VN,regexptranslate('wildcard','mgrad*'));
			ix=cellfun(@any,ix);
			VNoi1=VN(ix);
			Cat_Names=cellfun(@(x) strrep(x,'mgrad_',''),VNoi1,'UniformOutput',false);
			axis1='平均坡度（';
		elseif strcmp(cm1,'rlf')
			VN=T.Properties.VariableNames;
			srch_strng=['mr' num2str(rr) '*'];
			ix=regexp(VN,regexptranslate('wildcard',srch_strng));
			ix=cellfun(@any,ix);
			VNoi1=VN(ix);
			if isempty(VNoi1)
				if isdeployed
					errordlg('参数 "cat_mean1" 无法识别，请检查地形起伏半径')
				end
				error('参数 "cat_mean1" 无法识别，请检查地形起伏半径');
			end
			srch_strng=['mr' num2str(rr) '_'];
			Cat_Names=cellfun(@(x) strrep(x,srch_strng,''),VNoi1,'UniformOutput',false);	
			axis1=['平均' num2str(rr) 'm²地形起伏（'];	
		else
						VN=T.Properties.VariableNames;
			srch_strng=['m' cm1 '*'];
			ix=regexp(VN,regexptranslate('wildcard',srch_strng));
			ix=cellfun(@any,ix);
			VNoi1=VN(ix);
			if isempty(VNoi1)
				if isdeployed
					errordlg('参数 "cat_mean1" 无法识别')
				end
				error('参数 "cat_mean1" 无法识别');
			end	
			srch_strng=['m' cm1 '_'];			
			Cat_Names=cellfun(@(x) strrep(x,srch_strng,''),VNoi1,'UniformOutput',false);	
			axis1=['平均' cm1 '（'];
		end			

		if isempty(cm2)
			if isdeployed
				errordlg('对于绘图选项 "category_mean_compare"，必须为 "cat_mean2" 提供输入')
			end
			error('对于绘图选项 "category_mean_compare"，必须为 "cat_mean2" 提供输入');
		elseif strcmp(cm2,'ksn')
			VN=T.Properties.VariableNames;
			ix=regexp(VN,regexptranslate('wildcard','mksn*'));
			ix=cellfun(@any,ix);
			VNoi2=VN(ix);
			Cat_Names=cellfun(@(x) strrep(x,'mksn_',''),VNoi2,'UniformOutput',false);
			axis2='平均k_{sn}）';
		elseif strcmp(cm2,'gradient')
			VN=T.Properties.VariableNames;
			ix=regexp(VN,regexptranslate('wildcard','mgrad*'));
			ix=cellfun(@any,ix);
			VNoi2=VN(ix);
			Cat_Names=cellfun(@(x) strrep(x,'mgrad_',''),VNoi2,'UniformOutput',false);
			axis2='平均坡度）';
		elseif strcmp(cm2,'rlf')
			VN=T.Properties.VariableNames;
			srch_strng=['mr' num2str(rr) '*'];
			ix=regexp(VN,regexptranslate('wildcard',srch_strng));
			ix=cellfun(@any,ix);
			VNoi2=VN(ix);
			if isempty(VNoi2)
				if isdeployed
					errordlg('参数 "cat_mean2" 无法识别，请检查地形起伏半径')
				end
				error('参数 "cat_mean2" 无法识别，请检查地形起伏半径');
			end
			srch_strng=['mr' num2str(rr) '_'];
			Cat_Names=cellfun(@(x) strrep(x,srch_strng,''),VNoi2,'UniformOutput',false);	
			axis2=[num2str(rr) 'm²地形起伏）'];	
		else
			VN=T.Properties.VariableNames;
			srch_strng=['m' cm2 '*'];
			ix=regexp(VN,regexptranslate('wildcard',srch_strng));
			ix=cellfun(@any,ix);
			VNoi2=VN(ix);
			if isempty(VNoi2)
				if isdeployed
					errordlg('参数 "cat_mean2" 无法识别')
				end
				error('参数 "cat_mean2" 无法识别');
			end	
			srch_strng=['m' cm2 '_'];			
			Cat_Names=cellfun(@(x) strrep(x,srch_strng,''),VNoi2,'UniformOutput',false);	
			axis2=['平均' cm2 '）'];
		end	

		rng_v1=zeros(numel(VNoi1),2);
		rng_v2=zeros(numel(VNoi1),2);	

		for ii=1:numel(VNoi1)
			vals1=T.(VNoi1{ii});
			vals2=T.(VNoi2{ii});

			idx=~isnan(vals1) & ~isnan(vals2);
			vals1=vals1(idx);
			vals2=vals2(idx);

			if ~isempty(vals1) & ~isempty(vals2)
				rng_v1(ii,1)=min(vals1,[],'omitnan');
				rng_v2(ii,1)=min(vals2,[],'omitnan');
				rng_v1(ii,2)=max(vals1,[],'omitnan');
				rng_v2(ii,2)=max(vals2,[],'omitnan');				
			end
		end

		rng_v1=[min(rng_v1(:,1),[],'omitnan') max(rng_v1(:,2),[],'omitnan')];
		rng_v2=[min(rng_v2(:,1),[],'omitnan') max(rng_v2(:,2),[],'omitnan')];

		for ii=1:numel(VNoi1)
			vals1=T.(VNoi1{ii});
			vals2=T.(VNoi2{ii});

			idx=~isnan(vals1) & ~isnan(vals2);
			vals1=vals1(idx);
			vals2=vals2(idx);

			if ~isempty(vals1) & ~isempty(vals2)
				f(ii)=figure(ii);
				set(gcf,'Units','normalized','Position',[0.05 0.1 0.5 0.5],'renderer','painters');
				clf
				hold on
				scatter(vals1,vals2,30,'k','filled');
				xlim(rng_v1);
				ylim(rng_v2);
				xlabel([axis1 Cat_Names{ii}])
				ylabel([axis2 Cat_Names{ii}])
				if ~verLessThan('matlab','9.5')
					disableDefaultInteractivity(gca);
				end	
				hold off
			end
		end

	case 'compare_mean_and_dist'

		% 验证输入参数
		if strcmp(stOI,'ksn')
			val=T.mean_ksn;
			m_valN='mean_ksn';
			load_flag=1;
			title_str='河道陡度';
		elseif strcmp(stOI,'gradient')
			val=T.mean_gradient;
			m_valN='mean_gradient';
			load_flag=2;
			title_str='坡度';
		elseif strcmp(stOI,'elevation')
			val=T.mean_el;
			m_valN='mean_el';
			load_flag=3;
			title_str='高程';
		elseif strcmp(stOI,'relief')
			m_valN=['mean_rlf' num2str(rr)];
			if ismember(m_valN,T.Properties.VariableNames)
				val=T.(m_valN);	
				load_flag=4;
			else
				if isdeployed
					errordlg('地形起伏半径无法识别，请确认在运行"ProcessRiverBasins"时计算了该半径')
				end
				error('地形起伏半径无法识别，请确认在运行"ProcessRiverBasins"时计算了该半径')
			end
			title_str=[num2str(rr) 'm²地形起伏'];
		else
			m_valN=['mean_' stOI];
			if ismember(m_valN,T.Properties.VariableNames)
				val=T.(m_valN);
				load_flag=5;
			else
				if isdeployed
					errordlg(['提供的表格中未找到"mean_' stOI '"列，请确认网格名称'])
				end
				error(['提供的表格中未找到"mean_' stOI '"列，请确认网格名称'])
			end
			title_str=[upper(stOI(1)) stOI(2:end)];
		end			

		% 确定绘图类型
		if isempty(basin_num)
			out=cell(numel(val),1);

			w1=waitbar(0,'正在编译统计信息');
			for ii=1:numel(val)
				if load_flag==1
					load(T.file_path{ii,1},'MSNc');
					out{ii,1}=[MSNc.ksn]';
				elseif load_flag==2
					load(T.file_path{ii,1},'Goc');
					g=Goc.Z(:);
					g(isnan(g))=[];
					out{ii,1}=g;
				elseif load_flag==3
					load(T.file_path{ii,1},'DEMoc');
					d=DEMoc.Z(:);
					d(isnan(d))=[];
					out{ii,1}=d;
				elseif load_flag==4
					load(T.file_path{ii,1},'rlf');
					ix=find(cell2mat(rlf(:,2))==rr);
					r=rlf{ix,1}.Z(:);		
					r(isnan(r))=[];
					out{ii,1}=r;	
				elseif load_flag==5
					load(T.file_path{ii,1},'AGc');
					ix=find(strcmp(AGc(:,2),stOI));
					a=AGc{ix,1}.Z(:);
					a(isnan(a))=[];
					out{ii,1}=a;
				end				
				waitbar(ii/numel(val));
			end
			close(w1);

			out=vertcat(out{:});

			f=figure(1);
			set(f,'Units','normalized','Position',[0.05 0.1 0.5 0.5],'renderer','painters');
			clf

			% 数据筛选
			[N,ed]=histcounts(out,100);
			N=N/max(N);
			bin_list=[1:100]';
			[IX]=discretize(out,ed);
			idx=N>0.01;
			keep_bins=bin_list(idx);
			idx=ismember(IX,keep_bins);
			out_f=out(idx);

			sbplt1=subplot(1,2,1);
			hold on 
			[nF,edgesF]=histcounts(out_f,100);
			[nV,~]=histcounts(val,edgesF);
			nF=nF/max(nF);
			nV=nV/max(nV);
			histogram('BinEdges',edgesF,'BinCounts',nF);
			histogram('BinEdges',edgesF,'BinCounts',nV);
			title('移除<1%异常值');
			legend('全流域','流域均值','location','best');
			xlabel(title_str);
			ylabel('归一化计数');
			if ~verLessThan('matlab','9.5')
				disableDefaultInteractivity(sbplt1);
			end	
			hold off

			sbplt2=subplot(1,2,2);
			hold on 
			[n,edges]=histcounts(out,100);
			[nV,~]=histcounts(val,edges);
			n=n/max(n);
			nV=nV/max(nV);
			histogram('BinEdges',edges,'BinCounts',n);
			histogram('BinEdges',edges,'BinCounts',nV);
			title('完整数据集')
			legend('全流域','流域均值','location','best');			
			xlabel(title_str);	
			if ~verLessThan('matlab','9.5')
				disableDefaultInteractivity(sbplt2);
			end			
			hold off

		else
			disp(basin_num);
			ii=find(T.ID==basin_num);
			if isempty(ii)
				if isdeployed
					errordlg('提供的表格中未找到该流域编号')
				end
				error('提供的表格中未找到该流域编号')
			end

			if load_flag==1
				load(T.file_path{ii,1},'MSNc');
				out=[MSNc.ksn]';
			elseif load_flag==2
				load(T.file_path{ii,1},'Goc');
				g=Goc.Z(:);
				g(isnan(g))=[];
				out=g;
			elseif load_flag==3
				load(T.file_path{ii,1},'DEMoc');
				d=DEMoc.Z(:);
				d(isnan(d))=[];
				out=d;
			elseif load_flag==4
				load(T.file_path{ii,1},'rlf');
				ix=find(cell2mat(rlf(:,2))==rr);
				r=rlf{ix,1}.Z(:);		
				r(isnan(r))=[];
				out=r;	
			elseif load_flag==5
				load(T.file_path{ii,1},'AGc');
				ix=find(strcmp(AGc(:,2),stOI));
				a=AGc{ix,1}.Z(:);
				a(isnan(a))=[];
				out=a;
			end	

			[N,~]=histcounts(out,100);

			f=figure(1);
			set(f,'Units','normalized','Position',[0.05 0.1 0.5 0.5],'renderer','painters');
			clf

			hold on 
			histogram(out,100);
			plot([T.(m_valN)(ii,1),T.(m_valN)(ii,1)],[0,max(N)],'-k','LineWidth',2);
			title(title_str);
			xlabel(['流域 ' num2str(basin_num) ' 的数值']);
			if ~verLessThan('matlab','9.5')
				disableDefaultInteractivity(gca);
			end	
			hold off
		end		

	case 'xy'

		if isempty(xval) | isempty(yval)
			if isdeployed
				errordlg('使用"xy"绘图必须提供"xval"和"yval"参数')
			end
			error('使用"xy"绘图必须提供"xval"和"yval"参数');
		end

		VN=T.Properties.VariableNames;

		if ~any(strcmp(VN,xval))
			if isdeployed
				errordlg('输入的"xval"参数在表格中不存在')
			end
			error('输入的"xval"参数在表格中不存在');
		elseif ~any(strcmp(VN,yval))
			if isdeployed
				errordlg('输入的"yval"参数在表格中不存在')
			end
			error('输入的"yval"参数在表格中不存在');
		end

		x=T.(xval);
		y=T.(yval);

		if ~isnumeric(x) || ~isnumeric(y)
			if isdeployed
				errordlg('x和y值必须为数值型数据')
			end
			error('x和y值必须为数值型数据')
		end

		sxN=strrep(xval,'mean_',''); 
		syN=strrep(yval,'mean_','');

		if ismember(['std_' sxN],T.Properties.VariableNames) && strcmp(uncertainty,'std')
			sx=T.(['std_' sxN]);
		elseif ismember(['se_' sxN],T.Properties.VariableNames) && strcmp(uncertainty,'se')
			sx=T.(['se_' sxN]);
		elseif ismember(['std_' sxN],T.Properties.VariableNames) && ~ismember(['se_' sxN],T.Properties.VariableNames)
			sx=T.(['std_' sxN]);
		elseif ismember(['se_' sxN],T.Properties.VariableNames) && ~ismember(['std_' sxN],T.Properties.VariableNames)
			sx=T.(['se_' sxN]);
		else
			sx=[];
		end

		if ismember(['std_' syN],T.Properties.VariableNames) && strcmp(uncertainty,'std')
			sy=T.(['std_' syN]);
		elseif ismember(['se_' syN],T.Properties.VariableNames) && strcmp(uncertainty,'se')
			sy=T.(['se_' syN]);
		elseif ismember(['std_' syN],T.Properties.VariableNames) && ~ismember(['se_' syN],T.Properties.VariableNames)
			sy=T.(['std_' syN]);
		elseif ismember(['se_' syN],T.Properties.VariableNames) && ~ismember(['std_' syN],T.Properties.VariableNames)
			sy=T.(['se_' syN]);
		else
			sy=[];
		end

		f=figure(1);
		set(f,'Units','normalized','Position',[0.05 0.1 0.5 0.5],'renderer','painters');
		clf
		hold on 

		if ~strcmp(uncertainty,'none') && ~isempty(sx) && isempty(sy)
			errorbar(x,y,sx,'horizontal','.k','CapSize',0);
		elseif ~strcmp(uncertainty,'none') && isempty(sy) && ~isempty(sx)
			errorbar(x,y,sy,'vertical','.k','CapSize',0);
		elseif ~strcmp(uncertainty,'none') && ~isempty(sx) && ~isempty(sy)
			errorbar(x,y,sy,sy,sx,sx,'.k','CapSize',0);
		end

		if ~isempty(color_by) && isnumeric(T.(color_by))
			colormap(cmap);
			scatter(x,y,30,T.(color_by),'filled');
			cb=colorbar;
			color_by_label=strrep(color_by,'_',' ');
			ylabel(cb,color_by_label);
		elseif ~isempty(color_by) && isa(T.(color_by),'cell')
			colormap(cmap);
			scatter(x,y,30,categorical(T.(color_by)),'filled');
			cb=colorbar('Ticks',[1:numel(unique(T.(color_by)))],'YTickLabel',unique(T.(color_by)));
			color_by_label=strrep(color_by,'_',' ');
			ylabel(cb,color_by_label);
		else
			scatter(x,y,30,'k','filled');
		end

		xlabel(['平均' sxN])
		ylabel(['平均' syN]);

		if ~verLessThan('matlab','9.5')
			disableDefaultInteractivity(gca);
		end	
		hold off
	end


	if save_figure

		num_figs=numel(f);

		if num_figs>1
			for ii=1:num_figs
				orient(f(ii),'landscape');
				print(f(ii),'-dpdf','-fillpage',fullfile(out_dir,['Figure_' num2str(ii) '.pdf']));
			end
		else
			current=gcf;
			orient 'landscape';
			print(current,'-dpdf','-fillpage',fullfile(out_dir,['Figure_1.pdf']));
		end
	end
end

function [param_fit]=KGF(x0,y0,D_in,K_in,s_in,n)
	% n在此版本中不是自由参数

    model=@RC;
    stp=[D_in,K_in,s_in];
    lb=[0,0,0];
    ub=[1,1,2];
	opt=optimset('Display','iter'); % 保留迭代显示

	disp('正在进行拟合...')
    est=fmincon(model,stp,[],[],[],[],lb,ub);

    param_fit.D=est(1);
    param_fit.K=est(2);
    param_fit.Sc=est(3);
    param_fit.n=n;

    function [lad]=RC(params)
    	D=params(1);
    	K=params(2);
    	s=params(3);

    	er=ErKsn(x0,K,n);
    	er=er.';

		for ii=1:numel(er)
			l(1,ii)=CharLh(K,D,er(ii),n,s);
		end

		[mean_gradient]=NonLinearHillslope(er,s,D,l);
		res=(y0.')-mean_gradient;
		lad=sum(abs(res));
	end
end

function [me_grd,me_ksn,er,l]=KGG(D,K,s,n,ksn)
	% 定义侵蚀率范围
	ksn_vec=linspace(min(ksn),max(ksn)+50,100);
	er=ErKsn(ksn_vec,K,n);

	for ii=1:numel(er)
		l(1,ii)=CharLh(K,D,er(ii),n,s);
	end

	[me_grd]=NonLinearHillslope(er,s,D,l);
	[me_ksn]=KsnEr(er,K,n);
end

function [lh]=CharLh(K,D,er,n,s)
	b=2;
	m=n/2;

	fun=@(lh) (((D*(s^2))/(b*er*lh))*(sqrt(1+((b*er*lh)/(D*s))^2)-1))-((er/(2*K*(lh^(2*m))))^(1/n));

	lh=fzero(fun,100);
end

function [ksn]=KsnEr(er,K,n)
	ksn=(er./K).^(1/n);
end

function [er]=ErKsn(ksn,K,n)
	er=K.*(ksn.^n);
end

function [mean_gradient_r]=NonLinearHillslope(er,s,D,l)
	b=2;

	[z_0]=el_func_x(er,s,b,D,0);
	[z_l]=el_func_x(er,s,b,D,l);

	mean_gradient_r=(z_0-z_l)./l;
end

function [z]=el_func_x(er,s,b,D,xx)
	out_term=(-(s^2)./(2*b.*er));
	rep_term=(2*b.*er)./s;
	first_inner=sqrt((D^2)+(rep_term.*xx).^2);
	second_inner=D * log((first_inner+D)./rep_term);
	z=out_term.*(first_inner-second_inner);
end