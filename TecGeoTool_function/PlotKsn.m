function PlotKsn(DEM,FD,ksn,varargin)
	%
% 用法：
%	PlotKsn(DEM, FD, ksn);
%	PlotKsn(DEM, FD, ksn, bnd_list);
%	PlotKsn(DEM, FD, ksn, 'knicks.shp');
%
% 描述：
% 	该函数用于绘制规范化的河道陡峭度图，并以高程为背景着色的阴影图。
%
% 必需的输入：
%	DEM - 作为 GRIDobj 的数字高程，用于生成提供的 ksn 数据。
%	FD - 流向作为 FLOWobj，用于生成提供的 ksn 数据。
%	ksn - ksn 数据，可以是 shapefile（来自 KsnProfiler、ProcessRiverBasins 或 KsnChiBatch 的输出），
%		ASCII 文件（作为来自 KsnChiBatch 的连续 Ksn 值输出），
%		GRIDobj（作为来自 KsnChiBatch 的连续 ksn 输出），或 mapstructure 
%		（来自 ProcessRiverBasins 或 KsnChiBatch 的输出）。
% 
% 可选输入：
%	knicks [] - 河道拐点位置，可以是数组（来自 FindBasinKnicks 或 KsnProfiler 的 'bnd_list' 输出），
%		或 shapefile（由 FindBasinKnicks 或 KsnProfiler 输出）。
%	ksn_lim [] - 1 x n 向量，设置 ksn 颜色映射的最小值和最大值。若留空，
%		默认使用0到数据集最大值范围。
%
% 示例：
%	PlotKsn(DEMoc, FDc, MSNc); % 从 ProcessRiverBasins 绘制流域的 ksn 图
%	PlotKsn(DEM, FD, 'ksn.shp'); % 从 shapefile 绘制 ksn 图
%	PlotKsn(DEM, FD, 'ksn.shp', 'knicks', 'ksn_knicks.shp'); % 叠加 KsnProfiler 的河道拐点
%	PlotKsn(DEMoc, FDc, MSNc, 'knicks', KnickPoints); % 叠加 FindBasinKnicks 的拐点输出
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数编写者：Adam M. Forte - 更新日期：06/18/18 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if ischar(ksn)
		[~,~,ext]=fileparts(ksn);
	else
		ext=' ';
	end

    % 解析输入参数
    p = inputParser;
    p.FunctionName = 'PlotKsn';
    addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
    addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
    addRequired(p,'ksn',@(x) isstruct(x) || strcmpi(ext,'.shp') || strcmpi(ext,'.txt') || isa(x,'GRIDobj'));

    addParameter(p,'knicks',[],@(x) isnumeric(x) || regexp(x,regexptranslate('wildcard','*.shp')) || istable(x));
    addParameter(p,'ksn_lim',[],@(x) isnumeric(x) && numel(x)==2);

    parse(p,DEM,FD,ksn,varargin{:});
    DEM=p.Results.DEM;
    FD=p.Results.FD;
    ksn=p.Results.ksn;

    knks=p.Results.knicks;
    ksn_lim=p.Results.ksn_lim;

%%
	if ischar(ksn) & logical(regexp(ksn,regexptranslate('wildcard','*.shp')))
		ksn=shaperead(ksn);
		grid_flag=false;
	elseif ischar(ksn) && logical(regexp(ksn,regexptranslate('wildcard','*.txt')))
		ksn=GRIDobj(ksn);
		if ~validatealignment(DEM,ksn)
			ksn=resample(ksn,DEM);
		end
		grid_flag=true;
	elseif isstruct(ksn)
		ksn=ksn;
		grid_flag=false;
	elseif isa(ksn,'GRIDobj')
		if ~validatealignment(DEM,ksn)
			ksn=resample(ksn,DEM);
		end
		grid_flag=true;
	else
		if isdeployed
			errordlg('输入的"ksn"无法识别为shapefile或地图结构体')
		end
		error('输入的"ksn"无法识别为shapefile或地图结构体');
	end
	
	if grid_flag
		f1=figure(1);
		set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.8],'renderer','painters');
		hold on
		if isempty(ksn_lim)
			imageschs(DEM,ksn,'colormap','ksncolor');
		else
			imageschs(DEM,ksn,'colormap','ksncolor','caxis',ksn_lim);
		end

        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(gca);
        end 
		hold off
	else	
		num_seg=numel(ksn);

		sx=cell(num_seg,1);
		sy=cell(num_seg,1);
		sk=cell(num_seg,1);
		for ii=1:num_seg
			sx{ii,1}=ksn(ii,1).X(:);
			sy{ii,1}=ksn(ii,1).Y(:);
			if isfield(ksn,'ksn')
				sk{ii,1}=ones(numel(sx{ii,1}),1)*ksn(ii,1).ksn;
			elseif isfield(ksn,'fit_ksn')
				sk{ii,1}=ones(numel(sx{ii,1}),1)*ksn(ii,1).fit_ksn;
			else
				if isdeployed
					errordlg('提供的shapefile或地图结构体中不存在有效的"ksn"字段')
				end
				error('提供的shapefile或地图结构体中不存在有效的"ksn"字段')
			end
		end

		sx=vertcat(sx{:});
		sy=vertcat(sy{:});
		sk=vertcat(sk{:});

		ix=coord2ind(DEM,sx,sy);
		idx=isnan(ix);

		ix(idx)=[];
		sk(idx)=[];

		W=GRIDobj(DEM,'logical');
		W.Z(ix)=true;
		S=STREAMobj(FD,W);

		[~,loc,~]=unique(ix);
		sk=sk(loc);

		f1=figure(1);
		set(f1,'Visible','off');

		[RGB]=imageschs(DEM,DEM,'colormap','gray');
		[~,R]=GRIDobj2im(DEM);

		imshow(flipud(RGB),R);
		axis xy
		hold on
		colormap(ksncolor(20));
		plotc(S,sk);
		if isempty(ksn_lim)
			caxis([0 max(sk)]);
		else
			caxis([min(ksn_lim) max(ksn_lim)]);
		end
		c1=colorbar;
		ylabel(c1,'归一化河道陡峭度')
%%
		if ~isempty(knks)
			if ischar(knks) && logical(regexp(knks,regexptranslate('wildcard','*.shp')))
				knk=shaperead(knks);
				knkx=[knk.X];
				knky=[knk.Y];
				scatter(knkx,knky,100,'w','p','filled','MarkerEdgeColor','k');
			elseif istable(knks)
				knkx=knks.x_coord;
				knky=knks.y_coord;
				scatter(knkx,knky,100,'w','p','filled','MarkerEdgeColor','k');
			else
				scatter(knks(:,1),knks(:,2),100,'w','p','filled','MarkerEdgeColor','k');
			end
		end
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(gca);
        end 
		hold off
		set(f1,'Visible','on','Units','normalized','Position',[0.05 0.1 0.8 0.8],'renderer','painters');
	end
end