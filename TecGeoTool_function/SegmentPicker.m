function [Sc]=SegmentPicker(DEM,FD,A,S,basin_num,varargin)
%
% 用法:
% [Sc]=SegmentPicker(DEM,FD,A,S,basin_num);
% [Sc]=SegmentPicker(DEM,FD,A,S,basin_num,'参数名',参数值,...);
%
% 描述:
% 该函数用于从河流网络的顶部选择一段，并绘制该段的纵剖面和chi-Z关系图，同时输出提取的河段和chi结构（来自'chiplot'）。
% 允许用户交互式选择不同的河段并显示。保留所有已选择和确认的河段数据记录。
%
% 必需输入:
% DEM - 数字高程模型（GRIDobj格式），假定为未经过水文修正的原始DEM（例如来自ProcessRiverBasins的DEMoc）
% FD - 流向图（FLOWobj格式）
% A - 汇流累积量（GRIDobj格式）
% S - 河流网络（STREAMobj格式）
% basin_num - 用于输出命名的流域编号，或用于标识所选河段集合的识别号
%
% 可选参数:
% conditioned_DEM [] - 提供经过平滑的DEM用于高程提取（不要将此参数作为主DEM输入！）。参见'ConditionDEM'函数了解修正方法。
% 若未提供，默认使用mincosthydrocon函数生成修正DEM。
% direction ['down'] - 可选'up'或'down'，默认为'down'。选'down'时假设选择的是河道起点；选'up'时假设选择的是流域出口点。
% method ['new_picks'] - 可选'new_picks'或'prev_picks'。选'prev_picks'时必须提供'picks'参数。
% plot_style ['refresh'] - 可选'refresh'或'keep'。控制绘图是否重置或保留历史选择。
% plot_type ['vector'] - 可选'vector'或'grid'。控制河流网络以矢量线还是栅格形式绘制，后者适合大数据但精度较低。
% calc_full_slope_area [false] - 逻辑标志，控制是否计算整个网络的坡度-面积数据（true）或仅主干流（false）
% complete_networks_only [false] - 若为true，在选河前过滤不完整的河网部分
% picks - m×3矩阵（x,y,编号）或点shapefile名称。当direction为'down'时视为河道起点，为'up'时视为流域出口
% ref_concavity [0.50] - 计算Chi-Z的参考凹度值
% min_elev [] - 下游提取的最低高程截止点（仅direction为'down'时有效）
% max_area [] - 下游提取的最大集水面积截止点（仅direction为'down'时有效）
% recalc [false] - 是否重新计算chi值（当使用min_elev/max_area截断时）
% threshold_area [1e6] - 栅格绘图时的集水面积阈值
% interp_value [0.1] - mincosthydrocon的插值参数（未提供conditioned_DEM时使用）
% bin_size [500] - 坡度-面积数据的分箱尺寸（地图单位）
%
% 输出:
% Sc - 包含所有选定河段的STREAMobj
%
% 保存文件:
% 输出'PickedSegements_*.mat'文件包含：
% StreamSgmnts - 选定河段单元数组（STREAMobj格式）
% ChiSgmnts - chi结构单元数组
% SlpAreaSgmnts - 坡度-面积数据单元数组
% Sc - 合并后的河段STREAMobj
% 当direction为'down'时包含Heads矩阵（河道起点坐标和编号）
% 当direction为'up'时包含Outlets矩阵（流域出口坐标和编号，可作为ProcessRiverBasins的river_mouths参数）
%
% 示例:
% [Sc]=SegmentPicker(DEM,FD,A,S,2);
% [Sc]=SegmentPicker(DEM,FD,A,S,2,'direction','up','theta_ref',0.6,'min_grad',0.00001);
% [Sc]=SegmentPicker(DEM,FD,A,S,32,'direction','down','theta_ref',0.5,'min_elev',300,'recalc',true);
% [Sc]=SegmentPicker(DEM,FD,A,S,1,'method','prev_picks','picks',channel_heads); % 输入矩阵
% [Sc]=SegmentPicker(DEM,FD,A,S,1,'method','prev_picks','picks','channel_heads'); % 输入shapefile
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数由Adam M. Forte编写 - 最后更新: 2018年6月18日 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 解析输入参数
p = inputParser;
p.FunctionName = 'SegmentPicker';
addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
addRequired(p,'A',@(x) isa(x,'GRIDobj'));
addRequired(p,'S',@(x) isa(x,'STREAMobj'));
addRequired(p,'basin_num',@(x) isnumeric(x));

addParameter(p,'direction','down',@(x) ischar(validatestring(x,{'down','up'})));
addParameter(p,'method','new_picks',@(x) ischar(validatestring(x,{'new_picks','prev_picks'})));
addParameter(p,'plot_type','vector',@(x) ischar(validatestring(x,{'vector','grid'})));	
addParameter(p,'plot_style','refresh',@(x) ischar(validatestring(x,{'refresh','keep'})));
addParameter(p,'calc_full_slope_area',false,@(x) isscalar(x) && islogical(x));
addParameter(p,'complete_networks_only',false,@(x) isscalar(x) && islogical(x));
addParameter(p,'ref_concavity',0.50,@(x) isscalar(x) && isnumeric(x));
addParameter(p,'min_elev',[],@(x) isscalar(x) && isnumeric(x) || isempty(x));
addParameter(p,'max_area',[],@(x) isscalar(x) && isnumeric(x) || isempty(x));
addParameter(p,'recalc',false,@(x) isscalar(x));
addParameter(p,'picks',[],@(x) (isnumeric(x) && size(x,2)==3) | ischar(x));
addParameter(p,'threshold_area',1e6,@(x) isnumeric(x));
addParameter(p,'conditioned_DEM',[],@(x) isa(x,'GRIDobj') || isempty(x));
addParameter(p,'interp_value',0.1,@(x) isnumeric(x) && x>=0 && x<=1);
addParameter(p,'bin_size',500,@(x) isscalar(x) && isnumeric(x));
addParameter(p,'out_dir',[],@(x) isdir(x));

parse(p,DEM,FD,A,S,basin_num,varargin{:});
DEM=p.Results.DEM;
FD=p.Results.FD;
S=p.Results.S;
A=p.Results.A;
basin_num=p.Results.basin_num;

cno=p.Results.complete_networks_only;
csa=p.Results.calc_full_slope_area;
direction=p.Results.direction;
method=p.Results.method;
theta_ref=p.Results.ref_concavity;
plot_type=p.Results.plot_type;
plot_style=p.Results.plot_style;
threshold_area=p.Results.threshold_area;
points=p.Results.picks;
iv=p.Results.interp_value;
DEMc=p.Results.conditioned_DEM;
bin_size=p.Results.bin_size;
out_dir=p.Results.out_dir;

% 错误检查
if strcmp(direction,'up') && ~isempty(p.Results.min_elev)
	if isdeployed 
		warndlg('当从流域出口向上选择时，最小高程参数将被忽略')
	else
		warning('当从流域出口向上选择时，最小高程参数将被忽略')
	end
elseif strcmp(direction,'up') && ~isempty(p.Results.max_area)
	if isdeployed
		warndlg('当从流域出口向上选择时，最大集水面积参数将被忽略')
	else
		warning('当从流域出口向上选择时，最大集水面积参数将被忽略')
	end
elseif ~isempty(p.Results.min_elev) && ~isempty(p.Results.max_area)
	if isdeployed
		errordlg('不能同时指定最小高程和最大集水面积，请只提供一个参数')
	end
	error('不能同时指定最小高程和最大集水面积，请只提供一个参数')
elseif strcmp(method,'prev_picks') && isempty(p.Results.picks)
	if isdeployed
		errordlg('选择先前点方法时必须提供出口点或河道起点列表')
	end		
	error('选择先前点方法时必须提供出口点或河道起点列表') 
end

% 解析不同输入
if strcmp(direction,'down') && strcmp(plot_style,'keep')
	plot_switch='down_keep';
elseif strcmp(direction,'down') && strcmp(plot_style,'refresh')
	plot_switch='down_ref';
elseif strcmp(direction,'up') && strcmp(plot_style,'keep')
	plot_switch='up_keep';
elseif strcmp(direction,'up') && strcmp(plot_style,'refresh')
	plot_switch='up_ref';
end

if strcmp(method,'prev_picks') && isempty(points)
	if isdeployed
		errordlg('请提供m×3的点矩阵或有效的点shapefile路径');
	end
    error('请提供m×3的点矩阵或有效的点shapefile路径');
elseif strcmp(method,'prev_picks') && ischar(points)
    disp('正在读取shapefile')
    try
        pt_name=[points '.shp'];
        pt_shp=shaperead(pt_name);
        fn=fieldnames(pt_shp);
        pt=horzcat([pt_shp.X]',[pt_shp.Y]',[pt_shp.(fn{4})]');
    catch
    	if isdeployed
			errordlg('读取shapefile错误，请确保：1.文件名不含.shp扩展名 2.文件为点要素 3.包含单值列 4.已安装Mapping Toolbox')
    	end
        error('读取shapefile错误，请确保：1.文件名不含.shp扩展名 2.文件为点要素 3.包含单值列 4.已安装Mapping Toolbox');
    end
elseif strcmp(method,'prev_picks')
	pt=points;
end

% 移除边缘效应（当标记设置时）
if cno
	S=removeedgeeffects(S,FD,DEM);
end

if isempty(out_dir)
	out_dir=pwd;
end

% 水文修正DEM
if isempty(DEMc)
	zc=mincosthydrocon(S,DEM,'interp',iv);
	DEMc=GRIDobj(DEM);
	DEMc.Z(DEMc.Z==0)=NaN;
	DEMc.Z(S.IXgrid)=zc;
end

switch plot_type
case 'grid'
	DA=A.*(A.cellsize^2);
	DA.Z(DA.Z<threshold_area)=0;
	LA=log10(DA);
end

colcol=colorcube(25);
% 去除灰白色系
colcol=colcol(1:20,:);

switch method
case 'new_picks'
	switch plot_switch
	% 从河道起点向下游提取
	case 'down_keep'

		str1='N';
		str2='Y';

		ii=1;
		close all
		f1=figure(1);
		set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');
		clf
		switch plot_type
		case 'grid'
			hold on
			imageschs(DEM,LA,'colormap','parula','colorbar',false);
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(gca);
		    end 
			hold off
		case 'vector'
			hold on
			imageschs(DEM,DEM,'colormap','parula','colorbar',false);
			plot(S,'-w');
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(gca);
		    end				
			hold off
		end

		while strcmpi(str2,'Y')         
			while strcmpi(str1,'N')    	
				
				% 重置短路开关
				short_circ=0;

				figure(f1)
				hold on
				title('缩放平移至感兴趣区域，按回车键开始选择')
				hold off
				pause()

				hold on
				title('选择河道起点附近的点')
				hold off
				[x,y]=ginput(1);
				pOI=[x y];

				% 寻找最近河道起点
				[ch]=streampoi(S,'channelheads','xy');
				distance=sqrt(sum(bsxfun(@minus, ch, pOI).^2,2));
				chOI=ch(distance==min(distance),:);

				% 创建逻辑栅格
				ix=coord2ind(DEM,chOI(:,1),chOI(:,2));
				IX=GRIDobj(DEM);
				IX.Z(ix)=1;
				[ixmat,X,Y]=GRIDobj2mat(IX);
				ixmat=logical(ixmat);
				IX=GRIDobj(X,Y,ixmat);

				% 从河道起点提取至出口
				Sn_t=modify(S,'downstreamto',IX);

				% 检查额外约束条件
				if ~isempty(p.Results.max_area)
					AR=A.*(A.cellsize^2);
					ar=getnal(Sn_t,AR);

					sx=Sn_t.x; sy=Sn_t.y;
					d=Sn_t.distance;
					[d_s,d_ix]=sort(d,'ascend');
					ar_s=ar(d_ix);
					sx=sx(d_ix); sy=sy(d_ix);

					ma=p.Results.max_area;

					ix2=find(ar_s>=ma,1,'last');
					if isempty(ix2)
						if isdeployed 
							warndlg('最大集水面积过大，提取完整河段')
						else
							warning('最大集水面积过大，提取完整河段')
						end
						Sn=Sn_t;
						short_circ=1;
					elseif ix2==numel(ar)
						if isdeployed
							errordlg('最大集水面积过小，未选择到河段')
						end
						error('最大集水面积过小，未选择到河段')
					else
						xn=sx(ix2);
						yn=sy(ix2);

						ix2=coord2ind(DEM,xn,yn);
						IX2=GRIDobj(DEM);
						IX2.Z(ix2)=1;
						[ix2mat,X,Y]=GRIDobj2mat(IX2);
						ix2mat=logical(ix2mat);
						IX2=GRIDobj(X,Y,ix2mat);

						Sn=modify(Sn_t,'upstreamto',IX2);
					end

				elseif ~isempty(p.Results.min_elev);
					el=getnal(Sn_t,DEMc);

					sx=Sn_t.x; sy=Sn_t.y;
					d=Sn_t.distance;
					[d_s,d_ix]=sort(d,'ascend');
					el_s=el(d_ix);
					sx=sx(d_ix); sy=sy(d_ix);

					me=p.Results.min_elev;

					ix2=find(el_s>=me,1,'first');
					if ix2==1
						if isdeployed
							warndlg('最小高程过低，提取完整河段')
						else
							warning('最小高程过低，提取完整河段')
						end
						Sn=Sn_t;
						short_circ=1;
					elseif isempty(ix2)
						if isdeployed
							errordlg('最小高程过高，未选择到河段')
						end
						error('最小高程过高，未选择到河段')
					else
						xn=sx(ix2);
						yn=sy(ix2);

						ix2=coord2ind(DEM,xn,yn);
						IX2=GRIDobj(DEM);
						IX2.Z(ix2)=1;
						[ix2mat,X,Y]=GRIDobj2mat(IX2);
						ix2mat=logical(ix2mat);
						IX2=GRIDobj(X,Y,ix2mat);

						Sn=modify(Sn_t,'upstreamto',IX2);
					end
				else
					Sn=Sn_t;
				end

				hold on
				SP=plot(Sn,'-r','LineWidth',2);
				hold off

	            qa=questdlg('这是否是您要选择的河段？','河段选择','否','是','是');

	            switch qa
	            case '是'
	                str1 = 'Y';
	            case '否'
	                str1 = 'N';
	                delete(SP);
	            end
			end

			if isempty(p.Results.min_elev) && isempty(p.Results.max_area) 
				C=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);
			elseif ~isempty(p.Results.min_elev) || ~isempty(p.Results.max_area) && p.Results.recalc
				C=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);
			elseif short_circ==1
				C=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);
			elseif ~isempty(p.Results.min_elev) || ~isempty(p.Results.max_area) && p.Results.recalc==0
				C_t=chiplot(Sn_t,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);
				C_n=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);

				% 查找匹配部分
				txyz=[C_t.x C_t.y C_t.elev];
				nxyz=[C_n.x C_n.y C_n.elev];
				ix3=ismember(txyz,nxyz,'rows');
				% 重建chi结构
				C=struct;
				C.mn=C_n.mn;
				C.beta=C_n.beta;
				C.betase=C_n.betase;
				C.a0=C_n.a0;
				C.ks=C_n.ks;
				C.R2=C_n.R2;
				C.chi=C_t.chi(ix3);
				C.x=C_t.x(ix3); C.y=C_t.y(ix3);
				C.elev=C_t.elev(ix3); C.elevbl=C_t.elevbl(ix3);
				C.distance=C_t.distance(ix3); C.pred=C_t.pred(ix3);
				C.area=C_t.area(ix3); C.res=C_t.res(ix3);
			end

			f2=figure(2);
			set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');

			sbplt1=subplot(3,1,1);
			hold on
			plot(C.chi,C.elev,'Color',colcol(mod(ii,20)+1,:));
			xlabel('\chi')
			ylabel('高程（米）')
			title('\chi - Z 关系图')
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(sbplt1);
		    end				
			hold off

			sbplt2=subplot(3,1,2);
			hold on
			if isempty(p.Results.min_elev) && isempty(p.Results.max_area) 
				plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
				plotdz(Sn,DEMc,'dunit','km','Color',colcol(mod(ii,20)+1,:));
			elseif ~isempty(p.Results.min_elev) || ~isempty(p.Results.max_area) && p.Results.recalc
				plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
				plotdz(Sn,DEMc,'dunit','km','Color',colcol(mod(ii,20)+1,:));
			elseif short_circ==1
				plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
				plotdz(Sn,DEMc,'dunit','km','Color',colcol(mod(ii,20)+1,:));
			elseif ~isempty(p.Results.min_elev) || ~isempty(p.Results.max_area) && p.Results.recalc==0
				Cu=chiplot(Sn_t,DEM,A,'a0',1,'mn',theta_ref,'plot',false);
				plot((Cu.distance(ix3))./1000,Cu.elev(ix3),'Color',[0.5 0.5 0.5]);
				plot(C.distance./1000,C.elev,'Color',colcol(mod(ii,20)+1,:));
			end
			xlabel('距河口距离（公里）')
			ylabel('高程（米）')
			title('纵剖面图')
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(sbplt2);
		    end					
			hold off

			saax=subplot(3,1,3);
			hold on
			if isempty(p.Results.min_elev) && isempty(p.Results.max_area) 
				if csa
					[bs,ba,aa,ag]=sa(DEMc,Sn,A,bin_size);
				else
					[bs,ba,aa,ag]=sa(DEMc,trunk(Sn),A,bin_size);
				end
			elseif ~isempty(p.Results.min_elev) || ~isempty(p.Results.max_area) && p.Results.recalc								
				if csa
					[bs,ba,aa,ag]=sa(DEMc,Sn,A,bin_size);
				else
					[bs,ba,aa,ag]=sa(DEMc,trunk(Sn),A,bin_size);
				end
			elseif short_circ==1
				if csa
					[bs,ba,aa,ag]=sa(DEMc,Sn,A,bin_size);
				else
					[bs,ba,aa,ag]=sa(DEMc,trunk(Sn),A,bin_size);
				end
			elseif ~isempty(p.Results.min_elev) || ~isempty(p.Results.max_area) && p.Results.recalc==0
				if csa
					[bs,ba,aa,ag]=sa(DEMc,Sn_t,A,bin_size);
				else
					[bs,ba,aa,ag]=sa(DEMc,trunk(Sn_t),A,bin_size);
				end
			end

			scatter(aa,ag,5,[0.5 0.5 0.5],'+');
			scatter(ba,bs,20,colcol(mod(ii,20)+1,:),'filled');	
			set(saax,'Xscale','log','Yscale','log','XDir','reverse');
			xlabel('集水面积对数');
			ylabel('坡度对数');	
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(saax);
		    end					
			hold off

			StreamSgmnts{ii}=Sn;
			ChiSgmnts{ii}=C;
			SlpAreaSgmnts{ii,1}=[bs ba];
			SlpAreaSgmnts{ii,2}=[ag aa];				
			Heads(ii,1)=chOI(:,1);
			Heads(ii,2)=chOI(:,2);
			Heads(ii,3)=ii;

			ii=ii+1;

            qa2=questdlg('继续选择其他河段？','河段选择','否','是','是');
            switch qa2
            case '是'
            	str2 = 'Y';
                str1 = 'N';
            case '否'
                str2 = 'N';
            end				
		end

	% 从流域出口向上游提取
	case 'up_keep'

		str1='N';
		str2='Y';

		ii=1;

			f1=figure(1);
set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');
clf
switch plot_type
case 'grid'
    hold on
    imageschs(DEM,LA,'colormap','parula','colorbar',false);
    if ~verLessThan('matlab','9.5')
        disableDefaultInteractivity(gca);
    end    
    hold off
case 'vector'
    hold on
    imageschs(DEM,DEM,'colormap','parula','colorbar',false);
    plot(S,'-w');
    if ~verLessThan('matlab','9.5')
        disableDefaultInteractivity(gca);
    end                
    hold off
end

while strcmpi(str2,'Y')       
    while strcmpi(str1,'N')           

        figure(f1)
        hold on
        title('缩放和平移到感兴趣的区域，准备好选择时按 "回车" 键')
        hold off
        pause()

        hold on
        title('选择下游的一个点，以计算 Chi-Z')
        hold off
        [x,y]=ginput(1);

        % 构建逻辑栅格
        [xn,yn]=snap2stream(S,x,y);
        ix=coord2ind(DEM,xn,yn);
        IX=GRIDobj(DEM);
        IX.Z(ix)=1;
        [ixmat,X,Y]=GRIDobj2mat(IX);
        ixmat=logical(ixmat);
        IX=GRIDobj(X,Y,ixmat);

        Sn=modify(S,'upstreamto',IX);

        hold on
        SP=plot(Sn,'-r','LineWidth',2);
        hold off

        qa=questdlg('这是您想要的河流段吗？','河流选择','否','是','是');
        switch qa
        case '是'
            str1 = 'Y';
        case '否'
            str1 = 'N';
            delete(SP);
        end
    end

    C=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);

    f2=figure(2);
    set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');

    sbplt1=subplot(3,1,1);
    hold on
    plot(C.chi,C.elev,'Color',colcol(mod(ii,20)+1,:));
    xlabel('\chi')
    ylabel('高程 (m)')
    title('\chi - Z 关系图')
    if ~verLessThan('matlab','9.5')
        disableDefaultInteractivity(sbplt1);
    end                
    hold off

    sbplt2=subplot(3,1,2);
    hold on
    plotdz(Sn,DEMc,'dunit','km','Color',colcol(mod(ii,20)+1,:));
    xlabel('距河口距离 (km)')
    ylabel('高程 (m)')
    title('纵剖面图')
    if ~verLessThan('matlab','9.5')
        disableDefaultInteractivity(sbplt2);
    end                    
    hold off

    saax=subplot(3,1,3);
    hold on
    if csa
        [bs,ba,aa,ag]=sa(DEMc,Sn,A,bin_size);
    else
        [bs,ba,aa,ag]=sa(DEMc,trunk(Sn),A,bin_size);
    end
    scatter(aa,ag,5,[0.5 0.5 0.5],'+');
    scatter(ba,bs,20,colcol(mod(ii,20)+1,:),'filled');
    set(saax,'Xscale','log','Yscale','log','XDir','reverse');
    xlabel('对数集水面积');
    ylabel('对数坡度');    
    if ~verLessThan('matlab','9.5')
        disableDefaultInteractivity(saax);
    end                    
    hold off

    StreamSgmnts{ii}=Sn;
    ChiSgmnts{ii}=C;
    SlpAreaSgmnts{ii,1}=[bs ba];
    SlpAreaSgmnts{ii,2}=[ag aa];    
    Outlets(ii,1)=un;
    Outlets(ii,2)=yn;
    Outlets(ii,3)=ii;

    ii=ii+1;

    qa2=questdlg('继续选择河流？','河流选择','否','是','是');
    switch qa2
    case '是'
        str2 = 'Y';
        str1 = 'N';
    case '否'
        str2 = 'N';
    end
end

case 'down_ref'
    str1='N';
    str2='Y';
    ii=1;
    f1=figure(1);
    set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');
    clf
    switch plot_type
    case 'grid'
        hold on
        imageschs(DEM,LA,'colormap','parula','colorbar',false);
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(gca);
        end    
        hold off
    case 'vector'
        hold on
        imageschs(DEM,DEM,'colormap','parula','colorbar',false);
        plot(S,'-w');
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(gca);
        end    
        hold off
    end

    while strcmpi(str2,'Y')     
        while strcmpi(str1,'N')            
            % 重置短路开关
            short_circ=0;

            figure(1)
            hold on
            title('缩放和平移到感兴趣的区域，准备好选择时按 "回车" 键')
            hold off
            pause()

            hold on
            title('选择靠近感兴趣的河道源头的点')
            hold off
            [x,y]=ginput(1);
            pOI=[x y];

            % 寻找最近的河道源头
            [ch]=streampoi(S,'channelheads','xy');
            distance=sqrt(sum(bsxfun(@minus, ch, pOI).^2,2));
            chOI=ch(distance==min(distance),:);

            % 构建逻辑栅格
            ix=coord2ind(DEM,chOI(:,1),chOI(:,2));
            IX=GRIDobj(DEM);
            IX.Z(ix)=1;
            [ixmat,X,Y]=GRIDobj2mat(IX);
            ixmat=logical(ixmat);
            IX=GRIDobj(X,Y,ixmat);

            % 从河道源头到出口提取河流
            Sn_t=modify(S,'downstreamto',IX);

            % 检查是否指定了额外约束条件
            if ~isempty(p.Results.max_area)
                AR=A.*(A.cellsize^2);
                ar=getnal(Sn_t,AR);

                sx=Sn_t.x; sy=Sn_t.y;
                d=Sn_t.distance;
                [d_s,d_ix]=sort(d,'ascend');
                ar_s=ar(d_ix);
                sx=sx(d_ix); sy=sy(d_ix);

                ma=p.Results.max_area;

                ix2=find(ar_s>=ma,1,'last');
                if isempty(ix2)
                    if isdeployed
                        warndlg('输入的最大流域面积过大，已超出完整河段范围');
                    else
                        warning('输入的最大流域面积过大，已超出完整河段范围');
                    end
                    Sn=Sn_t;
                    short_circ=1;
                elseif ix2==numel(ar)
                    if isdeployed
                        errordlg('输入的最大流域面积过小，未选中任何河段');
                    end
                    error('输入的最大流域面积过小，未选中任何河段');
                else
                    xn=sx(ix2);
                    yn=sy(ix2);

                    ix2=coord2ind(DEM,xn,yn);
                    IX2=GRIDobj(DEM);
                    IX2.Z(ix2)=1;
                    [ix2mat,X,Y]=GRIDobj2mat(IX2);
                    ix2mat=logical(ix2mat);
                    IX2=GRIDobj(X,Y,ix2mat);

                    Sn=modify(Sn_t,'upstreamto',IX2);
                end

            elseif ~isempty(p.Results.min_elev)
                el=getnal(Sn_t,DEMc);

                sx=Sn_t.x; sy=Sn_t.y;
                d=Sn_t.distance;
                [d_s,d_ix]=sort(d,'ascend');
                el_s=el(d_ix);
                sx=sx(d_ix); sy=sy(d_ix);

                me=p.Results.min_elev;

                ix2=find(el_s>=me,1,'first');
                if ix2==1
                    if isdeployed
                        warndlg('输入的最大流域面积过大，已超出完整河段范围');
                    else
                        warning('输入的最大流域面积过大，已超出完整河段范围');
                    end
                    Sn=Sn_t;
                    short_circ=1;
                elseif isempty(ix2)
                    if isdeployed
                        errordlg('输入的最小海拔过高，未选中任何河段');
                    end
                    error('输入的最小海拔过高，未选中任何河段');
                else
                    xn=sx(ix2);
                    yn=sy(ix2);

                    ix2=coord2ind(DEM,xn,yn);
                    IX2=GRIDobj(DEM);
                    IX2.Z(ix2)=1;
                    [ix2mat,X,Y]=GRIDobj2mat(IX2);
                    ix2mat=logical(ix2mat);
                    IX2=GRIDobj(X,Y,ix2mat);

                    Sn=modify(Sn_t,'upstreamto',IX2);
                end
            else
                Sn=Sn_t;
            end

            hold on
            SP=plot(Sn,'-r','LineWidth',2);
            hold off

            qa=questdlg('这是您想要的河流段吗？','河流选择','否','是','是');
            switch qa
            case '是'
                str1 = 'Y';
            case '否'
                str1 = 'N';
                delete(SP);
            end
        end

        if isempty(p.Results.min_elev) && isempty(p.Results.max_area) 
            C=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);
        elseif ~isempty(p.Results.min_elev) || ~isempty(p.Results.max_area) && p.Results.recalc
            C=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);
        elseif short_circ==1
            C=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);
        elseif ~isempty(p.Results.min_elev) || ~isempty(p.Results.max_area) && p.Results.recalc==0
            C_t=chiplot(Sn_t,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);
            C_n=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);

            % 查找范围
            txyz=[C_t.x C_t.y C_t.elev];
            nxyz=[C_n.x C_n.y C_n.elev];
            ix3=ismember(txyz,nxyz,'rows');
            % 重建chi结构体
            C=struct;
            C.mn=C_n.mn;
            C.beta=C_n.beta;
            C.betase=C_n.betase;
            C.a0=C_n.a0;
            C.ks=C_n.ks;
            C.R2=C_n.R2;
            C.chi=C_t.chi(ix3);
            C.x=C_t.x(ix3); C.y=C_t.y(ix3);
            C.elev=C_t.elev(ix3); C.elevbl=C_t.elevbl(ix3);
            C.distance=C_t.distance(ix3); C.pred=C_t.pred(ix3);
            C.area=C_t.area(ix3); C.res=C_t.res(ix3);
        end

        f2=figure(2);
        set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');

        sbplt1=subplot(3,1,1);
        hold on
        plot(C.chi,C.elev,'-k');
        xlabel('\chi')
        ylabel('高程 (m)')
        title('\chi - Z 关系图')
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(sbplt1);
        end                    
        hold off

        sbplt2=subplot(3,1,2);
        hold on
        if isempty(p.Results.min_elev) && isempty(p.Results.max_area) 
            plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
            plotdz(Sn,DEMc,'dunit','km','Color','k');
        elseif ~isempty(p.Results.min_elev) || ~isempty(p.Results.max_area) && p.Results.recalc
            plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
            plotdz(Sn,DEMc,'dunit','km','Color','k');
        elseif short_circ==1
            plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
            plotdz(Sn,DEMc,'dunit','km','Color','k');
        elseif ~isempty(p.Results.min_elev) || ~isempty(p.Results.max_area) && p.Results.recalc==0
            Cu=chiplot(Sn_t,DEM,A,'a0',1,'mn',theta_ref,'plot',false);
            plot((Cu.distance(ix3))./1000,Cu.elev(ix3),'Color',[0.5 0.5 0.5]);
            plot(C.distance./1000,C.elev,'-k');
        end
        xlabel('距河口距离 (km)')
        ylabel('高程 (m)')
        legend('未平滑DEM','平滑后DEM','location','best');
        title('纵剖面图')
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(sbplt2);
        end                    
        hold off

        saax=subplot(3,1,3);
        hold on
        if isempty(p.Results.min_elev) && isempty(p.Results.max_area) 
            if csa
                [bs,ba,aa,ag]=sa(DEMc,Sn,A,bin_size);
            else
                [bs,ba,aa,ag]=sa(DEMc,trunk(Sn),A,bin_size);
            end
        elseif ~isempty(p.Results.min_elev) || ~isempty(p.Results.max_area) && p.Results.recalc                                
            if csa
                [bs,ba,aa,ag]=sa(DEMc,Sn,A,bin_size);
            else
                [bs,ba,aa,ag]=sa(DEMc,trunk(Sn),A,bin_size);
            end
        elseif short_circ==1
            if csa
                [bs,ba,aa,ag]=sa(DEMc,Sn,A,bin_size);
            else
                [bs,ba,aa,ag]=sa(DEMc,trunk(Sn),A,bin_size);
            end
        elseif ~isempty(p.Results.min_elev) || ~isempty(p.Results.max_area) && p.Results.recalc==0
            if csa
                [bs,ba,aa,ag]=sa(DEMc,Sn_t,A,bin_size);
            else
                [bs,ba,aa,ag]=sa(DEMc,trunk(Sn_t),A,bin_size);
            end
        end        
        scatter(aa,ag,5,[0.5 0.5 0.5],'+');
        scatter(ba,bs,20,'k','filled');
        set(saax,'Xscale','log','Yscale','log','XDir','reverse');
        xlabel('对数集水面积');
        ylabel('对数坡度');
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(saax);
        end                        
        hold off                                

        StreamSgmnts{ii}=Sn;
        ChiSgmnts{ii}=C;
        SlpAreaSgmnts{ii,1}=[bs ba];
        SlpAreaSgmnts{ii,2}=[ag aa];    
        Heads(ii,1)=chOI(:,1);
        Heads(ii,2)=chOI(:,2);
        Heads(ii,3)=ii;

        ii=ii+1;

        qa2=questdlg('继续选择河流？','河流选择','否','是','是');
        switch qa2
        case '是'
            str2 = 'Y';
            str1 = 'N';
            close figure 2
        case '否'
            str2 = 'N';
        end    
    end

case 'up_ref'
    str1='N';
    str2='Y';
    ii=1;
    f1=figure(1);
    set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');
    clf
    switch plot_type
    case 'grid'
        hold on
        imageschs(DEM,LA,'colormap','parula','colorbar',false);
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(gca);
        end                    
        hold off
    case 'vector'
        hold on
        imageschs(DEM,DEM,'colormap','parula','colorbar',false);
        plot(S,'-w');
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(gca);
        end                    
        hold off
    end

    while strcmpi(str2,'Y')        
        while strcmpi(str1,'N')               
            figure(1);
            hold on
            title('缩放和平移到感兴趣的区域，准备好选择时按 "回车"键')
            hold off
            pause()

            hold on
            title('选择下游的一个点，以计算 Chi-Z')
            hold off
            [x,y]=ginput(1);

            % 构建逻辑栅格
            [xn,yn]=snap2stream(S,x,y);
            ix=coord2ind(DEM,xn,yn);
            IX=GRIDobj(DEM);
            IX.Z(ix)=1;
            [ixmat,X,Y]=GRIDobj2mat(IX);
            ixmat=logical(ixmat);
            IX=GRIDobj(X,Y,ixmat);

            Sn=modify(S,'upstreamto',IX);

            hold on
            SP=plot(Sn,'-r','LineWidth',2);
            hold off

            qa=questdlg('这是您想要的河流段吗？','河流选择','否','是','是');
            switch qa
            case '是'
                str1 = 'Y';
            case '否'
                str1 = 'N';
                delete(SP);
            end
        end

        C=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);

        f2=figure(2);
        set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');
        sbplt1=subplot(3,1,1);
        hold on
        plot(C.chi,C.elev,'-k');
        xlabel('\chi')
        ylabel('高程 (m)')
        title('\chi - Z 关系图')
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(sbplt1);
        end                    
        hold off

        sbplt2=subplot(3,1,2);
        hold on
        plotdz(Sn,DEMc,'dunit','km','Color','k');
        xlabel('距河口距离 (km)')
        ylabel('高程 (m)')
        title('纵剖面图')
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(sbplt2);
        end                    
        hold off

        saax=subplot(3,1,3);
        hold on
        if csa
            [bs,ba,aa,ag]=sa(DEMc,Sn,A,bin_size);
        else
            [bs,ba,aa,ag]=sa(DEMc,trunk(Sn),A,bin_size);
        end
        scatter(aa,ag,5,[0.5 0.5 0.5],'+');
        scatter(ba,bs,20,'k','filled');
        set(saax,'Xscale','log','Yscale','log','XDir','reverse');
        xlabel('对数集水面积');
        ylabel('对数坡度');
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(saax);
        end                        
        hold off

        StreamSgmnts{ii}=Sn;
        ChiSgmnts{ii}=C;
        SlpAreaSgmnts{ii,1}=[bs ba];
        SlpAreaSgmnts{ii,2}=[ag aa];    
        Outlets(ii,1)=xn;
        Outlets(ii,2)=yn;
        Outlets(ii,3)=ii;

        ii=ii+1;

        qa2=questdlg('继续选择河流？','河流选择','否','是','是');
        switch qa2
        case '是'
            str2 = 'Y';
            str1 = 'N';
            close figure 2
        case '否'
            str2 = 'N';
        end    
    end
end

case 'prev_picks'
    switch direction
    case 'down'
        heads=pt;
        [num_heads,~]=size(heads);

        w1=waitbar(0,'提取河段中...');
        for ii=1:num_heads
            short_circ=0;
            x=heads(ii,1); y=heads(ii,2);
            pOI=[x y];

            % 寻找最近的河道源头
            [ch]=streampoi(S,'channelheads','xy');
            distance=sqrt(sum(bsxfun(@minus, ch, pOI).^2,2));
            chOI=ch(distance==min(distance),:);

            % 构建逻辑栅格
            ix=coord2ind(DEM,chOI(:,1),chOI(:,2));
            IX=GRIDobj(DEM);
            IX.Z(ix)=1;
            [ixmat,X,Y]=GRIDobj2mat(IX);
            ixmat=logical(ixmat);
            IX=GRIDobj(X,Y,ixmat);

            % 从河道源头到出口提取河流
            Sn_t=modify(S,'downstreamto',IX);

            % 检查是否指定了额外约束条件
            if ~isempty(p.Results.max_area)
                AR=A.*(A.cellsize^2);

                sx=Sn_t.x; sy=Sn_t.y;
                d=Sn_t.distance;
                [d_s,d_ix]=sort(d,'ascend');
                ar_s=ar(d_ix);
                sx=sx(d_ix); sy=sy(d_ix);

                ma=p.Results.max_area;

                ix2=find(ar_s>=ma,1,'last');
                if isempty(ix2)
                    if isdeployed
                        warndlg('输入的最大流域面积过大，已超出完整河段范围');
                    else
                        warning('输入的最大流域面积过大，已超出完整河段范围');
                    end
                    Sn=Sn_t;
                    short_circ=1;
                elseif ix2==numel(ar)
                    if isdeployed
                        errordlg('输入的最大流域面积过小，未选中任何河段');
                    end
                    error('输入的最大流域面积过小，未选中任何河段');
                else
                    xn=sx(ix2);
                    yn=sy(ix2);

                    ix2=coord2ind(DEM,xn,yn);
                    IX2=GRIDobj(DEM);
                    IX2.Z(ix2)=1;
                    [ix2mat,X,Y]=GRIDobj2mat(IX2);
                    ix2mat=logical(ix2mat);
                    IX2=GRIDobj(X,Y,ix2mat);

                    Sn=modify(Sn_t,'upstreamto',IX2);
                end

            elseif ~isempty(p.Results.min_elev)
                el=getnal(Sn_t,DEMc);

                sx=Sn_t.x; sy=Sn_t.y;
                d=Sn_t.distance;
                [d_s,d_ix]=sort(d,'ascend');
                el_s=el(d_ix);
                sx=sx(d_ix); sy=sy(d_ix);

                me=p.Results.min_elev;

                ix2=find(el_s>=me,1,'first');
                if ix2==1
                    if isdeployed
                        warndlg('输入的最大流域面积过大，已超出完整河段范围');
                    else
                        warning('输入的最大流域面积过大，已超出完整河段范围');
                    end
                    Sn=Sn_t;
                elseif isempty(ix2)
                    if isdeployed
                        errordlg('输入的最大流域面积过小，未选中任何河段');
                    end
                    error('输入的最大流域面积过小，未选中任何河段');
                else
                    xn=sx(ix2);
                    yn=sy(ix2);

                    ix2=coord2ind(DEM,xn,yn);
                    IX2=GRIDobj(DEM);
                    IX2.Z(ix2)=1;
                    [ix2mat,X,Y]=GRIDobj2mat(IX2);
                    ix2mat=logical(ix2mat);
                    IX2=GRIDobj(X,Y,ix2mat);

                    Sn=modify(Sn_t,'upstreamto',IX2);
                end
            else
                Sn=Sn_t;
            end

            if csa
                [bs,ba,aa,ag]=sa(DEMc,Sn,A,bin_size);
            else
                [bs,ba,aa,ag]=sa(DEMc,trunk(Sn),A,bin_size);
            end

            if isempty(p.Results.min_elev) && isempty(p.Results.max_area) 
                C=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);
            elseif ~isempty(p.Results.min_elev) | ~isempty(p.Results.max_area) && p.Results.recalc
                C=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);
            elseif short_circ==1;
                C=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);
            elseif ~isempty(p.Results.min_elev) | ~isempty(p.Results.max_area) && p.Results.recalc==0
                C_t=chiplot(Sn_t,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);
                C_n=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);

                % 查找范围
                txyz=[C_t.x C_t.y C_t.elev];
                nxyz=[C_n.x C_n.y C_n.elev];
                ix3=ismember(txyz,nxyz,'rows');

                % 重建chi结构体
                C=C_t;
                C.chi=C_t.chi(ix3);
                C.x=C_t.x(ix3); C.y=C_t.y(ix3);
                C.elev=C_t.elev(ix3); C.elevbl=C_t.elevbl(ix3);
                C.distance=C_t.distance(ix3); C.pred=C_t.pred(ix3);
                C.area=C_t.area(ix3); C.res=C_t.res(ix3);
            end

            StreamSgmnts{ii}=Sn;
            ChiSgmnts{ii}=C;
            SlpAreaSgmnts{ii,1}=[bs ba];
            SlpAreaSgmnts{ii,2}=[ag aa];    
            Heads(ii,1)=chOI(:,1);
            Heads(ii,2)=chOI(:,2);
            Heads(ii,3)=ii;
            waitbar(ii/num_heads);
        end
        close(w1);

    case 'up'
        outlets=pt;
        [num_outs,~]=size(outlets);

        w1=waitbar(0,'提取河段中...');
        for ii=1:num_outs
            x=outlets(ii,1); y=outlets(ii,2);
            % 构建逻辑栅格
            [xn,yn]=snap2stream(S,x,y);
            ix=coord2ind(DEM,xn,yn);
            IX=GRIDobj(DEM);
            IX.Z(ix)=1;
            [ixmat,X,Y]=GRIDobj2mat(IX);
            ixmat=logical(ixmat);
            IX=GRIDobj(X,Y,ixmat);

            Sn=modify(S,'upstreamto',IX);
            Sn=klargestconncomps(Sn,1);
            C=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);
            if csa
                [bs,ba,aa,ag]=sa(DEMc,Sn,A,bin_size);
            else
                [bs,ba,aa,ag]=sa(DEMc,trunk(Sn),A,bin_size);
            end

            StreamSgmnts{ii}=Sn;
            ChiSgmnts{ii}=C;
            SlpAreaSgmnts{ii,1}=[bs ba];
            SlpAreaSgmnts{ii,2}=[ag aa];    
            Outlets(ii,1)=xn;
            Outlets(ii,2)=yn;
            Outlets(ii,3)=ii;
            waitbar(ii/num_outs);
        end
        close(w1);
    end
end

% 清理并生成输出
num_picks=numel(StreamSgmnts);
if num_picks==1
    Sc=StreamSgmnts{1};
else
    Sc=StreamSgmnts{1};
    for ii=2:num_picks
        Sc=union(Sc,StreamSgmnts{ii});
    end
end

fileOut=fullfile(out_dir,['PickedSegments_' num2str(basin_num) '.mat']);
switch direction
case 'up'
    save(fileOut,'StreamSgmnts','ChiSgmnts','SlpAreaSgmnts','Outlets','Sc','-v7.3');
case 'down'
    save(fileOut,'StreamSgmnts','ChiSgmnts','SlpAreaSgmnts','Heads','Sc','-v7.3');    
end
end

function [bs,ba,a,g]=sa(DEM,S,A,bin_size)
    % 改进的坡度-面积函数，使用平滑长度确定分箱数量，并使用相同分箱计算chi和距离的平均值用于绘图
    
    minX=min(S.distance);
    maxX=max(S.distance);
    b=[minX:bin_size:maxX+bin_size];
    
    numbins=round(max([numel(b) numel(S.IXgrid)/10]));
    
    an=getnal(S,A.*A.cellsize^2);
    z=getnal(S,DEM);
    gn=gradient(S,z,'unit','tangent','method','robust','drop',20);
    gn=smooth(gn,3);
    
    % 通过STREAMobj2XY处理确保chi和其他参数尺寸一致
    [~,~,a,g]=STREAMobj2XY(S,an,gn);
    % 移除NaN值
    a(isnan(a))=[];
    g(isnan(g))=[];
    
    mina=min(a);
    maxa=max(a);
    
    edges = logspace(log10(mina-0.1),log10(maxa+1),numbins+1);
    try
        % histc已弃用
        [ix]=discretize(a,edges);
    catch
        [~,ix] = histc(a,edges);
    end
    
    ba=accumarray(ix,a,[numbins 1],@median,nan);
    bs=accumarray(ix,g,[numbins 1],@(x) mean(x(~isnan(x))),nan);
    
    % 过滤负值
    idx=bs>=0 & ba>=0;
    bs=bs(idx);
    ba=ba(idx);
    
    idx=a>=0 & g>=0;
    a=a(idx);
    g=g(idx);
end

