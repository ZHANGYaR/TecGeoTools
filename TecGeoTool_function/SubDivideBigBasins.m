function SubDivideBigBasins(basin_dir,max_basin_size,divide_method,varargin)
	%
% 用法：
%   SubDivideBigBasins(basin_dir,max_basin_size,divide_method);
%   SubDivideBigBasins(basin_dir,max_basin_size,divide_method,'name',value,...);
%
% 描述：
%   该函数处理由 'ProcessRiverBasins' 函数生成的输出，将任何汇水面积超过指定大小的流域细分，
%   并输出裁剪后的DEM、河网、各种地形指标以及河流值（ks、ksn、chi）。
%
% 必需输入：
%     basin_dir - 包含由 'ProcessRiverBasins' 函数生成的mat文件的文件夹完整路径。
%     max_basin_size - 超过该面积的流域将被细分，单位为平方公里。
%     divide_method - 用于细分流域的方法，可选项如下（'confluences' 和 'up_confluences' 不推荐用于大数据集）：
%         'order' - 使用用户通过可选参数 's_order' 指定的河流阶数的出水口进行细分。
%         'confluences' - 使用河流汇合处的位置（将生成大量子流域）。包含内部参数用于移除
%             极短河流以防止代码出错。
%         'up_confluences' - 使用汇合点上游的位置（将生成大量子流域）。包含内部参数
%             用于移除极短河流以防止代码出错。
%         'filtered_confluences' - 使用汇合点的位置，仅当汇合点以上的流域面积达到用户通过
%             可选参数 'min_basin_size' 指定的大小。
%         'p_filtered_confluences' - 类似于 'filtered_confluences'，但用户通过可选参数 'min_basin_size'
%             以主流域面积的百分比形式指定阈值。
%         'trunk' - 使用主流域内与干流交汇的支流交汇点作为子流域的出水口。包含内部参数用于
%             移除极短河流以防止代码出错。
%         'filtered_trunk' - 与 'trunk' 类似，但仅包含面积大于 'min_basin_size' 的流域。
%         'p_filtered_trunk' - 与 'filtered_trunk' 类似，但 'min_basin_size' 被解释为主流域面积的百分比。
%
% 可选输入：
%     SBFiles_Dir ['SubBasins'] - 子流域文件存储的文件夹名称（位于主流域文件夹中）。子流域文件
%         现在存储在单独的文件夹中，方便基于不同需求创建不同的子流域集合。
%     recursive [true] - 逻辑标志，确保输出中没有子流域面积超过提供的 'max_basin_size'。
%         如果 'divide_method' 是干流类型之一，代码将继续重新定义干流并进一步分割子流域，
%         直到没有提取的流域面积超过 'max_basin_size'。如果 'divide_method' 是汇合类型之一，
%         则面积超过 'max_basin_size' 的子流域将不会包含在输出中。递归检查未对 'order' 方法实现。
%     threshold_area [1e6] - 定义河流的最小积累面积，单位为平方米。
%     segment_length [1000] - 计算ksn时的平滑距离，单位为米，建议值为1000米。
%     ref_concavity [0.5] - 计算ksn的参考凹度。
%     write_arc_files [false] - 设置为 true 输出多种网格的ASCII文件和ksn的shapefile，设置为 false 不输出。
%     s_order [3] - 如果 'divide_method' 是 'order'，用于定义细分河流出水口的河流阶数
%         （较低的阶数会产生更多子流域）。
%     min_basin_size [10] - 自动选择子流域的最小流域面积。如果 'divide_method' 是 'filtered_confluences'，
%         此值被解释为最小汇水面积，单位为平方公里。如果 'divide_method' 是 'p_filtered_confluences'，
%         此值被解释为输入流域面积的百分比，应输入0到100之间的值。
%     no_nested [false] - 逻辑标志，与 'filtered_confluences' 或 'p_filtered_confluences' 结合使用时，
%         仅提取满足汇水面积要求的最低阶流域（避免生成嵌套流域）。
%
% 示例：
%     SubdivideBigBasins('/Users/JoeBlow/Project',100,'confluences');
%     SubdivideBigBasins('/Users/JoeBlow/Project',100,'order','s_order',2,'threshold_area',1e5,'write_arc_files',true);
%
% 注意：
%  - 只有 'order'、'trunk'、'filtered_trunk' 和 'p_filtered_trunk' 方法不会生成嵌套子流域。
%  - 生成连续ksn网格所需的插值会在极小流域上失败。这不会导致代码失败，但会导致
%     没有 'KsnOBJc' 被保存到这些流域。
%  - 方法 'confluences'、'up_confluences' 和 'trunk' 可能尝试提取非常小的流域。内部检查会尝试
%     移除这些非常小的流域，但不总是有效，有时可能导致错误。如果遇到错误，
%     请尝试使用基于汇水面积过滤的方法。
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数作者：Adam M. Forte - 更新日期：2018年6月18日 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	% 解析输入参数
	p = inputParser;
	p.FunctionName = 'SubDivideBigBasins';
	addRequired(p,'basin_dir',@(x) isdir(x));
	addRequired(p,'max_basin_size',@(x) isnumeric(x));
	addRequired(p,'divide_method',@(x) ischar(validatestring(x,{'order','confluences','up_confluences','filtered_confluences','p_filtered_confluences','trunk','filtered_trunk','p_filtered_trunk'})));

	addParameter(p,'SBFiles_Dir','SubBasins',@(x) ischar(x));
	addParameter(p,'recursive',true,@(x) isscalar(x) && islogical(x));
	addParameter(p,'ref_concavity',0.5,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'threshold_area',1e6,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'segment_length',1000,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'write_arc_files',false,@(x) isscalar(x));
	addParameter(p,'s_order',[3],@(x) isscalar(x));
	addParameter(p,'min_basin_size',[10],@(x) isnumeric(x) & isscalar(x));
	addParameter(p,'no_nested',false,@(x) isscalar(x) && islogical(x));

	parse(p,basin_dir,max_basin_size,divide_method,varargin{:});
	location_of_data_files=p.Results.basin_dir;
	max_basin_size=p.Results.max_basin_size;
	divide_method=p.Results.divide_method;

	SBFiles_Dir=p.Results.SBFiles_Dir;
	recursive=p.Results.recursive;
	theta_ref=p.Results.ref_concavity;
	threshold_area=p.Results.threshold_area;
	segment_length=p.Results.segment_length;
	write_arc_files=p.Results.write_arc_files;
	s_order=p.Results.s_order;
	min_basin_size=p.Results.min_basin_size;
	no_nested=p.Results.no_nested;


	FileList=dir(fullfile(location_of_data_files,'*Data.mat'));
	num_files=numel(FileList);

	% 创建子流域目录（如果不存在）
	sb_path=fullfile(location_of_data_files,SBFiles_Dir);
	if ~isdir(sb_path)
		mkdir(sb_path);
	end

	% 参数合法性检查
	if strcmp(divide_method,'p_filtered_confluences') || strcmp(divide_method,'p_filtered_trunk') && (min_basin_size>100 || min_basin_size<=0)
		min_basin_size
		if isdeployed
			errordlg('当分割方法为"p_filtered_confluences"时，"min_basin_size"参数必须介于0到100之间')
		end
		error('当分割方法为"p_filtered_confluences"时，"min_basin_size"参数必须介于0到100之间')
	end

	% 主文件循环开始
	w1=waitbar(0,'正在分割流域');
	for ii=1:num_files
		FileName=fullfile(FileList(ii,1).folder,FileList(ii,1).name);

		% 载入流域面积进行初步检查
		load(FileName,'drainage_area');
		DA=drainage_area;

		% 判断当前流域是否需要处理
		if DA>=max_basin_size
			
			% 载入必要的基础数据文件并重命名
			load(FileName,'RiverMouth','DEMoc','DEMcc','FDc','Ac','Sc','ksn_method','gradient_method','radius');
			DEM=DEMoc;
			DEMhc=DEMcc;
			S=Sc;
			FD=FDc;
			A=Ac;
			RM=RiverMouth;	
			basin_num=RM(:,3);

			if strcmp(ksn_method,'trunk')
				load(FileName,'min_order');
			end

			% 检查段长度参数合理性
			if (DEM.cellsize*3)>segment_length
				segment_length=DEM.cellsize*3;
				disp('警告：段长度参数过小，已自动调整为DEM像元尺寸的三倍');
			end

			waitbar(ii/num_files,w1,['正在处理流域编号 ' num2str(basin_num) ' - 正在确定分割数量']);

			% 计算网格化流域面积（平方公里）
			DAG=(A.*(A.cellsize^2))/1e6;

			% 根据分割方法选择不同的处理逻辑
			switch divide_method
			case 'order' % 基于河流阶数的方法
				so=streamorder(S);
				if s_order<max(so) 
					Se=modify(S,'streamorder',s_order);
					outs=streampoi(Se,'outlets','xy');
					x=outs(:,1);
					y=outs(:,2);
					num_new_basins=numel(x);
				elseif s_order>=max(so) && max(so)>1
					s_order=s_order-1;
					Se=modify(S,'streamorder',s_order);
					outs=streampoi(Se,'outlets','xy');
					x=outs(:,1);
					y=outs(:,2);
					num_new_basins=numel(x);	
				else
					s_order=max(so);
					Se=modify(S,'streamorder',s_order);
					outs=streampoi(Se,'outlets','xy');
					x=outs(:,1);
					y=outs(:,2);
					num_new_basins=numel(x);
				end			
			case 'confluences' % 基于所有汇合点的方法
				S=removeshortstreams(S,DEM.cellsize*10);	
				cons=streampoi(S,'confluences','xy');
				if recursive % 递归检查流域面积
					cons_ix=streampoi(S,'confluences','ix');
					idx=DAG.Z(cons_ix)<max_basin_size;
					x=cons(idx,1);
					y=cons(idx,2);
				else
					x=cons(:,1);
					y=cons(:,2);
				end
				num_new_basins=numel(x);
			case 'up_confluences' % 基于上游汇合点的方法
				S=removeshortstreams(S,DEM.cellsize*10);
				cons=streampoi(S,'bconfluences','xy');
				if recursive
					cons_ix=streampoi(S,'bconfluences','ix');
					idx=DAG.Z(cons_ix)<max_basin_size;
					x=cons(idx,1);
					y=cons(idx,2);
				else
					x=cons(:,1);
					y=cons(:,2);
				end
				num_new_basins=numel(x);
			case 'filtered_confluences' % 基于面积过滤的汇合点方法
				if no_nested % 排除嵌套流域模式
					cons_ix=streampoi(S,'bconfluences','ix');
					if recursive
						da_idx=DAG.Z(cons_ix)>=min_basin_size & DAG.Z(cons_ix)<max_basin_size;
						cons_ix=cons_ix(da_idx);
					else
						da_idx=DAG.Z(cons_ix)>=min_basin_size;
						cons_ix=cons_ix(da_idx);
					end
					[x,y]=CheckUpstream(DEM,FD,cons_ix);
					num_new_basins=numel(x);
				else % 允许嵌套流域模式
					cons_ix=streampoi(S,'confluences','ix');
					cons=streampoi(S,'confluences','xy');
					if recursive
						da_idx=DAG.Z(cons_ix)>=min_basin_size & DAG.Z(cons_ix)<max_basin_size;
					else
						da_idx=DAG.Z(cons_ix)>=min_basin_size;
					end
					cons=cons(da_idx,:);
					x=cons(:,1);
					y=cons(:,2);
					num_new_basins=numel(x);
				end
			case 'p_filtered_confluences' % 基于百分比面积过滤的汇合点方法
				if no_nested
					cons_ix=streampoi(S,'bconfluences','ix');
					da_cons=DAG.Z(cons_ix);
					mbz=DA*(min_basin_size/100);
					if recursive
						da_idx=da_cons>=mbz & da_cons<max_basin_size;
					else
						da_idx=da_cons>=mbz;
					end
					[x,y]=CheckUpstream(DEM,FD,cons_ix(da_idx));
					num_new_basins=numel(x);
				else
					cons_ix=streampoi(S,'confluences','ix');
					cons=streampoi(S,'confluences','xy');
					da_cons=DAG.Z(cons_ix);
					mbz=DA*(min_basin_size/100);
					if recursive
						da_idx=da_cons>=mbz & da_cons<max_basin_size;
					else
						da_idx=da_cons>=mbz;
					end
					cons=cons(da_idx,:);
					x=cons(:,1);
					y=cons(:,2);
					num_new_basins=numel(x);
				end
			case 'trunk' % 基于主干河道的方法
				ST=trunk(klargestconncomps(S,1));
				S=removeshortstreams(S,DEM.cellsize*10);
				tix=streampoi(S,'bconfluences','ix');
				tix=ismember(ST.IXgrid,tix);
				ds=ST.distance;
				ds(~tix)=NaN;
				[~,tix]=max(ds);
				SupT=modify(S,'tributaryto',ST);
				cons=streampoi(SupT,'outlets','xy');
				cons_ix=streampoi(SupT,'outlets','ix');
				cons_ix=vertcat(cons_ix,ST.IXgrid(tix));
				x=cons(:,1); x=vertcat(x,ST.x(tix));
				y=cons(:,2); y=vertcat(y,ST.y(tix));
				num_new_basins=numel(x);

				if recursive % 递归处理大流域
					try
						rec_count=1;
						while any(DAG.Z(cons_ix)>=max_basin_size) && rec_count<=10
							nidx=DAG.Z(cons_ix)>=max_basin_size;
							if any(nidx)
								x(nidx)=[];
								y(nidx)=[];
								ixs=cons_ix(nidx);
								for jj=1:numel(ixs)
									TIX=GRIDobj(DEM,'logical');
									TIX.Z(ixs(jj))=true;
									S_sub=modify(S,'upstreamto',TIX);
									S_sub=removeshortstreams(S_sub,DEM.cellsize*10);
									ST_sub=trunk(S_sub);
									tix=streampoi(S_sub,'bconfluences','ix');
									tix=ismember(ST_sub.IXgrid,tix);
									ds=ST_sub.distance;
									ds(~tix)=NaN;
									[~,tix]=max(ds);
									SupT_sub=modify(S_sub,'tributaryto',ST_sub);
									cons=streampoi(SupT_sub,'outlets','xy');
									cons_ix=streampoi(SupT_sub,'outlets','ix');
									cons_ix=vertcat(cons_ix,ST_sub.IXgrid(tix));
									xx=cons(:,1); xx=vertcat(xx,ST_sub.x(tix));
									yy=cons(:,2); yy=vertcat(yy,ST_sub.y(tix));
									x=vertcat(x,xx);
									y=vertcat(y,yy);
								end
							end
							rec_count=rec_count+1;
							if rec_count>10
								if isdeployed
									warndlg(['流域编号 ' num2str(basin_num) ' 的分割提前终止以防止无限循环'])
								end
								warning(['流域编号 ' num2str(basin_num) ' 的分割提前终止以防止无限循环']);
							end
						end
						num_new_basins=numel(x);
					catch
						if isdeployed
							warndlg(['流域编号 ' num2str(basin_num) ' 的递归分割失败，将进行常规分割'])
						end
						warning(['流域编号 ' num2str(basin_num) ' 的递归分割失败，将进行常规分割']);
					end
				end

			case 'filtered_trunk' % 基于面积过滤的主干河道方法
				ST=trunk(klargestconncomps(S,1));
				S=removeshortstreams(S,DEM.cellsize*10);
				tix=streampoi(S,'bconfluences','ix');
				tix=ismember(ST.IXgrid,tix);
				ds=ST.distance;
				ds(~tix)=NaN;
				[~,tix]=max(ds);
				SupT=modify(S,'tributaryto',ST);				
				cons_ix=streampoi(SupT,'outlets','ix');
				cons_ix=vertcat(cons_ix,ST.IXgrid(tix));
				cons=streampoi(SupT,'outlets','xy');
				cons=vertcat(cons,[ST.x(tix) ST.y(tix)]);
				da_cons=DAG.Z(cons_ix);
				da_idx=da_cons>=min_basin_size;
				cons=cons(da_idx,:);
				cons_ix=cons_ix(da_idx);
				x=cons(:,1);
				y=cons(:,2);
				num_new_basins=numel(x);

				if recursive
					try
						rec_count=1;
						while any(DAG.Z(cons_ix)>=max_basin_size) && rec_count<=10
							nidx=DAG.Z(cons_ix)>=max_basin_size;
							if any(nidx)
								x(nidx)=[];
								y(nidx)=[];
								ixs=cons_ix(nidx);
								for jj=1:numel(ixs)
									TIX=GRIDobj(DEM,'logical');
									TIX.Z(ixs(jj))=true;
									S_sub=modify(S,'upstreamto',TIX);
									S_sub=removeshortstreams(S_sub,DEM.cellsize*10);
									ST_sub=trunk(S_sub);
									tix=streampoi(S_sub,'bconfluences','ix');
									tix=ismember(ST_sub.IXgrid,tix);
									ds=ST_sub.distance;
									ds(~tix)=NaN;
									[~,tix]=max(ds);
									SupT_sub=modify(S_sub,'tributaryto',ST_sub);
									cons=streampoi(SupT_sub,'outlets','xy');
									cons_ix=streampoi(SupT_sub,'outlets','ix');
									cons_ix=vertcat(cons_ix,ST_sub.IXgrid(tix));
									xx=cons(:,1); xx=vertcat(xx,ST_sub.x(tix));
									yy=cons(:,2); yy=vertcat(yy,ST_sub.y(tix));
									da_cons=DAG.Z(cons_ix);
									da_idx=da_cons>=min_basin_size;
									x=vertcat(x,xx(da_idx));
									y=vertcat(y,yy(da_idx));
								end
							end
							rec_count=rec_count+1;
							if rec_count>10
								if isdeployed
									warndlg(['流域编号 ' num2str(basin_num) ' 的分割提前终止以防止无限循环'])
								end								
								warning(['流域编号 ' num2str(basin_num) ' 的分割提前终止以防止无限循环']);
							end						
						end
						num_new_basins=numel(x);
					catch
						if isdeployed
							warndlg(['流域编号 ' num2str(basin_num) ' 的递归分割失败，将进行常规分割'])
						end
						warning(['流域编号 ' num2str(basin_num) ' 的递归分割失败，将进行常规分割']);
					end
				end

			case 'p_filtered_trunk' % 基于百分比面积过滤的主干河道方法
				ST=trunk(klargestconncomps(S,1));
				S=removeshortstreams(S,DEM.cellsize*10);
				tix=streampoi(S,'bconfluences','ix');
				tix=ismember(ST.IXgrid,tix);
				ds=ST.distance;
				ds(~tix)=NaN;
				[~,tix]=max(ds);
				SupT=modify(S,'tributaryto',ST);
				cons_ix=streampoi(SupT,'confluences','ix');
				cons_ix=vertcat(cons_ix,ST.IXgrid(tix));
				cons=streampoi(SupT,'confluences','xy');
				cons=vertcat(cons,[ST.x(tix) ST.y(tix)]);
				da_cons=DAG.Z(cons_ix);
				mbz=DA*(min_basin_size/100);
				da_idx=da_cons>=mbz;
				cons=cons(da_idx,:);
				cons_ix=cons_ix(da_idx);
				x=cons(:,1);
				y=cons(:,2);
				num_new_basins=numel(x);

				if recursive
					try
						rec_count=1;
						while any(DAG.Z(cons_ix)>=max_basin_size) && rec_count<=10
							nidx=DAG.Z(cons_ix)>=max_basin_size;
							if any(nidx)
								x(nidx)=[];
								y(nidx)=[];
								ixs=cons_ix(nidx);
								for jj=1:numel(ixs)
									TIX=GRIDobj(DEM,'logical');
									TIX.Z(ixs(jj))=true;
									S_sub=modify(S,'upstreamto',TIX);
									S_sub=removeshortstreams(S_sub,DEM.cellsize*10);
									ST_sub=trunk(S_sub);
									tix=streampoi(S_sub,'bconfluences','ix');
									tix=ismember(ST_sub.IXgrid,tix);
									ds=ST_sub.distance;
									ds(~tix)=NaN;
									[~,tix]=max(ds);
									SupT_sub=modify(S_sub,'tributaryto',ST_sub);
									cons=streampoi(SupT_sub,'outlets','xy');
									cons_ix=streampoi(SupT_sub,'outlets','ix');
									cons_ix=vertcat(cons_ix,ST_sub.IXgrid(tix));
									xx=cons(:,1); xx=vertcat(xx,ST_sub.x(tix));
									yy=cons(:,2); yy=vertcat(yy,ST_sub.y(tix));
									da_cons=DAG.Z(cons_ix);
									da_idx=da_cons>=mbz;
									x=vertcat(x,xx(da_idx));
									y=vertcat(y,yy(da_idx));
								end
							end
							rec_count=rec_count+1;
							if rec_count>10
								if isdeployed
									warndlg(['流域编号 ' num2str(basin_num) ' 的分割提前终止以防止无限循环'])
								end
								warning(['流域编号 ' num2str(basin_num) ' 的分割提前终止以防止无限循环']);
							end
						end
						num_new_basins=numel(x);
					catch
						if isdeployed
							warndlg(['流域编号 ' num2str(basin_num) ' 的递归分割失败，将进行常规分割'])
						end
						warning(['流域编号 ' num2str(basin_num) ' 的递归分割失败，将进行常规分割']);
					end
				end

			end

			% 嵌套进度条初始化
			w2=waitbar(0,['正在处理 ' num2str(num_new_basins) ' 个新流域']);
			pos_w1=get(w1,'position');
			pos_w2=[pos_w1(1) pos_w1(2)-pos_w1(4) pos_w1(3) pos_w1(4)];
			set(w2,'position',pos_w2,'doublebuffer','on');

			waitbar(ii/num_files,w1,['正在分割流域编号 ' num2str(basin_num)]);

			% 逐个处理新生成的子流域
			for jj=1:num_new_basins
				waitbar(jj/num_new_basins,w2,['处理中：第 ' num2str(jj) ' 个/共 ' num2str(num_new_basins) ' 个子流域']);

				% 获取当前出水口坐标
				xx=x(jj);
				yy=y(jj);
				% 生成唯一流域编号
				basin_string=sprintf([num2str(basin_num) '%03d'],jj);
				RiverMouth=[xx yy str2num(basin_string)];

				% 构建依赖关系图并裁剪流域
				I=dependencemap(FD,xx,yy);
				DEMoc=crop(DEM,I,nan);
				DEMcc=crop(DEMhc,I,nan);
				FDc=crop(FD,I);
				Ac=crop(A,I,nan);

				% 计算流域面积
				dep_map=GRIDobj2mat(I);
				num_pix=sum(sum(dep_map));
				drainage_area=(num_pix*DEMoc.cellsize*DEMoc.cellsize)/(1e6);

				% 计算高程面积曲线
				[rb,eb]=hypscurve(DEMoc,100);
				hyps=[rb eb];

				% 计算流域质心
				[Cx,Cy]=FindCentroid(DEMoc);
				Centroid=[Cx Cy];

				% 生成新河网
				Sc=STREAMobj(FDc,'minarea',threshold_area,'unit','mapunits');

				% 检查河网是否为空
				if isempty(Sc.x)
					if isdeployed
						warndlg(['流域 ' num2str(RiverMouth(:,3)) ' 的阈值面积过大，正在自动调整'])
					end
					warning(['流域 ' num2str(RiverMouth(:,3)) ' 的阈值面积过大，正在自动调整']);
					new_thresh=threshold_area;
					while isempty(Sc.x)
						new_thresh=new_thresh/2;
						Sc=STREAMobj(FDc,'minarea',new_thresh,'unit','mapunits');
					end
				end

				% 计算Chi值并创建Chi网格
				Cc=chitransform(Sc,Ac,'a0',1,'mn',theta_ref);
				ChiOBJc=GRIDobj(DEMoc);
				ChiOBJc.Z(Sc.IXgrid)=Cc;

				% 计算坡度
				switch gradient_method
				case 'gradient8'
					Goc=gradient8(DEMoc);
				case 'arcslope'
					Goc=arcslope(DEMoc);
                end
				% 寻找最佳拟合凹度
				SLc=klargestconncomps(Sc,1);
				Chic=chiplot(SLc,DEMcc,Ac,'a0',1,'plot',false);

				% 计算ksn
				switch ksn_method
				case 'quick'
					[MSc]=KSN_Quick(DEMoc,DEMcc,Ac,Sc,Chic.mn,segment_length);
					[MSNc]=KSN_Quick(DEMoc,DEMcc,Ac,Sc,theta_ref,segment_length);
				case 'trunk'
					[MSc]=KSN_Trunk(DEMoc,DEMcc,Ac,Sc,Chic.mn,segment_length,min_order);
					[MSNc]=KSN_Trunk(DEMoc,DEMcc,Ac,Sc,theta_ref,segment_length,min_order);
				case 'trib'
					% 对于非常小的流域，覆盖选择，因为KSN_Trib在小流域会失败
					if drainage_area>2.5
						[MSc]=KSN_Trib(DEMoc,DEMcc,FDc,Ac,Sc,Chic.mn,segment_length);
						[MSNc]=KSN_Trib(DEMoc,DEMcc,FDc,Ac,Sc,theta_ref,segment_length);
					else
						[MSc]=KSN_Quick(DEMoc,DEMcc,Ac,Sc,Chic.mn,segment_length);
						[MSNc]=KSN_Quick(DEMoc,DEMcc,Ac,Sc,theta_ref,segment_length);
					end
				end

				% 计算流域范围的ksn统计量
				min_ksn=min([MSNc.ksn],[],'omitnan');
				mean_ksn=mean([MSNc.ksn],'omitnan');
				max_ksn=max([MSNc.ksn],[],'omitnan');
				std_ksn=std([MSNc.ksn],'omitnan');
				se_ksn=std_ksn/sqrt(numel(MSNc)); % 标准误差

				% 计算流域范围的坡度统计量
				min_grad=min(Goc.Z(:),[],'omitnan');
				mean_grad=mean(Goc.Z(:),'omitnan');
				max_grad=max(Goc.Z(:),[],'omitnan');
				std_grad=std(Goc.Z(:),'omitnan');
				se_grad=std_grad/sqrt(sum(~isnan(Goc.Z(:)))); % 标准误差

				% 计算流域范围的高程统计量
				min_z=min(DEMoc.Z(:),[],'omitnan');
				mean_z=mean(DEMoc.Z(:),'omitnan');
				max_z=max(DEMoc.Z(:),[],'omitnan');
				std_z=std(DEMoc.Z(:),'omitnan');
				se_z=std_z/sqrt(sum(~isnan(DEMoc.Z(:)))); % 标准误差

				KSNc_stats=[mean_ksn se_ksn std_ksn min_ksn max_ksn];
				Gc_stats=double([mean_grad se_grad std_grad min_grad max_grad]);
				Zc_stats=double([mean_z se_z std_z min_z max_z]);

				% 寻找出口高程
				out_ix=coord2ind(DEMoc,xx,yy);
				out_el=double(DEMoc.Z(out_ix));

				SubFileName=fullfile(sb_path,['Basin_' num2str(basin_num) '_DataSubset_' num2str(jj) '.mat']);

				save(SubFileName,'RiverMouth','DEMcc','DEMoc','out_el','drainage_area','hyps','FDc','Ac','Sc','SLc','Chic','Goc','MSc','MSNc','KSNc_stats','Gc_stats','Zc_stats','Centroid','ChiOBJc','ksn_method','gradient_method','theta_ref','-v7.3');

				if strcmp(ksn_method,'trunk')
					save(SubFileName,'min_order','-append');
				end
				
				% 创建插值的ksn网格
				if ~isempty(radius)
					try 
						[KsnOBJc] = KsnAvg(DEMoc,MSNc,radius);
						save(SubFileName,'KsnOBJc','radius','-append');
					catch
						if isdeployed
							warndlg(['KSN网格插值失败于盆地 ' num2str(RiverMouth(:,3))])
						end
						warning(['KSN网格插值失败于盆地 ' num2str(RiverMouth(:,3))]);
						save(SubFilename,'radius','-append');
					end
				else
					save(SubFileName,'radius','-append');
				end

				VarList=whos('-file',FileName);

				VarInd=find(strcmp(cellstr(char(VarList.name)),'KSNQc_stats'));
				if ~isempty(VarInd)
					% 提取降水加权的汇流累积量
					load(FileName,'WAc');
					WA=WAc;
					WAc=crop(WA,I,nan);
					switch ksn_method
					case 'quick'
						[WMSc]=KSN_Quick(DEMoc,DEMcc,WAc,Sc,Chic.mn,segment_length);
						[WMSNc]=KSN_Quick(DEMoc,DEMcc,WAc,Sc,theta_ref,segment_length);
					case 'trunk'
						[WMSc]=KSN_Trunk(DEMoc,DEMcc,WAc,Sc,Chic.mn,segment_length,min_order);
						[WMSNc]=KSN_Trunk(DEMoc,DEMcc,WAc,Sc,theta_ref,segment_length,min_order);			
					case 'trib'
						% 小流域覆盖选择
						if drainage_area>2.5
							[WMSc]=KSN_Trib(DEMoc,DEMcc,FDc,WAc,Sc,Chic.mn,segment_length);
							[WMSNc]=KSN_Trib(DEMoc,DEMcc,FDc,WAc,Sc,theta_ref,segment_length);
						else
							[WMSc]=KSN_Quick(DEMoc,DEMcc,WAc,Sc,Chic.mn,segment_length);
							[WMSNc]=KSN_Quick(DEMoc,DEMcc,WAc,Sc,theta_ref,segment_length);
						end
					end

					% 计算流域范围的水文加权ksn统计量
					min_ksnq=min([WMSNc.ksn],[],'omitnan');
					mean_ksnq=mean([WMSNc.ksn],'omitnan');
					max_ksnq=max([WMSNc.ksn],[],'omitnan');
					std_ksnq=std([WMSNc.ksn],'omitnan');
					se_ksnq=std_ksnq/sqrt(numel(WMSNc)); % 标准误差

					KSNQc_stats=[mean_ksnq se_ksnq std_ksnq min_ksnq max_ksnq];
					save(SubFileName,'KSNQc_stats','-append');
				end	
				
				VarInd=find(strcmp(cellstr(char(VarList.name)),'AGc'), 1);
				if ~isempty(VarInd)
					load(FileName,'AGc');
					AG=AGc;
					num_grids=size(AG,1);
					AGc=cell(size(AG));
					for kk=1:num_grids
						AGcOI=crop(AG{kk,1},I,nan);
						AGc{kk,1}=AGcOI;
						AGc{kk,2}=AG{kk,2};
						mean_AGc=mean(AGcOI.Z(:),'omitnan');
						min_AGc=min(AGcOI.Z(:),[],'omitnan');
						max_AGc=max(AGcOI.Z(:),[],'omitnan');
						std_AGc=std(AGcOI.Z(:),'omitnan');
						se_AGc=std_AGc/sqrt(sum(~isnan(AGcOI.Z(:))));
						AGc_stats(kk,:)=[mean_AGc se_AGc std_AGc min_AGc max_AGc];
					end
					save(SubFileName,'AGc','AGc_stats','-append');
				end

				VarInd=find(strcmp(cellstr(char(VarList.name)),'ACGc'), 1);
				if ~isempty(VarInd)
					load(FileName,'ACGc');
					ACG=ACGc;
					num_grids=size(ACG,1);
					ACGc=cell(size(ACG));
					for kk=1:num_grids
						ACGcOI=crop(ACG{kk,1},I,nan);
						ACGc{kk,1}=ACGcOI;
						ACGc{kk,3}=ACG{kk,3};
						edg=ACG{kk,2}.Numbers;
						edg=edg+0.5;
						edg=vertcat(0.5,edg);
						[N,~]=histcounts(ACGcOI.Z(:),edg);
						T=ACG{kk,2};
						T.Counts=N';
						ACGc{kk,2}=T;
						ACGc_stats(kk,1)=[mode(ACGcOI.Z(:))];
					end
					save(SubFileName,'ACGc','ACGc_stats','-append');	
				end	

				VarInd=find(strcmp(cellstr(char(VarList.name)),'rlf'), 1);
				if ~isempty(VarInd)
					load(FileName,'rlf');
					rlf_full=rlf; 
					num_rlf=size(rlf_full,1);
					rlf=cell(size(rlf_full));
					rlf_stats=zeros(num_rlf,6);
					for kk=1:num_rlf
						% 计算地形起伏
						radOI=rlf_full{kk,2};
						rlf{kk,2}=radOI;
						rlfOI=localtopography(DEMoc,radOI);
						rlf{kk,1}=rlfOI;
						% 计算统计量
						mean_rlf=mean(rlfOI.Z(:),'omitnan');
						min_rlf=min(rlfOI.Z(:),[],'omitnan');
						max_rlf=max(rlfOI.Z(:),[],'omitnan');
						std_rlf=std(rlfOI.Z(:),'omitnan');
						se_rlf=std_rlf/sqrt(sum(~isnan(rlfOI.Z(:))));
						rlf_stats(kk,:)=[mean_rlf se_rlf std_rlf min_rlf max_rlf radOI];
					end
					save(SubFileName,'rlf','rlf_stats','-append');
				end					

				if write_arc_files
					% 将DEM中的NaN替换为-9999
					Didx=isnan(DEMoc.Z);
					DEMoc_temp=DEMoc;
					DEMoc_temp.Z(Didx)=-9999;

					DEMFileName=fullfile(sb_path,['Basin_' num2str(basin_num) '_DataSubset_' num2str(jj) '_DEM.txt']);
					GRIDobj2ascii(DEMoc_temp,DEMFileName);
					CHIFileName=fullfile(sb_path,['Basin_' num2str(basin_num) '_DataSubset_' num2str(jj) '_CHI.txt']);
					GRIDobj2ascii(ChiOBJc,CHIFileName);
					KSNFileName=fullfile(sb_path,['Basin_' num2str(basin_num) '_DataSubset_' num2str(jj) '_KSN.shp']);
					shapewrite(MSNc,KSNFileName);

						for kk=1:num_rlf
							RLFFileName=fullfile(sb_path,['Basin_' num2str(basin_num) '_DataSubset_' num2str(jj) '_RLF_' num2str(rlf{kk,2}) '.txt']);
							GRIDobj2ascii(rlf{kk,1},RLFFileName);
						end


				 if exist('AG', 'var') && ~isempty(AG)
                        for kk = 1:num_grids
                          AGcFileName = fullfile(['Basin_' num2str(basin_num) '_DataSubset_' num2str(jj) '_' AGc{kk,2} '.txt']);
                          GRIDobj2ascii(AGc{kk,1}, AGcFileName);
                        end
                    end

                if exist('ACG', 'var') && ~isempty(ACG)
                    for jj = 1:num_grids
                        ACGcFileName = fullfile(sb_path, ['Basin_' num2str(basin_num) '_' ACGc{jj,3} '.txt']);
                        GRIDobj2ascii(ACGc{jj,1}, ACGcFileName);
                    end
                end

                
				end
			end % 新盆地循环结束
			close(w2);
		end % 流域面积检查结束
	end % 主循环结束
	close(w1);
end % 主函数结束

function [ksn_ms]=KSN_Quick(DEM,DEMc,A,S,theta_ref,segment_length)
	g=gradient(S,DEMc);
	G=GRIDobj(DEM);
	G.Z(S.IXgrid)=g;

	Z_RES=DEMc-DEM;

	ksn=G./(A.*(A.cellsize^2)).^(-theta_ref);

	SD=GRIDobj(DEM);
	SD.Z(S.IXgrid)=S.distance;
	
	ksn_ms=STREAMobj2mapstruct(S,'seglength',segment_length,'attributes',...
		{'ksn' ksn @mean 'uparea' (A.*(A.cellsize^2)) @mean 'gradient' G @mean 'cut_fill' Z_RES @mean...
		'min_dist' SD @min 'max_dist' SD @max});

	seg_dist=[ksn_ms.max_dist]-[ksn_ms.min_dist];
	distcell=num2cell(seg_dist');
	[ksn_ms(1:end).seg_dist]=distcell{:};
	ksn_ms=rmfield(ksn_ms,{'min_dist','max_dist'});
end

function [ksn_ms]=KSN_Trunk(DEM,DEMc,A,S,theta_ref,segment_length,min_order)

	order_exp=['>=' num2str(min_order)];

    Smax=modify(S,'streamorder',order_exp);
	Smin=modify(S,'rmnodes',Smax);

	g=gradient(S,DEMc);
	G=GRIDobj(DEM);
	G.Z(S.IXgrid)=g;

	Z_RES=DEMc-DEM;

	ksn=G./(A.*(A.cellsize^2)).^(-theta_ref);

	SDmax=GRIDobj(DEM);
	SDmin=GRIDobj(DEM);
	SDmax.Z(Smax.IXgrid)=Smax.distance;
	SDmin.Z(Smin.IXgrid)=Smin.distance;

	ksn_ms_min=STREAMobj2mapstruct(Smin,'seglength',segment_length,'attributes',...
		{'ksn' ksn @mean 'uparea' (A.*(A.cellsize^2)) @mean 'gradient' G @mean 'cut_fill' Z_RES @mean...
		'min_dist' SDmin @min 'max_dist' SDmin @max});

	ksn_ms_max=STREAMobj2mapstruct(Smax,'seglength',segment_length,'attributes',...
		{'ksn' ksn @mean 'uparea' (A.*(A.cellsize^2)) @mean 'gradient' G @mean 'cut_fill' Z_RES @mean...
		'min_dist' SDmax @min 'max_dist' SDmax @max});

	ksn_ms=vertcat(ksn_ms_min,ksn_ms_max);
	seg_dist=[ksn_ms.max_dist]-[ksn_ms.min_dist];
	distcell=num2cell(seg_dist');
	[ksn_ms(1:end).seg_dist]=distcell{:};
	ksn_ms=rmfield(ksn_ms,{'min_dist','max_dist'});
end

function [ksn_ms]=KSN_Trib(DEM,DEMc,FD,A,S,theta_ref,segment_length)

	% 定义不交叉的河段
	[as]=networksegment_slim(DEM,FD,S);
	seg_bnd_ix=as.ix;
	% 预计算后续需要的值
	z=getnal(S,DEMc);
	zu=getnal(S,DEM);
	z_res=z-zu;
	g=gradient(S,DEMc);
	c=chitransform(S,A,'a0',1,'mn',theta_ref);
	d=S.distance;
	da=getnal(S,A.*(A.cellsize^2));
	ixgrid=S.IXgrid;
	% 提取有序的河道节点列表
	s_node_list=S.orderednanlist;
	streams_ix=find(isnan(s_node_list));
	streams_ix=vertcat(1,streams_ix);
	% 生成空节点属性列表
	ksn_nal=zeros(size(d));
	% 主循环遍历河道
	num_streams=numel(streams_ix)-1;
	seg_count=1;
	for ii=1:num_streams
		if ii==1
			snlOI=s_node_list(streams_ix(ii):streams_ix(ii+1)-1);
		else
			snlOI=s_node_list(streams_ix(ii)+1:streams_ix(ii+1)-1);
		end

		[~,~,dn]=intersect(snlOI,seg_bnd_ix(:,1));
		[~,~,up]=intersect(snlOI,seg_bnd_ix(:,2));
		seg_ix=intersect(up,dn);

		num_segs=numel(seg_ix);
		dn_up=seg_bnd_ix(seg_ix,:);
		for jj=1:num_segs
			dnix=find(snlOI==dn_up(jj,1));
			upix=find(snlOI==dn_up(jj,2));
			seg_ix_oi=snlOI(upix:dnix);
			dOI=d(seg_ix_oi);
			dnOI=dOI-min(dOI);
			num_bins=ceil(max(dnOI)/segment_length);
			bin_edges=[0:segment_length:num_bins*segment_length];
			for kk=1:num_bins
				idx=dnOI>bin_edges(kk) & dnOI<=bin_edges(kk+1);
				bin_ix=seg_ix_oi(idx);
				cOI=c(bin_ix);
				zOI=z(bin_ix);
				if numel(cOI)>2
					[ksn_val,r2]=Chi_Z_Spline(cOI,zOI);
					ksn_nal(bin_ix)=ksn_val;

					% 构建地图结构
					ksn_ms(seg_count).Geometry='Line';
					ksm_ms(seg_count).BoundingBox=[min(S.x(bin_ix)),min(S.y(bin_ix));max(S.x(bin_ix)),max(S.y(bin_ix))];
					ksn_ms(seg_count).X=S.x(bin_ix);
					ksn_ms(seg_count).Y=S.y(bin_ix);
					ksn_ms(seg_count).ksn=ksn_val;
					ksn_ms(seg_count).uparea=mean(da(bin_ix));
					ksn_ms(seg_count).gradient=mean(g(bin_ix));
					ksn_ms(seg_count).cut_fill=mean(z_res(bin_ix));
					ksn_ms(seg_count).seg_dist=max(S.distance(bin_ix))-min(S.distance(bin_ix));
					ksn_ms(seg_count).chi_r2=r2;
					
					seg_count=seg_count+1;
				end
			end
		end
	end
end

function seg = networksegment_slim(DEM,FD,S)
	% 'networksegment'函数的精简版本，移除了零长度和单节点长度的河段

	%% 识别河道起点、汇合点、分支点和出口
	Vhead = streampoi(S,'channelheads','logical');  ihead=find(Vhead==1);  IXhead=S.IXgrid(ihead);
	Vconf = streampoi(S,'confluences','logical');   iconf=find(Vconf==1);  IXconf=S.IXgrid(iconf);
	Vout = streampoi(S,'outlets','logical');        iout=find(Vout==1);    IXout=S.IXgrid(iout);
	Vbconf = streampoi(S,'bconfluences','logical'); ibconf=find(Vbconf==1);IXbconf=S.IXgrid(ibconf);

	%% 识别关联流域
	DB   = drainagebasins(FD,vertcat(IXbconf,IXout));DBhead=DB.Z(IXhead); DBbconf=DB.Z(IXbconf); DBconf=DB.Z(IXconf); DBout=DB.Z(IXout);

	%% 计算流程距离
	D = flowdistance(FD);

	%% 识别河段
	[~,ind11,ind12]=intersect(DBbconf,DBhead);
	[~,ind21,ind22]=intersect(DBbconf,DBconf);
	[~,ind31,ind32]=intersect(DBout,DBhead);
	[~,ind41,ind42]=intersect(DBout,DBconf);
	IX(:,1) = [ IXbconf(ind11)' IXbconf(ind21)' IXout(ind31)'  IXout(ind41)'  ];   ix(:,1)= [ ibconf(ind11)' ibconf(ind21)' iout(ind31)'  iout(ind41)'  ];
	IX(:,2) = [ IXhead(ind12)'  IXconf(ind22)'  IXhead(ind32)' IXconf(ind42)' ];   ix(:,2)= [ ihead(ind12)'  iconf(ind22)'  ihead(ind32)' iconf(ind42)' ];

	% 计算河段长度
	flength=double(abs(D.Z(IX(:,1))-D.Z(IX(:,2))));

	% 移除零长度和单节点河段
	idx=flength>=2*DEM.cellsize;
	seg.IX=IX(idx,:);
	seg.ix=ix(idx,:);
	seg.flength=flength(idx);

	seg.n=numel(IX(:,1));
	end

function [KSN,R2] = Chi_Z_Spline(c,z)

	% 使用三次样条插值重采样Chi-高程关系
	[~,minIX]=min(c);
	zb=z(minIX);
	chiF=c-min(c);
	zabsF=z-min(z);
	chiS=linspace(0,max(chiF),numel(chiF)).';
	zS=spline(chiF,zabsF,chiS);

	% 通过斜率计算ksn
	KSN= chiS\(zS); % 因a0固定为1，无需mn参数

	% 计算R平方
	z_pred=chiF.*KSN;
	sstot=sum((zabsF-mean(zabsF)).^2);
	ssres=sum((zabsF-z_pred).^2);
	R2=1-(ssres/sstot);

end

function [KSNGrid] = KsnAvg(DEM,ksn_ms,radius)

	% 计算像素半径
	radiuspx = ceil(radius/DEM.cellsize);

	% 记录当前NaN掩膜
	MASK=isnan(DEM.Z);

	% 创建河道值网格
	KSNGrid=GRIDobj(DEM);
	KSNGrid.Z(:,:)=NaN;
	for ii=1:numel(ksn_ms)
		ix=coord2ind(DEM,ksn_ms(ii).X,ksn_ms(ii).Y);
		KSNGrid.Z(ix)=ksn_ms(ii).ksn;
	end

	% 基于半径的局部均值
	ISNAN=isnan(KSNGrid.Z);
    [~,L] = bwdist(~ISNAN,'e');
    ksng = KSNGrid.Z(L);           
    FLT   = fspecial('disk',radiuspx);
    ksng   = imfilter(ksng,FLT,'symmetric','same','conv');

    % 恢复原始NaN区域
    ksng(MASK)=NaN;

    KSNGrid.Z=ksng;
end

function [x,y] = CheckUpstream(DEM,FD,ix)
	% 构建影响区域单元列表
	inflcs=cell(numel(ix),1);
	for ii=1:numel(ix)
	    IX=influencemap(FD,ix(ii));
	    inflcs{ii}=find(IX.Z);
	end
	    
	% 构建索引
	idx=zeros(numel(ix),1);
	idx=logical(idx);
	for ii=1:numel(ix)
	    inflcs_temp=inflcs;
	    inflcs_temp{ii}=[0];
	    up_member=cellfun(@(x) ismember(ix(ii),x),inflcs_temp);
	    if any(up_member)
	        idx(ii)=false;
	    else
	        idx(ii)=true;
	    end
	end

	[x,y]=ind2coord(DEM,ix(idx));
end