function [links]=JunctionLinks(FD,S,IX,junctions,varargin)
% 用法：
%   [links] = JunctionLinks(FD, S, IX, junctions);
%   [links] = JunctionLinks(FD, S, IX, junctions, 'name', value);
%
% 描述：
%   该函数用于识别流域网络中的链接。结果与networksegments函数相似（但不完全相同）。
%   链接按照James & Krumbein（1969）的方法进行分类，具体参考"Frequency distributions of stream link lengths"，
%   将链接分为内链接和外链接（外链接的上游端为水头），以及顺向（cis）和反向（trans）链接（顺向链接的两端具有相同的手性，
%   反向链接的两端手性相反）。外链接、连接到出口的内链接、一个或多个端点的交汇点有超过两个上游链接的链接，以及手性不明确的内链接（即进入流的Shreve等级相同的链接）将在cis或trans分类上未定义。
%
% 必需的输入：
%   FD - FLOWobj
%   S - STREAMobj
%   IX - JunctionAngle的IX输出
%   junctions - JunctionAngle的交汇点表输出，如果您为JunctionAngle提供了多个fit_distance，
%                则输出的'links'将是与输入'junctions'相似维度的单元数组。
%
% 可选输入：
%   make_shape [false] - 逻辑标志，用于生成一个包含由交汇点定义的流域网络分段的shapefile，并根据link_type和link_side分类（见输出）。
%                        生成用于创建shapefile的mapstructure可能非常耗时。
%   par [true] - 逻辑标志，用于并行处理创建mapstructure，因为这是一个耗时的过程。此标志仅在'make_shape'设置为true时有效。
%                如果您没有并行计算工具箱的许可，函数将自动将此标志设置为false。
%
% 输出：
%   links - 表格，包含所有流域链接的信息：
%     link_number - 链接的ID号
%     downstream_x - 链接下游端的x坐标
%     downstream_y - 链接下游端的y坐标
%     downstream_IX - 链接下游端在GRIDobj中的索引
%     upstream_x - 链接上游端的x坐标
%     upstream_y - 链接上游端的y坐标
%     upstream_IX - 链接上游端在GRIDobj中的索引
%     link_type - 链接分类，可能为Exterior（链接的上游端为水头）或Interior
%     downstream_junc_num - 链接下游端的交汇点编号，引用提供给J的交汇点表，NaN表示下游端不是交汇点
%     upstream_junc_num - 链接上游端的交汇点编号，引用提供给J的交汇点表，NaN表示上游端不是交汇点
%     link_side - 链接分类，可能为Cis（链接两端具有相同的手性）、Trans（链接两端的手性相反）或Undefined（如果链接的两端没有明确的手性分类）
%
% 示例：
%   [links] = JunctionLinks(FD, S, IX, junctions);
%   [links] = JunctionLinks(FD, S, IX, junctions, 'make_shape', true);
%
% 相关函数：
%   JunctionAngle, InspectJunction, networksegment
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数编写者：Adam M. Forte - 更新日期：06/27/19 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'JunctionLinks';
	addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
	addRequired(p,'S',@(x) isa(x,'STREAMobj'));
	addRequired(p,'IX',@(x) iscell(x));
	addRequired(p,'junctions',@(x) istable(x) | iscell(x));

	addParameter(p,'make_shape',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'par',true,@(x) isscalar(x) && islogical(x));
	addParameter(p,'out_dir',[],@(x) isdir(x));

	parse(p,FD,S,IX,junctions,varargin{:});
	FD=p.Results.FD;
	S=p.Results.S;
	IX=p.Results.IX;
	JIN=p.Results.junctions;

	make_shape=p.Results.make_shape;
	par=p.Results.par;
	out_dir=p.Results.out_dir;

	if make_shape
		if isempty(out_dir)
			out_dir=pwd;
		end
	end

	% Determine type of J input
	if istable(JIN)
		J=JIN;

		% Grab handedness
		handedness=J.handedness;	

		% Find number of juctions (valid and invalid)
		num_con=size(IX,1);
		con=IX(:,1);
		% Determine junctions with invalid downstream nodes
		vld=~cellfun(@isempty,IX(:,2));

		% Generate original junction position list
		olist=1:num_con;
		olist=olist(vld)';

		% Populate lists of down- and up-stream indices
		dsix=zeros(num_con,1);
		us1ix=zeros(num_con,1);
		us2ix=zeros(num_con,1);

		for ii=1:nnz(vld)
			% Find nodes up and downstream of junctions
			JOI=IX(olist(ii),:);
			dsix(olist(ii),1)=JOI{1,2}(1);
			us1ix(olist(ii),1)=JOI{1,3}(1);
			us2ix(olist(ii),1)=JOI{1,4}(1);
		end

		%% Build links - Modified from networksegment.m
		% Build logical of head, up-stream confluences, and outlet
		% position in nal. 
		headIDX = streampoi(S,'channelheads','logical');      
		outIDX = streampoi(S,'outlets','logical');            
		bconfIDX = streampoi(S,'bconfluences','logical'); 

		% Find positions in nal of heads, up-stream confluences, down-stream confluences
		% and outlets
		ihead=find(headIDX==1);
		iout=find(outIDX==1);
		ibconf=find(bconfIDX==1);
		iconf=dsix(dsix>0);

		% Find GRID indices of heads, up-stream confluences, down-stream confluences and
		% outlets
		IXhead=S.IXgrid(ihead);
		IXout=S.IXgrid(iout);
		IXbconf=S.IXgrid(ibconf);
		IXconf=S.IXgrid(iconf);

		% Build drainage basin object establishing links
		DB=drainagebasins(FD,vertcat(IXbconf,IXout));

		% Extract link numbers for heads, up-stream confluences, down-stream confluences
		% and outlets
		DBhead=DB.Z(IXhead); 
		DBout=DB.Z(IXout);
		DBbconf=DB.Z(IXbconf); 
		DBconf=DB.Z(IXconf); 

		% Find links between channel heads and up-stream confluences
		[~,ind11,ind12]=intersect(DBbconf,DBhead);
		% Find links between down-stream confluences and up-stream confluencres
		[~,ind21,ind22]=intersect(DBbconf,DBconf);
		% Find inks between channel heads and outlets
		[~,ind31,ind32]=intersect(DBout,DBhead);
		% Find links between down-stream confluences and outlets
		[~,ind41,ind42]=intersect(DBout,DBconf);

		% Connect nodes into links
		ix(:,1)= [ ibconf(ind11)' ibconf(ind21)' iout(ind31)'  iout(ind41)'  ]; % Downstream end of link
		ix(:,2)= [ ihead(ind12)'  iconf(ind22)'  ihead(ind32)' iconf(ind42)' ]; % Upstream end of link

		% Build link table
		num_links=numel(ix(:,1));
		link_number=1:num_links;
		link_number=link_number(:);

		downstream_x=S.x(ix(:,1));
		downstream_y=S.y(ix(:,1));

		upstream_x=S.x(ix(:,2));
		upstream_y=S.y(ix(:,2));

		downstream_IX=S.IXgrid(ix(:,1));
		upstream_IX=S.IXgrid(ix(:,2));

		% Allocate arrays for junction numbers
		% downstream_junction=zeros(num_links,1);
		upstream_junc_num=zeros(num_links,1);
		downstream_junc_num=zeros(num_links,1);
		% Allocate link_side
		link_side=categorical(zeros(num_links,1));
		% Precompute data for determining interior vs exterior links
		inex=ismember(ix(:,2),ihead);
		link_type=categorical(zeros(num_links,1));
		for ii=1:num_links

			% Reset logical flags
			down1=false;
			down2=false;
			up=false;

			% Extract indices defining link
			down_ix=ix(ii,1);
			up_ix=ix(ii,2);

			% Label interior vs exterior links
			if inex(ii)
				link_type(ii,1)='Exterior';
				upstream_junc_num(ii,1)=NaN;
				link_side(ii,1)='Undefined';
			else
				link_type(ii,1)='Interior';
				% Find upstream junction number
				pos_in_ds=find(dsix==up_ix);
				if isempty(pos_in_ds)
					upstream_junc_num(ii,1)=NaN;
					link_side(ii,1)='Undefined';
				else
					upstream_junc_num(ii,1)=pos_in_ds;
					up=true;
				end

			end

			% Find downstream junction number
			pos_in_us1=find(us1ix==down_ix);
			pos_in_us2=find(us2ix==down_ix);
			if isempty(pos_in_us1) & isempty(pos_in_us2)
				downstream_junc_num(ii,1)=NaN;
				link_side(ii,1)='Undefined';
			elseif isempty(pos_in_us1)
				downstream_junc_num(ii,1)=pos_in_us2;
				down2=true;
			elseif isempty(pos_in_us2)
				downstream_junc_num(ii,1)=pos_in_us1;
				down1=true;
			else
				downstream_junc_num(ii,1)=NaN;
				link_side(ii,1)='Undefined';
			end

			% Determine link_type if applicable
			if up & down1
				hu=handedness(pos_in_ds);
				hd=handedness(pos_in_us1);
				if hu==hd & hu~='Undefined' & hd~='Undefined'
					link_side(ii,1)='Cis';
				elseif hu~=hd & hu~='Undefined' & hd~='Undefined'
					link_side(ii,1)='Trans';
				else 
					link_side(ii,1)='Undefined';
				end

			elseif up & down2
				hu=handedness(pos_in_ds);
				hd=handedness(pos_in_us2);

				if hu==hd & hu~='Undefined' & hd~='Undefined'
					link_side(ii,1)='Cis';
				elseif hu~=hd & hu~='Undefined' & hd~='Undefined'
					link_side(ii,1)='Trans';
				else
					link_side(ii,1)='Undefined';
				end

			end			
		end

		% Compile into link table
		links=table(link_number,downstream_x,downstream_y,downstream_IX,...
			upstream_x,upstream_y,upstream_IX,...
			link_type,downstream_junc_num,upstream_junc_num,link_side);

		if make_shape
			ms=makelinkshape(S,links,par);
			disp('Writing shapefile');
			shapewrite(ms,fullfile(out_dir,'junction_links.shp'));
		end

	else
		
		% Create links output cell
		links=cell(size(JIN));
		num_lengths=size(JIN,2);
		for jj=1:num_lengths;
			links{1,jj}=JIN{1,jj};
		end

		for jj=1:num_lengths

			% Grab junction table
			J=JIN{2,jj};

			% Grab handedness
			handedness=J.handedness;	

			% Find number of juctions (valid and invalid)
			num_con=size(IX,1);
			con=IX(:,1);
			% Determine junctions with invalid downstream nodes
			vld=~cellfun(@isempty,IX(:,2));

			% Generate original junction position list
			olist=1:num_con;
			olist=olist(vld)';

			% Populate lists of down- and up-stream indices
			dsix=zeros(num_con,1);
			us1ix=zeros(num_con,1);
			us2ix=zeros(num_con,1);

			for ii=1:nnz(vld)
				% Find nodes up and downstream of junctions
				JOI=IX(olist(ii),:);
				dsix(olist(ii),1)=JOI{1,2}(1);
				us1ix(olist(ii),1)=JOI{1,3}(1);
				us2ix(olist(ii),1)=JOI{1,4}(1);
			end

			%% Build links - Modified from networksegment.m
			% Build logical of head, up-stream confluences, and outlet
			% position in nal. 
			headIDX = streampoi(S,'channelheads','logical');      
			outIDX = streampoi(S,'outlets','logical');            
			bconfIDX = streampoi(S,'bconfluences','logical'); 

			% Find positions in nal of heads, up-stream confluences, down-stream confluences
			% and outlets
			ihead=find(headIDX==1);
			iout=find(outIDX==1);
			ibconf=find(bconfIDX==1);
			iconf=dsix(dsix>0);

			% Find GRID indices of heads, up-stream confluences, down-stream confluences and
			% outlets
			IXhead=S.IXgrid(ihead);
			IXout=S.IXgrid(iout);
			IXbconf=S.IXgrid(ibconf);
			IXconf=S.IXgrid(iconf);

			% Build drainage basin object establishing links
			DB=drainagebasins(FD,vertcat(IXbconf,IXout));

			% Extract link numbers for heads, up-stream confluences, down-stream confluences
			% and outlets
			DBhead=DB.Z(IXhead); 
			DBout=DB.Z(IXout);
			DBbconf=DB.Z(IXbconf); 
			DBconf=DB.Z(IXconf); 

			% Find links between channel heads and up-stream confluences
			[~,ind11,ind12]=intersect(DBbconf,DBhead);
			% Find links between down-stream confluences and up-stream confluencres
			[~,ind21,ind22]=intersect(DBbconf,DBconf);
			% Find inks between channel heads and outlets
			[~,ind31,ind32]=intersect(DBout,DBhead);
			% Find links between down-stream confluences and outlets
			[~,ind41,ind42]=intersect(DBout,DBconf);

			% Connect nodes into links
			ix(:,1)= [ ibconf(ind11)' ibconf(ind21)' iout(ind31)'  iout(ind41)'  ]; % Downstream end of link
			ix(:,2)= [ ihead(ind12)'  iconf(ind22)'  ihead(ind32)' iconf(ind42)' ]; % Upstream end of link

			% Build link table
			num_links=numel(ix(:,1));
			link_number=1:num_links;
			link_number=link_number(:);

			downstream_x=S.x(ix(:,1));
			downstream_y=S.y(ix(:,1));

			upstream_x=S.x(ix(:,2));
			upstream_y=S.y(ix(:,2));

			downstream_IX=S.IXgrid(ix(:,1));
			upstream_IX=S.IXgrid(ix(:,2));

			% Allocate arrays for junction numbers
			% downstream_junction=zeros(num_links,1);
			upstream_junc_num=zeros(num_links,1);
			downstream_junc_num=zeros(num_links,1);
			% Allocate link_side
			link_side=categorical(zeros(num_links,1));
			% Precompute data for determining interior vs exterior links
			inex=ismember(ix(:,2),ihead);
			link_type=categorical(zeros(num_links,1));
			for ii=1:num_links

				% Reset logical flags
				down1=false;
				down2=false;
				up=false;

				% Extract indices defining link
				down_ix=ix(ii,1);
				up_ix=ix(ii,2);

				% Label interior vs exterior links
				if inex(ii)
					link_type(ii,1)='Exterior';
					upstream_junc_num(ii,1)=NaN;
					link_side(ii,1)='Undefined';
				else
					link_type(ii,1)='Interior';
					% Find upstream junction number
					pos_in_ds=find(dsix==up_ix);
					if isempty(pos_in_ds)
						upstream_junc_num(ii,1)=NaN;
						link_side(ii,1)='Undefined';
					else
						upstream_junc_num(ii,1)=pos_in_ds;
						up=true;
					end

				end

				% Find downstream junction number
				pos_in_us1=find(us1ix==down_ix);
				pos_in_us2=find(us2ix==down_ix);
				if isempty(pos_in_us1) & isempty(pos_in_us2)
					downstream_junc_num(ii,1)=NaN;
					link_side(ii,1)='Undefined';
				elseif isempty(pos_in_us1)
					downstream_junc_num(ii,1)=pos_in_us2;
					down2=true;
				elseif isempty(pos_in_us2)
					downstream_junc_num(ii,1)=pos_in_us1;
					down1=true;
				else
					downstream_junc_num(ii,1)=NaN;
					link_side(ii,1)='Undefined';
				end

				% Determine link_type if applicable
				if up & down1
					hu=handedness(pos_in_ds);
					hd=handedness(pos_in_us1);
					if hu==hd & hu~='Undefined' & hd~='Undefined'
						link_side(ii,1)='Cis';
					elseif hu~=hd & hu~='Undefined' & hd~='Undefined'
						link_side(ii,1)='Trans';
					else 
						link_side(ii,1)='Undefined';
					end

				elseif up & down2
					hu=handedness(pos_in_ds);
					hd=handedness(pos_in_us2);

					if hu==hd & hu~='Undefined' & hd~='Undefined'
						link_side(ii,1)='Cis';
					elseif hu~=hd & hu~='Undefined' & hd~='Undefined'
						link_side(ii,1)='Trans';
					else
						link_side(ii,1)='Undefined';
					end

				end			
			end

			% Compile into link table
			l=table(link_number,downstream_x,downstream_y,downstream_IX,...
				upstream_x,upstream_y,upstream_IX,...
				link_type,downstream_junc_num,upstream_junc_num,link_side);
			links{2,jj}=l;

			if make_shape
				ms=makelinkshape(S,l,par);
				shp_name=fullfile(out_dir,['junction_links_' num2str(jj) '.shp']);
				disp('Writing shapefile');
				shapewrite(ms,shp_name);
			end
		end
	end
% Function End	
end

function [ms]=makelinkshape(S,L,par)
	disp('Starting shapefile construction')
	% Split stream network
	SS=split(S);

	% Make mapstructure based on split & modify
	disp('Building original mapstructure')
	ms=STREAMobj2mapstruct(SS);
	ms=rmfield(ms,{'streamorder','IX','tribtoIX'});
	ms(1,1).link_type=[]; ms(1,1).link_side=[]; ms(1,1).link_num=[];

	%Grab out data of interest from link table
	dx=L.downstream_x; dy=L.downstream_y; ux=L.upstream_x; uy=L.upstream_y;
	lty=L.link_type; lsd=L.link_side; ln=L.link_number;

	if par & license('test','distrib_computing_toolbox')

		[ms,w1]=par_loop_proc(ms,dx,dy,ux,uy,lty,lsd,ln);
		close(w1);
	else

		w1=waitbar(0,'Populating link information into mapstructure');
		for ii=1:numel(ms)
			xl=ms(1,ii).X;
			yl=ms(1,ii).Y;

			fun=@(DXL,UXL,DYL,UYL) any(ismember(DXL,xl)) & any(ismember(UXL,xl)) & any(ismember(DYL,yl)) & any(ismember(UYL,yl)); 
			idx=arrayfun(fun,dx,ux,dy,uy);
			ix=find(idx);

			if numel(ix)==1
				ms(1,ii).link_type=char(lty(ix));
				ms(1,ii).link_side=char(lsd(ix));
				ms(1,ii).link_num=ln(ix);
			elseif numel(ix>1)
				ms(1,ii).link_type=char(lty(ix(1)));
				ms(1,ii).link_side=char(lsd(ix(1)));
				ms(1,ii).link_num=ln(ix(1));
			end

			waitbar(ii/numel(ms));
		end
		close(w1);

	end

end

function [ms,w1]=par_loop_proc(ms,dx,dy,ux,uy,lty,lsd,ln)

	DQ=parallel.pool.DataQueue;
	w1=waitbar(0,'Populating link information into mapstructure');
	afterEach(DQ,@updateBar);
	num_loop=numel(ms);
	cnt=1;

	parfor ii=1:num_loop
		xl=ms(1,ii).X;
		yl=ms(1,ii).Y;

		fun=@(DXL,UXL,DYL,UYL) any(ismember(DXL,xl)) & any(ismember(UXL,xl)) & any(ismember(DYL,yl)) & any(ismember(UYL,yl)); 
		idx=arrayfun(fun,dx,ux,dy,uy);
		ix=find(idx);

		if numel(ix)==1
			ms(1,ii).link_type=char(lty(ix));
			ms(1,ii).link_side=char(lsd(ix));
			ms(1,ii).link_num=ln(ix);
		elseif numel(ix>1)
			ms(1,ii).link_type=char(lty(ix(1)));
			ms(1,ii).link_side=char(lsd(ix(1)));
			ms(1,ii).link_num=ln(ix(1));
		end

		send(DQ,ii);
	end

		function updateBar(~)
			waitbar(cnt/num_loop,w1);
			cnt=cnt+1;
		end
end