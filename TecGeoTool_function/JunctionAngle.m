function [junctions,IX,varargout]=JunctionAngle(S,A,DEM,fit_distance,varargin);
	% 用法：
% [junctions,IX]=JunctionAngle(S,A,DEM,fit_distance);
% [junctions,IX]=JunctionAngle(S,A,DEM,fit_distance,'name',value);
% [junctions,IX,mean_junctions]=JunctionAngle(S,A,DEM,[# # #]);
%
% 描述：
% 该函数用于计算河流交汇处的角度。本函数的基本策略和术语参考自Howard（1971）的《Optimal angles of stream junction:
% Geometric, stability to capture, and minimum power criteria》。代码首先对交汇点的上下游河道进行线性拟合，
% 拟合所用的河段长度由'fit_distance'参数控制。交汇角度通过计算各支流（支流1和支流2）与下游河道上游投影之间的两个角度（e1和e2）来确定。
% 总交汇角度为两个上游河道之间的夹角，在理想情况下应为e1和e2之和，但并非总是如此（例如当两个上游河道位于下游河道上游投影的同一侧时，
% 此时交汇角度并非e1与e2之和，而是支流1与支流2之间的实际夹角，而e1和e2仍从上游投影面测量）。本函数还基于Howard（1971）的
% 几何准则和James & Krumbein（1969）的《Frequency distributions of stream link lengths》中定义的方位特性计算预测交汇角度。
% 需注意，若交汇点有超过两个上游河道汇入或没有下游节点（即出口点），相关角度将被设为NaN（因这些情况在Howard准则中未定义）。
% 若提供了以立方米/秒为单位的流量GRIDobj，则还会使用Howard（1971）的最小功率准则进行预测。
%
% 必需输入：
% S - 用于计算交汇角度的STREAMobj
% A - 流量累积量GRIDobj
% DEM - 高程GRIDobj，用于计算坡度
% fit_distance - 河道距离（地图单位）的拟合长度，
% 若提供多个值，则表示使用不同拟合长度计算交汇角度。若设置'use_n_nodes'为true，
% 该值可解释为用于拟合的节点数量
%
% 可选输入：
% use_n_nodes [false] - 逻辑标志，指示将fit_distance解释为节点数量而非河道距离
% previous_IX [] - 若提供先前计算的IX结果可加速重复计算，需确保S、A、DEM输入一致
% discharge [] - 流量GRIDobj（立方米/秒），用于最小功率准则计算
% predict_angle_method ['area'] - 预测方法选择：'slope','area','minimum_power'或'all'
% ref_concavity [0.5] - 面积法预测时使用的水系凹度参考值
% make_shape [false] - 是否生成包含交汇点信息的shapefile
% verbose [false] - 是否显示计算进度信息
%
% 输出：
% junctions - 结果表格（单拟合长度）或元胞数组（多拟合长度），包含：
% junction_number 交汇点编号
% junction_x/y 交汇点坐标
% junction_angle 总交汇角度
% handedness 方位特性（左/右/未定义）
% split 是否分叉标志
% e1/e2观测角及预测角
% e1/e2旋转方向、Shreve级数、流向等
% IX - 节点索引信息
% mean_junctions - 多拟合长度时的统计结果（均值/标准差）
%
% 注意：
% 支流1始终为较大支流（按Shreve级数或集水面积判断）
%
% 示例：
% [J,IX]=JunctionAngle(S,A,DEM,1000);
% [J,IX]=JunctionAngle(S,A,DEM,1000,'verbose',true,'make_shape',true,'ref_concavity',0.45);
% [J,IX,MJ]=JunctionAngle(S,A,[250 500 1000 2500 5000]);
% [J,IX]=JunctionAngle(S,A,DEM,1000,'predict_angle_method','slope');
% % 最小功率准则示例
% PRECIPr=resample(PRECIP_M_S,DEM,'nearest');
% Q=flowacc(FD,PRECIPr);
% [J,IX]=JunctionAngle(S,A,DEM,1000,'discharge',Q,'predict_angle_method','minimum_power');
%
% 相关函数：
% InspectJunction, JunctionLinks, JunctionSegments
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数作者：Adam M. Forte - 最后更新日期：2021/09/27 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'JunctionAngle';
	addRequired(p,'S',@(x) isa(x,'STREAMobj'));
	addRequired(p,'A',@(x) isa(x,'GRIDobj'));
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'fit_distance',@(x) isnumeric(x));

	addParameter(p,'use_n_nodes',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'previous_IX',[],@(x) iscell(x) || isempty(x));
	addParameter(p,'discharge',[],@(x) isa(x,'GRIDobj'));
	addParameter(p,'predict_angle_method','area',@(x) ischar(validatestring(x,{'slope','area','minimum_power','all'})));
	addParameter(p,'ref_concavity',0.50,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'make_shape',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'verbose',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'out_dir',[],@(x) isdir(x));
	addParameter(p,'out_name',[],@(x) ischar(x));

	parse(p,S,A,DEM,fit_distance,varargin{:});
	S=p.Results.S;
	A=p.Results.A;
	DEM=p.Results.DEM;
	fit_distance=p.Results.fit_distance;

	use_n_nodes=p.Results.use_n_nodes;
	PIX=p.Results.previous_IX;
	Q=p.Results.discharge;
	method=p.Results.predict_angle_method;
	ref_theta=p.Results.ref_concavity;
	make_shape=p.Results.make_shape;
	verbose=p.Results.verbose;
	out_dir=p.Results.out_dir;
	out_name=p.Results.out_name;

	if isempty(out_dir)
		out_dir=pwd;
	end
	 
	if use_n_nodes
		if numel(fit_distance)==1
			multi=false;
			n=fit_distance;
			if n<1
				n=1;
			end
			fit_distance=n*(hypot(S.cellsize,S.cellsize));
		else
			multi=true;
			n=max(fit_distance);
			if n<1
				n=1;
			end

			nlist=fit_distance;
			nlist(nlist<1)=1;
			nlist=unique(nlist);
			if numel(nlist)==1
				n=nlist;
				multi=false;
			end			

			fit_distance=nlist.*(hypot(S.cellsize,S.cellsize));
			fit_distance=sort(fit_distance);
		end
	else
		% Convert distance to number of nodes and check if multiple distances were provided.
		if numel(fit_distance)==1
			multi=false;
			n=floor(fit_distance/hypot(S.cellsize,S.cellsize)); % 计算节点数
			if n<1
				if isdeployed
					warndlg('输入的拟合距离过短，将使用最小节点数1','参数警告') % GUI警告
				end
				warning('输入的拟合距离过短，将使用最小节点数1'); % 命令行警告
				n=1;
			end
		else
			n=floor(max(fit_distance)/hypot(S.cellsize,S.cellsize));
			multi=true;
			if n<1
				if isdeployed
					warndlg('所有输入拟合距离过短，将使用单节点计算','参数警告')
				end
				warning('所有输入拟合距离过短，将使用单节点计算');
				n=1;
				mutli=false;
			end

			nlist=floor(fit_distance/hypot(S.cellsize,S.cellsize));
			if any(nlist)<1
				if isdeployed
					warndlg('部分输入距离过短，已调整为最小节点数','参数调整')
				end
				warning('部分输入距离过短，已调整为最小节点数');
				nlist(nlist<1)=1;
				nlist=unique(nlist);
				if numel(nlist)==1
					n=nlist;
					multi=false;
				end
            end 

			fit_distance=sort(fit_distance); 
		end
	end

	% Extract stream distances
	dst=S.distance;

	% Calculate Shreve stream orders
	so=streamorder(S,'shreve');

	% 检查必要参数
	if strcmp(method,'minimum_power') & isempty(Q)
		error('使用最小功耗法必须提供流量GRIDobj参数');
	end

	if ~isempty(Q) & ~validatealignment(S,Q)
		error('流量GRIDobj与STREAMobj空间参考不一致');
	end

	% 计算坡度（坡度法/最小功耗法需要）
	switch method
	case {'slope','minimum_power','all'}
		if verbose
			disp('正在计算河流梯度...')
		end
		z=mincosthydrocon(S,DEM,'interp',0.1); % 水文校正高程
		sg=gradient(S,z); % 计算坡度
		G=GRIDobj(DEM);
		G.Z(S.IXgrid)=sg; % 生成坡度栅格
	end

	% Find indices of nodes up and downstream
	if isempty(PIX)
		IX=nnodesupdown(S,A,n,so,verbose);
	else
		% Check that max distance is consistent 
		UPIX=PIX(:,2);
		nn=cellfun(@numel,UPIX);
		maxn=max(nn);
		if maxn<n
			if isdeployed
				errordlg('提供的previous_IX参数与最大拟合距离不兼容','输入错误')
			end
			error('提供的previous_IX参数与最大拟合距离不兼容');
		else
			IX=PIX;
		end
	end

	% 坐标转换处理
	if verbose
		disp('正在转换为地理坐标系...')
	end

	if isempty(S.georef)
		xl=S.x; yl=S.y; % 无投影时使用原始坐标
		if isdeployed
			warndlg('未找到投影信息，将基于投影坐标系计算角度','投影警告')
		end
		warning('未找到投影信息，将基于投影坐标系计算角度');
		projcoord=false;
	else
		try
			[yl,xl]=projinv(S.georef.mstruct,S.x,S.y);
			projcoord=true;
		catch
			xl=S.x; yl=S.y;
			if isdeployed
				warndlg('Projection was not recognized, angles will be calculated based on projected coordinates')
			end
			warning('Projection was not recognized, angles will be calculated based on projected coordinates')
			projcoord=false;
		end
	end

	% Break up output of node indicees
	num_con=size(IX,1);
	con=IX(:,1);
	ds=IX(:,2);
	us=IX(:,3:end);
	% Check and mark if there any streams with more than 2 tribs and if any
	% downstream segments of confluences don't exist to populate logical index
	if size(us,2)>2
		usvld=cellfun(@isempty,us(:,3));
	else
		usvld=true(size(con));
	end
	dsvld=~cellfun(@isempty,ds);
	vld=dsvld & usvld;

	if verbose
		disp('Finding junction angles...')
	end

	if ~multi
		% Preallocate arrays
		link_angle=zeros(num_con,11);
		strm_dir=zeros(num_con,3);
		fit_streams=zeros(num_con,9);

		switch method
		case {'area','slope'}
			pred_angle=zeros(num_con,2);
		case 'both'
			pred_angle=zeros(num_con,4);
		end

		for ii=1:num_con
			if vld(ii)
				% Extract x y node lists
				% Confluence point
				xc=xl(con{ii}); yc=yl(con{ii});
				% Upstream links
				xu1=xl(us{ii,1}); yu1=yl(us{ii,1});
				xu2=xl(us{ii,2}); yu2=yl(us{ii,2});
				% Downstream link
				xd=xl(ds{ii,1}); yd=yl(ds{ii,1});

				% Translate so all links originate from confluence 
				xu1=xu1-xc; xu2=xu2-xc;
				yu1=yu1-yc; yu2=yu2-yc;
				xd=xd-xc; yd=yd-yc;

				%% Depending on orientation of links, doing a linear fit
				%  on the link may produce spurious results. To improve 
				%  the result, each link is rotated towards horizontal, 
				%  using the mean orientation of the link to find the angle
				%  of rotation

				% Find mean angle of stream points from stream orientation
				theta1ap=median(atan2(yu1,xu1));
				theta2ap=median(atan2(yu2,xu2));
				theta3ap=median(atan2(yd,xd));

				% Extract and calculate stream distances
				e1_dst=max(dst(us{ii,1}))-dst(con{ii});
				e2_dst=max(dst(us{ii,2}))-dst(con{ii});
				ds_dst=dst(con{ii})-min(dst(ds{ii,1}));

				% Calculate predicted angles
				switch method
				case 'area'
					% Howard 1971 area method
					ix1=S.IXgrid(us{ii,1});
					ix2=S.IXgrid(us{ii,2});
					if isempty(Q)
						da1=max(A.Z(ix1).*A.cellsize^2);
						da2=max(A.Z(ix2).*A.cellsize^2);
					else
						da1=max(Q.Z(ix1));
						da2=max(Q.Z(ix2));
					end
					e1p=acos(((da1+da2)/da1)^-ref_theta);
					e2p=acos(((da1+da2)/da2)^-ref_theta);
					% Convert to degrees
					e1p=rad2deg(e1p);
					e2p=rad2deg(e2p);
					% Store out
					pred_angle(ii,:)=[e1p e2p];
				case 'slope'
					% Howard 1971 slope method
					ix1=S.IXgrid(us{ii,1});
					ix2=S.IXgrid(us{ii,2});
					ix3=S.IXgrid(ds{ii,1});
					sl1=mean(G.Z(ix1));
					sl2=mean(G.Z(ix2));
					sl3=mean(G.Z(ix3));
					r1=sl3/sl1;
					r2=sl3/sl2;
					if r1>1
						e1p=0;
					else
						e1p=acos(r1);
					end

					if r2>1
						e2p=0;
					else
						e2p=acos(r2);
					end
					% Convert to degrees
					e1p=rad2deg(e1p);
					e2p=rad2deg(e2p);
					% Store out
					pred_angle(ii,:)=[e1p e2p];
				case 'minimum_power'
					% Howard 1971 minimum power method
					ix1=S.IXgrid(us{ii,1});
					ix2=S.IXgrid(us{ii,2});
					ix3=S.IXgrid(ds{ii,1});
					sl1=mean(G.Z(ix1));
					sl2=mean(G.Z(ix2));
					sl3=mean(G.Z(ix3));
					q1=max(Q.Z(ix1));
					q2=max(Q.Z(ix2));
					q3=q1+q2;
					C1=q1*sl1;
					C2=q2*sl2;
					C3=q3*sl3;
					inner1=((C3^2)+(C1^2)-(C2^2))/(2*C1*C3);
					inner2=((C3^2)+(C2^2)-(C1^2))/(2*C2*C3);
					if inner1>1
						e1p=0;
					else
						e1p=acos(inner1);
					end

					if inner>1
						e2p=0;
					else
						e2p=acos(inner2);
					end
					% Convert to degrees
					e1p=rad2deg(e1p);
					e2p=rad2deg(e2p);				
					% Store out
					pred_angle(ii,:)=[e1p e2p];
				case 'all'
					% Howard 1971 area method
					ix1=S.IXgrid(us{ii,1});
					ix2=S.IXgrid(us{ii,2});
					if isempty(Q)
						da1=max(A.Z(ix1).*A.cellsize^2);
						da2=max(A.Z(ix2).*A.cellsize^2);
					else
						da1=max(Q.Z(ix1));
						da2=max(Q.Z(ix2));
					end
					e1Ap=acos(((da1+da2)/da1)^-ref_theta);
					e2Ap=acos(((da1+da2)/da2)^-ref_theta);
					% Convert to degrees
					e1Ap=rad2deg(e1Ap);
					e2Ap=rad2deg(e2Ap);

					% Howard 1971 slope method
					ix3=S.IXgrid(ds{ii,1});
					sl1=mean(G.Z(ix1));
					sl2=mean(G.Z(ix2));
					sl3=mean(G.Z(ix3));
					r1=sl3/sl1;
					r2=sl3/sl2;
					if r1>1
						e1Sp=0;
					else
						e1Sp=acos(r1);
					end

					if r2>1
						e2Sp=0;
					else
						e2Sp=acos(r2);
					end
					% Convert to degrees
					e1Sp=rad2deg(e1Sp);
					e2Sp=rad2deg(e2Sp);

					% Howard 1971 Minimum Power
					if ~isempty(Q)
						C1=da1*sl1;
						C2=da2*sl2;
						C3=(da1+da2)*sl3;
						inner1=((C3^2)+(C1^2)-(C2^2))/(2*C1*C3);
						inner2=((C3^2)+(C2^2)-(C1^2))/(2*C2*C3);
						if inner1>1
							e1Mp=0;
						else
							e1Mp=acos(inner1);
						end

						if inner>1
							e2Mp=0;
						else
							e2Mp=acos(inner2);
						end
						% Convert to degrees
						e1Mp=rad2deg(e1Mp);
						e2Mp=rad2deg(e2Mp);	
						% Store out
						pred_angle(ii,:)=[e1Ap e2Ap e1Sp e2Sp e1Mp e2Mp];					
					else
						% Store out
						pred_angle(ii,:)=[e1Ap e2Ap e1Sp e2Sp];
					end
				end

				% Rotate all points to horizontal
				[xu1r,yu1r]=rotcoord(xu1,yu1,-theta1ap,0,0); 
				[xu2r,yu2r]=rotcoord(xu2,yu2,-theta2ap,0,0);	
				[xdr,ydr]=rotcoord(xd,yd,-theta3ap,0,0);

				% Simple approximation of linear segment on
				% rotated positions
				a1=xu1r\yu1r;
				a2=xu2r\yu2r; 
				a3=xdr\ydr;

				% Find projected y positions
				yp1r=xu1r*a1;
				yp2r=xu2r*a2;
				yp3r=xdr*a3;

				% Rotate back
				[xp1,yp1]=rotcoord(xu1r,yp1r,theta1ap,0,0);
				[xp2,yp2]=rotcoord(xu2r,yp2r,theta2ap,0,0);
				[xp3,yp3]=rotcoord(xdr,yp3r,theta3ap,0,0);

				% Calculate r-squared
				r21=rsquared(yu1,yp1);
				r22=rsquared(yu2,yp2);
				r23=rsquared(yd,yp3);

				% Convert to polar
				[theta1,rho1]=cart2pol(xp1,yp1);
				[theta2,rho2]=cart2pol(xp2,yp2);
				[theta3,rho3]=cart2pol(xp3,yp3);

				% Find maximum radii and thetas for those radii
				[mrho1,ix1]=max(rho1);
				[mrho2,ix2]=max(rho2);
				[mrho3,ix3]=max(rho3);
				theta1=theta1(ix1);
				theta2=theta2(ix2);
				theta3=theta3(ix3);

				% Calculate interlink angles
				[e1,e2,ba,split,d1,d2]=interangle(theta1,theta2,theta3);
				e1=rad2deg(e1);
				e2=rad2deg(e2);
				ba=rad2deg(ba);

				% Determine order of streams and define handedness
				% sensu James & Krumbein, 1969
				so1=max(so(us{ii,1}));
				so2=max(so(us{ii,2}));
				% so1 should be greater than or equal to so2 as output
				% from nnodesupdown function
				if so1>so2 
					if d1==1 & d2==2
						hnd=1;
					elseif d1==2 & d2==1
						hnd=2;
					elseif d1==1 & d2==1 & e1>e2 
						hnd=1;
					elseif d1==1 & d2==1 & e1<e2 
						hnd=2;
					elseif d1==2 & d2==2 & e1>e2
						hnd=2;
					elseif d1==2 & d2==2 & e1<e2
						hnd=1;
					else
						hnd=3;
					end
				else
					hnd=3;
				end

				% Convert direction angles to cardinal
				card1=mod(-90-rad2deg(theta1),360);
				card2=mod(-90-rad2deg(theta2),360);
				card3=mod(-90-rad2deg(theta3),360);

				link_angle(ii,:)=[S.x(con{ii}) S.y(con{ii}) ba e1 e2 split d1 d2 hnd so1 so2];
				strm_dir(ii,:)=[card1 card2 card3];
				fit_streams(ii,:)=[e1_dst numel(xu1) r21 e2_dst numel(xu2) r22 ds_dst numel(xd) r23];
			else
				link_angle(ii,:)=[S.x(con{ii}) S.y(con{ii}) NaN NaN NaN NaN NaN NaN NaN NaN NaN];
				strm_dir(ii,:)=[NaN NaN NaN];
				fit_streams(ii,:)=[NaN NaN NaN NaN NaN NaN NaN NaN NaN];

				switch method
				case {'area','slope','minimum_power'}
					pred_angle(ii,:)=[NaN NaN];
				case 'all'
					if isempty(Q)
						pred_angle(ii,:)=[NaN NaN NaN NaN];
					else
						pred_angle(ii,:)=[NaN NaN NaN NaN NaN NaN];
					end
				end
			end
		end

		if make_shape
			makejunctionshape(link_angle,pred_angle,strm_dir,fit_streams,method,1,out_dir,out_name);
		end

		junctions=makejunctiontable(link_angle,pred_angle,strm_dir,fit_streams,method);

	else
		num_n=numel(nlist);

		% Preallocate cell array
		junctions=cell(2,num_n);

		for jj=1:num_n
			nOI=nlist(jj);

			% Preallocate arrays
			la=zeros(num_con,11);
			sd=zeros(num_con,3);
			fs=zeros(num_con,9);

			switch method
			case {'area','slope','minimum_power'}
				pa=zeros(num_con,2);
			case 'all'
				if isempty(Q)
					pa=zeros(num_con,4);
				else
					pa=zeros(num_con,6);
				end
			end			

			for ii=1:num_con
				if vld(ii)
					% Extract ix lists and clip
					conix=con{ii};
					us1ix=us{ii,1};
					us2ix=us{ii,2};
					dsix=ds{ii,1};

					if numel(us1ix)>nOI
						us1ix=us1ix(1:nOI);
					end

					if numel(us2ix)>nOI
						us2ix=us2ix(1:nOI);
					end

					if numel(dsix)>nOI
						dsix=dsix(1:nOI);
					end

					% Extract x y node lists
					% Confluence point
					xc=xl(conix); yc=yl(conix);
					% Upstream links
					xu1=xl(us1ix); yu1=yl(us1ix);
					xu2=xl(us2ix); yu2=yl(us2ix);
					% Downstream link
					xd=xl(dsix); yd=yl(dsix);

					% Translate so all links originate from confluence 
					xu1=xu1-xc; xu2=xu2-xc;
					yu1=yu1-yc; yu2=yu2-yc;
					xd=xd-xc; yd=yd-yc;

					%% Depending on orientation of links, doing a linear fit
					%  on the link may produce spurious results. To improve 
					%  the result, each link is rotated towards horizontal, 
					%  using the mean orientation of the link to find the angle
					%  of rotation

					% Find mean angle of stream points from stream orientation
					theta1ap=median(atan2(yu1,xu1));
					theta2ap=median(atan2(yu2,xu2));
					theta3ap=median(atan2(yd,xd));

					% Extract and calculate stream distances
					e1_dst=max(dst(us1ix))-dst(conix);
					e2_dst=max(dst(us2ix))-dst(conix);
					ds_dst=dst(conix)-min(dst(dsix));

					% Calculate predicted angles 
					switch method
					case 'area'
						% Howard 1971 area method
						ix1=S.IXgrid(us1ix);
						ix2=S.IXgrid(us2ix);
						if isempty(Q)
							da1=max(A.Z(ix1).*A.cellsize^2);
							da2=max(A.Z(ix2).*A.cellsize^2);
						else
							da1=max(Q.Z(ix1));
							da2=max(Q.Z(ix2));
						end
						e1p=acos(((da1+da2)/da1)^-ref_theta);
						e2p=acos(((da1+da2)/da2)^-ref_theta);
						% Convert to degrees
						e1p=rad2deg(e1p);
						e2p=rad2deg(e2p);
						% Store out
						pa(ii,:)=[e1p e2p];
					case 'slope'
						% Howard 1971 slope method
						ix1=S.IXgrid(us1ix);
						ix2=S.IXgrid(us2ix);
						ix3=S.IXgrid(dsix);
						sl1=mean(G.Z(ix1));
						sl2=mean(G.Z(ix2));
						sl3=mean(G.Z(ix3));
						r1=sl3/sl1;
						r2=sl3/sl2;
						if r1>1
							e1p=0;
						else
							e1p=acos(r1);
						end

						if r2>1
							e2p=0;
						else
							e2p=acos(r2);
						end
						% Convert to degrees
						e1p=rad2deg(e1p);
						e2p=rad2deg(e2p);
						% Store out
						pa(ii,:)=[e1p e2p];
					case 'minimum_power'
						% Howard 1971 minimum power method
						ix1=S.IXgrid(us1ix);
						ix2=S.IXgrid(us2ix);
						ix3=S.IXgrid(dsix);
						sl1=mean(G.Z(ix1));
						sl2=mean(G.Z(ix2));
						sl3=mean(G.Z(ix3));
						q1=max(Q.Z(ix1));
						q2=max(Q.Z(ix2));
						q3=q1+q2;
						C1=q1*sl1;
						C2=q2*sl2;
						C3=q3*sl3;
						inner1=((C3^2)+(C1^2)-(C2^2))/(2*C1*C3);
						inner2=((C3^2)+(C2^2)-(C1^2))/(2*C2*C3);
						if inner1>1
							e1p=0;
						else
							e1p=acos(inner1);
						end

						if inner>1
							e2p=0;
						else
							e2p=acos(inner2);
						end
						% Convert to degrees
						e1p=rad2deg(e1p);
						e2p=rad2deg(e2p);				
						% Store out
						pa(ii,:)=[e1p e2p];
					case 'all'
						% Howard 1971 area method
						ix1=S.IXgrid(us1ix);
						ix2=S.IXgrid(us2ix);
						if isempty(Q)
							da1=max(A.Z(ix1).*A.cellsize^2);
							da2=max(A.Z(ix2).*A.cellsize^2);
						else
							da1=max(Q.Z(ix1));
							da2=max(Q.Z(ix2));
						end
						e1Ap=acos(((da1+da2)/da1)^-ref_theta);
						e2Ap=acos(((da1+da2)/da2)^-ref_theta);
						% Convert to degrees
						e1Ap=rad2deg(e1Ap);
						e2Ap=rad2deg(e2Ap);	

						% Howard 1971 slope method
						ix3=S.IXgrid(dsix);
						sl1=mean(G.Z(ix1));
						sl2=mean(G.Z(ix2));
						sl3=mean(G.Z(ix3));
						r1=sl3/sl1;
						r2=sl3/sl2;
						if r1>1
							e1Sp=0;
						else
							e1Sp=acos(r1);
						end

						if r2>1
							e2Sp=0;
						else
							e2Sp=acos(r2);
						end
						% Convert to degrees
						e1Sp=rad2deg(e1Sp);
						e2Sp=rad2deg(e2Sp);

						% Howard 1971 Minimum Power
						if ~isempty(Q)
							C1=da1*sl1;
							C2=da2*sl2;
							C3=(da1+da2)*sl3;
							inner1=((C3^2)+(C1^2)-(C2^2))/(2*C1*C3);
							inner2=((C3^2)+(C2^2)-(C1^2))/(2*C2*C3);
							if inner1>1
								e1Mp=0;
							else
								e1Mp=acos(inner1);
							end

							if inner>1
								e2Mp=0;
							else
								e2Mp=acos(inner2);
							end
							% Convert to degrees
							e1Mp=rad2deg(e1Mp);
							e2Mp=rad2deg(e2Mp);	
							% Store out
							pa(ii,:)=[e1Ap e2Ap e1Sp e2Sp e1Mp e2Mp];					
						else
							% Store out
							pa(ii,:)=[e1Ap e2Ap e1Sp e2Sp];
						end
					end

					% Rotate all points to horizontal
					[xu1r,yu1r]=rotcoord(xu1,yu1,-theta1ap,0,0); 
					[xu2r,yu2r]=rotcoord(xu2,yu2,-theta2ap,0,0);	
					[xdr,ydr]=rotcoord(xd,yd,-theta3ap,0,0);

					% Simple approximation of linear segment on
					% rotated positions
					a1=xu1r\yu1r;
					a2=xu2r\yu2r; 
					a3=xdr\ydr;

					% Find projected y positions
					yp1r=xu1r*a1;
					yp2r=xu2r*a2;
					yp3r=xdr*a3;

					% Rotate back
					[xp1,yp1]=rotcoord(xu1r,yp1r,theta1ap,0,0);
					[xp2,yp2]=rotcoord(xu2r,yp2r,theta2ap,0,0);
					[xp3,yp3]=rotcoord(xdr,yp3r,theta3ap,0,0);

					% Calculate r-squared
					r21=rsquared(yu1,yp1);
					r22=rsquared(yu2,yp2);
					r23=rsquared(yd,yp3);

					% Convert to polar
					[theta1,rho1]=cart2pol(xp1,yp1);
					[theta2,rho2]=cart2pol(xp2,yp2);
					[theta3,rho3]=cart2pol(xp3,yp3);

					% Find maximum radii and thetas for those radii
					[mrho1,ix1]=max(rho1);
					[mrho2,ix2]=max(rho2);
					[mrho3,ix3]=max(rho3);
					theta1=theta1(ix1);
					theta2=theta2(ix2);
					theta3=theta3(ix3);

					% Calculate interlink angles
					[e1,e2,ba,split,d1,d2]=interangle(theta1,theta2,theta3);
					e1=rad2deg(e1);
					e2=rad2deg(e2);
					ba=rad2deg(ba);

					% Determine order of streams and define handedness
					% sensu James & Krumbein, 1969
					so1=max(so(us1ix));
					so2=max(so(us2ix));
					% so1 should be greater than or equal to so2 as output
					% from nnodesupdown function
					if so1>so2 
						if d1==1 & d2==2
							hnd=1;
						elseif d1==2 & d2==1
							hnd=2;
						elseif d1==1 & d2==1 & e1>e2 
							hnd=1;
						elseif d1==1 & d2==1 & e1<e2 
							hnd=2;
						elseif d1==2 & d2==2 & e1>e2
							hnd=2;
						elseif d1==2 & d2==2 & e1<e2
							hnd=1;
						else
							hnd=3;
						end
					else
						hnd=3;
					end

					% Convert direction angles to cardinal
					card1=mod(-90-rad2deg(theta1),360);
					card2=mod(-90-rad2deg(theta2),360);
					card3=mod(-90-rad2deg(theta3),360);

					la(ii,:)=[S.x(conix) S.y(conix) ba e1 e2 split d1 d2 hnd so1 so2];
					sd(ii,:)=[card1 card2 card3];
					fs(ii,:)=[e1_dst numel(xu1) r21 e2_dst numel(xu2) r22 ds_dst numel(xd) r23];
				else
					la(ii,:)=[S.x(con{ii}) S.y(con{ii}) NaN NaN NaN NaN NaN NaN NaN NaN NaN];
					sd(ii,:)=[NaN NaN NaN];
					fs(ii,:)=[NaN NaN NaN NaN NaN NaN NaN NaN NaN];

					switch method
					case {'slope','area','minimum_power'}
						pa(ii,:)=[NaN NaN];
					case 'both'
						if isempty(Q)
							pa(ii,:)=[NaN NaN NaN NaN];
						else
							pa(ii,:)=[NaN NaN NaN NaN NaN NaN];
						end
					end
				end
			end

			if make_shape
				makejunctionshape(la,pa,sd,fs,method,jj,out_dir,out_name);
			end

			junctions{1,jj}=fit_distance(jj);
			junctions{2,jj}=makejunctiontable(la,pa,sd,fs,method);

		end	% end nodes for

		% Calculate sensitivity to fit distance
		jamat=zeros(num_con,num_n);
		e1mat=zeros(size(jamat));
		e2mat=zeros(size(jamat));
		for jj=1:num_n
			J=junctions{2,jj};
			jamat(:,jj)=J.junction_angle;
			e1mat(:,jj)=J.e1_obs_angle;
			e2mat(:,jj)=J.e2_obs_angle;
		end

		% Calculate mean and standard deviation
		mean_junction_angle=mean(jamat,2);
		mean_e1_angle=mean(e1mat,2);
		mean_e2_angle=mean(e2mat,2);
		std_junction_angle=std(jamat,1,2);
		std_e1_angle=std(e1mat,1,2);
		std_e2_angle=std(e2mat,1,2);

		if strcmp(method,'slope') | strcmp(method,'all')
			e1pmat=zeros(size(jamat));
			e2pmat=zeros(size(jamat));
			for jj=1:num_n
				J=junctions{2,jj};
				e1pmat(:,jj)=J.e1_Spred_angle;
				e2pmat(:,jj)=J.e2_Spred_angle;
			end

			mean_e1Sp_angle=mean(e1pmat,2);
			std_e1Sp_angle=std(e1pmat,1,2);
			mean_e2Sp_angle=mean(e2pmat,2);
			std_e2Sp_angle=std(e2pmat,1,2);	

			MJ=table(mean_junction_angle,std_junction_angle,...
				mean_e1_angle,std_e1_angle,mean_e1Sp_angle,std_e1Sp_angle,...
				mean_e2_angle,std_e2_angle,mean_e2Sp_angle,std_e2Sp_angle);						
		else
			MJ=table(mean_junction_angle,std_junction_angle,...
				mean_e1_angle,std_e1_angle,...
				mean_e2_angle,std_e2_angle);	
		end

		varargout{1}=MJ;
	end % end multi if
end

function [r2]=rsquared(d,pred);
	sstot=sum((d-mean(d)).^2);
	ssres=sum((d-pred).^2);
	r2=1-(ssres/sstot);
end

function [IX]=nnodesupdown(S,A,n,so,verbose)
	% Find indices of a specified number of nodes up and downstream of a confluence
	% Output IX will be a cell array containing indices. All indices refer to positions
	% in the node attributed list of the provided STREAMobj. First column is index of confluence
	% Second column are indices of n nodes downstream of confluence (not including the confluence).
	% Additional columns are n nodes upstream of confluence. Nominally, there should be two columns 
	% beyond the second column, but if there any confluences where more than 2 streams meet, then
	% there will be additional columns (which will be empty for all confluences that only have
	% two links upstream).
	%
	% In the downstream direction, if the search reaches the end of a stream (outlet) or 
	% a b-confluence point, then there will be less than n nodes downstream of this particular confluence.
	% If a confluence coincides with an outlet, then the downstream nodes cell will be an empty array.
	%
	% In the upstream direction, if the search reaches the end of a stream (channel head) or another
	% confluence, there will be less than n nodes upstream of this particular confluence (in the
	% relevant link).
	%
	% Upstream links are sorted so that tributary 1 has higher shreve order or higher drainage area
	% if shreve order of tributary 1 and 2 is the same.

	if verbose
		disp('Finding confluence points...')
	end

	% Calculate drainage area and shreve order
	da=A.Z(S.IXgrid).*A.cellsize^2;

	% Extract confluences
	cix=streampoi(S,'confluences','ix');
	clo=streampoi(S,'bconfluences','logical');

	% Find index of confluences within nal
	cixnal=find(ismember(S.IXgrid,cix));

	%% Find downstream nodes
	IXDS=cell(numel(cixnal),2);
	% Loop through confluences and find n number of nodes downstream
	if verbose
		w1=waitbar(0,'Finding nodes downstream of confluences...');
	end	

	for ii=1:numel(cixnal)
		% Position of confluence in giver list
		ix=find(S.ix==cixnal(ii));

		ixl=zeros(n,1);
		jj=1;
		while jj<=n & ~isempty(ix)
			ixl(jj)=S.ixc(ix);

			% Stop if point is b-confluence
			if clo(ixl(jj))
				break 
			end

			ix=find(S.ix==ixl(jj));
			jj=jj+1;
		end
		ixl(ixl==0)=[];

		IXDS{ii,1}=cixnal(ii);
		IXDS{ii,2}=ixl;

		if verbose
			waitbar(ii/numel(cixnal));
		end
	end

	if verbose
		close(w1);
	end

	%% Find upstream nodes

	% Preallocate cell array of most likely size (may expand columns if there are more
	%	than one tributary at a confluence)
	IXUS=cell(numel(cixnal),2);
	% Loop through confluences and find n number of nodes downstream
	if verbose
		w1=waitbar(0,'Finding nodes upstream of confluences...');
	end	

	for ii=1:numel(cixnal)
		% Position of confluence in receiver list
		lix=find(S.ixc==cixnal(ii));
		% Determine number of links
		num_links=numel(lix);
		sort_mat=zeros(num_links,2);
		for jj=1:num_links
			ixl=zeros(n,1);
			ix=lix(jj);
			kk=1;
			while kk<=n & ~isempty(ix) & numel(ix)==1;
				ixl(kk)=S.ix(ix);
				ix=find(S.ixc==ixl(kk));
				kk=kk+1;
			end
			ixl(ixl==0)=[];
			IXUS{ii,jj}=ixl;
			sort_mat(jj,1)=max(so(ixl));
			sort_mat(jj,2)=max(da(ixl));
		end
		% Sort links
		[~,six]=sortrows(sort_mat,'descend');
		% IXUS(ii,:)=IXUS(ii,six);	
		IXUS(ii,1:numel(six))=IXUS(ii,six);	

		if verbose
			waitbar(ii/numel(cixnal));
		end		
	end	

	if verbose
		close(w1);
	end

	% Join into master cell array
	IX=horzcat(IXDS,IXUS);
end

function [e1,e2,ba,split,d1,d2]=interangle(t1,t2,t3)
	% Assumes all angles are as returned from
	% cart2pol between -pi and pi
	%
	% t1 is tributary 1 theta
	% t2 is tributary 2 theta
	% t3 is downstream theta
	% 
	% e1 is angle between upstream projection and t1
	% e2 is angle between upstream projection and t2
	% ba is angle between t1 and t2
	% split is a logical indicating whether the upstream
	%	projection is between t1 and t2 (true) or not (false)
	% d1 is direction from t1 to upstream projection, CW is 1
	%	CCW is 2
	% d2 is direction from t2 to upstream direction, CW is 1
	%	CCW is 2

	% Rotate entire system so that downstream link always
	% points toward 0 and upstream projection is pi
	t1r=t1-t3;
	t2r=t2-t3;
	t3r=t3-t3;

	% More efficient way to remap rotate angles into range of
	% -pi to pi
	t1r=atan2(sin(t1r),cos(t1r));
	t2r=atan2(sin(t2r),cos(t2r));

	% Calculate angles with respect upstream projection of downstream
	% reach and interlink angle. Also determine whether the upstream
	% projection is between the two tributaries
	if t1r<0 && t2r>=0
		e1=pi-abs(t1r);
		e2=pi-t2r;
		ba=e1+e2;
		split=true;
		d1=2;
		d2=1;
	elseif t1r>=0 && t2r<0
		e1=pi-t1r;
		e2=pi-abs(t2r);
		ba=e1+e2;
		split=true;
		d1=1;
		d2=2;
	elseif t1r>=0 && t2r>=0
		e1=pi-t1r;
		e2=pi-t2r;
		ba=abs(t1r-t2r);
		split=false;
		d1=1;
		d2=1;
	elseif t1r<0 && t2r<0
		e1=pi-abs(t1r);
		e2=pi-abs(t2r);
		ba=abs(abs(t1r)-abs(t2r));
		split=false;
		d1=2;
		d2=2;
	end
end

function [n_x,n_y]=rotcoord(x,y,theta,x0,y0)
	n_x=(x-x0).*cos(theta)-(y-y0).*sin(theta);
	n_y=(x-x0).*sin(theta)+(y-y0).*cos(theta);
end

function [T]=makejunctiontable(la,pa,cr,fs,method)

	num_juncs=numel(la(:,1));

	junction_number=1:num_juncs;
	junction_number=junction_number(:);

	junction_x=la(:,1);
	junction_y=la(:,2);

	junction_angle=la(:,3);
	e1_obs_angle=la(:,4);
	e2_obs_angle=la(:,5);

	e1_shreve=la(:,10);
	e2_shreve=la(:,11);

	e1_direction=cr(:,1);
	e2_direction=cr(:,2);
	e3_direction=cr(:,3);

	e1_distance=fs(:,1);
	e2_distance=fs(:,4);
	e3_distance=fs(:,7);

	e1_num=fs(:,2);
	e2_num=fs(:,5);
	e3_num=fs(:,8);

	split=logical(zeros(num_juncs,1));
	e1_rotation=categorical(zeros(num_juncs,1));
	e2_rotation=categorical(zeros(num_juncs,1));
	handedness=categorical(zeros(num_juncs,1));
	e1_R2=zeros(num_juncs,1);
	e2_R2=zeros(num_juncs,1);
	e3_R2=zeros(num_juncs,1);
	for ii=1:num_juncs

		if la(ii,6)==1
			split(ii,1)=true;
		else
			split(ii,1)=false;
		end

		if la(ii,7)==1
			e1_rotation(ii)='CW';
		elseif la(ii,7)==2
			e1_rotation(ii)='CCW';
		else
			e1_rotation(ii)='Undefined';
		end

		if la(ii,8)==1
			e2_rotation(ii)='CW';
		elseif la(ii,8)==2
			e2_rotation(ii)='CCW';
		else
			e2_rotation(ii)='Undefined';
		end

		if la(ii,9)==1
			handedness(ii)='Right';
		elseif la(ii,9)==2
			handedness(ii)='Left';
		else
			handedness(ii)='Undefined';
		end

		E1_r2=fs(ii,3);
		E2_r2=fs(ii,6);
		E3_r2=fs(ii,9);

		if ~isnan(E1_r2) & ~isinf(E1_r2)
			e1_R2(ii,1)=E1_r2;
		else
			e1_R2(ii,1)=0;
		end

		if ~isnan(E2_r2) & ~isinf(E2_r2)
			e2_R2(ii,1)=E2_r2;
		else
			e2_R2(ii,1)=0;
		end

		if ~isnan(E3_r2) & ~isinf(E3_r2)
			e3_R2(ii,1)=E3_r2;
		else
			e3_R2(ii,1)=0;
		end
	end

	switch method
	case 'area'
		e1_Apred_angle=pa(:,1);	
		e2_Apred_angle=pa(:,2);

		T=table(junction_number,junction_x,junction_y,junction_angle,handedness,split,...
			e1_obs_angle,e1_Apred_angle,e1_rotation,e1_shreve,e1_direction,e1_distance,e1_num,e1_R2,...
			e2_obs_angle,e2_Apred_angle,e2_rotation,e2_shreve,e2_direction,e2_distance,e2_num,e2_R2,...
			e3_direction,e3_distance,e3_num,e3_R2);

	case 'slope'
		e1_Spred_angle=pa(:,1);	
		e2_Spred_angle=pa(:,2);

		T=table(junction_number,junction_x,junction_y,junction_angle,handedness,split,...
			e1_obs_angle,e1_Spred_angle,e1_rotation,e1_shreve,e1_direction,e1_distance,e1_num,e1_R2,...
			e2_obs_angle,e2_Spred_angle,e2_rotation,e2_shreve,e2_direction,e2_distance,e2_num,e2_R2,...
			e3_direction,e3_distance,e3_num,e3_R2);		

	case 'minimum_power'
		e1_MPpred_angle=pa(:,1);	
		e2_MPpred_angle=pa(:,2);

		T=table(junction_number,junction_x,junction_y,junction_angle,handedness,split,...
			e1_obs_angle,e1_MPpred_angle,e1_rotation,e1_shreve,e1_direction,e1_distance,e1_num,e1_R2,...
			e2_obs_angle,e2_MPpred_angle,e2_rotation,e2_shreve,e2_direction,e2_distance,e2_num,e2_R2,...
			e3_direction,e3_distance,e3_num,e3_R2);	

	case 'all'
		if size(pa,2)==4
			e1_Apred_angle=pa(:,1);	
			e2_Apred_angle=pa(:,2);
			e1_Spred_angle=pa(:,3);	
			e2_Spred_angle=pa(:,4);


			T=table(junction_number,junction_x,junction_y,junction_angle,handedness,split,...
				e1_obs_angle,e1_Apred_angle,e1_Spred_angle,e1_rotation,e1_shreve,e1_direction,e1_distance,e1_num,e1_R2,...
				e2_obs_angle,e2_Apred_angle,e2_Spred_angle,e2_rotation,e2_shreve,e2_direction,e2_distance,e2_num,e2_R2,...
				e3_direction,e3_distance,e3_num,e3_R2);	
		else
			e1_Apred_angle=pa(:,1);	
			e2_Apred_angle=pa(:,2);
			e1_Spred_angle=pa(:,3);	
			e2_Spred_angle=pa(:,4);
			e1_MPpred_angle=pa(:,5);
			e2_MPpred_angle=pa(:,6);			

			T=table(junction_number,junction_x,junction_y,junction_angle,handedness,split,...
				e1_obs_angle,e1_Apred_angle,e1_Spred_angle,e1_MPpred_angle,e1_rotation,e1_shreve,e1_direction,e1_distance,e1_num,e1_R2,...
				e2_obs_angle,e2_Apred_angle,e2_Spred_angle,e2_MPpred_angle,e2_rotation,e2_shreve,e2_direction,e2_distance,e2_num,e2_R2,...
				e3_direction,e3_distance,e3_num,e3_R2);	
		end
	end

end

function makejunctionshape(la,pa,cr,fs,method,num_shp,out_dir,out_name)
	ms=struct;

	all_nums=1:numel(pa(:,1));

	idx=~isnan(pa(:,1));


	la=la(idx,:);
	pa=pa(idx,:);
	cr=cr(idx,:);
	fs=fs(idx,:);
	all_nums=all_nums(idx);

	split=logical(la(:,6));

	num=size(la,1);
	for ii=1:num
		ms(ii,1).Geometry='Point';
		ms(ii,1).X=double(la(ii,1));
		ms(ii,1).Y=double(la(ii,2));
		ms(ii,1).Confl_ID=double(all_nums(ii));
		ms(ii,1).Junc_Ang=double(la(ii,3));
		if split(ii)
			ms(ii,1).Split='true';
		else
			ms(ii,1).Split='false';
		end		

		ms(ii,1).E1_Ang=double(la(ii,4));

		switch method 
		case 'area'
			ms(ii,1).E1_APred=double(pa(ii,1));
		case 'slope'
			ms(ii,1).E1_SPred=double(pa(ii,1));
		case 'minimum_power'
			ms(ii,1).E1_MPpred=double(pa(ii,1));
		case 'all'
			if size(pa,2)==4
				ms(ii,1).E1_APred=double(pa(ii,1));
				ms(ii,1).E1_SPred=double(pa(ii,3));
			else
				ms(ii,1).E1_APred=double(pa(ii,1));
				ms(ii,1).E1_SPred=double(pa(ii,3));
				ms(ii,1).E1_MPpred=double(pa(ii,5));
			end
		end			

		ms(ii,1).E2_Ang=double(la(ii,5));

		switch method 
		case 'area'
			ms(ii,1).E2_APred=double(pa(ii,2));
		case 'slope'
			ms(ii,1).E2_SPred=double(pa(ii,2));
		case 'minimum_power'
			ms(ii,1).E2_MPpred=double(pa(ii,2));			
		case 'all'
			if size(pa,2)==4
				ms(ii,1).E2_APred=double(pa(ii,2));
				ms(ii,1).E2_SPred=double(pa(ii,4));
			else
				ms(ii,1).E2_APred=double(pa(ii,2));
				ms(ii,1).E2_SPred=double(pa(ii,4));
				ms(ii,1).E2_MPpred=double(pa(ii,6));
			end
		end	

		if la(ii,7)==1
			ms(ii,1).E1_Rot='CW';
		else
			ms(ii,1).E1_Rot='CCW';
		end

		if la(ii,8)==1
			ms(ii,1).E2_Rot='CW';
		else
			ms(ii,1).E2_Rot='CCW';
		end

		if la(ii,9)==1
			ms(ii,1).Handedness='R';
		elseif la(ii,9)==2
			ms(ii,1).Handedness='L';
		else
			ms(ii,1).Handedness='NA';
		end

		ms(ii,1).E1_Dir=double(cr(ii,1));
		ms(ii,1).E2_Dir=double(cr(ii,2));
		ms(ii,1).DS_Dir=double(cr(ii,3));

		ms(ii,1).E1_Dst=double(fs(ii,1));
		ms(ii,1).E2_Dst=double(fs(ii,4));
		ms(ii,1).DS_Dst=double(fs(ii,7));

		e1_r2=fs(ii,3);
		e2_r2=fs(ii,6);
		e3_r2=fs(ii,9);

		if ~isnan(e1_r2) & ~isinf(e1_r2)
			ms(ii,1).E1_R2=double(e1_r2);
		else
			ms(ii,1).E1_R2=double(0);
		end

		if ~isnan(e2_r2) & ~isinf(e2_r2)
			ms(ii,1).E2_R2=double(e2_r2);
		else
			ms(ii,1).E2_R2=double(0);
		end

		if ~isnan(e3_r2) & ~isinf(e3_r2)
			ms(ii,1).DS_R2=double(e3_r2);
		else
			ms(ii,1).DS_R2=double(0);
		end

	end

	if isempty(out_name)
		shp_name=fullfile(out_dir,['junction_angle_' num2str(num_shp) '.shp']);
	else
		shp_name=fullfile(out_dir,[out_name '_junction_angle_' num2str(num_shp) '.shp']);
	end
	shapewrite(ms,shp_name);
end