function [OUT]=SegmentProjector(DEM,FD,A,S,basin_num,varargin)
    %
    % 用法：
    %     [OUT]=SegmentProjector(DEM,FD,A,S,basin_num);
    %     [OUT]=SegmentProjector(DEM,FD,A,S,basin_num,'name',value,...);
    %
    % 描述：
    %     该函数用于交互式选择您想要投影的河道剖面段（例如投影具有不同ksn的剖面部分）。
    %     您可以使用 'SegmentPicker' 函数交互式选择要提供给 StreamProjector 函数的河道。
    %     为此，加载 PickedSegments_*.mat 文件，并在调用 SegmentProjector 时用 'Sc' STREAMobj 替换 'S'。
    %     如果 STREAMobj 包含多个河道源，该代码将遍历所有河道（即确保只提供您想要投影的河道，而不是整个网络！）。
    %     该函数计算并显示拟合的95%置信区间。
    %
    % 必需输入参数：
    %     DEM - 数字高程模型（GRIDobj类型），假定未进行平滑的DEM（例如来自 ProcessRiverBasins 的 DEMoc）
    %     FD - 流向（FLOWobj类型）
    %     A - 流积累量（GRIDobj类型）
    %     S - 要投影的河道（STREAMobj类型）
    %     basin_num - 盆地编号或用于标识要选择的河道集的其他编号
    %
    % 可选输入参数：
    %     conditioned_DEM [] - 可选的平滑DEM用于此函数（不要将平滑后的DEM提供给主要的DEM输入参数！）。
    %         用于提取高程。请参阅 'ConditionDEM' 函数以获取制作平滑DEM的选项。
    %         如果没有提供输入，代码默认为使用 mincosthydrocon 函数。
    %     concavity_method ['ref'] - 凹度选项
    %         'ref' - 使用参考凹度，用户可以使用参考凹度选项指定该值（见下文）
    %         'auto' - 函数会寻找最佳拟合凹度用于所提供的河道
    %     pick_method ['chi'] - 选择要投影的河道段的方法：
    %         'chi' - 在 chi - z 图上选择段
    %         'stream' - 在纵剖面上选择段
    %     ref_concavity [0.50] - 如果 'theta_method' 设为 'auto'，则使用参考凹度
    %     refit_streams [false] - 选项用于根据所选段的凹度重新计算chi值（true），如果您想精确匹配所选剖面段的形状，该选项很有用。
    %         仅在 'theta_method' 设置为 'auto' 时使用
    %     save_figures [false] - 选项用于在投影过程结束时保存（如果设置为true）图像
    %     interp_value [0.1] - 插值参数（0到1之间的值），用于 mincosthydrocon 中的插值（如果用户提供了平滑DEM，则不使用此参数）
    %
    % 输出：
    %     生成一个2 x n的单元数组，每个河道段（或网络中每个河道源）对应一列。第一行是该河道源的x-y坐标。
    %     第二行是包含以下信息的数组：x坐标、y坐标、距出口的距离、集水面积、chi值、凹度、真实高程、投影高程、投影高程的上限和下限。
    %     该输出还会保存为名为 'ProjectedSegments.mat' 的mat文件。
    %
    % 示例：
    %     [OUT]=StreamProjector(DEM,FD,A,S)
    %     [OUT]=StreamProjector(DEM,FD,A,S,'ref_concavity',0.55);
    %     [OUT]=StreamProjector(DEM,FD,A,S,'theta_method','auto','pick_method','stream','refit_streams',true);
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 函数由Adam M. Forte编写 - 更新日期：2018年6月18日 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	% 解析输入参数
	p = inputParser;
	p.FunctionName = 'SegmentProjector';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
	addRequired(p,'A',@(x) isa(x,'GRIDobj'));
	addRequired(p,'S',@(x) isa(x,'STREAMobj'));
	addRequired(p,'basin_num',@(x) isnumeric(x));	

	addParameter(p,'concavity_method','ref',@(x) ischar(validatestring(x,{'ref','auto'})));
	addParameter(p,'pick_method','chi',@(x) ischar(validatestring(x,{'chi','stream'})));
	addParameter(p,'smooth_distance',1000,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'ref_concavity',0.50,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'refit_streams',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'save_figures',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'conditioned_DEM',[],@(x) isa(x,'GRIDobj') || isempty(x));
	addParameter(p,'interp_value',0.1,@(x) isnumeric(x) && x>=0 && x<=1);
	addParameter(p,'out_dir',[],@(x) isdir(x));

	parse(p,DEM,FD,A,S,basin_num,varargin{:});
	DEM=p.Results.DEM;
	FD=p.Results.FD;
	S=p.Results.S;
	A=p.Results.A;
	basin_num=p.Results.basin_num;

	smooth_distance=p.Results.smooth_distance;
	theta_method=p.Results.concavity_method;
	ref_theta=p.Results.ref_concavity;
	refit_streams=p.Results.refit_streams;
	pick_method=p.Results.pick_method;
	save_figures=p.Results.save_figures;
	iv=p.Results.interp_value;
	DEMc=p.Results.conditioned_DEM;
	out_dir=p.Results.out_dir;

	if isempty(out_dir)
		out_dir=pwd;
	end

	% 查找河道源头
	ST=S;
	chix=streampoi(ST,'channelheads','ix');
	num_ch=numel(chix);

	% 平滑DEM
	if isempty(DEMc)
		zc=mincosthydrocon(ST,DEM,'interp',iv);
		DEMc=GRIDobj(DEM);
		DEMc.Z(DEMc.Z==0)=NaN;
		DEMc.Z(ST.IXgrid)=zc;
	end

	% 解析选项开关
	if strcmp(theta_method,'ref') && strcmp(pick_method,'chi')
		method=1;
	elseif strcmp(theta_method,'auto') && strcmp(pick_method,'chi')
		method=2;
	elseif strcmp(theta_method,'ref') && strcmp(pick_method,'stream')
		method=3;
	elseif strcmp(theta_method,'auto') && strcmp(pick_method,'stream')
		method=4;
	end

	OUT=cell(2,num_ch);

	% 设置字符串输入
	str1='N';

	switch method
	case 1
		% 自动计算ksn用于比较
		[auto_ksn]=KSN_Quick(DEM,A,ST,ref_theta);

		for ii=1:num_ch
			CHIX=GRIDobj(DEM);
			CHIX.Z(chix(ii))=1; CHIX.Z=logical(CHIX.Z);
			S=modify(ST,'downstreamto',CHIX);

			C=ChiCalc(S,DEMc,A,1,ref_theta);

			ak=getnal(S,auto_ksn);
			[DAvg,KsnAvg]=BinAverage(S.distance,ak,smooth_distance);
			[~,CAvg]=BinAverage(C.distance,C.chi,smooth_distance);

			f1=figure(1);
			set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');

			clf
			ax3=subplot(3,1,3);
			hold on
			pl1=plotdz(S,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
			pl2=plotdz(S,DEMc,'dunit','km','Color','k');
			xlabel('距河口距离 (km)')
			ylabel('高程 (m)')
			legend([pl1 pl2],'未平滑DEM','平滑后DEM','location','best');
			title('纵剖面图')
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(ax3);
		    end				
			hold off

			ax2=subplot(3,1,2);
			hold on
			scatter(CAvg,KsnAvg,20,'k','filled');
			xlabel('Chi')
			ylabel('k_{sn}');
			title('Chi - k_{sn}');
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(ax2);
		    end				
			hold off

			ax1=subplot(3,1,1);
			hold on
			plot(C.chi,C.elev,'-k');
			scatter(C.chi,C.elev,10,'k');
			xlabel('\chi','Color','r')
			ylabel('高程 (m)','Color','r')
			title(['\chi - Z : \theta = ' num2str(C.mn) ' : 请选择要投影的河段边界'],'Color','r')
			ax1.XColor='Red';
			ax1.YColor='Red';
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(ax1);
		    end				
			hold off

			linkaxes([ax1,ax2],'x');

			while strcmpi(str1,'N')
				[cv,e]=ginput(2);

				% 排序裂点列表并构建边界列表
				cvs=sortrows(cv);

				rc=C.chi;
				rx=C.x;
				ry=C.y;
				re=C.elev;

				lb=cvs(1);
				rb=cvs(2);

				% 裁剪出河段
				lb_chidist=sqrt(sum(bsxfun(@minus, rc, lb).^2,2));
				rb_chidist=sqrt(sum(bsxfun(@minus, rc, rb).^2,2));

				[~,lbix]=min(lb_chidist);
				[~,rbix]=min(rb_chidist);

				rcC=rc(rbix:lbix);
				reC=re(rbix:lbix);

				hold on
				p1=scatter(ax1,rcC,reC,20,'r','filled');
				hold off


	            qa=questdlg('这是您想要投影的河段吗？','古河道投影','否','是','是');
	            switch qa
	            case '是'
	                str1 = 'Y';
	            case '否'
	                str1 = 'N';
	            end

				delete(p1);
			end

			f=fit(rcC,reC,'poly1');
			cf=coeffvalues(f);
			ci=confint(f);
			ksn=cf(1);
			eint=cf(2);
			ksnl=ci(1,1);
			ksnu=ci(2,1);
			eintl=ci(1,2);
			eintu=ci(2,2);

			pred_el=(rc.*ksn)+eint;
			pred_el_u=(rc.*ksnl)+eintu;
			pred_el_l=(rc.*ksnu)+eintl;

			subplot(3,1,1)
			hold on
			plot(C.chi,pred_el,'-r','LineWidth',2);
			plot(C.chi,pred_el_u,'--r');
			plot(C.chi,pred_el_l,'--r');
			hold off

			subplot(3,1,3)
			hold on
			pl3=plot(C.distance./1000,pred_el,'-r','LineWidth',2);
			pl4=plot(C.distance./1000,pred_el_u,'--r');
			plot(C.distance./1000,pred_el_l,'--r');
			legend([pl1 pl2 pl3 pl4],'未平滑DEM','平滑后的DEM','投影河道','不确定性','location','best');
			hold off

			f2=figure(2);
			set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');
			sbplt1=subplot(2,1,1);
			hold on
			plot([0,max(C.chi)],[0,0],'--k');
			scatter(C.chi,pred_el-C.elev,10,'k')
			xlabel('\chi');
			ylabel('投影与真实剖面差值 (m)')
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(sbplt1);
		    end				
			hold off

			sbplt2=subplot(2,1,2);
			hold on
			plot([0,max(C.distance)/1000],[0,0],'--k');
			scatter(C.distance./1000,pred_el-C.elev,10,'k')
			xlabel('距离 (km)');
			ylabel('投影与真实剖面差值 (m)')
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(sbplt2);
		    end				
			hold off

			[chx,chy]=ind2coord(DEM,chix(ii));
			OUT{1,ii}=[chx chy];
			OUT{2,ii}=[C.x C.y C.distance C.area C.chi ones(size(C.chi)).*C.mn C.elev pred_el pred_el_u pred_el_l];

			if save_figures
				print(f1,'-dpdf',fullfile(out_dir,['ProjectedProfile_' num2str(basin_num) '_' num2str(ii) '.pdf']),'-bestfit');
				print(f2,'-dpdf',fullfile(out_dir,['Residual_' num2str(basin_num) '_' num2str(ii) '.pdf']),'-bestfit');
			else
				if ii<num_ch
					uiwait(msgbox('准备好继续时点击确定'))
				end
			end

			if ii<num_ch
				close(f1);
				close(f2);
			end

			% 重置字符串
			str1='N';
		end
	case 2
		for ii=1:num_ch
			CHIX=GRIDobj(DEM);
			CHIX.Z(chix(ii))=1; CHIX.Z=logical(CHIX.Z);
			S=modify(ST,'downstreamto',CHIX);

			C=ChiCalc(S,DEMc,A,1);

			% 自动计算ksn用于比较
			[auto_ksn]=KSN_Quick(DEM,A,S,C.mn);
			ak=getnal(S,auto_ksn);
			[DAvg,KsnAvg]=BinAverage(S.distance,ak,smooth_distance);
			[~,CAvg]=BinAverage(C.distance,C.chi,smooth_distance);

			f1=figure(1);
			set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');

			clf
			ax3=subplot(3,1,3);
			hold on
			pl1=plotdz(S,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
			pl2=plotdz(S,DEMc,'dunit','km','Color','k');
			xlabel('距河口距离 (km)')
			ylabel('高程 (m)')
			legend([pl1 pl2],'未平滑DEM','平滑后DEM','location','best');
			title('纵剖面图')
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(ax3);
		    end				
			hold off

			ax2=subplot(3,1,2);
			hold on
			scatter(CAvg,KsnAvg,20,'k','filled');
			xlabel('\chi')
			ylabel('k_{sn}');
			title('\chi - k_{sn}');
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(ax2);
		    end				
			hold off

			ax1=subplot(3,1,1);
			hold on
			plot(C.chi,C.elev,'-k');
			scatter(C.chi,C.elev,10,'k');
			xlabel('\chi','Color','r')
			ylabel('高程 (m)','Color','r')
			title(['\chi - Z : \theta = ' num2str(C.mn) ' : 请选择要投影的河段边界'],'Color','r')
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(ax1);
		    end				
			hold off

			linkaxes([ax1,ax2],'x');

			switch refit_streams
			case false

				while strcmpi(str1,'N')
					[cv,e]=ginput(2);

					% 排序裂点列表并构建边界列表
					cvs=sortrows(cv);

					rc=C.chi;
					rx=C.x;
					ry=C.y;
					re=C.elev;

					lb=cvs(1);
					rb=cvs(2);

					% 裁剪出河段
					lb_chidist=sqrt(sum(bsxfun(@minus, rc, lb).^2,2));
					rb_chidist=sqrt(sum(bsxfun(@minus, rc, rb).^2,2));

					[~,lbix]=min(lb_chidist);
					[~,rbix]=min(rb_chidist);

					rcC=rc(rbix:lbix);
					reC=re(rbix:lbix);

					hold on
					p1=scatter(ax1,rcC,reC,20,'r','filled');
					hold off

		            qa=questdlg('这是您想要投影的河段吗？','古河道投影','否','是','是');
		            switch qa
		            case '是'
		                str1 = 'Y';
		            case '否'
		                str1 = 'N';
		            end

					delete(p1);
				end

				f=fit(rcC,reC,'poly1');
				cf=coeffvalues(f);
				ci=confint(f);
				ksn=cf(1);
				eint=cf(2);
				ksnl=ci(1,1);
				ksnu=ci(2,1);
				eintl=ci(1,2);
				eintu=ci(2,2);

				pred_el=(rc.*ksn)+eint;
				pred_el_u=(rc.*ksnl)+eintu;
				pred_el_l=(rc.*ksnu)+eintl;

				subplot(3,1,1)
				hold on
				pl3=plot(C.chi,pred_el,'-r','LineWidth',2);
				pl4=plot(C.chi,pred_el_u,'--r');
				plot(C.chi,pred_el_l,'--r');
				legend([pl1 pl2 pl3 pl4],'未平滑DEM','平滑后的DEM','投影河道','不确定性','location','best');
				hold off

				subplot(3,1,3)
				hold on
				plot(C.distance./1000,pred_el,'-r','LineWidth',2);
				plot(C.distance./1000,pred_el_u,'--r');
				plot(C.distance./1000,pred_el_l,'--r');
				hold off

				f2=figure(2);
				set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');
				sbplt1=subplot(2,1,1);
				hold on
				plot([0,max(C.chi)],[0,0],'--k');
				scatter(C.chi,pred_el-C.elev,10,'k')
				xlabel('\chi');
				ylabel('投影与真实剖面差值 (m)')
				if ~verLessThan('matlab','9.5')
			        disableDefaultInteractivity(sbplt1);
			    end					
				hold off

				sbplt2=subplot(2,1,2);
				hold on
				plot([0,max(C.distance)/1000],[0,0],'--k');
				scatter(C.distance./1000,pred_el-C.elev,10,'k')
				xlabel('距离 (km)');
				ylabel('投影与真实剖面差值 (m)')
				if ~verLessThan('matlab','9.5')
			        disableDefaultInteractivity(sbplt2);
			    end					
				hold off

				[chx,chy]=ind2coord(DEM,chix(ii));
				OUT{1,ii}=[chx chy];
				OUT{2,ii}=[C.x C.y C.distance C.area C.chi ones(size(C.chi)).*C.mn C.elev pred_el pred_el_u pred_el_l];

			case true

				while strcmpi(str1,'N')
					[cv,e]=ginput(2);

					% 排序裂点列表并构建边界列表
					cvs=sortrows(cv);

					rc0=C.chi;
					rx0=C.x;
					ry0=C.y;
					re0=C.elev;

					lb=cvs(1);
					rb=cvs(2);

					% 裁剪出河段
					lb_chidist=sqrt(sum(bsxfun(@minus, rc0, lb).^2,2));
					rb_chidist=sqrt(sum(bsxfun(@minus, rc0, rb).^2,2));

					[~,lbix]=min(lb_chidist);
					[~,rbix]=min(rb_chidist);

					% 裁剪出河段
					lbx=rx0(lb_chidist==min(lb_chidist));
					lby=ry0(lb_chidist==min(lb_chidist));

					rbx=rx0(rb_chidist==min(rb_chidist));
					rby=ry0(rb_chidist==min(rb_chidist));	

					lix=coord2ind(DEM,lbx,lby);
					LIX=GRIDobj(DEM);
					LIX.Z(lix)=1;
					[lixmat,X,Y]=GRIDobj2mat(LIX);
					lixmat=logical(lixmat);
					LIX=GRIDobj(X,Y,lixmat);	

					rix=coord2ind(DEM,rbx,rby);
					RIX=GRIDobj(DEM);
					RIX.Z(rix)=1;
					[rixmat,X,Y]=GRIDobj2mat(RIX);
					rixmat=logical(rixmat);
					RIX=GRIDobj(X,Y,rixmat);	

					Seg=modify(S,'downstreamto',RIX);
					Seg=modify(Seg,'upstreamto',LIX);

					% 查找河段凹度
					Csegrf=ChiCalc(Seg,DEMc,A,1);

					% 使用新凹度重新计算整个河道的chi值
					CN=ChiCalc(S,DEMc,A,1,Csegrf.mn);
					rc=CN.chi;
					rx=CN.x;
					ry=CN.y;
					re=CN.elev;

					rcC0=rc0(rbix:lbix);
					reC0=re0(rbix:lbix);

					rcC=rc(rbix:lbix);
					reC=re(rbix:lbix);

					hold on
					p1=scatter(ax1,rcC0,reC0,20,'r','filled');
					hold off

		            qa=questdlg('这是您想要投影的河段?','古河道投影','No','Yes','Yes');
		            switch qa
		            case 'Yes'
		                str1 = 'Y';
		            case 'No'
		                str1 = 'N';
		            end

					delete(p1);
				end

				% 自动计算ksn以供比较
				[auto_ksn]=KSN_Quick(DEM,A,S,Csegrf.mn);
				ak=getnal(S,auto_ksn);
				[DAvg,KsnAvg]=BinAverage(S.distance,ak,smooth_distance);
				[~,CAvg]=BinAverage(CN.distance,CN.chi,smooth_distance);

				f=fit(rcC,reC,'poly1');
				cf=coeffvalues(f);
				ci=confint(f);
				ksn=cf(1);
				eint=cf(2);
				ksnl=ci(1,1);
				ksnu=ci(2,1);
				eintl=ci(1,2);
				eintu=ci(2,2);

				pred_el=(rc.*ksn)+eint;
				pred_el_u=(rc.*ksnl)+eintu;
				pred_el_l=(rc.*ksnu)+eintl;

				f1=figure(1);
				clf; cla;
				set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');

				ax1=subplot(3,1,1);
				hold on
				plot(CN.chi,CN.elev,'-k');
				scatter(CN.chi,CN.elev,10,'k');
				plot(CN.chi,pred_el,'-r','LineWidth',2);
				plot(CN.chi,pred_el_u,'--r');
				plot(CN.chi,pred_el_l,'--r');
				xlabel('\chi')
				ylabel('高程（米）')
				title(['\chi-Z 关系 : \theta = ' num2str(CN.mn)],'Color','r')
				if ~verLessThan('matlab','9.5')
			        disableDefaultInteractivity(ax1);
			    end	
				hold off

				ax2=subplot(3,1,2);
				hold on
				scatter(CAvg,KsnAvg,20,'k','filled');
				xlabel('Chi值')
				ylabel('计算k_{sn}');
				title('Chi-k_{sn}关系');
				if ~verLessThan('matlab','9.5')
			        disableDefaultInteractivity(ax2);
			    end					
				hold off

				ax3=subplot(3,1,3);
				hold on
				pl1=plotdz(S,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
				pl2=plotdz(S,DEMc,'dunit','km','Color','k');
				pl3=plot(CN.distance./1000,pred_el,'-r','LineWidth',2);
				pl4=plot(CN.distance./1000,pred_el_u,'--r');
				plot(CN.distance./1000,pred_el_l,'--r');
				xlabel('距河口距离（km）')
				ylabel('高程（m）')
				legend([pl1 pl2 pl3 pl4],'未平滑DEM','平滑DEM','投影河道','不确定区间','location','best');
				title('纵剖面图')
				if ~verLessThan('matlab','9.5')
			        disableDefaultInteractivity(ax3);
			    end					
				hold off	

				linkaxes([ax1,ax2],'x');

				f2=figure(2);
				set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');
				sbplt1=subplot(2,1,1);
				hold on
				plot([0,max(CN.chi)],[0,0],'--k');
				scatter(CN.chi,pred_el-CN.elev,10,'k')
				xlabel('\chi值');
				ylabel('投影剖面与真实剖面高差（m）')
				if ~verLessThan('matlab','9.5')
			        disableDefaultInteractivity(sbplt1);
			    end	
				hold off

				sbplt2=subplot(2,1,2);
				hold on
				plot([0,max(CN.distance)/1000],[0,0],'--k');
				scatter(CN.distance./1000,pred_el-CN.elev,10,'k')
				xlabel('距离（km）');
				ylabel('投影剖面与真实剖面高差（m）')
				if ~verLessThan('matlab','9.5')
			        disableDefaultInteractivity(sbplt2);
			    end	
				hold off

				[chx,chy]=ind2coord(DEM,chix(ii));
				OUT{1,ii}=[chx chy];
				OUT{2,ii}=[CN.x CN.y CN.distance CN.area CN.chi ones(size(CN.chi)).*CN.mn CN.elev pred_el pred_el_u pred_el_l];

			end


			if save_figures
				print(f1,'-dpdf',fullfile(out_dir,['ProjectedProfile_' num2str(basin_num) '_' num2str(ii) '.pdf']),'-bestfit');
				print(f2,'-dpdf',fullfile(out_dir,['Residual_' num2str(basin_num) '_' num2str(ii) '.pdf']),'-bestfit');
			else
				if ii<num_ch
					uiwait(msgbox('准备好后点击确定继续'))
				end
			end

			if ii<num_ch
				close(f1);
				close(f2);
			end

			% 重置输出字符串
			str1='N';
		end
	case 3
		% 自动计算ksn用于比较
		[auto_ksn]=KSN_Quick(DEM,A,ST,ref_theta);

		for ii=1:num_ch
			CHIX=GRIDobj(DEM);
			CHIX.Z(chix(ii))=1; CHIX.Z=logical(CHIX.Z);
			S=modify(ST,'downstreamto',CHIX);

			C=ChiCalc(S,DEMc,A,1,ref_theta);

			ak=getnal(S,auto_ksn);
			[DAvg,KsnAvg]=BinAverage(S.distance,ak,smooth_distance);

			f1=figure(1);
			set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');

			clf

			ax3=subplot(3,1,1);
			hold on
			plot(C.chi,C.elev,'-k');
			scatter(C.chi,C.elev,10,'k');
			xlabel('\chi值')
			ylabel('高程（m）')
			title('\chi-Z关系图')
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(ax3);
		    end				
			hold off

			ax2=subplot(3,1,2);
			hold on
			scatter(DAvg./1000,KsnAvg,20,'k','filled');
			xlabel('距离（km）')
			ylabel('计算k_{sn}');
			title('距离-k_{sn}关系');
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(ax2);
		    end				
			hold off

			ax1=subplot(3,1,3);
			hold on
			pl1=plotdz(S,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
			pl2=plotdz(S,DEMc,'dunit','km','Color','k');
			xlabel('距河口距离（km）','Color','r')
			ylabel('高程（m）','Color','r')
			legend([pl1 pl2],'未平滑DEM','平滑DEM','location','best');
			title(['纵剖面图 : \theta = ' num2str(C.mn) ' : 请选择要投影的河段范围'],'Color','r')
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(ax1);
		    end				
			hold off

			linkaxes([ax1,ax2],'x');

			while strcmpi(str1,'N')
				[d,e]=ginput(2);
				d=d*1000;

				% 对裂点列表进行排序并构建边界列表
				ds=sortrows(d);

				rd=C.distance;
				rx=C.x;
				ry=C.y;
				rc=C.chi;
				re=C.elev;

				lb=ds(1);
				rb=ds(2);

				lb_dist=sqrt(sum(bsxfun(@minus, rd, lb).^2,2));
				rb_dist=sqrt(sum(bsxfun(@minus, rd, rb).^2,2));

				[~,lbix]=min(lb_dist);
				[~,rbix]=min(rb_dist);

				rcC=rc(rbix:lbix);
				reC=re(rbix:lbix);
				rdC=rd(rbix:lbix);

				hold on
				p1=scatter(ax1,rdC/1000,reC,20,'r','filled');
				hold off

	            qa=questdlg('确认要投影此河段吗？','古河道投影','否','是','是');
	            switch qa
	            case '是'
	                str1 = 'Y';
	            case '否'
	                str1 = 'N';
	            end

				delete(p1);
			end				

			f=fit(rcC,reC,'poly1');
			cf=coeffvalues(f);
			ci=confint(f);
			ksn=cf(1);
			eint=cf(2);
			ksnl=ci(1,1);
			ksnu=ci(2,1);
			eintl=ci(1,2);
			eintu=ci(2,2);

			pred_el=(rc.*ksn)+eint;
			pred_el_u=(rc.*ksnl)+eintu;
			pred_el_l=(rc.*ksnu)+eintl;

			subplot(3,1,1)
			hold on
			plot(C.chi,pred_el,'-r','LineWidth',2);
			plot(C.chi,pred_el_u,'--r');
			plot(C.chi,pred_el_l,'--r');
			hold off

			subplot(3,1,3)
			hold on
			pl3=plot(C.distance./1000,pred_el,'-r','LineWidth',2);
			pl4=plot(C.distance./1000,pred_el_u,'--r');
			plot(C.distance./1000,pred_el_l,'--r');
			legend([pl1 pl2 pl3 pl4],'未平滑DEM','平滑DEM','投影河道','不确定区间','location','best');
			hold off

			f2=figure(2);
			set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');
			sbplt1=subplot(2,1,1);
			hold on
			plot([0,max(C.chi)],[0,0],'--k');
			scatter(C.chi,pred_el-C.elev,10,'k')
			xlabel('\chi值');
			ylabel('投影与真实剖面高差（m）')
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(sbplt1);
		    end				
			hold off

			sbplt2=subplot(2,1,2);
			hold on
			plot([0,max(C.distance)/1000],[0,0],'--k');
			scatter(C.distance./1000,pred_el-C.elev,10,'k')
			xlabel('距离（km）');
			ylabel('投影与真实剖面高差（m）')
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(sbplt2);
		    end				
			hold off

			[chx,chy]=ind2coord(DEM,chix(ii));
			OUT{1,ii}=[chx chy];
			OUT{2,ii}=[C.x C.y C.distance C.area C.chi ones(size(C.chi)).*C.mn C.elev pred_el pred_el_u pred_el_l];

			if save_figures
				print(f1,'-dpdf',fullfile(out_dir,['ProjectedProfile_' num2str(basin_num) '_' num2str(ii) '.pdf']),'-bestfit');
				print(f2,'-dpdf',fullfile(out_dir,['Residual_' num2str(basin_num) '_' num2str(ii) '.pdf']),'-bestfit');
			else
				if ii<num_ch
					uiwait(msgbox('准备好后点击确定继续'))
				end
			end

			if ii<num_ch
				close(f1);
				close(f2);
			end

			% 重置字符串
			str1='N';
		end

	case 4
		for ii=1:num_ch
			CHIX=GRIDobj(DEM);
			CHIX.Z(chix(ii))=1; CHIX.Z=logical(CHIX.Z);
			S=modify(ST,'downstreamto',CHIX);

			C=ChiCalc(S,DEMc,A,1);

			% 自动计算ksn用于比较
			[auto_ksn]=KSN_Quick(DEM,A,S,C.mn);
			ak=getnal(S,auto_ksn);
			[DAvg,KsnAvg]=BinAverage(S.distance,ak,smooth_distance);

			f1=figure(1);
			set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');

			clf

			ax3=subplot(3,1,1);
			hold on
			plot(C.chi,C.elev,'-k');
			scatter(C.chi,C.elev,10,'k');
			xlabel('\chi值')
			ylabel('高程（m）')
			title('\chi-Z关系图')
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(ax3);
		    end				
			hold off

			ax2=subplot(3,1,2);
			hold on
			scatter(DAvg./1000,KsnAvg,20,'k','filled');
			xlabel('距离（km）')
			ylabel('计算k_{sn}');
			title('距离-k_{sn}关系');
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(ax2);
		    end				
			hold off

			ax1=subplot(3,1,3);
			hold on
			pl1=plotdz(S,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
			pl2=plotdz(S,DEMc,'dunit','km','Color','k');
			xlabel('距河口距离（km）','Color','r')
			ylabel('高程（m）','Color','r')
			legend([pl1 pl2],'未平滑DEM','平滑DEM','location','best');
			title(['纵剖面图 : \theta = ' num2str(C.mn) ' : 请选择要投影的河段范围'],'Color','r')
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(ax1);
		    end				
			hold off

			linkaxes([ax1,ax2],'x');

			switch refit_streams
			case false

				while strcmpi(str1,'N')
					[d,e]=ginput(2);
					d=d*1000;

					% 对裂点列表排序并构建边界
					ds=sortrows(d);

					rd=C.distance;
					rx=C.x;
					ry=C.y;
					rc=C.chi;
					re=C.elev;

					lb=ds(1);
					rb=ds(2);

					lb_dist=sqrt(sum(bsxfun(@minus, rd, lb).^2,2));
					rb_dist=sqrt(sum(bsxfun(@minus, rd, rb).^2,2));

					[~,lbix]=min(lb_dist);
					[~,rbix]=min(rb_dist);

					rcC=rc(rbix:lbix);
					reC=re(rbix:lbix);
					rdC=rd(rbix:lbix);

					hold on
					p1=scatter(ax1,rdC/1000,reC,20,'r','filled');
					hold off

		            qa=questdlg('确认要投影此河段吗？','古河道投影','否','是','是');
		            switch qa
		            case '是'
		                str1 = 'Y';
		            case '否'
		                str1 = 'N';
		            end

					delete(p1);
				end

				f=fit(rcC,reC,'poly1');
				cf=coeffvalues(f);
				ci=confint(f);
				ksn=cf(1);
				eint=cf(2);
				ksnl=ci(1,1);
				ksnu=ci(2,1);
				eintl=ci(1,2);
				eintu=ci(2,2);

				pred_el=(rc.*ksn)+eint;
				pred_el_u=(rc.*ksnl)+eintu;
				pred_el_l=(rc.*ksnu)+eintl;

				subplot(3,1,1)
				hold on
				plot(C.chi,pred_el,'-r','LineWidth',2);
				plot(C.chi,pred_el_u,'--r');
				plot(C.chi,pred_el_l,'--r');
				hold off

				subplot(3,1,3)
				hold on
				pl3=plot(C.distance./1000,pred_el,'-r','LineWidth',2);
				pl4=plot(C.distance./1000,pred_el_u,'--r');
				plot(C.distance./1000,pred_el_l,'--r');
				legend([pl1 pl2 pl3 pl4],'未平滑DEM','平滑DEM','投影河道','不确定区间','location','best');
				hold off

				f2=figure(2);
				set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');
				sbplt1=subplot(2,1,1);
				hold on
				plot([0,max(C.chi)],[0,0],'--k');
				scatter(C.chi,pred_el-C.elev,10,'k')
				xlabel('\chi值');
				ylabel('投影与真实剖面高差（m）')
				if ~verLessThan('matlab','9.5')
			        disableDefaultInteractivity(sbplt1);
			    end					
				hold off

				sbplt2=subplot(2,1,2);
				hold on
				plot([0,max(C.distance)/1000],[0,0],'--k');
				scatter(C.distance./1000,pred_el-C.elev,10,'k')
				xlabel('距离（km）');
				ylabel('投影与真实剖面高差（m）')
				if ~verLessThan('matlab','9.5')
			        disableDefaultInteractivity(sbplt2);
			    end					
				hold off

				[chx,chy]=ind2coord(DEM,chix(ii));
				OUT{1,ii}=[chx chy];
				OUT{2,ii}=[C.x C.y C.distance C.area C.chi ones(size(C.chi)).*C.mn C.elev pred_el pred_el_u pred_el_l];

			case true

				while strcmpi(str1,'N')
					[d,e]=ginput(2);
					d=d*1000;

					% 对裂点列表进行排序并构建边界
					ds=sortrows(d);

					rd=C.distance;
					rx=C.x;
					ry=C.y;
					rc=C.chi;
					re=C.elev;

					lb=ds(1);
					rb=ds(2);

					lb_dist=sqrt(sum(bsxfun(@minus, rd, lb).^2,2));
					rb_dist=sqrt(sum(bsxfun(@minus, rd, rb).^2,2));

					[~,lbix]=min(lb_dist);
					[~,rbix]=min(rb_dist);

					lbx=rx(lb_dist==min(lb_dist));
					lby=ry(lb_dist==min(lb_dist));

					rbx=rx(rb_dist==min(rb_dist));
					rby=ry(rb_dist==min(rb_dist));	

					lix=coord2ind(DEM,lbx,lby);
					LIX=GRIDobj(DEM);
					LIX.Z(lix)=1;
					[lixmat,X,Y]=GRIDobj2mat(LIX);
					lixmat=logical(lixmat);
					LIX=GRIDobj(X,Y,lixmat);	

					rix=coord2ind(DEM,rbx,rby);
					RIX=GRIDobj(DEM);
					RIX.Z(rix)=1;
					[rixmat,X,Y]=GRIDobj2mat(RIX);
					rixmat=logical(rixmat);
					RIX=GRIDobj(X,Y,rixmat);	

					Seg=modify(S,'downstreamto',RIX);
					Seg=modify(Seg,'upstreamto',LIX);

					Csegrf=ChiCalc(Seg,DEMc,A,1);

					CN=ChiCalc(S,DEMc,A,1,Csegrf.mn);

					rc=CN.chi;
					rx=CN.x;
					ry=CN.y;
					re=CN.elev;

					rcC=rc(rbix:lbix);
					reC=re(rbix:lbix);
					rdC=rd(rbix:lbix);

					hold on
					p1=scatter(ax1,rdC/1000,reC,20,'r','filled');
					hold off

		            qa=questdlg('确认要投影此河段吗？','古河道投影','否','是','是');
		            switch qa
		            case '是'
		                str1 = 'Y';
		            case '否'
		                str1 = 'N';
		            end

					delete(p1);
				end

				% 自动计算ksn用于比较
				[auto_ksn]=KSN_Quick(DEM,A,S,Csegrf.mn);
				ak=getnal(S,auto_ksn);
				[DAvg,KsnAvg]=BinAverage(S.distance,ak,smooth_distance);
				[~,CAvg]=BinAverage(CN.distance,CN.chi,smooth_distance);

				f=fit(rcC,reC,'poly1');
				cf=coeffvalues(f);
				ci=confint(f);
				ksn=cf(1);
				eint=cf(2);
				ksnl=ci(1,1);
				ksnu=ci(2,1);
				eintl=ci(1,2);
				eintu=ci(2,2);

				pred_el=(rc.*ksn)+eint;
				pred_el_u=(rc.*ksnl)+eintu;
				pred_el_l=(rc.*ksnu)+eintl;

				f1=figure(1);
				clf; cla;
				set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');

				ax3=subplot(3,1,1);
				hold on
				plot(CN.chi,CN.elev,'-k');
				scatter(CN.chi,CN.elev,10,'k');
				plot(CN.chi,pred_el,'-r','LineWidth',2);
				plot(CN.chi,pred_el_u,'--r');
				plot(CN.chi,pred_el_l,'--r');
				xlabel('\chi值')
				ylabel('高程（m）')
				title('\chi-Z关系图')
				if ~verLessThan('matlab','9.5')
			        disableDefaultInteractivity(ax3);
			    end					
				hold off

				ax2=subplot(3,1,2);
				hold on
				scatter(DAvg./1000,KsnAvg,20,'k','filled');
				xlabel('距离（km）')
				ylabel('计算k_{sn}');
				title('距离-k_{sn}关系');
				if ~verLessThan('matlab','9.5')
			        disableDefaultInteractivity(ax2);
			    end					
				hold off

				ax1=subplot(3,1,3);
				hold on
				pl1=plotdz(S,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
				pl2=plotdz(S,DEMc,'dunit','km','Color','k');
				pl3=plot(CN.distance./1000,pred_el,'-r','LineWidth',2);
				pl4=plot(CN.distance./1000,pred_el_u,'--r');
				plot(CN.distance./1000,pred_el_l,'--r');
				xlabel('距河口距离（km）')
				ylabel('高程（m）')
				legend([pl1 pl2 pl3 pl4],'未平滑DEM','平滑DEM','投影河道','不确定区间','location','best');
				title(['纵剖面图 : \theta = ' num2str(CN.mn)],'Color','r')
				if ~verLessThan('matlab','9.5')
			        disableDefaultInteractivity(ax1);
			    end					
				hold off	

				linkaxes([ax1,ax2],'x');

				f2=figure(2);
				set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');
				sbplt1=subplot(2,1,1);
				hold on
				plot([0,max(CN.chi)],[0,0],'--k');
				scatter(CN.chi,pred_el-CN.elev,10,'k')
				xlabel('\chi值');
				ylabel('投影与真实剖面高差（m）')
				if ~verLessThan('matlab','9.5')
			        disableDefaultInteractivity(sbplt1);
			    end					
				hold off

				sbplt2=subplot(2,1,2);
				hold on
				plot([0,max(CN.distance)/1000],[0,0],'--k');
				scatter(CN.distance./1000,pred_el-CN.elev,10,'k')
				xlabel('距离（km）');
				ylabel('投影与真实剖面高差（m）')
				if ~verLessThan('matlab','9.5')
			        disableDefaultInteractivity(sbplt2);
			    end					
				hold off

				[chx,chy]=ind2coord(DEM,chix(ii));
				OUT{1,ii}=[chx chy];
				OUT{2,ii}=[CN.x CN.y CN.distance CN.area CN.chi ones(size(CN.chi)).*CN.mn CN.elev pred_el pred_el_u pred_el_l];
			end

			if save_figures
				print(f1,'-dpdf',fullfile(out_dir,['ProjectedProfile_' num2str(basin_num) '_' num2str(ii) '.pdf']),'-bestfit');
				print(f2,'-dpdf',fullfile(out_dir,['Residual_' num2str(basin_num) '_' num2str(ii) '.pdf']),'-bestfit');
			else
				if ii<num_ch
					uiwait(msgbox('准备好后点击确定继续'))
				end
			end

			if ii<num_ch
				close(f1);
				close(f2);
			end
			str1='N';
		end
	end

	save(fullfile(out_dir,['ProjectedSegments_' num2str(basin_num) '.mat']),'OUT','-v7.3');
end

function [OUT]=ChiCalc(S,DEM,A,a0,varargin)
% 修改自Wolfgang Schwanghart的chiplot函数，移除非必要功能
% 通过等间距chi插值避免高集水区数据聚集对ksn拟合的影响

mnmethod='ls'; % 计算最佳拟合凹度的选项（若无参考值）

% 河流网络总节点数
nrc = numel(S.x);
M   = sparse(double(S.ix),double(S.ixc),true,nrc,nrc);
% 查找出口点
outlet = sum(M,2) == 0 & sum(M,1)'~=0;
if nnz(outlet)>1
    error('河流网络不能有多个出口');
end

% 节点高程值
zx   = double(DEM.Z(S.IXgrid));
% 出口点高程
zb   = double(DEM.Z(S.IXgrid(outlet)));
% 公式6b括号内的项
a    = double(a0./(A.Z(S.IXgrid)*(A.cellsize.^2)));
% 上游方向累计水平距离
x    = S.distance;
Lib = true(size(x));

% 若无凹度值则寻找最佳拟合值，否则使用给定值
if isempty(varargin)
    mn0  = 0.5; % 初始值
    mn   = fminsearch(@mnfit,mn0);
else
	mn=varargin{1};
end

% 计算chi值
chi = netcumtrapz(x,a.^mn,S.ix,S.ixc);

% 使用三次样条插值重采样chi-高程关系
chiF=chi(Lib);
zabsF=zx(Lib)-zb;

chiS=linspace(0,max(chiF),numel(chiF)).';
zS=spline(chiF,zabsF,chiS);

% 计算ksn
ksn = chiS\(zS);

OUT=struct;
OUT.ks   = ksn;
OUT.mn   = mn;
[OUT.x,...
 OUT.y,...
 OUT.chi,...
 OUT.elev,...
 OUT.elevbl,...
 OUT.distance,...
 OUT.pred,...
 OUT.area] = STREAMobj2XY(S,chi,DEM,zx-zb,S.distance,ksn*chi,A.*(A.cellsize^2));
 OUT.res = OUT.elevbl - OUT.pred;

	% 计算mn比值的嵌套函数
	function sqres = mnfit(mn)
		% 用给定mn计算chi值并向上游积分
		CHI = netcumtrapz(x(Lib),a(Lib).^mn,S.ix,S.ixc);%*ab.^mn
		% 标准化变量
		CHI = CHI ./ max(CHI);
		z   = zx(Lib)-zb;
		z   = z./max(z);
        switch mnmethod
            case 'ls' % 最小二乘法
                sqres = sum((CHI - z).^2);
            case 'lad' % 最小绝对偏差
                sqres = sum(sqrt(abs(CHI-z)));
        end
	end

end

function z = netcumtrapz(x,y,ix,ixc)
% 沿上游方向在定向树状网络中进行梯形积分

z = zeros(size(x));
for lp = numel(ix):-1:1;
    z(ix(lp)) = z(ixc(lp)) + (y(ixc(lp))+(y(ix(lp))-y(ixc(lp)))/2) *(abs(x(ixc(lp))-x(ix(lp))));
end
end

function [Xavg,Yavg]=BinAverage(X,Y,bin_size);

	ix=~isnan(X);
	X=X(ix); Y=Y(ix);

	minX=min(X);
	maxX=max(X);

	b=[minX:bin_size:maxX+bin_size];

	try
		[~,~,idx]=histcounts(X,b);
	catch
		[~,idx]=histc(X,b);
	end

	Xavg=accumarray(idx(:),X,[],@mean);
	Yavg=accumarray(idx(:),Y,[],@mean);
end

function [ksn]=KSN_Quick(DEM,A,S,theta_ref)

	zc=mincosthydrocon(S,DEM,'interp',0.1);
	DEMc=GRIDobj(DEM);
	DEMc.Z(DEMc.Z==0)=NaN;
	DEMc.Z(S.IXgrid)=zc;
	G=gradient8(DEMc);

	ksn=G./(A.*(A.cellsize^2)).^(-theta_ref);
	
end