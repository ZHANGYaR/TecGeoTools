function [nal,res,ms,brkPnts,param]=AutoKsnProfiler(DEM,FD,A,S,thresh_ratio,varargin)
    % 用法：
    %   [nal, res, ms, brkPnts, param] = AutoKsnProfiler(DEM, FD, A, S, thresh_ratio);
    %   [nal, res, ms, brkPnts, param] = AutoKsnProfiler(DEM, FD, A, S, thresh_ratio, 'name', value);
    %
    % 描述：
    %   这是一个自动化 ksn 拟合的原型函数，其功能类似于 KsnProfiler。该函数借助 MATLAB 的 'ischange' 函数来查找分析 ksn 值中的断点。
    %   函数的敏感度，即 ksn 值的变化如何被判定为裂点，由 threshold_ratio 参数掌控。用户能够通过可选的 'explore_param' 选项，
    %   探究不同 threshold_ratio 值、参考凹度和段长对结果的影响。
    %
    % 必需输入：
    %   DEM - DEM 网格对象（假定为未平滑处理的 DEM）
    %   A - 流量累积的 GRID 对象
    %   S - STREAM 对象
    %   thresh_ratio - 取值范围在 0 到 1 之间，此参数控制函数对分箱 ksn 相对于距离的运行均值变化的敏感度。
    %                 数值越接近 0，阈值越低，检测到的断点就越多。
    %
    % 可选输入：
    %   prev_param [] - 可选项，用于包含上一次运行 'AutoKsnProfiler' 时的 'param' 输出内容。
    %                   这会覆盖其他用户提供的或者默认的 'thresh_ratio'、'ref_concavity' 和 'segment_length' 值。
    %                   你仍可使用 'explore_param' 选项来修改存储在 'prev_param' 输入中的值，这样做的目的是便于在多个数据集上
    %                   以完全相同的选项运行 AutoKsnProfiler。
    %   ref_concavity [0.50] - 计算 ksn 时所使用的参考凹度（需为正值）
    %   calc_concavity [false] - 逻辑标志，若设置为 true，会为每个流段计算最佳拟合凹度。此操作会显著降低函数运行速度，
    %                             因此仅当你特别关注单个段的凹度时，才建议将其设置为 true。
    %   segment_length [1000] - 以地图单位衡量的用于对 ksn 和距离数据进行分箱的长度，用于拟合均值 ksn 值。
    %                            同时也作为生成输出地图结构时平均值的长度尺度。
    %   plot_example [false] - 逻辑标志，设置为 true 时，会开启单条流的轮廓视图图展示结果。默认情况下，
    %                           所选的流将是提供的网络中最长的干流。你可以通过为 'channeloi' 提供有效的输入来更改此选择。
    %                           这是一种静态可视化 'thresh_ratio'、'ref_concavity' 和 'segment_length' 对结果选择影响的方式。
    %                           若要更动态地探索这些值的影响，请参考 'explore_param' 选项。
    %   explore_param [false] - 逻辑标志，设置为 true 时，会开启单条流的轮廓视图图展示结果，并且可以交互式地
    %                           更改 'thresh_ratio'、'ref_concavity' 和 'segment_length' 的值。默认情况下，
    %                           所选的流将是提供的网络中最长的干流。你可以通过为 'channeloi' 提供有效的输入来更改此选择。
    %   channeloi [] - 一个 1 x 2 的数组，包含用于选择流的河道头的 x-y 坐标。当 'plot_example' 或 'explore_param' 设置为 true 时，
    %                  会使用该流。若这两个选项都为 false，则此输入对代码的行为无影响。
    %   conditioned_DEM [] - 可提供一个平滑处理的 DEM 供该函数使用（注意不要将条件化后的 DEM 作为主必需的 DEM 输入！），
    %                        此 DEM 将用于提取高程信息。有关制作水文条件化 DEM 的选项，请参考 'ConditionDEM' 函数。
    %                        若未提供输入，代码将默认使用 mincosthydrocon 函数。
    %   interp_value [0.1] - 用于 mincosthydrocon 的插值参数的值（范围在 0 到 1 之间）（若用户提供了条件化 DEM，则此参数不使用）。
    %   generate_shape [false] - 逻辑标志，用于控制是否创建流网络和断点的 shapefile 文件。
    %   shape_name ['auto'] - 导出的 shapefile 文件的名称前缀，名称中不能包含空格，以确保其为 ArcGIS 的有效名称，
    %                         并且不应包含 '.shp' 后缀。
    %
    % 输出：
    %   nal - 自动拟合 ksn 值的节点属性列表
    %   res - ksn 拟合的节点属性列表的残差（单位为米，即预测高程减去真实高程）
    %   ms - 自动拟合 ksn 值的地图结构（适合使用 shapewrite 函数输出 shapefile 文件）
    %   brkPnts - 一个 n x 4 的数组，包含函数识别出的断点的 x-y 位置、流距离和高程信息
    %             （即拟合 ksn 值之间的界限）
    %   param - 一个结构体，用于存储生成输出数据时所使用的 thresh_ratio、ref_concavity 和 segment_length 值。
    %           若你使用了可选的 'explore_param' 选项更改了这些参数，此结构体将非常有用。
    %
    % 注意：
    %   此函数依赖于 MATLAB 的 'ischange' 函数，该函数于 MATLAB 2017b 版本引入。若你使用的是较旧版本的 MATLAB，
    %   此函数将无法正常工作（会报错）。建议考虑使用 TAK 的编译版本。
    %
    % 示例用法和结果绘图：
    %   [nal, res, ms, brkPnts, param] = AutoKsnProfiler(DEM, A, S, 0.1);
    %   % 绘制带有断点的 ksn 彩色流图
    %   plotc(S, nal); hold on; colorbar; scatter(brkPnts(:,1), brkPnts(:,2), 10, 'k', 'filled'); hold off; 
    %   % 绘制残差的彩色流图
    %   plotc(S, res); colorbar;
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 函数由 Adam M. Forte 编写 - 更新日期：2019 年 04 月 02 日 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 解析输入参数
    p = inputParser;         
    p.FunctionName = 'AutoKsnProfiler';
    addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
    addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
    addRequired(p,'A', @(x) isa(x,'GRIDobj'));
    addRequired(p,'S',@(x) isa(x,'STREAMobj'));
    addRequired(p,'thresh_ratio',@(x) isscalar(x) && isnumeric(x));

    addParameter(p,'segment_length',1000,@(x) isscalar(x) && isnumeric(x));
    addParameter(p,'ref_concavity',0.50,@(x) isscalar (x) && isnumeric(x));
    addParameter(p,'calc_concavity',false,@(x) isscalar(x) && islogical(x));
    addParameter(p,'plot_example',false,@(x) isscalar(x) && islogical(x));
    addParameter(p,'interp_value',0.1,@(x) isnumeric(x) && x>=0 && x<=1);
    addParameter(p,'conditioned_DEM',[],@(x) isa(x,'GRIDobj') || isempty(x));
    addParameter(p,'explore_param',false,@(x) isscalar(x) && islogical(x));
    addParameter(p,'channeloi',[],@(x) isnumeric(x) && size(x,1)==1 && size(x,2)==2);
    addParameter(p,'shape_name','auto',@(x) ischar(x));
    addParameter(p,'generate_shape',false,@(x) isscalar(x) && islogical(x));
    addParameter(p,'prev_param',[],@(x) isstruct(x));

    parse(p,DEM,FD,A,S,thresh_ratio,varargin{:});
    DEM=p.Results.DEM;
    FD=p.Results.FD;
    A=p.Results.A;
    S=p.Results.S;
    thresh_ratio=p.Results.thresh_ratio;

    segment_length=p.Results.segment_length;
    theta=p.Results.ref_concavity;
    plot_example=p.Results.plot_example;
    iv=p.Results.interp_value;
    DEMc=p.Results.conditioned_DEM;
    explore_param=p.Results.explore_param;
    channeloi=p.Results.channeloi;
    shape_name=p.Results.shape_name;
    generate_shape=p.Results.generate_shape;
    prev_param=p.Results.prev_param;
    cc=p.Results.calc_concavity;

    % 检查 MATLAB 版本，因为使用了 ischange 函数
    if verLessThan('matlab','9.3')
        error('此函数需要 MATLAB 2017a 或更高版本')
    end

    % 若未提供条件化 DEM，则生成一个
    if isempty(DEMc)
        zc=mincosthydrocon(S,DEM,'interp',iv);
        DEMc=GRIDobj(DEM);
        DEMc.Z(DEMc.Z==0)=NaN;
        DEMc.Z(S.IXgrid)=zc;
    end

    % 若提供了 'prev_param' 输入，则加载并覆盖其他重要参数的选择
    if ~isempty(prev_param)
        thresh_ratio=prev_param.thresh_ratio;
        theta=prev_param.reference_concavity;
        segment_length=prev_param.segment_length;
    end

    % 启动对话框以更改参数值
    if explore_param
        % 提取示例河道
        if isempty(channeloi)
            ST=trunk(klargestconncomps(S,1));
        else
            chxy=streampoi(S,'channelheads','xy');
            [~,mix]=min(hypot(chxy(:,1)-channeloi(:,1),chxy(:,2)-channeloi(:,2)));
            chix=coord2ind(DEM,chxy(mix,1),chxy(mix,2));
            CHIX=GRIDobj(DEM);
            CHIX.Z(chix)=true;
            ST=modify(S,'downstreamto',CHIX);
        end

        % 使用初始值运行
        [fig_handle]=PlotFunc(DEMc,FD,ST,A,segment_length,theta,thresh_ratio);
        % 首次询问用户
        qa1=questdlg('接受这些参数吗？','拟合参数','更改值','接受值','接受值');

        switch qa1
            case '接受值'
                % 关闭图形窗口，继续使用原始值运行代码
                close(fig_handle);			
            case '更改值'
                % 初始化循环状态
                cont='y';
                % 开始循环
                while strcmp(cont,'y')
                    % 询问用户输入新值
                    prompt={'输入新的阈值:','输入新的参考凹度:','输入新的段长:'};
                    definput={num2str(thresh_ratio),num2str(theta),num2str(segment_length)};
                    vals=inputdlg(prompt,'定义新值',[1 50],definput);
                    % 检查输出并转换为数值
                    if ~isempty(vals)
                        if ~isempty(vals{1})
                            thresh_ratio=str2num(vals{1});
                        end

                        if ~isempty(vals{2})
                            theta=str2num(vals{2});
                        end

                        if ~isempty(vals{3})
                            segment_length_temp=str2num(vals{3});
                            if segment_length_temp<DEM.cellsize*3;
                                if isdeployed
                                    warndlg('输入的段长小于 DEM 单元格大小的 3 倍，将恢复为先前的值')
                                else
                                    warning('输入的段长小于 DEM 单元格大小的 3 倍，将恢复为先前的值');
                                end
                            else
                                segment_length=segment_length_temp;
                            end
                        end
                    end

                    % 使用新值重新绘图
                    [fig_handle]=PlotFunc(DEMc,FD,ST,A,segment_length,theta,thresh_ratio);
                    % 再次询问用户
                    qa2=questdlg('接受这些参数吗？','拟合参数','更改值','接受值','接受值');

                    switch qa2
                        case '接受值'
                            close(fig_handle);
                            cont='n';
                        case '更改值'
                            cont='y';
                    end % 子开关结束
                end % 循环结束
        end % 主开关结束
    end 
    
    % 绘制网络中最长干流的代表性示例
    if plot_example
        if isempty(channeloi)
            ST=trunk(klargestconncomps(S,1));
        else
            chxy=streampoi(S,'channelheads','xy');
            [~,mix]=min(hypot(chxy(:,1)-channeloi(:,1),chxy(:,2)-channeloi(:,2)));
            chix=coord2ind(DEM,chxy(mix,1),chxy(mix,2));
            CHIX=GRIDobj(DEM);
            CHIX.Z(chix)=true;
            ST=modify(S,'downstreamto',CHIX);
        end

        [fig_handle]=PlotFunc(DEMc,FD,ST,A,segment_length,theta,thresh_ratio);
    end

    % 存储参数以便输出，以防参数被更改
    param=struct; 
    param.thresh_ratio=thresh_ratio;
    param.reference_concavity=theta;
    param.segment_length=segment_length;


    w1=waitbar(0,'正在计算ksn值...');
    % 逐个节点计算 ksn
    g=gradient(S,DEMc);
    a=getnal(S,(A.*A.cellsize^2));
    knal=g./a.^(-theta);
    
    waitbar(0.25/5,w1,'正在分割河流网络...');
    % 按河道头分割流网络
    [SC,LOCS]=STREAMobj2cell(S,'channelheads');

    % 分箱并查找 ksn 段和断点
    nalC=cell(size(SC));
    nalR=cell(size(SC));
    nalP=cell(size(SC));
    nalN=cell(size(SC));
    nalCo=cell(size(SC));
    brkX=cell(size(SC));
    brkY=cell(size(SC));
    brkD=cell(size(SC));
    brkZ=cell(size(SC));
    for ii=1:numel(SC)
        prog=(ii/numel(SC))*(3/5)+(1/5);
        waitbar(prog,w1,'正在计算自动 ksn 值...');
        % 分箱平均
        [dav,kav,idx]=BinAverage(SC{ii}.distance,knal(LOCS{ii}),segment_length);
        % 查找变化点
        T=mean(kav)*thresh_ratio;
        [bp,mk]=ischange(kav,'variance','Threshold',T);
        % 将分箱值转换回节点属性列表
        [nal_temp,brk]=UnBinNal(idx,bp,mk);
        % 计算残差
        [nalC{ii},nalR{ii},brkX{ii},brkY{ii},brkD{ii},brkZ{ii},~,~,nalP{ii},nalN{ii},nalCO{ii}]=ResidNal(SC{ii},DEMc,FD,A,theta,nal_temp,brk,cc);
    end

    % 累积值
    waitbar(4/5,w1,'正在累积值...');
    locl=vertcat(LOCS{:});
    nalcl=vertcat(nalC{:});
    nalrl=vertcat(nalR{:});
    nalpl=vertcat(nalP{:});
    nalnl=vertcat(nalN{:});
    if cc
        nalcol=vertcat(nalCO{:});
    end
    brkX=vertcat(brkX{:});
    brkY=vertcat(brkY{:});
    brkD=vertcat(brkD{:});
    brkZ=vertcat(brkZ{:});
    brkPnts=[brkX brkY brkD brkZ];
    brkPnts=unique(brkPnts,'rows');
    nal=accumarray(locl,nalcl,[],@mean);
    res=accumarray(locl,nalrl,[],@mean);
    pos=accumarray(locl,nalpl,[],@mean);
    neg=accumarray(locl,nalnl,[],@mean);
    if cc
        con=accumarray(locl,nalcol,[],@mean);
    end
    tta=ones(size(nal))*theta;

    waitbar(4.5/5,w1,'正在生成地图结构和 shapefile 输出...');
    % 将节点属性列表转换为地图结构
    if cc
        ms=STREAMobj2mapstruct(S,'seglength',segment_length,'attributes',...
            {'fit_ksn' nal @mean 'ksn_neg' neg @mean 'ksn_pos' pos @mean 'resid' res @mean 'theta' tta @mean 'seg_theta' con @mean...
             'rough_ksn' knal @mean 'up_area' a @mean 'gradient' g @mean 'cut_fill' DEMc-DEM @mean});
    else
        ms=STREAMobj2mapstruct(S,'seglength',segment_length,'attributes',...
            {'fit_ksn' nal @mean 'ksn_neg' neg @mean 'ksn_pos' pos @mean 'resid' res @mean 'theta' tta @mean...
             'rough_ksn' knal @mean 'up_area' a @mean 'gradient' g @mean 'cut_fill' DEMc-DEM @mean});
    end

    % 若需要，输出 shapefile
    if generate_shape
        ksn_name=[shape_name '_ksn.shp'];
        bnd_name=[shape_name '_bounds.shp'];

        shapewrite(ms,ksn_name);

        if ~isempty(brkPnts)
            bnds=struct;
            for jj=1:numel(brkPnts(:,1))
                bnds(jj,1).Geometry='Point';
                bnds(jj,1).X=double(brkPnts(jj,1));
                bnds(jj,1).Y=double(brkPnts(jj,2));
                bnds(jj,1).Dist=double(brkPnts(jj,3));
                bnds(jj,1).Elev=double(brkPnts(jj,4));
            end
            shapewrite(bnds,bnd_name);
        end
    end

    waitbar(1,w1,'完成');
    close(w1);
    plotc(S, nal); hold on; colorbar; 
    if ~isempty(brkPnts) && size(brkPnts,2)>=2
    scatter(brkPnts(:,1), brkPnts(:,2), 10, 'k', 'filled'); 
else
    warning('未找到有效裂点，跳过裂点绘制');
    end
end

function [f]=PlotFunc(DEMc,FD,S,A,segment_length,theta,thresh_ratio)
    gT=gradient(S,DEMc);
    aT=getnal(S,(A.*A.cellsize^2));
    knalT=gT./aT.^(-theta);

    [davT,kavT,idxT]=BinAverage(S.distance,knalT,segment_length);
    TT=mean(kavT)*thresh_ratio;
    [bpT,mkT]=ischange(kavT,'variance','Threshold',TT);
    [nalT,brkT]=UnBinNal(idxT,bpT,mkT);
    [rcnalT,resT,brkTX,brkTY,brkTD,~,zp,z,nalPT,nalNT,~]=ResidNal(S,DEMc,FD,A,theta,nalT,brkT,false);

    f=figure(1);
    clf 
    set(f,'Units','normalized','Position',[0.5 0.1 0.5 0.8],'renderer','painters');	

    sbplt1=subplot(3,1,1);
    hold on 
    xlim([0 max(S.distance)]);
    p1=scatter(davT,kavT,10,'k','filled');
    p2=plotdz(S,rcnalT,'color','r');
    p3=plotdz(S,rcnalT+nalPT,'color','r');
    p4=plotdz(S,rcnalT-nalNT,'color','r');
    set(p3,'LineStyle',':');
    set(p4,'LineStyle',':');
    for ii=1:numel(brkTD)
        plot([brkTD(ii) brkTD(ii)],ylim,':k','LineWidth',0.5);
    end
    xlabel('河流距离 (m)');
    ylabel('分段后的 k_{sn}');
    title(['阈值比率 = ' num2str(thresh_ratio) '; 参考凹度 = ' num2str(theta) '; 段长 = ' num2str(segment_length) ' m']);
    legend([p1 p2 p3],{'分段后的值','自动拟合值','不确定性'},'location','最佳');

    if ~verLessThan('matlab','9.5')
        disableDefaultInteractivity(sbplt1);
    end
    hold off

    sbplt2=subplot(3,1,2);
    hold on
    xlim([0 max(S.distance)]);
    p1=scatter(S.distance,zp,5,'r');
    p2=plotdz(S,z,'color','k');
    for ii=1:numel(brkTD)
        plot([brkTD(ii) brkTD(ii)],ylim,':k','LineWidth',0.5);
    end
    xlabel('河流距离 (m)');
    ylabel('高程 (m)');
    legend([p1 p2],{'预测高程','平滑后的高程'},'location','最佳');
    if ~verLessThan('matlab','9.5')
        disableDefaultInteractivity(sbplt2);
    end		
    hold off

    sbplt3=subplot(3,1,3);
    hold on
    xlim([0 max(S.distance)]);
    plot(S.distance,zeros(size(S.distance)),'-k');
    p1=scatter(S.distance,resT,5,'r','filled');
    for ii=1:numel(brkTD)
        plot([brkTD(ii) brkTD(ii)],ylim,':k','LineWidth',0.5);
    end
    xlabel('河流距离 (m)');
    ylabel('残差 (m)');
    title(['平均 (绝对值) 残差 = ' num2str(mean(abs(resT))) ' m']);
    legend(p1,{'预测值 - 平滑后的值'},'location','最佳');
    if ~verLessThan('matlab','9.5')
        disableDefaultInteractivity(sbplt2);
    end			
    hold off

    drawnow
end

function [Xavg,Yavg,idx]=BinAverage(X,Y,bin_size)

    ix=~isnan(X);
    X=X(ix); Y=Y(ix);

    minX=min(X);
    maxX=max(X);

    b=[minX:bin_size:maxX+bin_size];

    try
        [idx]=discretize(X,b);
    catch
        [~,idx]=histc(X,b);
    end

    Xavg=accumarray(idx(:),X,[],@mean);
    Yavg=accumarray(idx(:),Y,[],@mean);
end

function [nal,brk]=UnBinNal(idx,bp,mk)

    num_bins=numel(unique(idx));
    nal=zeros(size(idx));
    brk=zeros(size(idx));

    brk_num=1;

    for ii=1:num_bins
        bidx=idx==ii;
        % nal(bidx)=mk(ii);

        if bp(ii)
            brk(bidx)=brk_num;
            brk_num=brk_num+1;
            nal(bidx)=mk(ii);
        else
            brk(bidx)=0;
            nal(bidx)=mk(ii);
        end
    end
end

function [rcnal,res,brk_x,brk_y,brk_d,brk_z,pred_elev,elev,nalpos,nalneg,nalcon]=ResidNal(S,DEMc,FD,A,theta,nal,brk,cc_flag)
    
    if nnz(brk)>0
        brk_ix=zeros(size(S.x));
        brk_nal_ix(S.ixc)=(nal(S.ix)-nal(S.ixc))~=0;
        brk_d=S.distance(brk_nal_ix);
        brk_ix=S.IXgrid(logical(brk_nal_ix));
        [brk_x,brk_y]=ind2coord(DEMc,brk_ix);
        brk_z=DEMc.Z(brk_ix);

        SP=split(S,brk_ix);

        elev=getnal(S,DEMc);
        c=chitransform(SP,A,'a0',1,'mn',theta);
        [L,nc]=conncomps(SP);

        rcnal=zeros(size(nal));
        nalpos=zeros(size(nal));
        nalneg=zeros(size(nal));
        pred_elev=zeros(size(nal));
        nalcon=zeros(size(nal));

        for jj=1:nc
            idx=L==jj;
            out=ChiSpline(elev(idx),c(idx));
            rcnal(idx)=out.ks;
            nalpos(idx)=out.ks_pos;
            nalneg(idx)=out.ks_neg;
            pred_elev(idx)=c(idx).*rcnal(idx);
            pred_elev(idx)=pred_elev(idx)+min(elev(idx));

            if cc_flag
                W=GRIDobj(DEMc);
                W.Z(S.IXgrid(idx))=true;
                SPT=STREAMobj(FD,W);
                nalcon(idx)=BestFitTheta(S,A,DEMc);
            else 
                nalcon=[];
            end
        end

        res=pred_elev-elev;
    else
        brk_x=[]; brk_y=[]; brk_ix=[]; brk_d=[]; brk_z=[];	
        elev=getnal(S,DEMc);
        c=chitransform(S,A,'a0',1,'mn',theta);

        out=ChiSpline(elev,c);
        rcnal=ones(size(nal))*out.ks;
        nalpos=ones(size(nal))*out.ks_pos;
        nalneg=ones(size(nal))*out.ks_neg;

        pred_elev=(c.*rcnal)+min(elev);

        res=pred_elev-elev;

        if cc_flag
            nalcon=ones(size(nal));
            mnF=BestFitTheta(S,A,DEMc);
            nalcon=nalcon*mnF;
        else
            nalcon=[];
        end
    end
end

function [OUT]=ChiSpline(z,c)
    zabsF=z-min(z);
    chiF=c;

    % 样条插值对于小线段会产生很多警告
    warning off
    chiS=linspace(0,max(chiF),numel(chiF)).';
    try
        zS=spline(chiF,zabsF,chiS);
    catch
        % 若线段几乎是直线，三次样条插值会失败，此时跳过样条拟合以避免出错
        zS=zabsF;
        chiS=chiF;
    end

    OUT=struct;
    try
        ft=fittype('a*x');
        fobj=fit(chiS,zS,ft,'StartPoint',chiS\zS);
        BETA=coeffvalues(fobj);
        BETA_UNC=confint(fobj);
        OUT.ks   = BETA;
        OUT.ks_neg = (BETA)-min(BETA_UNC);
        OUT.ks_pos = max(BETA_UNC)-(BETA);
    catch
        BETA = chiS\(zS);
        OUT.ks   = BETA;
        OUT.ks_neg = 0;
        OUT.ks_pos = 0;
    end

    warning on

end

function [mnF]=BestFitTheta(S,A,DEM)

    z=getnal(S,DEM);
    mnF=fminsearch(@mnfit,0.5);

    function [sqres]=mnfit(mn)
        c=chitransform(S,A,'a0',1,'mn',mn);
        c=c./max(c);
        z=z-min(z);
        z=z./max(z);
        sqres=sum((c-z).^2);
    end

end