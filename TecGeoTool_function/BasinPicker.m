function [Outlets]=BasinPicker(DEM,FD,A,S,varargin)
 %
% 用法:
%   [Outlets]=BasinPicker(DEM,FD,A,S);
%   [Outlets]=BasinPicker(DEM,FD,A,S,'属性名',值,...);
%
% 描述:
%   此函数使用生成的河流数据，允许用户交互式选择流域（集水区）。该函数最初设计用于选择适合碎屑分析的流域
% （例如 Be-10 同位素分析）。显示一个包含高程着色地形图和局部起伏图的双面板图形，以供选择流域。图形显示后，
% 会等待用户按下回车键以开始流域选择过程。这一设计使用户可以缩放、平移等以找到感兴趣的河流。当按下回车键后，
% 十字光标会出现在高程图上以便选择一个出水点。一旦选择了出水点，新的图形会显示此流域及其河流，以确认是否为
% 所需的流域（同时会显示其流域面积）。用户可以接受此流域或拒绝（若为误点击）。如果接受，接下来会显示包含此流域
% 的 χ-z 和纵向剖面的新图形。然后会提示用户保存此选择或丢弃它。最后询问用户是否继续选择河流，若选择是（默认），
% 则重新开始选择过程。注意，所有已选择（并保存）的出水点将显示在主图中。用户在选择流域时，函数会保存一个名为
% 'Outlets.mat' 的文件，其中包含到目前为止所选择的出水点。如果退出函数并稍后重新启动，它会在当前工作目录中查找
% 此 Outlets 文件，以便从中断处继续。
%
% 必需输入:
%       DEM - 数字高程模型的 GRIDobj 对象
%       FD - 从提供的 DEM 创建的 FLOWobj 对象
%       A - 流量累积栅格 (GRIDobj)
%       S - 从 DEM 导出的 STREAMobj 对象
%
% 可选输入:
%       ref_concavity [0.50]- χ-z 图的参考凹度
%       rlf_radius [2500] - 用于计算局部起伏的半径（地图单位）或
%       rlf_grid [] - 如果已有计算好的局部起伏栅格，可以通过 'rlf_grid' 提供，必须是一个 GRIDobj，且与提供的 DEM 具有相同的坐标和维度。
%       extra_grid [] - 有时可能需要查看附加的栅格（例如地理参考的道路图、降水栅格等）与 DEM 和起伏图一起显示。
%           此栅格可以是不同的尺寸或像元大小（但必须与提供的 DEM 具有相同的投影和坐标系统），它将被重新采样以匹配提供的 DEM。
%       cmap ['landcolor'] - 用于显示地图的颜色图。输入可以是标准颜色图的名称或一个 nx3 的 RGB 值数组作为颜色图。
%       conditioned_DEM [] - 提供用于此函数的平滑 DEM 的选项（请勿将平滑DEM 用作主要必需的 DEM 输入），此 DEM 将用于提取高程值。有关生成水文条件化 DEM 的选项，请参阅 'ConditionDEM' 函数。如果未提供输入，则代码默认为使用 mincosthydrocon 函数。
%       interp_value [0.1] - mincosthydrocon 插值参数的值（介于 0 和 1 之间）（如果用户提供了平滑 DEM，则不会使用）。
%       plot_type ['vector'] - 接受 'vector' 或 'grid'，默认值为 'vector'。控制所有河流是以独立线段 ('vector') 绘制，还是以栅格形式绘制并下采样 ('grid')。对于大型数据集，'grid' 选项更快，但可能导致站点选择不准确；'vector' 选项更易于观察，但加载和交互速度较慢。
%       threshold_area [1e6] - 如果 'plot_type' 设置为 'grid'，则用于重新绘制下采样的河流网络。
%       refine_positions [] - 期望为 m x 2 的矩阵，包含接近河流网络的 x y 位置，用户希望手动将其对齐到适当的河流网络。例如，一组 GPS 记录的河流样本位置未完全位于由流向计算确定的河流网络上。如果为 'refine_positions' 提供了输入，代码将依次处理这些点，每次显示其位置，并附带一个放大的嵌套窗口（大小由 'window_size' 控制），用户可以在其中精确定位河流出口位置。
%       window_size [1] - 如果提供了 'refine_positions'，则嵌套窗口的大小（单位：公里）。
%
% 输出:
%       Outlets - n x 3 矩阵，包含样本位置的流域编号、x 坐标和 y 坐标列（可作为 'ProcessRiverBasins' 的 'river_mouths' 参数的有效输入）。
%
% 示例:
%       [Outs]=BasinPicker(DEM,FD,A,S);
%       [Outs]=BasinPicker(DEM,FD,A,S,'rlf_radius',5000);
%       [Outs]=BasinPicker(DEM,FD,A,S,'rlf_grid',RLF);
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数作者：Adam M. Forte - 更新时间：2018 年 6 月 18 日 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % 解析输入参数
    p = inputParser;
    p.FunctionName = 'BasinPicker';
    addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
    addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
    addRequired(p,'A',@(x) isa(x,'GRIDobj'));
    addRequired(p,'S',@(x) isa(x,'STREAMobj'));

    addParameter(p,'ref_concavity',0.5,@(x) isscalar(x) && isnumeric(x));
    addParameter(p,'rlf_radius',2500,@(x) isscalar(x) && isnumeric(x));
    addParameter(p,'plot_type','vector',@(x) ischar(validatestring(x,{'vector','grid'})));
    addParameter(p,'cmap','landcolor',@(x) ischar(x) || isnumeric(x) & size(x,2)==3);
    addParameter(p,'threshold_area',1e6,@(x) isnumeric(x));   
    addParameter(p,'rlf_grid',[],@(x) isa(x,'GRIDobj') || isempty(x));
    addParameter(p,'extra_grid',[],@(x) isa(x,'GRIDobj') || isempty(x));
    addParameter(p,'conditioned_DEM',[],@(x) isa(x,'GRIDobj') || isempty(x));
    addParameter(p,'interp_value',0.1,@(x) isnumeric(x) && x>=0 && x<=1);
    addParameter(p,'refine_positions',[],@(x) isnumeric(x) && size(x,2)==2 || isempty(x));
    addParameter(p,'window_size',1,@(x) isnumeric(x) && isscalar(x));
    addParameter(p,'out_dir',[],@(x) ischar(x)); % 用于GUI版本的隐藏选项

    parse(p,DEM,FD,A,S,varargin{:});
    DEM=p.Results.DEM;
    FD=p.Results.FD;
    S=p.Results.S;
    A=p.Results.A;

    theta_ref=p.Results.ref_concavity;
    rlf_radius=p.Results.rlf_radius;
    plot_type=p.Results.plot_type;
    threshold_area=p.Results.threshold_area;
    cmap=p.Results.cmap;
    RLF=p.Results.rlf_grid;
    EG=p.Results.extra_grid;
    DEMf=p.Results.conditioned_DEM;
    iv=p.Results.interp_value;
    rp=p.Results.refine_positions;
    ws=p.Results.window_size;
    out_dir=p.Results.out_dir;

    % 检查当前目录是否存在之前运行的出口文件
    if isempty(out_dir)
        current=pwd;
    else
        current=out_dir;
    end

    out_file=fullfile(current,'Outlets.mat');

    if exist(out_file,'file')==2
        an=questdlg('您想从当前目录下的 "Outlets.mat" 文件中加载之前选择的出口吗？',...
            '之前的出口','否','是','是');
        switch an
        case '是'           
            load(out_file,'Outlets');
            ii=max(Outlets(:,3))+1;
        case '否'
            ii=1; 
            Outlets=zeros(0,3); 
        end
    else 
       ii=1; 
       Outlets=zeros(0,3);
    end
             
    % 平滑DEM处理
    if isempty(DEMf)
        zc=mincosthydrocon(S,DEM,'interp',iv);
        DEMf=GRIDobj(DEM);
        DEMf.Z(DEMf.Z==0)=NaN;
        DEMf.Z(S.IXgrid)=zc;
    end

    % 计算地形起伏
    if isempty(RLF)
        disp('正在计算局部地形起伏');
        RLF=localtopography(DEM,rlf_radius);
    end

    % 设置操作模式标志
    if ~isempty(rp)
        ws=ws*1000;
        num_points=size(rp,1);
        op_mode='refine';
    else
        op_mode='normal';
    end

    % 计算坡度和简单KSN
    G=gradient8(DEMf);
    KSN=G./(A.*(A.cellsize^2)).^(-theta_ref);

    switch plot_type
    case 'vector'

        % 绘制主图
        if isempty(EG)
            f1=figure(1);
            clf
            set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8]);        

            ax(2)=subplot(2,1,2);
            hold on
            imageschs(DEM,RLF,'colormap',cmap);
            plot(S,'-k','LineWidth',1.5);
            scatter(Outlets(:,1),Outlets(:,2),20,'r','filled');
            title('局部地形起伏')

            if ~verLessThan('matlab','9.5')
                disableDefaultInteractivity(ax(2));
            end             
            hold off
                
            ax(1)=subplot(2,1,1);
            hold on
            imageschs(DEM,DEM,'colormap',cmap);
            plot(S,'-k','LineWidth',1.5);
            scatter(Outlets(:,1),Outlets(:,2),20,'r','filled');
            title('高程')
            hold off
            if ~verLessThan('matlab','9.5')
                disableDefaultInteractivity(ax(1));
            end             
            linkaxes(ax,'xy');
        else

            if ~validatealignment(EG,DEM)
                disp('正在重采样附加栅格');
                EG=resample(EG,DEM,'bicubic');
            end

            % 绘制主图
            f1=figure(1);
            clf
            set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8]);
            
            ax(3)=subplot(3,1,3);
            hold on
            imageschs(DEM,EG,'colormap',cmap)
            plot(S,'-k','LineWidth',1.5);
            scatter(Outlets(:,1),Outlets(:,2),20,'r','filled');
            title('用户提供的附加栅格')
            if ~verLessThan('matlab','9.5')
                disableDefaultInteractivity(ax(3));
            end 
            hold off

            ax(2)=subplot(3,1,2);
            hold on
            imageschs(DEM,RLF,'colormap',cmap);
            plot(S,'-k','LineWidth',1.5);
            scatter(Outlets(:,1),Outlets(:,2),20,'r','filled');
            title('局部地形起伏')
            if ~verLessThan('matlab','9.5')
                disableDefaultInteractivity(ax(2));
            end 
            hold off
                
            ax(1)=subplot(3,1,1);
            hold on
            imageschs(DEM,DEM,'colormap',cmap);
            plot(S,'-k','LineWidth',1.5);
            scatter(Outlets(:,1),Outlets(:,2),20,'r','filled');
            title('高程')
            if ~verLessThan('matlab','9.5')
                disableDefaultInteractivity(ax(1));
            end 
            hold off
            
            linkaxes(ax,'xy');
        end

    case 'grid'
        disp('正在下采样数据集以用于显示')  
        % 重新计算流向  
        DEMr=resample(DEM,DEM.cellsize*4,'bicubic');
        FDr=FLOWobj(DEMr,'preprocess','carve');
        % 真实出口点
        out_T_xy=streampoi(S,'outlets','xy');
        % 下采样后的总河网
        Sr_temp=STREAMobj(FDr,'minarea',threshold_area,'unit','mapunits');
        out_D_xy=streampoi(Sr_temp,'outlets','xy');
        out_D_ix=streampoi(Sr_temp,'outlets','ix');
        % 检查出口列表是否不同
        dists=pdist2(out_T_xy,out_D_xy);
        [~,s_out_ix]=min(dists,[],2);
        out_D_ix=out_D_ix(s_out_ix);
        % 重建下采样河网
        Sr=STREAMobj(FDr,'minarea',threshold_area,'unit','mapunits','outlets',out_D_ix);
        % 转换为栅格
        SG=STREAMobj2GRIDobj(Sr);
        % 重新计算局部起伏
        RLFr=resample(RLF,DEMr,'bicubic');

        if isempty(EG)
            % 绘制主图
            f1=figure(1);
            clf
            set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8]);
            
            ax(2)=subplot(2,1,2);
            hold on
            imageschs(DEMr,RLFr,'colormap',cmap);
            plot(Sr,'-k','LineWidth',1.5);
            scatter(Outlets(:,1),Outlets(:,2),20,'r','filled');
            title('局部地形起伏')
            if ~verLessThan('matlab','9.5')
                disableDefaultInteractivity(ax(2));
            end 
            hold off
                
            ax(1)=subplot(2,1,1);
            hold on
            imageschs(DEMr,DEMr,'colormap',cmap);
            plot(Sr,'-k','LineWidth',1.5);
            scatter(Outlets(:,1),Outlets(:,2),20,'r','filled');
            title('高程')
            if ~verLessThan('matlab','9.5')
                disableDefaultInteractivity(ax(1));
            end 
            hold off
            
            linkaxes(ax,'xy');
        else

            EGr=resample(EG,DEMr,'bicubic');

            % 绘制主图
            f1=figure(1);
            clf
            set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8]);
            
            ax(3)=subplot(3,1,3);
            hold on
            imageschs(DEMr,EGr,'colormap',cmap)
            plot(Sr,'-k','LineWidth',1.5);
            scatter(Outlets(:,1),Outlets(:,2),20,'r','filled');
            title('用户提供的附加栅格')
            if ~verLessThan('matlab','9.5')
                disableDefaultInteractivity(ax(3));
            end 
            hold off

            ax(2)=subplot(3,1,2);
            hold on
            imageschs(DEMr,RLFr,'colormap',cmap);
            plot(Sr,'-k','LineWidth',1.5);
            scatter(Outlets(:,1),Outlets(:,2),20,'r','filled');
            title('局部地形起伏')
            if ~verLessThan('matlab','9.5')
                disableDefaultInteractivity(ax(2));
            end 
            hold off
                
            ax(1)=subplot(3,1,1);
            hold on
            imageschs(DEMr,DEMr,'colormap',cmap);
            plot(Sr,'-k','LineWidth',1.5);
            scatter(Outlets(:,1),Outlets(:,2),20,'r','filled');
            title('高程')
            if ~verLessThan('matlab','9.5')
                disableDefaultInteractivity(ax(1));
            end 
            hold off
            
            linkaxes(ax,'xy');
        end

    end
    
    switch op_mode
    case 'normal'

        % 开始样本选择流程
        str1='N'; % 流域选择标志
        str2='Y'; % 继续选择标志

        while strcmpi(str2,'Y')
            while strcmpi(str1,'N')	
                
                if isempty(EG)
                    subplot(2,1,1);
                    hold on
                    title('缩放并平移到感兴趣的区域，准备选取时按 "Enter" 键')
                    hold off
                    pause()

                    subplot(2,1,1);
                    hold on
                    title('在高程图上选择样本点')
                    hold off
                else
                    subplot(3,1,1);
                    hold on
                    title('缩放并平移到感兴趣的区域，准备选取时按 "Enter" 键')
                    hold off
                    pause()

                    subplot(3,1,1);
                    hold on
                    title('在高程图上选择样本点')
                    hold off
                end

                [x,y]=ginput(1);

                % 构建逻辑栅格
                [xn,yn]=snap2stream(S,x,y);
                ix=coord2ind(DEM,xn,yn);
                IX=GRIDobj(DEM,'logical');
                IX.Z(ix)=true;

                % 裁剪出河网
                Sn=modify(S,'upstreamto',IX);
                
                % 裁剪栅格数据
                I=dependencemap(FD,xn,yn);
                DEMc=crop(DEM,I,nan);
                RLFc=crop(RLF,I,nan);
                if ~isempty(EG)
                    EGc=crop(EG,I,nan);
                end
                            
                % 计算流域面积
                dep_map=GRIDobj2mat(I);
                num_pix=sum(sum(dep_map));
                drainage_area=(num_pix*DEMc.cellsize*DEMc.cellsize)/(1e6);
            
                f2=figure(2);
                clf
                set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.45])
                hold on
                title(['流域面积 = ' num2str(round(drainage_area)) '平方公里']);
                colormap(cmap);
                imagesc(DEMc)
                plot(Sn,'-r','LineWidth',2);
                scatter(xn,yn,20,'w','filled');
                if ~verLessThan('matlab','9.5')
                    disableDefaultInteractivity(gca);
                end 
                hold off

                if isempty(EG)
                    figure(1);

                    subplot(2,1,2);
                    hold on
                    pl(1)=plot(Sn,'-r','LineWidth',2);
                    sc(1)=scatter(xn,yn,20,'w','filled');
                    hold off
                              
                    subplot(2,1,1);
                    hold on
                    pl(2)=plot(Sn,'-r','LineWidth',2);        
                    sc(2)=scatter(xn,yn,20,'w','filled');      
                    hold off
                else
                    figure(1);

                    subplot(3,1,3);
                    hold on
                    pl(1)=plot(Sn,'-r','LineWidth',2);
                    sc(1)=scatter(xn,yn,20,'w','filled');
                    hold off

                    subplot(3,1,2);
                    hold on
                    pl(2)=plot(Sn,'-r','LineWidth',2);
                    sc(2)=scatter(xn,yn,20,'w','filled');
                    hold off
                              
                    subplot(3,1,1);
                    hold on
                    pl(3)=plot(Sn,'-r','LineWidth',2);
                    sc(3)=scatter(xn,yn,20,'w','filled');
                    hold off
                end

                qa=questdlg('这是您想要的流域吗？','流域选择','否','是','是');

                switch qa
                case '是'
                    str1 = 'Y';
                case '否'
                    str1 = 'N';
                end

                delete(pl);
                delete(sc);
                close figure 2
            end
            
            Sn=klargestconncomps(Sn,1);
            C=chiplot(Sn,DEMf,A,'a0',1,'mn',theta_ref,'plot',false);

            ksn=getnal(Sn,KSN);
            mksn=mean(ksn,'omitnan');
            mrlf=mean(RLFc.Z(:),'omitnan');
            if ~isempty(EG)
                meg=mean(EGc.Z(:),'omitnan');
            end

            f2=figure(2);
            clf
            set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');
            sbplt1=subplot(2,1,1);
            hold on
            plot(C.chi,C.elev);
            xlabel('\chi')
            ylabel('高程 (m)')
            if isempty(EG)
                title(['\chi-Z 图：平均 k_{sn} = ' num2str(round(mksn)) ' : 平均起伏 = ' num2str(round(mrlf)) ' : 流域面积 = ' num2str(round(drainage_area)) '平方公里'])
            else
                title(['\chi-Z 图：平均 k_{sn} = ' num2str(round(mksn)) ' : 平均起伏 = ' num2str(round(mrlf)) ' : 附加栅格均值 = ' num2str(round(meg)) ' : 流域面积 = ' num2str(round(drainage_area)) '平方公里'])
            end

            if ~verLessThan('matlab','9.5')
                disableDefaultInteractivity(sbplt1);
            end 
            hold off

            sbplt2=subplot(2,1,2);
            hold on
            plotdz(Sn,DEMf,'dunit','km','Color','k');
            xlabel('溯源距离(km)')
            ylabel('高程 (m)')
            title('纵剖面图')

            if ~verLessThan('matlab','9.5')
                disableDefaultInteractivity(sbplt2);
            end 
            hold off

            qa2=questdlg('确定选择该流域吗？','流域选择','否','是','是');       

            switch qa2
            case '是'
                Outlets(ii,1)=xn;
                Outlets(ii,2)=yn;
                Outlets(ii,3)=ii;
                ii=ii+1;
                
                % 在主图上标记选中点

                if isempty(EG)
                    figure(1);
                    subplot(2,1,2);
                    hold on
                    scatter(xn,yn,20,'r','filled');
                    hold off
                              
                    subplot(2,1,1);
                    hold on
                    scatter(xn,yn,20,'r','filled');
                    hold off
                else
                    figure(1);

                    subplot(3,1,3);
                    hold on
                    scatter(xn,yn,20,'r','filled');
                    hold off

                    subplot(3,1,2);
                    hold on
                    scatter(xn,yn,20,'r','filled');
                    hold off
                              
                    subplot(3,1,1);
                    hold on
                    scatter(xn,yn,20,'r','filled');
                    hold off
                end
                
                save(out_file,'Outlets','-v7.3');
            end
            
            qa3=questdlg('是否继续选择流域？','流域选择','否','是','是'); 
            switch qa3
            case '是'
                str2='Y';
                str1='N';
            case '否'
                str2='N';
            end  
         
            %close figure 2;
        end

    case 'refine'

        for jj=1:num_points
            pnt=rp(jj,:);

            box_x=[pnt(1)-ws/2 pnt(1)+ws/2 pnt(1)+ws/2 pnt(1)-ws/2 pnt(1)-ws/2];
            box_y=[pnt(2)-ws/2 pnt(2)-ws/2 pnt(2)+ws/2 pnt(2)+ws/2 pnt(2)-ws/2];


            % 开始样本选择流程
            str1='N'; % 流域选择标志
            str2='Y'; % 继续选择标志

            while strcmpi(str2,'Y')
                while strcmpi(str1,'N')    
                    
                    if isempty(EG)
                        figure(1);
                        subplot(2,1,1);
                        hold on
                        bb(1)=plot(box_x,box_y,'-r','LineWidth',2);
                        hold off

                        subplot(2,1,2);
                        hold on
                        bb(2)=plot(box_x,box_y,'-r','LineWidth',2);
                        hold off

                    else
                        figure(1);
                        subplot(3,1,1);
                        hold on
                        bb(1)=plot(box_x,box_y,'-r','LineWidth',2);
                        hold off

                        subplot(3,1,2);
                        hold on
                        bb(2)=plot(box_x,box_y,'-r','LineWidth',2);
                        hold off

                        subplot(3,1,3);
                        hold on
                        bb(3)=plot(box_x,box_y,'-r','LineWidth',2);
                        hold off
                    end

                    f2=figure(2);
                    clf
                    set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.45])
                    hold on
                    title('请在高程图上选择样本点');
                    xlim([pnt(1)-ws/2 pnt(1)+ws/2]);
                    ylim([pnt(2)-ws/2 pnt(2)+ws/2]);
                    colormap(cmap);
                    switch plot_type
                    case 'vector'
                        imageschs(DEM,RLF,'colormap',cmap);
                        plot(S,'-k','LineWidth',1.5);
                    case 'grid'
                        imageschs(DEMr,EGr,'colormap',cmap)
                        plot(Sr,'-k','LineWidth',1.5);
                    end
                    scatter(pnt(1),pnt(2),20,'r','filled');
                    xlim([pnt(1)-ws/2 pnt(1)+ws/2]);
                    ylim([pnt(2)-ws/2 pnt(2)+ws/2]);

                    if ~verLessThan('matlab','9.5')
                        disableDefaultInteractivity(gca);
                    end 
                    hold off

                    figure(2);
                    [x,y]=ginput(1);
                    close figure 2

                    % 构建逻辑栅格
                    [xn,yn]=snap2stream(S,x,y);
                    ix=coord2ind(DEM,xn,yn);
                    IX=GRIDobj(DEM,'logical');
                    IX.Z(ix)=true;

                    % 裁剪出河网
                    Sn=modify(S,'upstreamto',IX);
                    
                    % 裁剪栅格数据
                    I=dependencemap(FD,xn,yn);
                    DEMc=crop(DEM,I,nan);
                    RLFc=crop(RLF,I,nan);
                    if ~isempty(EG)
                        EGc=crop(EG,I,nan);
                    end
                                
                    % 计算流域面积
                    dep_map=GRIDobj2mat(I);
                    num_pix=sum(sum(dep_map));
                    drainage_area=(num_pix*DEMc.cellsize*DEMc.cellsize)/(1e6);
                
                    f2=figure(2);
                    clf
                    set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.45])
                    hold on
                    title(['流域面积 = ' num2str(round(drainage_area)) '平方公里']);
                    colormap(cmap);
                    imagesc(DEMc)
                    plot(Sn,'-r','LineWidth',2);
                    scatter(xn,yn,20,'w','filled');

                    if ~verLessThan('matlab','9.5')
                        disableDefaultInteractivity(gca);
                    end 
                    hold off

                    if isempty(EG)
                        figure(1);

                        subplot(2,1,2);
                        hold on
                        pl(1)=plot(Sn,'-r','LineWidth',2);
                        sc(1)=scatter(xn,yn,20,'w','filled');
                        hold off
                                  
                        subplot(2,1,1);
                        hold on
                        pl(2)=plot(Sn,'-r','LineWidth',2);        
                        sc(2)=scatter(xn,yn,20,'w','filled');      
                        hold off
                    else
                        figure(1);

                        subplot(3,1,3);
                        hold on
                        pl(1)=plot(Sn,'-r','LineWidth',2);
                        sc(1)=scatter(xn,yn,20,'w','filled');
                        hold off

                        subplot(3,1,2);
                        hold on
                        pl(2)=plot(Sn,'-r','LineWidth',2);
                        sc(2)=scatter(xn,yn,20,'w','filled');
                        hold off
                                  
                        subplot(3,1,1);
                        hold on
                        pl(3)=plot(Sn,'-r','LineWidth',2);
                        sc(3)=scatter(xn,yn,20,'w','filled');
                        hold off
                    end

                    qa=questdlg('这是您想要的流域吗？','流域选择','否','是','是');

                    switch qa
                    case '是'
                        str1 = 'Y';
                    case '否'
                        str1 = 'N';
                    end

                    delete(pl);
                    delete(sc);
                    delete(bb);
                    close figure 2
                end
                
                Sn=klargestconncomps(Sn,1);
                C=chiplot(Sn,DEMf,A,'a0',1,'mn',theta_ref,'plot',false);

                ksn=getnal(Sn,KSN);
                mksn=mean(ksn,'omitnan');
                mrlf=mean(RLFc.Z(:),'omitnan');
                if ~isempty(EG)
                    meg=mean(EGc.Z(:),'omitnan');
                end

                f2=figure(2);
                clf
                set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');
                sbplt1=subplot(2,1,1);
                hold on
                plot(C.chi,C.elev);
                xlabel('\chi')
                ylabel('高程 (m)')
                if isempty(EG)
                    title(['\chi-Z 图：平均 k_{sn} = ' num2str(round(mksn)) ' : 平均起伏 = ' num2str(round(mrlf)) ' : 流域面积 = ' num2str(round(drainage_area)) '平方公里'])
                else
                    title(['\chi-Z 图：平均 k_{sn} = ' num2str(round(mksn)) ' : 平均起伏 = ' num2str(round(mrlf)) ' : 附加栅格均值 = ' num2str(round(meg)) ' : 流域面积 = ' num2str(round(drainage_area)) '平方公里'])
                end

                if ~verLessThan('matlab','9.5')
                    disableDefaultInteractivity(sbplt1);
                end 
                hold off

                sbplt2=subplot(2,1,2);
                hold on
                plotdz(Sn,DEMf,'dunit','km','Color','k');
                xlabel('溯源距离(km)')
                ylabel('高程 (m)')
                title('纵剖面图')
                if ~verLessThan('matlab','9.5')
                    disableDefaultInteractivity(sbplt2);
                end 
                hold off

                qa2=questdlg('确定选择该流域吗？','流域选择','否','是','是');       

                switch qa2
                case '是'
                    Outlets(ii,1)=xn;
                    Outlets(ii,2)=yn;
                    Outlets(ii,3)=ii;
                    ii=ii+1;
                    
                    % 在主图上标记选中点

                    if isempty(EG)
                        figure(1);
                        subplot(2,1,2);
                        hold on
                        scatter(xn,yn,20,'r','filled');
                        hold off
                                  
                        subplot(2,1,1);
                        hold on
                        scatter(xn,yn,20,'r','filled');
                        hold off
                    else
                        figure(1);

                        subplot(3,1,3);
                        hold on
                        scatter(xn,yn,20,'r','filled');
                        hold off

                        subplot(3,1,2);
                        hold on
                        scatter(xn,yn,20,'r','filled');
                        hold off
                                  
                        subplot(3,1,1);
                        hold on
                        scatter(xn,yn,20,'r','filled');
                        hold off
                    end
                    
                    save(out_file,'Outlets','-v7.3');
                end
                
                % 退出循环继续处理下一个点
                str2='N';
             
                close figure 2;
            % While循环结束
            end
        % For循环结束
        end
    % Switch结束
    end
    %close figure 1;
% 函数结束
end   