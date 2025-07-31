function PlotChi(DEM,S,chi,chi_type,varargin)
%
% 使用方法：
% PlotChi(DEM, S, chi, 'chimap');
% PlotChi(DEM, S, chi, 'chigrid');
% PlotChi(DEM, S, chi, 'chimap', '属性名', '值', ...);
%
% 功能描述：
% 此函数用于绘制标准化河道陡度图（Chi图），并以地形阴影图为背景展示高程信息。
%
% 必需输入参数：
% DEM - 数字高程模型，以GRIDobj格式存储，用于生成chi数据。
% S - 生成chi值时所用的STREAMobj类型的河流网络对象。
% chi - chi数据，支持ASCII文件或GRIDobj格式（例如KsnChiBatch的输出）。
% chi_type - 指定chi数据的类型：'chimap'（仅沿河道）或'chigrid'（连续栅格）。
%
% 可选输入参数：
% chi_lim [] - 一个1×n向量，用于设置颜色映射的范围[最小值 最大值]。若留空，默认范围为[0 数据集最大值]。
%
% 示例：
% PlotChi(DEM, S, chimap, 'chimap'); % 绘制GRIDobj格式的河道chi图
% PlotChi(DEM, S, 'Topo_chigrid.txt', 'chigrid', 'chi_lim', [0 10]);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 编写者：Adam M. Forte - 最后更新日期：2019年9月30日 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 解析输入参数
p = inputParser;
p.FunctionName = 'PlotChi';
addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
addRequired(p,'S',@(x) isa(x,'STREAMobj'));
addRequired(p,'chi',@(x) isa(x,'GRIDobj') || regexp(x,regexptranslate('wildcard','*.txt')));
addRequired(p,'chi_type',@(x) ischar(validatestring(x,{'chimap','chigrid'})));

addParameter(p,'chi_lim',[],@(x) isnumeric(x) && numel(x)==2);
addParameter(p,'override_resample',false,@(x) isscalar(x) && islogical(x)); % 图形用户界面的隐藏选项

parse(p,DEM,S,chi,chi_type,varargin{:});
DEM=p.Results.DEM;
S=p.Results.S;
chi=p.Results.chi;
chi_type=p.Results.chi_type;

chi_lim=p.Results.chi_lim;
os=p.Results.override_resample;

%%
% 若chi为字符串且符合文本文件格式
if ischar(chi) & logical(regexp(chi,regexptranslate('wildcard','*.txt')))
    chi=GRIDobj(chi);
    % 若DEM与chi未对齐且不覆盖重采样
    if ~validatealignment(DEM,chi) && ~os
        chi=resample(chi,DEM);
    % 若DEM与chi未对齐且覆盖重采样
    elseif ~validatealignment(DEM,chi) && os
        chi.refmat=DEM.refmat;
        chi.georef=DEM.georef;
    % 若chi为GRIDobj对象
    elseif isa(chi,'GRIDobj');
        % 若DEM与chi未对齐
        if ~validatealignment(DEM,chi)
            chi=resample(chi,DEM);
        end
    else
        % 若代码为部署模式
        if isdeployed
            errordlg('输入的 "chi" 未被识别为有效的ASCII文件或GRIDobj对象')
        end
        error('输入的 "chi" 未被识别为有效的ASCII文件或GRIDobj对象');
    end
end

% 根据chi_type的值进行不同处理
switch chi_type
    % 若为'chigrid'
    case 'chigrid'
        f1=figure(1);
        set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.8],'renderer','painters');
        hold on
        % 若chi_lim为空
        if isempty(chi_lim)
            imageschs(DEM,chi,'colormap','jet');
        else
            imageschs(DEM,chi,'colormap','jet','caxis',chi_lim);
        end
        % 若MATLAB版本不低于9.5
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(gca);
        end 
        hold off
    % 若为'chimap'
    case 'chimap'
        nal=getnal(S,chi);

        f1=figure(1);
        set(f1,'Visible','off');

        [RGB]=imageschs(DEM,DEM,'colormap','gray');
        [~,R]=GRIDobj2im(DEM);

        imshow(flipud(RGB),R);
        axis xy
        hold on
        colormap(jet);
        plotc(S,nal);
        % 若chi_lim为空
        if isempty(chi_lim)
            caxis([0 max(nal)]);
        else
            caxis([min(chi_lim) max(chi_lim)])
        end
        c1=colorbar;
        ylabel(c1,'\chi');

        % 若MATLAB版本不低于9.5
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(gca);
        end 
        hold off
        set(f1,'Visible','on','Units','normalized','Position',[0.05 0.1 0.8 0.8],'renderer','painters');
end
end