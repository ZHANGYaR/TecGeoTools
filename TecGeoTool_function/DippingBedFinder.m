function [BED]=DippingBedFinder(DEM,xy,hght_abv_base,thickness,strike,dip)
    %
    % 用法：
    %   [GRIDobj]=DippingBedFinder(DEM,xy,hght_abv_base,thickness,strike,dip); 
    %
    % 功能描述：
    %   本函数根据输入坐标点预测倾斜岩层在数字地形中的空间展布，生成岩层出露区域的二值栅格
    %
    % 必需输入参数：
    % 	DEM - 数字高程模型 GRIDobj 对象
    %   xy - 1x2向量，表示目标位置的x,y坐标（东坐标和北坐标）。若输入空数组[]，
    %        可通过交互方式在DEM图上选取位置
    % 	hght_abv_base - 目标露头在岩层基底以上的高度（即露头在岩层剖面中的垂直位置）
    %   thickness - 岩层总厚度（hght_abv_base必须小于总厚度）
    % 	strike - 岩层走向，按右手法则报告
    % 	dip - 岩层倾角
    %
    % 输出：
    % 	生成展示岩层预测位置的图形，并返回二值GRIDobj（1表示岩层应出露区域，0表示非出露区）
    %
    % 输入示例：
    %   创建DEM GRIDobj：
    %   [DEM]=GRIDobj('tif文件名.tif'); 或 [DEM]=GRIDobj('文本格式DEM.txt');
    %   使用示例：
    % 	[BED]=DippingBedFinder(DEM,[45325.23 1024567.2],10,50,270,40); % 指定坐标点
    %   [BED]=DippingBedFinder(DEM,[],10,50,270,40); % 交互式选取坐标点
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 函数作者：Adam M. Forte - 最近更新：2018/06/18        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isempty(xy)
        f1=figure(1);
        set(f1,'Units','normalized','Position',[0.1 0.1 0.8 0.8],'renderer','painters');
        imageschs(DEM,DEM);
        [x_coord,y_coord]=ginput(1); % 交互式获取坐标点
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(gca); % 禁用新版MATLAB默认交互功能
        end 
        close(f1)
    else
        x_coord=xy(1);
        y_coord=xy(2);
    end

    % 获取DEM空间信息
    gs=DEM.cellsize; % 栅格分辨率
    [xvec,yvec]=getcoordinates(DEM);
    dim=DEM.size; % DEM维度

    % 定位输入坐标对应的行列号及高程值
    x_dif=abs(xvec-x_coord);
    y_dif=abs(yvec-y_coord);

    [~,col]=min(x_dif);
    [~,row]=min(y_dif);

    sam_el=DEM.Z(row,col); % 输入点处的高程值

    % 生成网格系统
    xmax=dim(2)*gs;
    ymax=dim(1)*gs;

    x_spacing=(xmax)/(gs);
    y_spacing=(ymax)/(gs);

    x=linspace(0,xmax,x_spacing);
    y=linspace(0,ymax,y_spacing);

    [xi,yi]=meshgrid(x,y);

    % 计算岩层几何参数
    de=PlaneEq(strike,dip,xi,yi); % 调用平面方程计算

    % 基于输入点高程校正平面高程
    orig_plane_el=de(row,col);
    el_dif=sam_el-orig_plane_el;
    bed_el=de+el_dif; % 校正后的岩层高程面

    % 计算岩层厚度相关参数
    theta=90-dip; % 岩层法线与垂直方向的夹角

    if hght_abv_base>thickness
        error('露头基底高度不可超过岩层总厚度') % 参数校验错误提示
    end

    tu=thickness-hght_abv_base; % 岩层上覆厚度
    tb=thickness-tu; % 岩层下伏厚度
    atu=tu/sind(theta); % 沿倾向的上覆距离
    atb=tb/sind(theta); % 沿倾向的下伏距离

    % 确定岩层出露区间
    up_int=bed_el+atu; % 岩层上界面高程
    lb_int=bed_el-atb; % 岩层下界面高程

    z=DEM.Z;
    inter=z<=up_int & z>=lb_int; % 高程区间判断
    inter=double(inter); % 转换为双精度型

    BED=GRIDobj(xvec,yvec,inter); % 生成结果栅格
    BED.georef=DEM.georef; % 继承空间参考

    % 可视化结果
    f1=figure(1);
    set(f1,'Units','normalized','Position',[0.1 0.1 0.8 0.8],'renderer','painters');
    hold on
    imageschs(DEM,BED,'colormap','parula','colorbar',false); % 叠加显示预测结果
    scatter(x_coord,y_coord,20,'w','filled'); % 标记输入点位置
    title('岩层投影交汇区');
    if ~verLessThan('matlab','9.5')
        disableDefaultInteractivity(gca); % 保持图形稳定性
    end     
    hold off
end

function [de]=PlaneEq(strike,dip,xi,yi)
    % 法向量方程来自Charles Ammon：
    % http://eqseis.geosc.psu.edu/~cammon/HTML/UsingMATLAB/PDF/ML1%20FaultNormals.pdf
    %
    % 计算法向量分量
    a=sind(dip)*cosd(strike);
    b=-sind(dip)*sind(strike);
    c=-cosd(dip);

    % 生成以原点为中心的平面方程
    de=(a.*xi + b.*yi)./-c;

    % 坐标系翻转校正
    de=fliplr(de);
end