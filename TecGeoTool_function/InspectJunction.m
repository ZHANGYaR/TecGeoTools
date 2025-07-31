function InspectJunction(S,IX,num,junctions,varargin)
    % 用法：
    %   InspectJunction(S, IX, num, junctions);
    %   InspectJunction(S, IX, num, junctions, 'fit_distance', val);
    %
    % 描述：
    %   该函数用于可视化上下游连接的拟合，以评估交汇角度的值。它将显示交汇点的放大视图（使用极坐标显示）
    %   以及交汇点在流网地图中的空间位置。
    %
    % 必需输入参数：
    %   S - STREAMobj对象
    %   IX - JunctionAngle函数输出的IX元胞数组
    %   num - JunctionAngle表中交汇点编号。若为空数组，将提示用户在流网地图上交互选择交汇点
    %   junctions - 包含交汇点信息的表格，必须包含预测角度列
    %
    % 可选输入参数：
    %   fit_distance - 上游拟合距离（地图单位），需小于JunctionAngle使用的最大距离
    %   num_nodes - 用于拟合的节点数，需小于JunctionAngle提取的节点数
    %                 注意：fit_distance和num_nodes不能同时指定
    %
    % 示例：
    %   InspectJunction(S, IX, 50, junctions);
    %   InspectJunction(S, IX, 50, junctions, 'fit_distance', 500);
    %   InspectJunction(S, IX, 50, junctions, 'num_nodes', 10);
    %
    % 相关函数：
    %   JunctionAngle, JunctionLinks
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 作者：Adam M. Forte - 更新于2019/06/26
    % 修改者：张亚荣       更新于2024/12/26
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 解析输入参数
    p = inputParser;
    p.FunctionName = 'InspectJunction';
    addRequired(p,'S',@(x) isa(x,'STREAMobj'));
    addRequired(p,'IX',@(x) iscell(x));
    addRequired(p,'num',@(x) isnumeric(x) && isscalar(x) || isempty(x));
    addRequired(p,'junctions',@(x) istable(x)); 

    addParameter(p,'fit_distance',[],@(x) isscalar(x) && isnumeric(x) || isempty(x));
    addParameter(p,'num_nodes',[],@(x) isscalar(x) && isnumeric(x) || isempty(x));

    parse(p,S,IX,num,junctions,varargin{:});
    S = p.Results.S;
    IX = p.Results.IX;
    num = p.Results.num;
    junctions = p.Results.junctions;

    fit_distance = p.Results.fit_distance;
    num_nodes = p.Results.num_nodes;

    % 参数冲突检查
    if ~isempty(fit_distance) && ~isempty(num_nodes)
        if isdeployed
            errordlg('请仅指定fit_distance或num_nodes中的一个参数')
        end
        error('参数冲突：不能同时指定fit_distance和num_nodes');
    end

    % 初始化继续选择标志
    continueSelect = true;

    % 创建交互式地图窗口
    f1 = figure(1);
    hold on
    plot(S,'-k');
    axis equal
    title('请点击目标交汇点附近');
    if ~verLessThan('matlab','9.5')
        disableDefaultInteractivity(gca);
    end 
    hold off

    while continueSelect
        if isempty(num)
            % 获取所有交汇点坐标
            cix = cell2mat(IX(:,1));
            axc = S.x(cix);
            ayc = S.y(cix);

            % 交互式点选择
            figure(f1);
            [xpick, ypick] = ginput(1);
            [~, num] = min(hypot(axc - xpick, ayc - ypick));
        end

        c = IX{num,1};
        ds = IX{num,2};
        us1 = IX{num,3};
        us2 = IX{num,4};

        % 处理拟合距离参数
        if ~isempty(fit_distance)
            n = floor(fit_distance/hypot(S.cellsize,S.cellsize));
            if n < 1
                if isdeployed
                    warndlg('输入的拟合距离过短，将使用最小节点距离')
                end
                warning('输入的拟合距离过短，将使用最小节点距离');
                n = 1;
            end

            if n < numel(ds)
                ds = ds(1:n);
            end

            if n < numel(us1)
                us1 = us1(1:n);
            end

            if n < numel(us2)
                us2 = us2(1:n);
            end
        elseif ~isempty(num_nodes)
            n = floor(num_nodes);
            if n < 1
                n = 1;
            end

            if n < numel(ds)
                ds = ds(1:n);
            end

            if n < numel(us1)
                us1 = us1(1:n);
            end

            if n < numel(us2)
                us2 = us2(1:n);
            end        
        end

        % 坐标转换处理
        if isempty(S.georef)
            xl = S.x; yl = S.y;
            projcoord = false;
        else
            try
                [yl, xl] = projinv(S.georef.mstruct,S.x,S.y);
                projcoord = true;
            catch
                xl = S.x; yl = S.y;
                projcoord = false;
            end
        end

        % 交汇点坐标
        xc = xl(c); yc = yl(c);
        % 上游连接点坐标
        xu1 = xl(us1); yu1 = yl(us1);
        xu2 = xl(us2); yu2 = yl(us2);
        % 下游连接点坐标
        xd = xl(ds); yd = yl(ds);

        % 坐标归一化处理
        xu1 = xu1 - xc; xu2 = xu2 - xc;
        yu1 = yu1 - yc; yu2 = yu2 - yc;
        xd = xd - xc; yd = yd - yc;

        % 角度计算与坐标旋转
        theta1ap = median(atan2(yu1,xu1)); 
        theta2ap = median(atan2(yu2,xu2));
        theta3ap = median(atan2(yd,xd));

        % 执行坐标旋转
        [xu1r, yu1r] = RotCoord(xu1, yu1, -theta1ap, 0, 0); 
        [xu2r, yu2r] = RotCoord(xu2, yu2, -theta2ap, 0, 0);    
        [xdr, ydr] = RotCoord(xd, yd, -theta3ap, 0, 0);

        % 线性拟合处理
        a1 = xu1r\yu1r;
        a2 = xu2r\yu2r; 
        a3 = xdr\ydr;

        % 计算投影坐标
        yp1r = xu1r*a1;
        yp2r = xu2r*a2;
        yp3r = xdr*a3;

        % 反向旋转坐标
        [xp1, yp1] = RotCoord(xu1r, yp1r, theta1ap, 0, 0);
        [xp2, yp2] = RotCoord(xu2r, yp2r, theta2ap, 0, 0);
        [xp3, yp3] = RotCoord(xdr, yp3r, theta3ap, 0, 0);    

        % 转换为极坐标系
        [theta1, rho1] = cart2pol(xp1, yp1);
        [theta2, rho2] = cart2pol(xp2, yp2);
        [theta3, rho3] = cart2pol(xp3, yp3);

        [theta1o, rho1o] = cart2pol(xu1, yu1);
        [theta2o, rho2o] = cart2pol(xu2, yu2);
        [theta3o, rho3o] = cart2pol(xd, yd);

        % 寻找最大半径对应角度
        [mrho1, ix1] = max(rho1);
        [mrho2, ix2] = max(rho2);
        [mrho3, ix3] = max(rho3);
        theta1 = theta1(ix1);
        theta2 = theta2(ix2);
        theta3 = theta3(ix3);

        % 在地图上绘制选择结果
        figure(f1);
        hold on
        scatter(S.x(us1), S.y(us1), 10, 'r');
        scatter(S.x(us2), S.y(us2), 10, 'b');
        scatter(S.x(ds), S.y(ds), 10, 'k');
        scatter(S.x(c), S.y(c), 20, 'g', 'filled');
        hold off

        % 创建极坐标可视化窗口
        f2 = figure;
        set(f2, 'unit', 'normalized', 'position', [0.1 0.1 0.4 0.4]);
        clf 
        polarplot([theta1 theta1], [0 mrho1], '-r', 'LineWidth', 2);
        hold on 
        polarplot([theta2 theta2], [0 mrho2], '-b', 'LineWidth', 2);
        polarplot([theta3 theta3], [0 mrho3], '-k', 'LineWidth', 2);
        if theta3 >= 0
            theta3proj = mod(theta3, -pi);
        else
            theta3proj = mod(theta3, pi);
        end

        polarplot([theta3proj theta3proj], [0 mrho3], ':k', 'LineWidth', 2);
        polarscatter(theta1o, rho1o, 20, 'r');
        polarscatter(theta2o, rho2o, 20, 'b');
        polarscatter(theta3o, rho3o, 20, 'k');

        % 获取并处理预测角度
        e1_pred_angle_deg = junctions.e1_Apred_angle(num);
        e2_pred_angle_deg = junctions.e2_Apred_angle(num);
        theta3proj_deg = rad2deg(theta3proj);
        
        % 角度规范化处理
        while theta3proj_deg < 0
            theta3proj_deg = theta3proj_deg + 360;
        end
        while theta3proj_deg >= 360
            theta3proj_deg = theta3proj_deg - 360;
        end

        % 根据流向调整预测角度
        handedness = junctions.handedness(num);
        handedness = string(handedness);
        if strcmp(handedness, 'Right')
            e1_pred_angle_deg = theta3proj_deg - e1_pred_angle_deg;
            e2_pred_angle_deg = theta3proj_deg + e2_pred_angle_deg;
        elseif strcmp(handedness, 'Left')
            e1_pred_angle_deg = theta3proj_deg + e1_pred_angle_deg;
            e2_pred_angle_deg = theta3proj_deg - e2_pred_angle_deg;
        else
            e1_pred_angle_deg = theta3proj_deg + e1_pred_angle_deg;
            e2_pred_angle_deg = theta3proj_deg + e2_pred_angle_deg;
        end

        % 转换回弧度制
        e1_pred_angle = deg2rad(e1_pred_angle_deg);
        e2_pred_angle = deg2rad(e2_pred_angle_deg);
        e1_pred_angle = atan2(sin(e1_pred_angle), cos(e1_pred_angle));
        e2_pred_angle = atan2(sin(e2_pred_angle), cos(e2_pred_angle));

        % 绘制预测角度线
        polarplot([e1_pred_angle e1_pred_angle], [0 mrho1], '--r', 'LineWidth', 2);
        polarplot([e2_pred_angle e2_pred_angle], [0 mrho2], '--b', 'LineWidth', 2);

        % 配置极坐标显示
        ax = gca;
        ax.ThetaTickLabel = {'90', '60', '30', '0', '330', '300', '270', '240', '210', '180', '150', '120'};
        ax.RTickLabel = {''};
        title(['交汇点 ' num2str(num)]);
        legend('支流1投影', '支流2投影', '下游投影', '上游平面',...
              '支流1实测', '支流2实测', '下游实测',...
              '支流1预测', '支流2预测', 'Location','best');
        if ~verLessThan('matlab', '9.5')
            disableDefaultInteractivity(ax);
        end 

        % 添加角度标注
        e1_obs_angle = junctions.e1_obs_angle(num);
        e2_obs_angle = junctions.e2_obs_angle(num);
        text(theta1 + 0.1, mrho1/2, sprintf('实测角度一: %.1f°', e1_obs_angle),...
             'Color', 'r', 'HorizontalAlignment', 'center');
        text(theta2 + 0.1, mrho2/2, sprintf('实测角度二: %.1f°', e2_obs_angle),...
             'Color', 'b', 'HorizontalAlignment', 'center');
        text(e1_pred_angle + 0.1, mrho1/2, sprintf('预测角度一: %.1f°', e1_pred_angle_deg),...
             'Color', 'r', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        text(e2_pred_angle + 0.1, mrho2/2, sprintf('预测角度二: %.1f°', e2_pred_angle_deg),...
             'Color', 'b', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

        hold off

        % 交互式继续选择对话框
        answer = questdlg('是否继续选择其他交汇点？', '继续选择', '是', '否', '否');
        if strcmp(answer, '是')
            continueSelect = true;
            num = []; 
        else
            continueSelect = false;
        end
    end
end

function [n_x, n_y] = RotCoord(x, y, theta, x0, y0)
    % 坐标旋转变换函数
    % 输入：原始坐标(x,y)，旋转角度theta，旋转中心(x0,y0)
    % 输出：旋转后的坐标(n_x,n_y)
    n_x = (x - x0).*cos(theta) - (y - y0).*sin(theta);
    n_y = (x - x0).*sin(theta) + (y - y0).*cos(theta);
end