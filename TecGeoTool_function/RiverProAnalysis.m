function RiverProAnalysis(river_i,str1,outpath)
% 描述:
% 该函数用于河流纵剖面分析，支持交互式裂点选择，计算河道凹度、ksn值等形态参数，
% 生成包含分段特征的矢量文件及统计分析结果。函数通过图形界面引导用户选择汇水面
% 积阈值和裂点位置，支持在logS-A图和chi-z图两种模式下进行裂点识别。
%
% 必要输入:
% river_i - 河流编号（数值型），用于结果文件标识
% str1 - 输入数据文件路径（字符串），应为包含以下列的文本文件：
% [Y坐标, X坐标, 溯源距离(m), 高程(m), 流向, 汇水面积(m²), Chi值]
% outpath - 输出文件保存路径（字符串），需包含完整目录路径及文件名前缀
%
% 输出文件:
% 1. [outpath]_Sfeature[river_i].shp - 矢量线文件，包含各河段形态参数：
% (ID, RiverNo, theta, std_theta, ksn, std_ksn, DWksn, DWstd_ksn等)
% 2. [outpath]_Pnt_divide.txt - 分界点特征文本文件，包含：
% [河流编号, X坐标, Y坐标, 凹度, 误差, ksn, 误差, 校正ksn, 误差,
% 汇水面积范围(km²), 高程极值(m), 最大溯源距(m), Chi值]
% 3. [outpath]_Pnt_analysis.txt - 河段分析结果文本文件，包含：
% [河流编号, 坐标, 溯源距, 凹度, log_ks]
% 4. [str1]New.txt - 增强版纵剖面数据文件，新增：
% [平滑高程, 局部ksn, 坡度]等字段
%
% 交互操作说明:
% 1. 程序运行后将弹出组合分析图窗，包含：
% - 汇水面积-距离半对数图
% - 坡度-距离半对数图
% - 面积-坡度双对数图
% - 高程-距离剖面图
% - Chi-高程关系图
% 2. 根据提示在图形界面完成以下操作：
% a) 在汇水面积图中点击选择面积阈值点
% b) 输入需要剔除的高程极值（可选）
% c) 选择裂点识别模式（logS-A图或chi-z图）
% d) 在选定图中点击确定裂点位置
%
% 注意事项:
% 1. 输入文本文件需确保列顺序与格式正确
% 2. 矢量文件生成需MATLAB Mapping Toolbox支持
% 3. 建议剖面数据点间距≤50m以保证计算精度
% 4. 图形界面操作时请按提示顺序点击，避免坐标轴切换
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数编写者: 王一舟 - 更新日期: 2025年2月4日
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 绘图时候的坐标显示约束
%minLim_Area=10^1;maxLim_Area=10^9;
% minLim_Len=0;maxLim_Len=5e4;
%minLim_Slope=10^(-3);maxLim_Slope=10^0;
%MaxChi=25; MinEle=1100; MaxEle=MinEle+1000; %(floor(max(Ele)/100)+1)*100; 区域相对高差最大5000，而河流最低处约为500

% ******
river_i=double(river_i);
str2='.txt'; %str_i=strcat(str1,str2);

riverPro_old=textread(str1);

% 计算logS A，新的riverPro 除了原先的，还包括平滑高程、局部ksn、坡度
SmoothLen=50; ResampleWin=40; %单位 m，平滑距离，高程重采样距离
riverPro=CalLogSA(riverPro_old,SmoothLen,ResampleWin); riverPro_old=[];
[m,n]=size(riverPro);

% 存储新的河流剖面: 纵剖面点的平面坐标Y、X、溯源距离(3)、高程(4)、流向、汇水面积(6)、Chi距离、分段ksn、平滑高程、局部ksn、坡度(10)
str_iNew=strcat(str1,'New',str2); fid1=fopen(str_iNew,'w');
for i=1:m
    for j=1:n-1; fprintf(fid1,'%f\t',riverPro(i,j)); end
    fprintf(fid1,'%f',riverPro(i,n));fprintf(fid1,'\n');
end
fclose(fid1);
%%

% 创建一个新图窗口
h_combined = figure; 
set(h_combined, 'units', 'normalized', 'position', [0.05 0.05 0.85 0.85]);

% 第一个图 (x-log(Area))
h1=subplot(3, 2, 1); 
set(gca, 'position', [0.1 0.71 0.35 0.25]);
semilogy(riverPro(:, 3), riverPro(:, 6));
xlabel('Length (m)'); ylabel('Area (m^2)'); hold on;

% 第二个图 (x-log(Slope))
h2=subplot(3, 2, 3); 
set(gca, 'position', [0.1 0.39 0.35 0.25]);
semilogy(riverPro(:, 3), riverPro(:, 10), '+');
xlabel('Length (m)'); ylabel('Slope'); hold on;

% 第三个图 (log(Area-Slope))
h3=subplot(3, 2, 5); 
set(gca, 'position', [0.1 0.07 0.35 0.25]);
loglog(riverPro(:, 6), riverPro(:, 10), '+');hold on;
xlabel('Area (m^2)'); ylabel('Slope'); 

% 第四个图 (x-Elevation)
h4=subplot(3, 2, 2); 
set(gca, 'position', [0.55 0.57 0.35 0.4]); % 确保位置不被覆盖
plot(riverPro(:, 3), riverPro(:, 4)); hold on;% 绘制完整的河道数据
xlabel('Length (m)'); ylabel('Elevation (m)');


% 弹出对话框提示用户
choice = questdlg('请开始选择汇水面积阈值', '选择提示', ...
    '确认', '取消', '确认');

% 根据用户选择执行操作
if strcmp(choice, '确认')
    % 选择汇水面积阈值
    [Ac, ~] = ginput(1); 

    % 后续计算逻辑
    i = m;
    while riverPro(i, 6) < Ac 
        i = i - 1; 
    end
else
    disp('操作已取消');
end
% 剔除汇水面积小于Ac的行
riverPro = riverPro(1:i, :); m = size(riverPro, 1);            % 新的数据量

% 第四个图 (x-Elevation)
subplot(h4);
plot(riverPro(i, 3), riverPro(i, 4), '*');hold on;% 标记河道顶点
xlabel('Length (m)'); ylabel('Elevation (m)');

% 在顶点处标记汇水面积阈值 Acr
Acr_km2 = riverPro(i, 6) / 1e6; % 转换为 km²
labelText = strcat('A_{cr} = ', num2str(Acr_km2, '%.2f'), ' km^2'); % 格式化字符串

% 在河道顶点位置添加标注
text(riverPro(i, 3), riverPro(i, 4), labelText, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
    'FontSize', 10, 'Color', 'red');

% 第五个图 (Chi-Elevation)
h5=subplot(3, 2, 6); 
set(gca, 'position', [0.55 0.05 0.35 0.4]);
plot(riverPro(:, 7), riverPro(:, 4), 'g', 'linewidth', 3);
xlabel('Chi'); ylabel('Elevation (m)');
hold on;

%% Fig3，用于DW校正后的chi-z
%h3=figure; 
% 根据高程剔除小的高程
minorEle = 0;
while true
    prompt = {'请输入剔除的小的高程值 (默认值-1，请大于0)：'};
    dlgtitle = '剔除小的高程';
    dims = [1 50];
    definput = {'-1'};
    answer = inputdlg(prompt, dlgtitle, dims, definput);
    
    if isempty(answer)
        break;
    end
    
    minorEle = str2double(answer{1});
    if minorEle <= 0
        break;
    end
    
    for iErase = 1:m
        if Ele(iErase) > minorEle
            riverPro = riverPro(iErase:m, :);
            break;
        end
    end
    
    % 更新高程数据
    Chi = riverPro(:, 7);  Ele = riverPro(:, 4);   m = length(Chi);
end

% 根据高程剔除大的高程
largeEle = 0;
while true
    prompt = {'请输入剔除的大的高程值 (默认值-1，请大于0)：'};
    dlgtitle = '剔除大的高程';
    dims = [1 50];
    definput = {'-1'};
    answer = inputdlg(prompt, dlgtitle, dims, definput);
    
    if isempty(answer)
        % 如果用户取消输入对话框，退出循环
        break;
    end
    
    largeEle = str2double(answer{1});
    if largeEle <= 0
        break;
    end
    
    for jErase = 1:m
        if Ele(jErase) > largeEle
            riverPro = riverPro(1:(jErase-1), :);
            break;
        end
    end
    
    % 更新高程数据
    Chi = riverPro(:, 7);    Ele = riverPro(:, 4);    m = length(Chi);
end


% 剔除设置的高程后新的河流纵剖面 % set(gca,'ytick',MinEle:100:MaxEle);% text(Chi(1:2:m),Ele(1:2:m),num2str(Ele(1:2:m)));
Chi=riverPro(:,7);Ele=riverPro(:,4);Area=riverPro(:,6);Slope=riverPro(:,10); m=length(Chi);
%figure(h2); subplot(h22); plot(Chi,Ele,'r','linewidth',3);xlabel('Chi (m)');ylabel('Elevation (m)');hold on;
% 第五个图 (Chi-Elevation)
subplot(h5); 
plot(Chi, Ele, 'r', 'linewidth', 3);
xlabel('Chi'); ylabel('Elevation (m)');hold on;
hold on;
%% 程序主干
PntKs = [];
Pnt_analysis = [];

% 1. 输入分割点数量
prompt = {'请输入裂点数量 ：'};
dlgtitle = '输入裂点数量';
dims = [1 50];
definput = {'0'};
answer = inputdlg(prompt, dlgtitle, dims, definput);

if isempty(answer)
    disp('操作取消');
    return;
end

PntSum = str2double(answer{1});

%% 2. 选择找寻分割点的方法：1—从logS-A；2—从chi-z图
if PntSum>0
choicePrompt = {'请选择找寻裂点的方法：1—从logS-A；2—从chi-z图'};
choiceDlgTitle = '选择方法';
choiceDefInput = {'1'};
choiceAnswer = inputdlg(choicePrompt, choiceDlgTitle, dims, choiceDefInput);

if isempty(choiceAnswer)
    %disp('操作取消');
    return;
end

choose_fig = str2double(choiceAnswer{1});
if isnan(choose_fig) || ~ismember(choose_fig, [1, 2])
    %disp('无效的选择，程序终止。');
    return;
end

if choose_fig == 1
    % 3. 从 logS-A 图选择分割点
    %disp(['从 logS-A 图选择 ', num2str(PntSum), ' 个分割点 (按汇水面积)。']);
    uiwait(msgbox('请从 logS-A 图中选择裂点。', '提示', 'modal'));
    figure(h_combined);  subplot(h3);  hold on; 
    [PntArea, ~] = ginput(PntSum);
    PntArea = sort(PntArea, 'descend');

    PntNum = [1]; % 初始化分界点序号，包含起点
    for di = 1:PntSum
        for k = 1:m-1
            if (PntArea(di) <= Area(k) && PntArea(di) >= Area(k+1))
                PntNum = [PntNum; k]; 
                %disp([num2str(Area(k) / 1e6), ' km²； ', num2str(Ele(k)), ' m']);
                break;
            end
        end
    end
    PntNum = [PntNum; m]; % 包含终点
    PntSum = length(PntNum) - 2; % 实际分割点数

elseif choose_fig == 2
    % 4. 从 chi-z 图选择分割点
    %disp(['从 chi-z 图选择 ', num2str(PntSum), ' 个裂点 (按 χ 值)。']);
    uiwait(msgbox('请从 chi-z 图中选择裂点。', '提示', 'modal'));
    figure(h_combined);subplot(h5);  hold on; % 保持原有图像
    [PntChi, ~] = ginput(PntSum);
    PntChi = sort(PntChi);

    PntNum = [1]; % 初始化分界点序号，包含起点
    for j = 1:PntSum
        for k = 1:m-1
            if (PntChi(j) >= Chi(k) && PntChi(j) <= Chi(k+1))
                PntNum = [PntNum; k];
                %disp([num2str(Area(k) / 1e6), ' km²； ', num2str(Ele(k)), ' m']);
                break;
            end
        end
    end
    PntNum = [PntNum; m]; % 包含终点
    PntSum = length(PntNum) - 2; % 实际分割点数
end
end
%disp(['最终，从大到小，实际分界点 ', num2str(PntSum), ' 个。']);

%%
if PntSum==0
    [R2,z0,ksn,std_ksn]=Bi_Regress(Ele,Chi); % 线性回归 Y=a0+a1*X
    figure(h_combined); subplot(h5);plot(Chi,Chi.*ksn+z0,'k');hold on;
    text(max(Chi)/2,max(Ele),['R=',num2str(R2),'z=Chi.*',num2str(ksn),'±',num2str(std_ksn),'+',num2str(z0)]); % 这个是在Fig2上面画的
    [DW_r,DW_Ele,DW_Chi]=DW_test(Ele,Chi,z0,ksn); %DW检验
    [R2,DWz0,DWksn,DWstd_ksn]=Bi_Regress(DW_Ele,DW_Chi); 
    figure(2); plot(DW_Chi,DW_Ele,'bo'); hold on; plot(DW_Chi,DW_Chi.*DWksn+DWz0,'k'); text(max(DW_Chi)/2,max(DW_Ele),['ro=',num2str(DW_r),'R=',num2str(R2),'ksn=',num2str(DWksn),'±',num2str(DWstd_ksn)]);
    [R2,log_ks,theta,std_theta]=Bi_Regress(log(Slope),log(Area));
    %[theta2,ksn2]=Nonlinear_Regress(Area,Slope);
    figure(h_combined); subplot(h3); loglog(Area,(Area.^(theta)).*(exp(log_ks)),'k');hold on; text(max(Area)/2,0.1,['R=',num2str(R2),'θ=',num2str(theta),'±',num2str(std_theta)]);
    %% 输出分界点信息：河流编号，X、Y坐标，河段凹度、误差，ksn、误差，校正后ksn、误差，折点汇水面积（Km2，min,max），折点高程(max,min)、最高折点的溯源距离、chi值
    PntKs=[PntKs;[river_i,riverPro(m,2),riverPro(m,1),theta,std_theta,ksn,std_ksn,DWksn,DWstd_ksn,(riverPro(m,6))./(1e6),(riverPro(1,6))./(1e6),riverPro(m,4),riverPro(1,4),riverPro(m,3),riverPro(m,7)]];
    Pnt_analysis=[Pnt_analysis;[river_i,riverPro(m,1),riverPro(m,2),riverPro(m,3),theta,log_ks]];    
    %disp(river_i);disp(riverPro(m,2));disp(riverPro(m,1));disp(theta);disp(std_theta);
    %% 生成矢量线,新的riverPro: 纵剖面点的平面坐标Y、X、溯源距离、高程、流向、汇水面积、Chi距离、平滑高程、局部ksn、坡度
    Sfea=struct();
    
    Sfea.Geometry='Line';
    Xcor=riverPro(:,2); Xcor=Xcor';Sfea.X=[Xcor,NaN];
    Ycor=riverPro(:,1); Ycor=Ycor';Sfea.Y=[Ycor,NaN];
    Sfea.ID=0; Sfea.RiverNo=river_i; 
    Sfea.theta=theta; Sfea.std_theta=std_theta;
    Sfea.ksn=ksn; Sfea.std_ksn=std_ksn;
    Sfea.DWksn=DWksn; Sfea.DWstd_ksn=DWstd_ksn;
    Sfea.Amin_km=(riverPro(m,6))./(1e6); Sfea.Amax_km=(riverPro(1,6))./(1e6);
    Sfea.Emax_m=riverPro(m,4); Sfea.Emin_m=riverPro(1,4);
    Sfea.Lenmax_m=riverPro(m,3); Sfea.Chimax=riverPro(m,7); 
    shapewrite(Sfea,strcat(outpath,'_Sfeature',num2str(river_i),'.shp'));
        
else % >=1个分界点 
    %% 填充 PntKs    % 第一段
    Ele_1=Ele(PntNum(1):PntNum(2)); Chi_1=Chi(PntNum(1):PntNum(2)); Area_1=Area(PntNum(1):PntNum(2)); Slope_1=Slope(PntNum(1):PntNum(2));
    figure(h_combined); subplot(h4); plot(riverPro(PntNum(2),3),riverPro(PntNum(2),4),'kx');hold on; %在UpLen-z上标记河道裂点 
    [R2,z0_1,ksn_1,std_ksn_1]=Bi_Regress(Ele_1,Chi_1); % 线性回归 Y=a0+a1*X
    %[theta2_1,ksn2_1]=Nonlinear_Regress(Area_1,Slope_1);
    figure(h_combined); subplot(h5); plot(Chi_1,Chi_1.*ksn_1+z0_1,'k');hold on; text(max(Chi_1)/2,max(Ele_1),['R1=',num2str(R2),'z1=Chi.*',num2str(ksn_1),'±',num2str(std_ksn_1),'+',num2str(z0_1)]);% 这个是在Fig2上面画的
    [DW_r_1,DW_Ele_1,DW_Chi_1]=DW_test(Ele_1,Chi_1,z0_1,ksn_1); %DW检验
    [R2,DWz0_1,DWksn_1,DWstd_ksn_1]=Bi_Regress(DW_Ele_1,DW_Chi_1); 
    figure(2); plot(DW_Chi_1,DW_Ele_1,'bo'); hold on; plot(DW_Chi_1,DW_Chi_1.*DWksn_1+DWz0_1,'k'); text(max(DW_Chi_1)/2,max(DW_Ele_1),['ro1=',num2str(DW_r_1),'R=',num2str(R2),'ksn=',num2str(DWksn_1),'±',num2str(DWstd_ksn_1)]);
    [R2,log_ks_1,theta_1,std_theta_1]=Bi_Regress(log(Slope_1),log(Area_1));
    figure(h_combined); subplot(h3); loglog(Area_1,(Area_1.^(theta_1)).*(exp(log_ks_1)),'k');hold on; text(max(Area_1)/2,max(Slope_1)/2,['R1=',num2str(R2),'θ=',num2str(theta_1),'±',num2str(std_theta_1)]);
    %% 输出分界点信息：河流编号，X、Y坐标，折点以下河段凹度、误差，ksn、误差，校正后ksn、误差，折点汇水面积（Km2，min,max），折点高程(max,min)、最高折点的溯源距离、chi值
    PntKs=[PntKs;[river_i,riverPro(PntNum(2),2),riverPro(PntNum(2),1),theta_1,std_theta_1,ksn_1,std_ksn_1,DWksn_1,DWstd_ksn_1,(riverPro(PntNum(2),6))./(1e6),(riverPro(PntNum(1),6))./(1e6),...
        riverPro(PntNum(2),4),riverPro(PntNum(1),4),riverPro(PntNum(2),3),riverPro(PntNum(2),7)]];
     Pnt_analysis=[Pnt_analysis;[river_i,riverPro(PntNum(2),1),riverPro(PntNum(2),2),riverPro(PntNum(2),3),theta_1,log_ks_1]];  
    %% 生成矢量线,新的riverPro: 纵剖面点的平面坐标Y、X、溯源距离、高程、流向、汇水面积、Chi距离、平滑高程、局部ksn、坡度
    Sfea=struct();
    
    Sfea(1).Geometry='Line';
    Xcor=riverPro(PntNum(1):PntNum(2),2); Xcor=Xcor'; Sfea(1).X=[Xcor,NaN];
    Ycor=riverPro(PntNum(1):PntNum(2),1); Ycor=Ycor'; Sfea(1).Y=[Ycor,NaN];
    Sfea(1).ID=0; Sfea(1).RiverNo=river_i; 
    Sfea(1).theta=theta_1; Sfea(1).std_theta=std_theta_1;  
    Sfea(1).ksn=ksn_1;     Sfea(1).std_ksn=std_ksn_1; 
    Sfea(1).DWksn=DWksn_1; Sfea(1).DWstd_ksn=DWstd_ksn_1;
    Sfea(1).Amin_km=(riverPro(PntNum(2),6))./(1e6); Sfea(1).Amax_km=(riverPro(PntNum(1),6))./(1e6);
    Sfea(1).Emax_m=riverPro(PntNum(2),4); Sfea(1).Emin_m=riverPro(PntNum(1),4);
    Sfea(1).Lenmax_m=riverPro(PntNum(2),3); Sfea(1).Chimax=riverPro(PntNum(2),7);
        
    %% ************剩余的河段*************************************************
    for j=2:(PntSum+1)
        Ele_j=Ele((PntNum(j)+1):PntNum(j+1)); Chi_j=Chi((PntNum(j)+1):PntNum(j+1)); Area_j=Area((PntNum(j)+1):PntNum(j+1)); Slope_j=Slope((PntNum(j)+1):PntNum(j+1));
        figure(h_combined); subplot(h4); plot(riverPro(PntNum(j+1),3),riverPro(PntNum(j+1),4),'kx');hold on; %在UpLen-z上标记河道裂点 
        [R2,z0_j,ksn_j,std_ksn_j]=Bi_Regress(Ele_j,Chi_j); % 线性回归 Y=a0+a1*X
        figure(h_combined); subplot(h5); plot(Chi_j,Chi_j.*ksn_j+z0_j,'k');hold on; text(max(Chi_j)/2,max(Ele_j),['R',num2str(j),'=',num2str(R2),'z=Chi.*',num2str(ksn_j),'±',num2str(std_ksn_j),'+',num2str(z0_j)]);% 这个是在Fig2上面画的       
        [DW_r_j,DW_Ele_j,DW_Chi_j]=DW_test(Ele_j,Chi_j,z0_j,ksn_j); %DW检验
        [R2,DWz0_j,DWksn_j,DWstd_ksn_j]=Bi_Regress(DW_Ele_j,DW_Chi_j);
        %[theta2_j,ksn2_j]=Nonlinear_Regress(Area_j,Slope_j);
        figure(2); plot(DW_Chi_j,DW_Ele_j,'bo'); hold on; plot(DW_Chi_j,DW_Chi_j.*DWksn_j+DWz0_j,'k'); text(max(DW_Chi_j)/2,max(DW_Ele_j),['ro',num2str(j),'=',num2str(DW_r_j),'R=',num2str(R2),'ksn=',num2str(DWksn_j),'±',num2str(DWstd_ksn_j)]);
 %      figure (3); plot(Chi_j,Chi_j.*DWksn_j+DWz0_j/(1-DW_r_j),'k');hold on; text(max(Chi_j)/2,max(Ele_j)/2,['zj=Chi.*',num2str(DWksn_j),'±',num2str(DWstd_ksn_j),'+',num2str(DWz0_j/(1-DW_r_j))]); % 这个是在Fig3上面画的
        [R2,log_ks_j,theta_j,std_theta_j]=Bi_Regress(log(Slope_j),log(Area_j));
        figure(h_combined); subplot(h3); loglog(Area_j,(Area_j.^(theta_j)).*(exp(log_ks_j)),'k');hold on; text(max(Area_j)/2,max(Slope_j)/2,['R',num2str(j),'=',num2str(R2),'θ=',num2str(theta_j),'±',num2str(std_theta_j)]);
        %% 输出分界点信息：河流编号，X、Y坐标，折点以下河段凹度、误差，ksn、误差，校正后ksn、误差，折点汇水面积（Km2，min,max），折点高程(max,min)、最高折点的溯源距离、chi值
        PntKs=[PntKs;[river_i,riverPro(PntNum(j+1),2),riverPro(PntNum(j+1),1),theta_j,std_theta_j,ksn_j,std_ksn_j,DWksn_j,DWstd_ksn_j,(riverPro(PntNum(j+1),6))./(1e6),(riverPro(PntNum(j),6))./(1e6),...
            riverPro(PntNum(j+1),4),riverPro(PntNum(j),4),riverPro(PntNum(j+1),3),riverPro(PntNum(j+1),7)]];
        Pnt_analysis=[Pnt_analysis;[river_i,riverPro(PntNum(j+1),1),riverPro(PntNum(j+1),2),riverPro(PntNum(j+1),3),theta_j,log_ks_j]];  
       %% 生成矢量线,新的riverPro: 纵剖面点的平面坐标Y、X、溯源距离、高程、流向、汇水面积、Chi距离、平滑高程、局部ksn、坡度
         Sfea(j).Geometry='Line';
        Xcor=riverPro((PntNum(j)+1):PntNum(j+1),2); Xcor=Xcor'; Sfea(j).X=[Xcor,NaN];
        Ycor=riverPro((PntNum(j)+1):PntNum(j+1),1); Ycor=Ycor'; Sfea(j).Y=[Ycor,NaN];
        Sfea(j).ID=0; Sfea(j).RiverNo=river_i; 
        Sfea(j).theta=theta_j; Sfea(j).std_theta=std_theta_j; 
        Sfea(j).ksn=ksn_j; Sfea(j).std_ksn=std_ksn_j; 
        Sfea(j).DWksn=DWksn_j; Sfea(j).DWstd_ksn=DWstd_ksn_j;
        Sfea(j).Amin_km=(riverPro(PntNum(j+1),6))./(1e6); Sfea(j).Amax_km=(riverPro(PntNum(j),6))./(1e6);
        Sfea(j).Emax_m=riverPro(PntNum(j+1),4); Sfea(j).Emin_m=riverPro(PntNum(j),4);
        Sfea(j).Lenmax_m=riverPro(PntNum(j+1),3); Sfea(j).Chimax=riverPro(PntNum(j+1),7);
    end
       
    shapewrite(Sfea,strcat(outpath,'_Sfeature',num2str(river_i),'.shp'));  
end
%% 
[Pnt_m,Pnt_n]=size(PntKs);
str_PntKs=strcat(outpath,'_Pnt_divide','.txt');

fid=fopen(str_PntKs,'a+');
removeLines(str_PntKs, river_i);%删除之前的记录
for i=1:Pnt_m
    for j=1:Pnt_n-1
        fprintf(fid,'%f\t',PntKs(i,j));
        %disp(PntKs(i,j));
    end
    fprintf(fid,'%f',PntKs(i,Pnt_n));fprintf(fid,'\n');
end
fclose(fid);
[analysis_m,analysis_n]=size(Pnt_analysis);

str_analysis=strcat(outpath,'_Pnt_analysis','.txt');
fid=fopen(str_analysis,'a+');
removeLines(str_analysis, river_i);%删除之前的记录
for i=1:analysis_m
    for j=1:analysis_n-1
        fprintf(fid,'%f\t',Pnt_analysis(i,j));
        %disp(Pnt_analysis(i,j))
    end
    fprintf(fid,'%f',Pnt_analysis(i,analysis_n));fprintf(fid,'\n');
end
fclose(fid);
end
%% 计算坡度
function riverPro_f=CalLogSA(riverPro_Old_f,SmoothLen_f,ResampleWin_f)
% riverPro_i：纵剖面点的平面坐标Y、X、溯源距离、高程、流向、汇水面积、Chi距离、粗略ksn
% 新的到的riverPro 除了原先的，还包括平滑高程、局部ksn、坡度
UpLen=riverPro_Old_f(:,3);Ele=riverPro_Old_f(:,4);Chi=riverPro_Old_f(:,7); m=length(UpLen); 

%计算平滑后的高程
SmoothedEle=Ele;Localksn=zeros(m,1);Slope=zeros(m,1);
for i=1:m
    if abs(UpLen(i)-UpLen(1))<=SmoothLen_f/2 % 当点i与出水口距离小于SmoothLen_f/2，则一直向后寻找直至达到阈值为止
        j=1; while i+j<=m && abs(UpLen(i+j)-UpLen(1))<=SmoothLen_f; j=j+1;end; j=j-1; 
        if (i+j)-1<2; j=3-i; end %当参与插值回归点数＜3，强制多选择点，确保插值、回归中有三个点。
        SmoothedEle(i)=interp1(UpLen(1:i+j),Ele(1:i+j),UpLen(i));
        b=regress(Ele(1:i+j),[ones(size(Chi(1:i+j))) Chi(1:i+j)]); %regress函数用法，y=a0+a1*x
        Localksn(i)=abs(b(2));
    elseif abs(UpLen(m)-UpLen(i))<=SmoothLen_f/2 % 当点i与水头距离小于SmoothLen_f/2，则一直向前寻找直至达到阈值为止
        j=1; while i-j>=1 && abs(UpLen(m)-UpLen(i-j))<=SmoothLen_f; j=j+1;end; j=j-1;
        if (i-j)>(m-2); j=i-(m-2); end %当参与插值回归点数＜3，强制多选择点，确保插值、回归中有三个点。
        SmoothedEle(i)=interp1(UpLen(i-j:m),Ele(i-j:m),UpLen(i));
        b=regress(Ele(i-j:m),[ones(size(Chi(i-j:m))) Chi(i-j:m)]); 
        Localksn(i)=abs(b(2));
    elseif abs(UpLen(i)-UpLen(1))>SmoothLen_f/2 && abs(UpLen(m)-UpLen(i))>SmoothLen_f/2 %点与出水口、水头距离都超过SmoothLen_f/2
        j=1; while i-j>=1 && i+j<=m && abs(UpLen(i+j)-UpLen(i-j))<=SmoothLen_f; j=j+1; end; j=j-1;
        if j==0; j=1; end %当搜索点的左右相邻点的距离＜SmoothLen_f，强制选择相邻的点，确保插值、回归中有三个点。
        SmoothedEle(i)=interp1(UpLen(i-j:i+j),Ele(i-j:i+j),UpLen(i));
        b=regress(Ele(i-j:i+j),[ones(size(Chi(i-j:i+j))) Chi(i-j:i+j)]);
        Localksn(i)=abs(b(2));
    end
end

%计算重采样后的Localksn和Slope
for i=1:m
    if abs(Ele(i)-Ele(1))<=ResampleWin_f/2 % 当点i与出水口高程小于ResampleWin_f/2，则一直向后寻找直至达到阈值为止
        j=1; while i+j<=m && abs(Ele(i+j)-Ele(1))<=ResampleWin_f; j=j+1;end; j=j-1; 
        if (i+j)-1<2; j=3-i; end %当参与插值回归点数＜3，强制多选择点，确保插值、回归中有三个点。
        b=regress(Ele(1:i+j),[ones(size(UpLen(1:i+j))) UpLen(1:i+j)]); %regress函数用法，y=a0+a1*x
        Slope(i)=abs(b(2));
    elseif abs(Ele(m)-Ele(i))<=ResampleWin_f/2 % 当点i与水头高程小于ResampleWin_f/2，则一直向前寻找直至达到阈值为止
        j=1; while i-j>=1 && abs(Ele(m)-Ele(i-j))<=ResampleWin_f; j=j+1;end; j=j-1;
        if (i-j)>(m-2); j=i-(m-2); end %当参与插值回归点数＜3，强制多选择点，确保插值、回归中有三个点。
        b=regress(Ele(i-j:m),[ones(size(UpLen(i-j:m))) UpLen(i-j:m)]); 
        Slope(i)=abs(b(2));
    elseif abs(Ele(i)-Ele(1))>ResampleWin_f/2 && abs(Ele(m)-Ele(i))>ResampleWin_f/2 %点与出水口、水头高程都超过ResampleWin_f
        j=1; while i-j>=1 && i+j<=m && abs(Ele(i+j)-Ele(i-j))<=ResampleWin_f; j=j+1; end; j=j-1;
        if j==0; j=1; end %当搜索点的左右相邻点的高差＜ResampleWin_f，强制选择相邻的点，确保插值、回归中有三个点。
        b=regress(Ele(i-j:i+j),[ones(size(UpLen(i-j:i+j))) UpLen(i-j:i+j)]);
        Slope(i)=abs(b(2));
       
    end
end

riverPro_f=[riverPro_Old_f SmoothedEle Localksn Slope];
end
%% 统计检验
function [DW_r_f,DW_Ele_f,DW_Chi_f]=DW_test(Ele_f,Chi_f,z0_f,ksn_f)
DW_e=Ele_f-(z0_f+ksn_f.*Chi_f); % 残差
m=length(DW_e);
DW_e1=DW_e(1:m-1); DW_e2=DW_e(2:m);
DW_r_f=sum(DW_e1.*DW_e2)/(sum(DW_e1.^2));
DW_Chi_f=Chi_f(2:m)-DW_r_f.*Chi_f(1:m-1);
DW_Ele_f=Ele_f(2:m)-DW_r_f.*Ele_f(1:m-1);
end
%% 线性回归
function [R2,a0,a1,std_a1]=Bi_Regress(Y,X)
% 线性回归 Y=a0+a1*X
Length_X=length(X);
if Length_X>=3
    R2=corrcoef(Y,X); if length(R2)>1;R2=R2(1,2); R2=R2^2;else R2=R2^2;end
    mean_X=mean(X); mean_Y=mean(Y);
    Sxx=sum((X-mean_X).^2);Syy=sum((Y-mean_Y).^2);
    Sxy=sum((X-mean_X).*(Y-mean_Y));
    a1=Sxy/Sxx; a0=mean_Y-a1*mean_X;
    sigma=sqrt((Syy-a1*Sxy)/(Length_X-2));
    std_a1=sigma/sqrt(Sxx);
elseif Length_X==2
    R2=1;
    a1=(Y(1)-Y(2))./(X(1)-X(2)); 
    std_a1=0;
    a0=Y(1)-a1*X(1);
else
    %回归数据不够，不可模拟
    R2=9999;a0=9999;a1=9999;std_a1=9999;    
end
end

%% 删除之前河流记录，避免多次选择造成导出数据出问题
function removeLines(file_path, starting_number)
    % 读取文件内容
    fileID = fopen(file_path, 'r');
    lines = textscan(fileID, '%s', 'Delimiter', '\n');
    fclose(fileID);

    % 找到以指定数开头的行并删除
    lines = lines{1};
    modified_lines = lines(~startsWith(lines, num2str(starting_number)));

    % 将修改后的内容写回文件
    fileID = fopen(file_path, 'w');
    fprintf(fileID, '%s\n', modified_lines{:});
    fclose(fileID);
end
