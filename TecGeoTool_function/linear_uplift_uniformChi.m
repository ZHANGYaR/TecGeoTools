function linear_uplift_uniformChi(str,startRiver,endRiver,numRange_1,K1,K2)
% 使用方法:
%   linear_uplift_uniformChi(str,startRiver,endRiver,numRange_1,K1,K2)
%
% 描述:
%   该函数用于模拟山体均匀隆升历史，通过合并多条河流纵剖面数据计算侵蚀系数。
%   基于chi空间分析方法，将多条河流的纵剖面标准化后进行分段线性拟合，
%   重建统一的山体隆升历史并生成可视化的速率-时间曲线。
%   输出包含：相对高程-chi关系图、标准化隆升速率曲线、多侵蚀系数对比曲线。
%
% 必要输入:
%   str        - 河流纵剖面文件路径前缀（文本文件命名规则示例：StudyArea_RiverPro1.txt）
%   startRiver - 起始河流编号（分析的河流范围起始编号）
%   endRiver   - 结束河流编号（分析的河流范围终止编号）
%   numRange_1 - chi值间隔（用于分段计算局部侵蚀系数的空间分段长度，单位：米）
%   K1         - 侵蚀系数1（基准侵蚀系数，单位：m/yr）
%   K2         - 侵蚀系数2（对比侵蚀系数，单位：m/yr）

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数编写者: 王一舟 - 更新日期: 2024年12月4日
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



format long; close all;
% 均匀间隔 线性模拟山体隆升历史，去掉最小、最大高程
% Line13 河道高程剖面存储路径

%% 计算侵蚀系数
%K1=1.5e-6; K2=1e-6; 
% K1=(1.08e-6); K2=(1.08e-6)-2.85e-7; K3=(1.08e-6)+2.85e-7;
% std=2.85e-7; 4.42e-7 侵蚀系数，所有的需要换算成国际单位制 m/a
% 1.47中,1.06东；1.10西，0.891中

%% 输入需要合并的河流个数
%所有的河流纵剖面都存在这个文件夹中
str2 = '_RiverPro';  Path3='.txt';
%startRiver=input('输入开始计算的河流编号：startRiver='); % 这里从5开始  注意 Inyo山从5开始
%endRiver=input('输入最终计算的河流编号：endRiver=');
sumRiver=endRiver-startRiver+1; %河流个数
AllRiverPro=cell(sumRiver,1); % 利用cell存储所有的河流纵剖面矩阵
lengthChi=[];
for i=startRiver:endRiver
    % numRiver_i=input('输入需要合并的河流编号：');
    Path2=num2str(i);
    Path_i=strcat(str,str2,Path2,Path3);
    RiverPro_i=textread(Path_i); % Y、X坐标，溯源距离、高程、流向、汇水面积、Chi距离
    chiLen_i=RiverPro_i(:,7); Ele_i=RiverPro_i(:,4);
    minEle_i=min(Ele_i); Ele_i=Ele_i-minEle_i;
    plot(chiLen_i,Ele_i,'Color', [0.5 0.5 0.5]);hold on; % 上述河流已经是相对高程
    AllRiverPro{i-startRiver+1,1}=[chiLen_i,Ele_i];% 利用cell存储所有的河流纵剖面矩阵
    lengthChi=[lengthChi;max(chiLen_i)]; % 记录每一河流长度 最大chhi
end

%% 将河流合并为一条
%numRange_1=input('输入最小的chi值间隔'); % 1.3;   %input('输入时间节点间隔：'); % 每段chi的长度

MaxChi=max(lengthChi); lengthChi=[]; % 寻找最长的河流长度
n_localksn=ceil(MaxChi/numRange_1); % 向上取整，chi空间分段数目；n_localksn=floor(MaxChi/numRange_1);
localksn=zeros(2,n_localksn); %Line1: 每一段的local_ksn值; Line2: 每一段计算了几个ksn

interChi=0:numRange_1:MaxChi; % 合并之后的chi空间，等分间隔
if max(interChi) < MaxChi; interChi=[interChi,MaxChi]; end  % interChi长度比localksn多1

interEle=zeros(1,length(interChi)); % interEle=zeros(1,n_localksn+1);

for j=1:sumRiver
    RiverProj=AllRiverPro{j,1}; 
    chiLen_j=RiverProj(:,1); Ele_j=RiverProj(:,2); % 都是列向量
    max_chiLen_j=max(chiLen_j);
    
    i=1; % 第1段，chi坐标（1,1+1）
    while i <= n_localksn % 对于每一条河chi-z，按照分段标准，判别其chi空间最大值 是否在小段内
        if max_chiLen_j >= interChi(i+1) % 该河道chi空间，完全包含 第(i,i+1)段
%             (interChi(i),interChi(i+1))段的local_ksn，编写function函数(chiLen_j,Ele_j,interChi(i),interChi(i+1))
            local_ksn_i=cal_local_ksn(chiLen_j,Ele_j,interChi(i),interChi(i+1));
            localksn(1,i)=localksn(1,i)+local_ksn_i; % 这一段累计ksn
            localksn(2,i)=localksn(2,i)+1; % 这一段累计次数
            i=i+1; % 下一段
        elseif max_chiLen_j > interChi(i) && max_chiLen_j <= interChi(i+1) % 第(i,i+1)段包含局部
%             (interChi(i),max_chiLen_j)段的local_ksn
            local_ksn_i=cal_local_ksn(chiLen_j,Ele_j,interChi(i),max_chiLen_j);
            localksn(1,i)=localksn(1,i)+local_ksn_i;
            localksn(2,i)=localksn(2,i)+1;
            break;
        elseif max_chiLen_j < interChi(i) % 第(i,i+1)段不包含
            break;
        end
    end
end
% figure;plot(localksn(1,:),localksn(2,:),'o')
localksn=localksn(1,:)./localksn(2,:); 

for i=1:n_localksn
    interEle(i+1)=interEle(i)+(interChi(i+1)-interChi(i)).*localksn(i);
end

%%

% interChi=0:numRange_1:MaxChi; 
% if max(interChi) < MaxChi; interChi=[interChi,MaxChi]; end
interChi=interChi';
interEle=interEle';

%% % 以上步骤得到合并后的河流纵剖面
plot(interChi,interEle,'linewidth',3); hold on;
max_xlabel=ceil(ceil(max(interChi))/2)*2; xlim([0,max_xlabel]);
xlabel('Chi (m)'); ylabel('Relative elevation (m)');

%% 计算无维度的隆升历史
len_path=length(interChi);
t_path=interChi; 
U_path=(interEle(2:len_path)-interEle(1:len_path-1))./(t_path(2:len_path)-t_path(1:len_path-1));

U_path=[U_path;U_path(length(U_path))];

figure;
subplot(2,1,1);stairs(t_path,U_path); 
% axis([0,20,0,300]); 
xlabel('t* (m)'); ylabel('U*');
% set(gca,'ytick',0:50:300);
% subplot(2,1,1);stairs(t_path(1:upLen-1),U_path);
%%
upRate=(U_path).*K1*1e3; % 显示速率转换为mm/a
% upTime=zeros(upLen-1,1); %upTime1=zeros(upLen,1);
upTime=t_path./K1./1e6;
subplot(2,1,2);stairs(upTime,upRate,'k','LineWidth',2.5); hold on;

upRate=(U_path).*K2*1e3; % 显示速率转换为mm/a
% upTime=zeros(upLen-1,1); %upTime1=zeros(upLen,1);
upTime=t_path./K2./1e6;
subplot(2,1,2);stairs(upTime,upRate,'b','LineWidth',1.5); hold on;
% 添加图例
legend('K1速率', 'K2速率');

xlim([0,20]);
% axis([0,26,0,0.4]); 
xlabel('t (Ma)'); ylabel('U (mm/a)');
set(gca,'xtick',0:2:20);
% set(gca,'ytick',0:0.05:0.4);

%%
function local_ksn=cal_local_ksn(chiLen_j,Ele_j,interChi_i,interChi_i1)
% 计算local_ksn，河道整体chi-z，待计算的chi小段(interChi_i,interChi_i1)
% 都是列向量，这个小段一定在河道整个chi空间之内

m_len=length(chiLen_j);
chi_cal=[];Ele_cal=[];
for i=1:m_len
    if chiLen_j(i) >= interChi_i && chiLen_j(i) <= interChi_i1
        chi_cal=[chi_cal;chiLen_j(i)];
        Ele_cal=[Ele_cal;Ele_j(i)];
    elseif chiLen_j(i) > interChi_i1
        break;
    end
end

if length(chi_cal) < 2
    chi_cal=[];Ele_cal=[];
    chi_cal=[interChi_i;interChi_i1];
    Ele_cal=interp1(chiLen_j,Ele_j,chi_cal);
end

chi_cal=[ones(length(chi_cal),1),chi_cal];
localksn_ij=regress(Ele_cal,chi_cal);
local_ksn=localksn_ij(2);

%%


















