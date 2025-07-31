function MisFit=LinearUp_MisFit(chiDivision,str,startRiver,endRiver)
% 使用方法:
%   MisFit = LinearUp_MisFit(chiDivision,str,startRiver,endRiver)
%
% 描述:
%   该函数用于评估河流纵剖面与线性隆升模型的拟合程度。通过合并多条河流的chi-z数据建立标准剖面，
%   计算各河道实际剖面与标准剖面的均方根误差(RMSE)，量化模型预测与实际观测的偏离程度。
%   适用于山体隆升模型验证和侵蚀过程敏感性分析。
%
% 必要输入:
%   chiDivision - chi空间分段间隔（单位：米），用于标准化处理的不同河段划分长度
%   str         - 河流数据文件路径前缀（示例："StudyArea_RiverPro"，对应文件命名规则为前缀+编号.txt）
%   startRiver  - 起始河流编号（分析的河流序列起始编号，整数）
%   endRiver    - 终止河流编号（分析的河流序列结束编号，整数）
%
% 输出:
%   MisFit      - 标准化平均拟合误差，反映所有河流剖面与模拟隆升剖面的整体偏离程度（无量纲）
%
% 可选输入:
%   无
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数编写者: 王一舟 - 更新日期: 2025年2月4日
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
format long; close all;

%% 输入需要合并的河流个数
%所有的河流纵剖面都存在这个文件夹中
Path3='.txt';str2 = '_RiverPro';
%startRiver=1; % input('输入开始计算的河流编号：startRiver='); % 这里从5开始  注意 Inyo山从5开始
%endRiver=622; % input('输入最终计算的河流编号：endRiver=');
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
%     plot(chiLen_i,Ele_i,'m');hold on; % 上述河流已经是相对高程
    AllRiverPro{i-startRiver+1,1}=[chiLen_i,Ele_i];% 利用cell存储所有的河流纵剖面矩阵
    lengthChi=[lengthChi;max(chiLen_i)]; % 记录每一河流长度 最大chhi
end

%% 将河流合并为一条
% chiDivision=input('输入时间节点间隔：'); % 每段chi的长度
MaxChi=max(lengthChi); lengthChi=[]; % 寻找最长的河流长度

n_localksn=ceil(MaxChi/chiDivision); % 向上取整，chi空间分段数目；n_localksn=floor(MaxChi/chiDivision);
localksn=zeros(2,n_localksn); %Line1: 每一段的local_ksn值; Line2: 每一段计算了几个ksn

localksn_mat=zeros(sumRiver,n_localksn); %% 每一行，表示一条河流；每一列，都是相同的Chi。

interChi=0:chiDivision:MaxChi; % 合并之后的chi空间，等分间隔
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
            localksn_mat(j,i)=local_ksn_i;
        elseif max_chiLen_j > interChi(i) && max_chiLen_j <= interChi(i+1) % 第(i,i+1)段包含局部
%             (interChi(i),max_chiLen_j)段的local_ksn
            local_ksn_i=cal_local_ksn(chiLen_j,Ele_j,interChi(i),max_chiLen_j);
            localksn(1,i)=localksn(1,i)+local_ksn_i;
            localksn(2,i)=localksn(2,i)+1;
            localksn_mat(j,i)=local_ksn_i;
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
interChi=interChi';
interEle=interEle';
% plot(interChi,interEle,'g','linewidth',3);

% 计算与平均chi-z的差MisFit
MisFit=0;
for j=1:sumRiver
    RiverProj=AllRiverPro{j,1}; 
    chiLen_j=RiverProj(:,1); Ele_j=RiverProj(:,2); % 都是列向量
    yi = interp1(interChi,interEle,chiLen_j);
    MisFit=MisFit+sqrt(mean((Ele_j-yi).^2));
end
MisFit=MisFit/sumRiver;




%% 记录每一条河流，每一个间隔的local_ksn
% % 行列localksn_mat(sumRiver,n_localksn)
% Chi_ksnmat=chiDivision:chiDivision:chiDivision*n_localksn;
% mean_localksn=zeros(1,n_localksn);
% std_localksn=zeros(1,n_localksn);
% for kk=1:n_localksn
%     a=localksn_mat(:,kk);
%     b=a(a~=0);
%     mean_localksn(kk)=mean(b);
%     if length(b)>1 % 仅有一个值，就不计算std。
%         std_localksn(kk)=std(b);
%     end
% end
% figure; plot(Chi_ksnmat,mean_localksn); hold on; title('mean_loacal_ksn');
% figure; plot(Chi_ksnmat,std_localksn); hold on; title('std_loacal_ksn')
% mean(std_localksn(std_localksn~=0)); % 对非零的std计算均值



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


















