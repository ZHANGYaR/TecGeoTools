function Theta_FitChiplot(Path1,startRiver,endRiver,str1)
% ʹ�÷���:
% Theta_FitChiplot(Path1, startRiver, endRiver, str1)
%
% ����:
% ����������Mudd et al (2018)������ͨ������������������Chi-�̹߳�ϵ��
% ����ʹ�߳�ɢ����С�������źӵ����Ȧ�ֵ��֧����������������ݣ����ɦ�-MisFit
% ��ϵ���߼���ͬ��ֵ��Chi-�̶߳Ա�ͼ��
%
% ��Ҫ����:
% Path1 - �������������ݴ洢·�����ַ�����ʽ��Ҫ��·���°���[N].txt��ʽ�ļ�
% startRiver - ��ʼ������ţ����������ļ�������ʼ���
% endRiver - ��ֹ������ţ����������ļ������������
% str1 - ����ļ�ǰ׺���ַ�����������������洢·��
%
% ����ļ�:
% [str1]_theta-min_Scatter_rather.txt - ���Ŧ�ֵ��¼�ļ���������
% ��СMisFitֵ(m) ��Ѧ�ֵ
%
% ͼ�����:
% 1) ��-MisFit��ϵ����ͼ
% 2) ��Ѧ�ֵ��Ӧ��Chi-�̷ֲ߳�ͼ
% 3) ��-0.2ʱ��Chi-�̶߳Ա�ͼ
% 4) ��+0.2ʱ��Chi-�̶߳Ա�ͼ
%
% �ؼ���������:
% 1. ���ݶ�ȡ����������ָ����ŵĺ��������ļ����ļ��������
% ���ξ���(m) | �߳�(m) | ��ˮ���(m�0�5) �ȹؼ�����
% 2. Chi���㣺���ݹ�ʽ �� = ��(A0/A)^�� dx ���ɲ�ͬ��ֵ��Chi����
% 3. ����ʶ��ѡ��Chi��ĺ�����Ϊ��׼����
% 4. �����㣺ͨ�����Բ�ֵ�Ƚϸ�֧�������ɵĸ߳�ƫ�������������(MisFit)
% 5. �����Ż���Ѱ��ʹ��MisFit��С�Ħ�ֵ��Ϊ���źӵ�����
%
% ע������:
% 1. �����ļ�Ҫ��
% - �ı���ʽ��ÿ�а�����X���� Y���� ���ξ��� �߳� ���� ��ˮ���
% - �ļ�������������ţ���1.txt, 2.txt,...n.txt
% 2. �������ã�
% - ��������Χ�̶�Ϊ0.1-1.0������0.05�����޸����������
% - �̱߳�׼��������Сֵ���㴦��
% 3. ͼ�����ԣ�
% - ��ͼ1��ʾ��-MisFit���ߣ��Զ��궨y�᷶Χ
% - ��ͼ2-4ʹ�û�ɫ���Ƹ�����Chi-�߳�����
% 4. �쳣����
% - �Զ����������ڻ��ʽ������ļ�
% - �ȳ�����Χʱ�Զ�ȡ�߽�ֵ
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ������д��: ��һ�� - ��������: 2025��2��4��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ����
format long; close all;
theta = 0.1:0.05:1;
len_theta = length(theta);

%% ��ȡ���к����������㲻ͬ�ȶ�Ӧ��chi-z
Path3 = '.txt';
sumRiver = endRiver - startRiver + 1; %��������
AllRiverPro = cell(sumRiver,1); % ����cell�洢���еĺ������������
lengthChi = [];

% ��ȡÿ������������
for i = startRiver:endRiver
    Path2 = num2str(i);
    Path_i = strcat(Path1, Path2, Path3);
    RiverPro_i = textread(Path_i); % ��ȡ����
    Ele_i = RiverPro_i(:,4); 
    Area_i = RiverPro_i(:,6); 
    upLen_i = RiverPro_i(:,3);

    chiLen_i = cal_chiLen_i(Area_i, upLen_i, theta); % ����chiLen����
    
    minEle_i = min(Ele_i); 
    Ele_i = Ele_i - minEle_i; % ��׼���߶�
    AllRiverPro{i - startRiver + 1, 1} = [chiLen_i, Ele_i]; % �洢���к�������
end

%% ����ÿ��theta��MisFitֵ��������ͼ��
MisFit_i = zeros(1, len_theta);
figure;

for i = 1:len_theta
    RiverPro_Chosen = cell(sumRiver, 1);
    for j = 1:sumRiver
        ChiZ_ij = AllRiverPro{j, 1};
        RiverPro_Chosen{j, 1} = [ChiZ_ij(:, i), ChiZ_ij(:, len_theta + 1)]; % ��ȡ��ǰ���µ�Chi-z
    end
    MisFit_i(i) = co_stem(RiverPro_Chosen); % ����MisFit
end

% ���Ʀ���MisFit�Ĺ�ϵͼ
subplot(2, 2, 1);
plot(theta, MisFit_i); 
xlabel('theta');
ylabel('Elevation scatter (m)');
[min_MisFit_i, num_MisFit_i] = min(MisFit_i);
min_theta = theta(num_MisFit_i);
title(strcat('Min MisFit: ', num2str(min_MisFit_i), ', theta = ', num2str(min_theta)));
max_MisFit_i = max(MisFit_i);
ylim_down = floor(min_MisFit_i / 10) * 10;
ylim_up = ceil(max_MisFit_i / 10) * 10;
ylim([ylim_down, ylim_up]);

%% ������Ѱ��������е�chi��߳�ͼ
subplot(2, 2, 2);
for j = 1:sumRiver
    RiverPro_ij = AllRiverPro{j, 1};
    plot(RiverPro_ij(:, num_MisFit_i), RiverPro_ij(:, len_theta + 1),'Color', [0.5 0.5 0.5]);
    hold on;
end
xlabel('Chi Length');
ylabel('Elevation');
title(strcat('Chi-Z for best theta  = ', num2str(min_theta)));

%% ������Ѱ��� - 0.2 ʱ��chi��߳�ͼ
theta_minus_02 = min_theta - 0.2;
subplot(2, 2, 3);
for j = 1:sumRiver
    RiverPro_ij = AllRiverPro{j, 1};
    [~, idx] = min(abs(theta - theta_minus_02)); % �ҵ���ӽ�min_theta-0.2��theta����
    plot(RiverPro_ij(:, idx), RiverPro_ij(:, len_theta + 1), 'Color', [0.5 0.5 0.5]);
    hold on;
end
xlabel('Chi Length');
ylabel('Elevation');
title(strcat('Chi-Z for theta  = ', num2str(theta_minus_02)));

%% ������Ѱ��� + 0.2 ʱ��chi��߳�ͼ
theta_plus_02 = min_theta + 0.2;
subplot(2, 2, 4);
for j = 1:sumRiver
    RiverPro_ij = AllRiverPro{j, 1};
    [~, idx] = min(abs(theta - theta_plus_02)); % �ҵ���ӽ�min_theta+0.2��theta����
    plot(RiverPro_ij(:, idx), RiverPro_ij(:, len_theta + 1), 'Color', [0.5 0.5 0.5]);
    hold on;
end
xlabel('Chi Length');
ylabel('Elevation');
title(strcat('Chi-Z for theta  = ', num2str(theta_plus_02)+0.2));

%% ������С���ֵ
str2 = 'theta-min_Scatter_rather.txt';
str = strcat(str1, str2);
fid = fopen(str, 'a+');
fprintf(fid, '%f\t', min_MisFit_i);       
fprintf(fid, '%f\n', min_theta);
fclose(fid);

end

%%  ���㲻ͬ��ֵ��Ӧ��Chi
function chiLen_i=cal_chiLen_i(Area_new,upLen_new,concavity)

chiLen_i=[];

len_con=length(concavity);      mLen_new=length(Area_new);
for i=1:len_con
    Chi=[]; sum_i=0; Chi=[Chi;sum_i];
    for j=2:mLen_new
        sum_i=sum_i+(1/Area_new(j-1)).^(concavity(i))*(upLen_new(j)-upLen_new(j-1));
        Chi=[Chi;sum_i]; % Chi��������
    end
    chiLen_i=[chiLen_i,Chi];
end
end
%% ����Chi-z��misFit,��ӽ����ɵ�
function MisFit_i=co_stem(RiverPro_Chosen)
%RiverPro_Chosen��cell��ÿ��Ԫ�ض���chi-z

sumRiver=length(RiverPro_Chosen);

% Ѱ��stem����Chi��ĺӶ�
max_Chi=0;
for i=1:sumRiver
    RiverPro_i=RiverPro_Chosen{i,1};    Chi_i=RiverPro_i(:,1);      max_Chi_i=max(Chi_i);
    if max_Chi_i>max_Chi
        max_Chi=max_Chi_i;      num_max_Chi=i;
    end
end
RiverPro_Stem=RiverPro_Chosen{num_max_Chi,1};    Stem_Chi=RiverPro_Stem(:,1);   Stem_z=RiverPro_Stem(:,2);

% ����MisFit_i
num_node=0;     MisFit_i=0;
for j=1:sumRiver
    if j ~=  num_max_Chi
        RiverProj=RiverPro_Chosen{j,1};    Chi_j=RiverProj(:,1);    z_j=RiverProj(:,2); % ����������
        z_j_recal=interp1(Stem_Chi,Stem_z,Chi_j);
        MisFit_i=MisFit_i+sum((z_j-z_j_recal).^2);      num_node=num_node+length(Chi_j);
    end
end
MisFit_i=sqrt(MisFit_i/num_node);
end

%%


