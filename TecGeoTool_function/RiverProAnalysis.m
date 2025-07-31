function RiverProAnalysis(river_i,str1,outpath)
% ����:
% �ú������ں��������������֧�ֽ���ʽ�ѵ�ѡ�񣬼���ӵ����ȡ�ksnֵ����̬������
% ���ɰ����ֶ�������ʸ���ļ���ͳ�Ʒ������������ͨ��ͼ�ν��������û�ѡ���ˮ��
% ����ֵ���ѵ�λ�ã�֧����logS-Aͼ��chi-zͼ����ģʽ�½����ѵ�ʶ��
%
% ��Ҫ����:
% river_i - ������ţ���ֵ�ͣ������ڽ���ļ���ʶ
% str1 - ���������ļ�·�����ַ�������ӦΪ���������е��ı��ļ���
% [Y����, X����, ��Դ����(m), �߳�(m), ����, ��ˮ���(m�0�5), Chiֵ]
% outpath - ����ļ�����·�����ַ����������������Ŀ¼·�����ļ���ǰ׺
%
% ����ļ�:
% 1. [outpath]_Sfeature[river_i].shp - ʸ�����ļ����������Ӷ���̬������
% (ID, RiverNo, theta, std_theta, ksn, std_ksn, DWksn, DWstd_ksn��)
% 2. [outpath]_Pnt_divide.txt - �ֽ�������ı��ļ���������
% [�������, X����, Y����, ����, ���, ksn, ���, У��ksn, ���,
% ��ˮ�����Χ(km�0�5), �̼߳�ֵ(m), �����Դ��(m), Chiֵ]
% 3. [outpath]_Pnt_analysis.txt - �Ӷη�������ı��ļ���������
% [�������, ����, ��Դ��, ����, log_ks]
% 4. [str1]New.txt - ��ǿ�������������ļ���������
% [ƽ���߳�, �ֲ�ksn, �¶�]���ֶ�
%
% ��������˵��:
% 1. �������к󽫵�����Ϸ���ͼ����������
% - ��ˮ���-��������ͼ
% - �¶�-��������ͼ
% - ���-�¶�˫����ͼ
% - �߳�-��������ͼ
% - Chi-�̹߳�ϵͼ
% 2. ������ʾ��ͼ�ν���������²�����
% a) �ڻ�ˮ���ͼ�е��ѡ�������ֵ��
% b) ������Ҫ�޳��ĸ̼߳�ֵ����ѡ��
% c) ѡ���ѵ�ʶ��ģʽ��logS-Aͼ��chi-zͼ��
% d) ��ѡ��ͼ�е��ȷ���ѵ�λ��
%
% ע������:
% 1. �����ı��ļ���ȷ����˳�����ʽ��ȷ
% 2. ʸ���ļ�������MATLAB Mapping Toolbox֧��
% 3. �����������ݵ����50m�Ա�֤���㾫��
% 4. ͼ�ν������ʱ�밴��ʾ˳�����������������л�
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ������д��: ��һ�� - ��������: 2025��2��4��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ��ͼʱ���������ʾԼ��
%minLim_Area=10^1;maxLim_Area=10^9;
% minLim_Len=0;maxLim_Len=5e4;
%minLim_Slope=10^(-3);maxLim_Slope=10^0;
%MaxChi=25; MinEle=1100; MaxEle=MinEle+1000; %(floor(max(Ele)/100)+1)*100; ������Ը߲����5000����������ʹ�ԼΪ500

% ******
river_i=double(river_i);
str2='.txt'; %str_i=strcat(str1,str2);

riverPro_old=textread(str1);

% ����logS A���µ�riverPro ����ԭ�ȵģ�������ƽ���̡߳��ֲ�ksn���¶�
SmoothLen=50; ResampleWin=40; %��λ m��ƽ�����룬�߳��ز�������
riverPro=CalLogSA(riverPro_old,SmoothLen,ResampleWin); riverPro_old=[];
[m,n]=size(riverPro);

% �洢�µĺ�������: ��������ƽ������Y��X����Դ����(3)���߳�(4)�����򡢻�ˮ���(6)��Chi���롢�ֶ�ksn��ƽ���̡߳��ֲ�ksn���¶�(10)
str_iNew=strcat(str1,'New',str2); fid1=fopen(str_iNew,'w');
for i=1:m
    for j=1:n-1; fprintf(fid1,'%f\t',riverPro(i,j)); end
    fprintf(fid1,'%f',riverPro(i,n));fprintf(fid1,'\n');
end
fclose(fid1);
%%

% ����һ����ͼ����
h_combined = figure; 
set(h_combined, 'units', 'normalized', 'position', [0.05 0.05 0.85 0.85]);

% ��һ��ͼ (x-log(Area))
h1=subplot(3, 2, 1); 
set(gca, 'position', [0.1 0.71 0.35 0.25]);
semilogy(riverPro(:, 3), riverPro(:, 6));
xlabel('Length (m)'); ylabel('Area (m^2)'); hold on;

% �ڶ���ͼ (x-log(Slope))
h2=subplot(3, 2, 3); 
set(gca, 'position', [0.1 0.39 0.35 0.25]);
semilogy(riverPro(:, 3), riverPro(:, 10), '+');
xlabel('Length (m)'); ylabel('Slope'); hold on;

% ������ͼ (log(Area-Slope))
h3=subplot(3, 2, 5); 
set(gca, 'position', [0.1 0.07 0.35 0.25]);
loglog(riverPro(:, 6), riverPro(:, 10), '+');hold on;
xlabel('Area (m^2)'); ylabel('Slope'); 

% ���ĸ�ͼ (x-Elevation)
h4=subplot(3, 2, 2); 
set(gca, 'position', [0.55 0.57 0.35 0.4]); % ȷ��λ�ò�������
plot(riverPro(:, 3), riverPro(:, 4)); hold on;% ���������ĺӵ�����
xlabel('Length (m)'); ylabel('Elevation (m)');


% �����Ի�����ʾ�û�
choice = questdlg('�뿪ʼѡ���ˮ�����ֵ', 'ѡ����ʾ', ...
    'ȷ��', 'ȡ��', 'ȷ��');

% �����û�ѡ��ִ�в���
if strcmp(choice, 'ȷ��')
    % ѡ���ˮ�����ֵ
    [Ac, ~] = ginput(1); 

    % ���������߼�
    i = m;
    while riverPro(i, 6) < Ac 
        i = i - 1; 
    end
else
    disp('������ȡ��');
end
% �޳���ˮ���С��Ac����
riverPro = riverPro(1:i, :); m = size(riverPro, 1);            % �µ�������

% ���ĸ�ͼ (x-Elevation)
subplot(h4);
plot(riverPro(i, 3), riverPro(i, 4), '*');hold on;% ��Ǻӵ�����
xlabel('Length (m)'); ylabel('Elevation (m)');

% �ڶ��㴦��ǻ�ˮ�����ֵ Acr
Acr_km2 = riverPro(i, 6) / 1e6; % ת��Ϊ km�0�5
labelText = strcat('A_{cr} = ', num2str(Acr_km2, '%.2f'), ' km^2'); % ��ʽ���ַ���

% �ںӵ�����λ����ӱ�ע
text(riverPro(i, 3), riverPro(i, 4), labelText, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', ...
    'FontSize', 10, 'Color', 'red');

% �����ͼ (Chi-Elevation)
h5=subplot(3, 2, 6); 
set(gca, 'position', [0.55 0.05 0.35 0.4]);
plot(riverPro(:, 7), riverPro(:, 4), 'g', 'linewidth', 3);
xlabel('Chi'); ylabel('Elevation (m)');
hold on;

%% Fig3������DWУ�����chi-z
%h3=figure; 
% ���ݸ߳��޳�С�ĸ߳�
minorEle = 0;
while true
    prompt = {'�������޳���С�ĸ߳�ֵ (Ĭ��ֵ-1�������0)��'};
    dlgtitle = '�޳�С�ĸ߳�';
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
    
    % ���¸߳�����
    Chi = riverPro(:, 7);  Ele = riverPro(:, 4);   m = length(Chi);
end

% ���ݸ߳��޳���ĸ߳�
largeEle = 0;
while true
    prompt = {'�������޳��Ĵ�ĸ߳�ֵ (Ĭ��ֵ-1�������0)��'};
    dlgtitle = '�޳���ĸ߳�';
    dims = [1 50];
    definput = {'-1'};
    answer = inputdlg(prompt, dlgtitle, dims, definput);
    
    if isempty(answer)
        % ����û�ȡ������Ի����˳�ѭ��
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
    
    % ���¸߳�����
    Chi = riverPro(:, 7);    Ele = riverPro(:, 4);    m = length(Chi);
end


% �޳����õĸ̺߳��µĺ��������� % set(gca,'ytick',MinEle:100:MaxEle);% text(Chi(1:2:m),Ele(1:2:m),num2str(Ele(1:2:m)));
Chi=riverPro(:,7);Ele=riverPro(:,4);Area=riverPro(:,6);Slope=riverPro(:,10); m=length(Chi);
%figure(h2); subplot(h22); plot(Chi,Ele,'r','linewidth',3);xlabel('Chi (m)');ylabel('Elevation (m)');hold on;
% �����ͼ (Chi-Elevation)
subplot(h5); 
plot(Chi, Ele, 'r', 'linewidth', 3);
xlabel('Chi'); ylabel('Elevation (m)');hold on;
hold on;
%% ��������
PntKs = [];
Pnt_analysis = [];

% 1. ����ָ������
prompt = {'�������ѵ����� ��'};
dlgtitle = '�����ѵ�����';
dims = [1 50];
definput = {'0'};
answer = inputdlg(prompt, dlgtitle, dims, definput);

if isempty(answer)
    disp('����ȡ��');
    return;
end

PntSum = str2double(answer{1});

%% 2. ѡ����Ѱ�ָ��ķ�����1����logS-A��2����chi-zͼ
if PntSum>0
choicePrompt = {'��ѡ����Ѱ�ѵ�ķ�����1����logS-A��2����chi-zͼ'};
choiceDlgTitle = 'ѡ�񷽷�';
choiceDefInput = {'1'};
choiceAnswer = inputdlg(choicePrompt, choiceDlgTitle, dims, choiceDefInput);

if isempty(choiceAnswer)
    %disp('����ȡ��');
    return;
end

choose_fig = str2double(choiceAnswer{1});
if isnan(choose_fig) || ~ismember(choose_fig, [1, 2])
    %disp('��Ч��ѡ�񣬳�����ֹ��');
    return;
end

if choose_fig == 1
    % 3. �� logS-A ͼѡ��ָ��
    %disp(['�� logS-A ͼѡ�� ', num2str(PntSum), ' ���ָ�� (����ˮ���)��']);
    uiwait(msgbox('��� logS-A ͼ��ѡ���ѵ㡣', '��ʾ', 'modal'));
    figure(h_combined);  subplot(h3);  hold on; 
    [PntArea, ~] = ginput(PntSum);
    PntArea = sort(PntArea, 'descend');

    PntNum = [1]; % ��ʼ���ֽ����ţ��������
    for di = 1:PntSum
        for k = 1:m-1
            if (PntArea(di) <= Area(k) && PntArea(di) >= Area(k+1))
                PntNum = [PntNum; k]; 
                %disp([num2str(Area(k) / 1e6), ' km�0�5�� ', num2str(Ele(k)), ' m']);
                break;
            end
        end
    end
    PntNum = [PntNum; m]; % �����յ�
    PntSum = length(PntNum) - 2; % ʵ�ʷָ����

elseif choose_fig == 2
    % 4. �� chi-z ͼѡ��ָ��
    %disp(['�� chi-z ͼѡ�� ', num2str(PntSum), ' ���ѵ� (�� �� ֵ)��']);
    uiwait(msgbox('��� chi-z ͼ��ѡ���ѵ㡣', '��ʾ', 'modal'));
    figure(h_combined);subplot(h5);  hold on; % ����ԭ��ͼ��
    [PntChi, ~] = ginput(PntSum);
    PntChi = sort(PntChi);

    PntNum = [1]; % ��ʼ���ֽ����ţ��������
    for j = 1:PntSum
        for k = 1:m-1
            if (PntChi(j) >= Chi(k) && PntChi(j) <= Chi(k+1))
                PntNum = [PntNum; k];
                %disp([num2str(Area(k) / 1e6), ' km�0�5�� ', num2str(Ele(k)), ' m']);
                break;
            end
        end
    end
    PntNum = [PntNum; m]; % �����յ�
    PntSum = length(PntNum) - 2; % ʵ�ʷָ����
end
end
%disp(['���գ��Ӵ�С��ʵ�ʷֽ�� ', num2str(PntSum), ' ����']);

%%
if PntSum==0
    [R2,z0,ksn,std_ksn]=Bi_Regress(Ele,Chi); % ���Իع� Y=a0+a1*X
    figure(h_combined); subplot(h5);plot(Chi,Chi.*ksn+z0,'k');hold on;
    text(max(Chi)/2,max(Ele),['R=',num2str(R2),'z=Chi.*',num2str(ksn),'��',num2str(std_ksn),'+',num2str(z0)]); % �������Fig2���滭��
    [DW_r,DW_Ele,DW_Chi]=DW_test(Ele,Chi,z0,ksn); %DW����
    [R2,DWz0,DWksn,DWstd_ksn]=Bi_Regress(DW_Ele,DW_Chi); 
    figure(2); plot(DW_Chi,DW_Ele,'bo'); hold on; plot(DW_Chi,DW_Chi.*DWksn+DWz0,'k'); text(max(DW_Chi)/2,max(DW_Ele),['ro=',num2str(DW_r),'R=',num2str(R2),'ksn=',num2str(DWksn),'��',num2str(DWstd_ksn)]);
    [R2,log_ks,theta,std_theta]=Bi_Regress(log(Slope),log(Area));
    %[theta2,ksn2]=Nonlinear_Regress(Area,Slope);
    figure(h_combined); subplot(h3); loglog(Area,(Area.^(theta)).*(exp(log_ks)),'k');hold on; text(max(Area)/2,0.1,['R=',num2str(R2),'��=',num2str(theta),'��',num2str(std_theta)]);
    %% ����ֽ����Ϣ��������ţ�X��Y���꣬�Ӷΰ��ȡ���ksn����У����ksn�����۵��ˮ�����Km2��min,max�����۵�߳�(max,min)������۵����Դ���롢chiֵ
    PntKs=[PntKs;[river_i,riverPro(m,2),riverPro(m,1),theta,std_theta,ksn,std_ksn,DWksn,DWstd_ksn,(riverPro(m,6))./(1e6),(riverPro(1,6))./(1e6),riverPro(m,4),riverPro(1,4),riverPro(m,3),riverPro(m,7)]];
    Pnt_analysis=[Pnt_analysis;[river_i,riverPro(m,1),riverPro(m,2),riverPro(m,3),theta,log_ks]];    
    %disp(river_i);disp(riverPro(m,2));disp(riverPro(m,1));disp(theta);disp(std_theta);
    %% ����ʸ����,�µ�riverPro: ��������ƽ������Y��X����Դ���롢�̡߳����򡢻�ˮ�����Chi���롢ƽ���̡߳��ֲ�ksn���¶�
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
        
else % >=1���ֽ�� 
    %% ��� PntKs    % ��һ��
    Ele_1=Ele(PntNum(1):PntNum(2)); Chi_1=Chi(PntNum(1):PntNum(2)); Area_1=Area(PntNum(1):PntNum(2)); Slope_1=Slope(PntNum(1):PntNum(2));
    figure(h_combined); subplot(h4); plot(riverPro(PntNum(2),3),riverPro(PntNum(2),4),'kx');hold on; %��UpLen-z�ϱ�Ǻӵ��ѵ� 
    [R2,z0_1,ksn_1,std_ksn_1]=Bi_Regress(Ele_1,Chi_1); % ���Իع� Y=a0+a1*X
    %[theta2_1,ksn2_1]=Nonlinear_Regress(Area_1,Slope_1);
    figure(h_combined); subplot(h5); plot(Chi_1,Chi_1.*ksn_1+z0_1,'k');hold on; text(max(Chi_1)/2,max(Ele_1),['R1=',num2str(R2),'z1=Chi.*',num2str(ksn_1),'��',num2str(std_ksn_1),'+',num2str(z0_1)]);% �������Fig2���滭��
    [DW_r_1,DW_Ele_1,DW_Chi_1]=DW_test(Ele_1,Chi_1,z0_1,ksn_1); %DW����
    [R2,DWz0_1,DWksn_1,DWstd_ksn_1]=Bi_Regress(DW_Ele_1,DW_Chi_1); 
    figure(2); plot(DW_Chi_1,DW_Ele_1,'bo'); hold on; plot(DW_Chi_1,DW_Chi_1.*DWksn_1+DWz0_1,'k'); text(max(DW_Chi_1)/2,max(DW_Ele_1),['ro1=',num2str(DW_r_1),'R=',num2str(R2),'ksn=',num2str(DWksn_1),'��',num2str(DWstd_ksn_1)]);
    [R2,log_ks_1,theta_1,std_theta_1]=Bi_Regress(log(Slope_1),log(Area_1));
    figure(h_combined); subplot(h3); loglog(Area_1,(Area_1.^(theta_1)).*(exp(log_ks_1)),'k');hold on; text(max(Area_1)/2,max(Slope_1)/2,['R1=',num2str(R2),'��=',num2str(theta_1),'��',num2str(std_theta_1)]);
    %% ����ֽ����Ϣ��������ţ�X��Y���꣬�۵����ºӶΰ��ȡ���ksn����У����ksn�����۵��ˮ�����Km2��min,max�����۵�߳�(max,min)������۵����Դ���롢chiֵ
    PntKs=[PntKs;[river_i,riverPro(PntNum(2),2),riverPro(PntNum(2),1),theta_1,std_theta_1,ksn_1,std_ksn_1,DWksn_1,DWstd_ksn_1,(riverPro(PntNum(2),6))./(1e6),(riverPro(PntNum(1),6))./(1e6),...
        riverPro(PntNum(2),4),riverPro(PntNum(1),4),riverPro(PntNum(2),3),riverPro(PntNum(2),7)]];
     Pnt_analysis=[Pnt_analysis;[river_i,riverPro(PntNum(2),1),riverPro(PntNum(2),2),riverPro(PntNum(2),3),theta_1,log_ks_1]];  
    %% ����ʸ����,�µ�riverPro: ��������ƽ������Y��X����Դ���롢�̡߳����򡢻�ˮ�����Chi���롢ƽ���̡߳��ֲ�ksn���¶�
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
        
    %% ************ʣ��ĺӶ�*************************************************
    for j=2:(PntSum+1)
        Ele_j=Ele((PntNum(j)+1):PntNum(j+1)); Chi_j=Chi((PntNum(j)+1):PntNum(j+1)); Area_j=Area((PntNum(j)+1):PntNum(j+1)); Slope_j=Slope((PntNum(j)+1):PntNum(j+1));
        figure(h_combined); subplot(h4); plot(riverPro(PntNum(j+1),3),riverPro(PntNum(j+1),4),'kx');hold on; %��UpLen-z�ϱ�Ǻӵ��ѵ� 
        [R2,z0_j,ksn_j,std_ksn_j]=Bi_Regress(Ele_j,Chi_j); % ���Իع� Y=a0+a1*X
        figure(h_combined); subplot(h5); plot(Chi_j,Chi_j.*ksn_j+z0_j,'k');hold on; text(max(Chi_j)/2,max(Ele_j),['R',num2str(j),'=',num2str(R2),'z=Chi.*',num2str(ksn_j),'��',num2str(std_ksn_j),'+',num2str(z0_j)]);% �������Fig2���滭��       
        [DW_r_j,DW_Ele_j,DW_Chi_j]=DW_test(Ele_j,Chi_j,z0_j,ksn_j); %DW����
        [R2,DWz0_j,DWksn_j,DWstd_ksn_j]=Bi_Regress(DW_Ele_j,DW_Chi_j);
        %[theta2_j,ksn2_j]=Nonlinear_Regress(Area_j,Slope_j);
        figure(2); plot(DW_Chi_j,DW_Ele_j,'bo'); hold on; plot(DW_Chi_j,DW_Chi_j.*DWksn_j+DWz0_j,'k'); text(max(DW_Chi_j)/2,max(DW_Ele_j),['ro',num2str(j),'=',num2str(DW_r_j),'R=',num2str(R2),'ksn=',num2str(DWksn_j),'��',num2str(DWstd_ksn_j)]);
 %      figure (3); plot(Chi_j,Chi_j.*DWksn_j+DWz0_j/(1-DW_r_j),'k');hold on; text(max(Chi_j)/2,max(Ele_j)/2,['zj=Chi.*',num2str(DWksn_j),'��',num2str(DWstd_ksn_j),'+',num2str(DWz0_j/(1-DW_r_j))]); % �������Fig3���滭��
        [R2,log_ks_j,theta_j,std_theta_j]=Bi_Regress(log(Slope_j),log(Area_j));
        figure(h_combined); subplot(h3); loglog(Area_j,(Area_j.^(theta_j)).*(exp(log_ks_j)),'k');hold on; text(max(Area_j)/2,max(Slope_j)/2,['R',num2str(j),'=',num2str(R2),'��=',num2str(theta_j),'��',num2str(std_theta_j)]);
        %% ����ֽ����Ϣ��������ţ�X��Y���꣬�۵����ºӶΰ��ȡ���ksn����У����ksn�����۵��ˮ�����Km2��min,max�����۵�߳�(max,min)������۵����Դ���롢chiֵ
        PntKs=[PntKs;[river_i,riverPro(PntNum(j+1),2),riverPro(PntNum(j+1),1),theta_j,std_theta_j,ksn_j,std_ksn_j,DWksn_j,DWstd_ksn_j,(riverPro(PntNum(j+1),6))./(1e6),(riverPro(PntNum(j),6))./(1e6),...
            riverPro(PntNum(j+1),4),riverPro(PntNum(j),4),riverPro(PntNum(j+1),3),riverPro(PntNum(j+1),7)]];
        Pnt_analysis=[Pnt_analysis;[river_i,riverPro(PntNum(j+1),1),riverPro(PntNum(j+1),2),riverPro(PntNum(j+1),3),theta_j,log_ks_j]];  
       %% ����ʸ����,�µ�riverPro: ��������ƽ������Y��X����Դ���롢�̡߳����򡢻�ˮ�����Chi���롢ƽ���̡߳��ֲ�ksn���¶�
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
removeLines(str_PntKs, river_i);%ɾ��֮ǰ�ļ�¼
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
removeLines(str_analysis, river_i);%ɾ��֮ǰ�ļ�¼
for i=1:analysis_m
    for j=1:analysis_n-1
        fprintf(fid,'%f\t',Pnt_analysis(i,j));
        %disp(Pnt_analysis(i,j))
    end
    fprintf(fid,'%f',Pnt_analysis(i,analysis_n));fprintf(fid,'\n');
end
fclose(fid);
end
%% �����¶�
function riverPro_f=CalLogSA(riverPro_Old_f,SmoothLen_f,ResampleWin_f)
% riverPro_i����������ƽ������Y��X����Դ���롢�̡߳����򡢻�ˮ�����Chi���롢����ksn
% �µĵ���riverPro ����ԭ�ȵģ�������ƽ���̡߳��ֲ�ksn���¶�
UpLen=riverPro_Old_f(:,3);Ele=riverPro_Old_f(:,4);Chi=riverPro_Old_f(:,7); m=length(UpLen); 

%����ƽ����ĸ߳�
SmoothedEle=Ele;Localksn=zeros(m,1);Slope=zeros(m,1);
for i=1:m
    if abs(UpLen(i)-UpLen(1))<=SmoothLen_f/2 % ����i���ˮ�ھ���С��SmoothLen_f/2����һֱ���Ѱ��ֱ���ﵽ��ֵΪֹ
        j=1; while i+j<=m && abs(UpLen(i+j)-UpLen(1))<=SmoothLen_f; j=j+1;end; j=j-1; 
        if (i+j)-1<2; j=3-i; end %�������ֵ�ع������3��ǿ�ƶ�ѡ��㣬ȷ����ֵ���ع����������㡣
        SmoothedEle(i)=interp1(UpLen(1:i+j),Ele(1:i+j),UpLen(i));
        b=regress(Ele(1:i+j),[ones(size(Chi(1:i+j))) Chi(1:i+j)]); %regress�����÷���y=a0+a1*x
        Localksn(i)=abs(b(2));
    elseif abs(UpLen(m)-UpLen(i))<=SmoothLen_f/2 % ����i��ˮͷ����С��SmoothLen_f/2����һֱ��ǰѰ��ֱ���ﵽ��ֵΪֹ
        j=1; while i-j>=1 && abs(UpLen(m)-UpLen(i-j))<=SmoothLen_f; j=j+1;end; j=j-1;
        if (i-j)>(m-2); j=i-(m-2); end %�������ֵ�ع������3��ǿ�ƶ�ѡ��㣬ȷ����ֵ���ع����������㡣
        SmoothedEle(i)=interp1(UpLen(i-j:m),Ele(i-j:m),UpLen(i));
        b=regress(Ele(i-j:m),[ones(size(Chi(i-j:m))) Chi(i-j:m)]); 
        Localksn(i)=abs(b(2));
    elseif abs(UpLen(i)-UpLen(1))>SmoothLen_f/2 && abs(UpLen(m)-UpLen(i))>SmoothLen_f/2 %�����ˮ�ڡ�ˮͷ���붼����SmoothLen_f/2
        j=1; while i-j>=1 && i+j<=m && abs(UpLen(i+j)-UpLen(i-j))<=SmoothLen_f; j=j+1; end; j=j-1;
        if j==0; j=1; end %����������������ڵ�ľ��룼SmoothLen_f��ǿ��ѡ�����ڵĵ㣬ȷ����ֵ���ع����������㡣
        SmoothedEle(i)=interp1(UpLen(i-j:i+j),Ele(i-j:i+j),UpLen(i));
        b=regress(Ele(i-j:i+j),[ones(size(Chi(i-j:i+j))) Chi(i-j:i+j)]);
        Localksn(i)=abs(b(2));
    end
end

%�����ز������Localksn��Slope
for i=1:m
    if abs(Ele(i)-Ele(1))<=ResampleWin_f/2 % ����i���ˮ�ڸ߳�С��ResampleWin_f/2����һֱ���Ѱ��ֱ���ﵽ��ֵΪֹ
        j=1; while i+j<=m && abs(Ele(i+j)-Ele(1))<=ResampleWin_f; j=j+1;end; j=j-1; 
        if (i+j)-1<2; j=3-i; end %�������ֵ�ع������3��ǿ�ƶ�ѡ��㣬ȷ����ֵ���ع����������㡣
        b=regress(Ele(1:i+j),[ones(size(UpLen(1:i+j))) UpLen(1:i+j)]); %regress�����÷���y=a0+a1*x
        Slope(i)=abs(b(2));
    elseif abs(Ele(m)-Ele(i))<=ResampleWin_f/2 % ����i��ˮͷ�߳�С��ResampleWin_f/2����һֱ��ǰѰ��ֱ���ﵽ��ֵΪֹ
        j=1; while i-j>=1 && abs(Ele(m)-Ele(i-j))<=ResampleWin_f; j=j+1;end; j=j-1;
        if (i-j)>(m-2); j=i-(m-2); end %�������ֵ�ع������3��ǿ�ƶ�ѡ��㣬ȷ����ֵ���ع����������㡣
        b=regress(Ele(i-j:m),[ones(size(UpLen(i-j:m))) UpLen(i-j:m)]); 
        Slope(i)=abs(b(2));
    elseif abs(Ele(i)-Ele(1))>ResampleWin_f/2 && abs(Ele(m)-Ele(i))>ResampleWin_f/2 %�����ˮ�ڡ�ˮͷ�̶߳�����ResampleWin_f
        j=1; while i-j>=1 && i+j<=m && abs(Ele(i+j)-Ele(i-j))<=ResampleWin_f; j=j+1; end; j=j-1;
        if j==0; j=1; end %����������������ڵ�ĸ߲ResampleWin_f��ǿ��ѡ�����ڵĵ㣬ȷ����ֵ���ع����������㡣
        b=regress(Ele(i-j:i+j),[ones(size(UpLen(i-j:i+j))) UpLen(i-j:i+j)]);
        Slope(i)=abs(b(2));
       
    end
end

riverPro_f=[riverPro_Old_f SmoothedEle Localksn Slope];
end
%% ͳ�Ƽ���
function [DW_r_f,DW_Ele_f,DW_Chi_f]=DW_test(Ele_f,Chi_f,z0_f,ksn_f)
DW_e=Ele_f-(z0_f+ksn_f.*Chi_f); % �в�
m=length(DW_e);
DW_e1=DW_e(1:m-1); DW_e2=DW_e(2:m);
DW_r_f=sum(DW_e1.*DW_e2)/(sum(DW_e1.^2));
DW_Chi_f=Chi_f(2:m)-DW_r_f.*Chi_f(1:m-1);
DW_Ele_f=Ele_f(2:m)-DW_r_f.*Ele_f(1:m-1);
end
%% ���Իع�
function [R2,a0,a1,std_a1]=Bi_Regress(Y,X)
% ���Իع� Y=a0+a1*X
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
    %�ع����ݲ���������ģ��
    R2=9999;a0=9999;a1=9999;std_a1=9999;    
end
end

%% ɾ��֮ǰ������¼��������ѡ����ɵ������ݳ�����
function removeLines(file_path, starting_number)
    % ��ȡ�ļ�����
    fileID = fopen(file_path, 'r');
    lines = textscan(fileID, '%s', 'Delimiter', '\n');
    fclose(fileID);

    % �ҵ���ָ������ͷ���в�ɾ��
    lines = lines{1};
    modified_lines = lines(~startsWith(lines, num2str(starting_number)));

    % ���޸ĺ������д���ļ�
    fileID = fopen(file_path, 'w');
    fprintf(fileID, '%s\n', modified_lines{:});
    fclose(fileID);
end
