function head = Search_WaterHead(FlowAcc,FlowDir,LeftX,DownY,cellsize,ouletX,ouletY,Ac,str1)
% ʹ�÷���:
% head = Search_WaterHead(FlowAcc, FlowDir, LeftX, DownY, cellsize, ouletX, ouletY, Ac, str1)
%
% ����:
% �ú���ͨ������׷�ٷ��Զ�ʶ�������ˮ�����ֵ�ĺӵ�Դͷ��ˮͷ��������D8�㷨����������
% ��ϵ��β�������ˮͷ�������ļ��������ں����Զ���ȡ�������ݻ��о���
%
% ��Ҫ����:
% FlowAcc - ��ˮ�������(m��n)����ֵΪ��դ����ˮ���դ��������Ԥ�ȳ���cellsize^2ת��Ϊ��ʵ�����
% FlowDir - D8�������(m��n)��ȡֵ��Χ{1,2,4,8,16,32,64,128}����ѭArcGIS�����׼
% LeftX - դ�����½�X���꣨��ͼ����ϵ��
% DownY - դ�����½�Y���꣨��ͼ����ϵ��
% cellsize - դ��Ԫ����ߴ磨�ף�
% ouletX - Ŀ��Ӷγ��ڵ�X�������飨��ͼ����ϵ��
% ouletY - Ŀ��Ӷγ��ڵ�Y�������飨��ͼ����ϵ��
% Ac - ��ˮ�����ֵ(m�0�5)������ʶ����ʼˮͷ����ٽ�ֵ
% str1 - ����ļ�·��ǰ׺��ʾ����'D:/data/output_'
%
% ������:
% head - ˮͷ���������(N��2)����ʽΪ��[X����, Y����]
% �����ļ� [str1]_Head.txt - ASCII�ı��ļ���ÿ�д洢ˮͷ�����꣨X,Y�����ʽ��
%
% �㷨����:
% 1. ����������ת��Ϊ�������к�
% 2. ��ʼ��ջ�ṹ��������׷�٣�
% - ��ÿ�����ڵ�������ط���D8������������֧��
% - ��̬��¼���ڵ�Ļ���֧��������Ч֧����
% 3. ˮͷ�ж�������
% - ��ĳ�����л���֧·�Ļ�ˮ��� �� Ac ʱ�����Ϊˮͷ
% 4. ����ת����
% - �����к�ת��Ϊ�������꣨��դ�����ĵ�У����
%
% ʾ������:
% % ������������
% cellsize = 30; % DEM�ֱ���30��
% oulet = [120.5, 25.3]; % ��ע�ӿ�����
% Ac = 1e6; % ��ˮ�����ֵ1km�0�5
% [heads] = Search_WaterHeader(FlowAcc, FlowDir, 120000, 2500000,...
% cellsize, oulet(1), oulet(2), Ac, 'output/watershed1');
%
% �ؼ���������:
% 1. Acѡ��
% - �о�������1e6 ~ 1e7 m�0�5��1~10 km�0�5��
% - С���������1e4 ~ 1e5 m�0�5��0.01~0.1 km�0�5��
% 2. ������
% - ȷ��FlowDir�б�������ֵΪ0�����Զ�����
% - ����ǰӦͨ������Ԥ����DEM
%
% ע������:
% 1. ��������ϵͳ��
% - Ҫ��Y�ᱱ��������ʺ�UTM����ϵ��
% - ����ת����ʽ��Y = UpY - row��cellsize
% 2. �������Ҫ��
% - FlowAcc��Ϊ����դ��������ת��Ϊ��ʵ����ۼӽ��
% - ����������ϸ���ѭD8����淶
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ������д��: ��һ�� - ��������: 2025��2��4��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long;
[m,n]=size(FlowAcc);
UpY=DownY+cellsize*m;
Ac=floor(Ac/(cellsize^2)); % Acλ��ˮ����ٽ�ֵ
FlowDir = double(FlowDir);
len_ouletX=length(ouletX);
stepPoint=[]; WaterHead=[];

num_ToPoint=0; % ����õ��դ����
num_slope=0; % ����õ��դ���У�<Ac��դ����Ŀ
% ��������ȣ�˵���õ�Ϊˮͷ

for ilen=1:len_ouletX
    % ��ʼ��������С��к�
    startRow=floor((UpY-ouletY(ilen))/cellsize)+1; 
    startCol=floor((ouletX(ilen)-LeftX)/cellsize)+1; 

    stepPoint=[stepPoint;[startRow,startCol,FlowAcc(startRow,startCol),0,0]];
    empty_stepPoint=1;
    % �������С��кţ��õ�Ļ�ˮ������õ����������С��кţ�0,0��ʾ�޻��룩
    while empty_stepPoint
        startRow=stepPoint(1,1);
        startCol=stepPoint(1,2);
        
        if FlowAcc(startRow,startCol)>=Ac
            % ȷ���õ��Ǻ�����
        if startCol<n && FlowDir(startRow,startCol+1)==2^4 % && FlowAcc(startRow,startCol)>Ac % && FlowAcc(startRow,startCol+1)>Ac
            num_ToPoint=num_ToPoint+1; % ����õ�
            % �õ��Ҳ�(>Ac)����ջ
            if FlowAcc(startRow,startCol+1)>Ac; stepPoint=[stepPoint;[startRow,startCol+1,FlowAcc(startRow,startCol+1),startRow,startCol]]; end
            % �õ��Ҳ�(<Ac)���õ�����Ϊˮͷ�Ŀ���
            if FlowAcc(startRow,startCol+1)<=Ac; num_slope=num_slope+1; end
        end  
        if startRow<m && startCol<n && FlowDir(startRow+1,startCol+1)==2^5 % && FlowAcc(startRow,startCol)>Ac % && FlowAcc(startRow+1,startCol+1)>Ac
             num_ToPoint=num_ToPoint+1; % ����õ�
             if FlowAcc(startRow+1,startCol+1)>Ac;stepPoint=[stepPoint;[startRow+1,startCol+1,FlowAcc(startRow+1,startCol+1),startRow,startCol]];end
             if FlowAcc(startRow+1,startCol+1)<=Ac;num_slope=num_slope+1;end
        end
        if startRow<m && FlowDir(startRow+1,startCol)==2^6 % && FlowAcc(startRow,startCol)>Ac % && FlowAcc(startRow+1,startCol)>Ac
            num_ToPoint=num_ToPoint+1; % ����õ�
            if FlowAcc(startRow+1,startCol)>Ac;stepPoint=[stepPoint;[startRow+1,startCol,FlowAcc(startRow+1,startCol),startRow,startCol]];end
            if FlowAcc(startRow+1,startCol)<=Ac;num_slope=num_slope+1;end
        end
        if startRow<m && startCol>1 && FlowDir(startRow+1,startCol-1)==2^7 % && FlowAcc(startRow,startCol)>Ac  % && FlowAcc(startRow+1,startCol-1)>Ac
            num_ToPoint=num_ToPoint+1; % ����õ�
            if FlowAcc(startRow+1,startCol-1)>Ac;stepPoint=[stepPoint;[startRow+1,startCol-1,FlowAcc(startRow+1,startCol-1),startRow,startCol]];end
            if FlowAcc(startRow+1,startCol-1)<=Ac;num_slope=num_slope+1;end
        end
        if startCol>1 && FlowDir(startRow,startCol-1)==2^0 % && FlowAcc(startRow,startCol)>Ac % && FlowAcc(startRow,startCol-1)>Ac
            num_ToPoint=num_ToPoint+1; % ����õ�
            if FlowAcc(startRow,startCol-1)>Ac;stepPoint=[stepPoint;[startRow,startCol-1,FlowAcc(startRow,startCol-1),startRow,startCol]];end
            if FlowAcc(startRow,startCol-1)<=Ac;num_slope=num_slope+1;end
        end
        if startRow>1 && startCol>1 && FlowDir(startRow-1,startCol-1)==2^1 % && FlowAcc(startRow,startCol)>Ac % && FlowAcc(startRow-1,startCol-1)>Ac
            num_ToPoint=num_ToPoint+1; % ����õ�
            if FlowAcc(startRow-1,startCol-1)>Ac;stepPoint=[stepPoint;[startRow-1,startCol-1,FlowAcc(startRow-1,startCol-1),startRow,startCol]];end
            if FlowAcc(startRow-1,startCol-1)<=Ac;num_slope=num_slope+1;end
        end
        if startRow>1 && FlowDir(startRow-1,startCol)==2^2 % && FlowAcc(startRow,startCol)>Ac % && FlowAcc(startRow-1,startCol)>Ac
            num_ToPoint=num_ToPoint+1; % ����õ�
            if FlowAcc(startRow-1,startCol)>Ac;stepPoint=[stepPoint;[startRow-1,startCol,FlowAcc(startRow-1,startCol),startRow,startCol]];end
            if FlowAcc(startRow-1,startCol)<=Ac;num_slope=num_slope+1;end
        end
        if startRow>1 && startCol<n && FlowDir(startRow-1,startCol+1)==2^3 % && FlowAcc(startRow,startCol)>Ac  % && FlowAcc(startRow-1,startCol+1)>Ac
            num_ToPoint=num_ToPoint+1; % ����õ�
            if FlowAcc(startRow-1,startCol+1)>Ac;stepPoint=[stepPoint;[startRow-1,startCol+1,FlowAcc(startRow-1,startCol+1),startRow,startCol]];end
            if FlowAcc(startRow-1,startCol+1)<=Ac;num_slope=num_slope+1;end
        end
        % �������Ŀ �� �����<Ac��Ŀ ��ȣ���õ�Ϊˮͷ
        if num_slope==num_ToPoint; WaterHead=[WaterHead;[startRow,startCol]]; end
        end
        num_slope=0; num_ToPoint=0;
        stepPoint(1,:)=[]; % % ÿһ��������ɺ�ԭ������ʼ�� ��ջ, ��һ����Ϊ��һ����������
        empty_stepPoint=~isempty(stepPoint);
    end
end
WaterHead(:,1)=UpY-WaterHead(:,1).*cellsize+cellsize/2; % row��Y����
WaterHead(:,2)=LeftX+WaterHead(:,2).*cellsize-cellsize/2; % Col��X����
length_WaterHead = length(WaterHead);

% �ļ�·��
str2 = '_Head';    str3 = '.txt';     
Pathi_Out = strcat(str1, str2, str3);
% ȷ���ļ��д���
folderPath = fileparts(Pathi_Out);  % ��ȡ�ļ���·��
if exist(folderPath, 'dir') == 0
    mkdir(folderPath);  % ����ļ��в����ڣ��򴴽�
end

% ���ļ�
fid = fopen(Pathi_Out, 'w');  
if fid == -1
    error('�޷����ļ�: %s', Pathi_Out);  % ����ļ���ʧ�ܣ����������Ϣ
end
for j = 1:length_WaterHead
    fprintf(fid, '%f\t', WaterHead(j, 2));  % д�� X ����
    fprintf(fid, '%f\n', WaterHead(j, 1));  % д�� Y ����
end
fclose(fid);
head = WaterHead(:, [2, 1]);
end























