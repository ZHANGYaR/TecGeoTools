function [outpath,ChiMap]=Auto_ChiMap(FlowAcc,FlowDir,LeftX,DownY,cellsize,startX,startY,Ac,mVSn,str1)
% ʹ�÷���:
%   [outpath,ChiMap]=Auto_ChiMap(FlowAcc,FlowDir,LeftX,DownY,cellsize,startX,startY,Ac,mVSn,str1);
%
% ����:
%   �ú������ڼ��� chi - map��������ʸ���ļ����ο� willitte 2014 Science����
%   ��ͨ�������Ļ�ˮ���������������Լ���صĵ���������Ϣ����ָ������ʼ�㿪ʼ������
%   ����ÿ������������� chi ֵ�������������Ϊʸ���ļ����ı��ļ���
%
% ��Ҫ����:
%   FlowAcc - ��ˮ������󣬰�����ÿ��դ��Ԫ�Ļ�ˮ�����Ϣ��
%   FlowDir - �������ָʾ��ÿ��դ��Ԫ��ˮ������
%   LeftX - դ��ͼ����ʼ���µ�� X ���ꡣ
%   DownY - դ��ͼ����ʼ���µ�� Y ���ꡣ
%   cellsize - դ��ͼ��Ĳ�������ÿ��դ��Ԫ�ı߳���
%   startX - ������ʼ��� X �������顣
%   startY - ������ʼ��� Y �������顣
%   Ac - ��ˮ�����ֵ����λΪƽ���ס�
%   mVSn - ���ڼ��� chi ֵ�Ĳ�����
%   str1 - ����ļ�������ǰ׺��
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ������д��: ��һ��- ��������: 2025��2��4��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


format long;
[m,n]=size(FlowAcc);
ChiMap=zeros(m,n);
% ChiMap=zeros(m,n);
UpY=DownY+cellsize*m;
Ac=floor(Ac/(cellsize^2)); % Acλ��ˮ����ٽ�ֵ

len_startX=length(startX);
stepPoint=[];

Sfea_ChiMap=struct(); %����polyline�Ľṹ��
Sfea_num=1; % chimapʸ���ļ����ۻ�����

for ilen=1:len_startX
    % ��ʼ��������С��к�
    startRow=floor((UpY-startY(ilen))/cellsize)+1; 
    startCol=floor((startX(ilen)-LeftX)/cellsize)+1; 
    ChiMap(startRow,startCol)=0;

    stepPoint=[stepPoint;[startRow,startCol,FlowAcc(startRow,startCol),0,0]];
    empty_stepPoint=1;
    
    % �������С��кţ��õ�Ļ�ˮ������õ����������С��кţ�0,0��ʾ�޻��룩
    while empty_stepPoint
        startRow=stepPoint(1,1);
        startCol=stepPoint(1,2);
        
        if startCol<n && FlowDir(startRow,startCol+1)==2^4 && FlowAcc(startRow,startCol+1)>Ac
            % �õ� �Ҳ� ����õ㣬���ڻ�ˮ������ޣ���ջ
            stepPoint=[stepPoint;[startRow,startCol+1,FlowAcc(startRow,startCol+1),startRow,startCol]];
            ChiMap(startRow,startCol+1)=ChiMap(startRow,startCol)+(1/(FlowAcc(startRow,startCol).*cellsize^2))^mVSn*cellsize;
            
            Sfea_ChiMap(Sfea_num).Geometry='Line';
            Sfea_ChiMap(Sfea_num).ID=Sfea_num;
            CorX=[startCol,startCol+1]; CorX=CorX.*cellsize-cellsize/2+LeftX; Sfea_ChiMap(Sfea_num).X=[CorX,NaN];
            CorY=[startRow,startRow];   CorY=UpY-CorY.*cellsize+cellsize/2; Sfea_ChiMap(Sfea_num).Y=[CorY,NaN];
            Sfea_ChiMap(Sfea_num).ChiValue=ChiMap(startRow,startCol+1);
            Sfea_num=Sfea_num+1;
        end  
        if startRow<m && startCol<n && FlowDir(startRow+1,startCol+1)==2^5 && FlowAcc(startRow+1,startCol+1)>Ac
             % �õ� ���·� ����õ㣬���ڻ�ˮ������ޣ���ջ
             stepPoint=[stepPoint;[startRow+1,startCol+1,FlowAcc(startRow+1,startCol+1),startRow,startCol]];
             ChiMap(startRow+1,startCol+1)=ChiMap(startRow,startCol)+(1/(FlowAcc(startRow,startCol).*cellsize^2))^mVSn*cellsize*sqrt(2);
             
            Sfea_ChiMap(Sfea_num).Geometry='Line';
            Sfea_ChiMap(Sfea_num).ID=Sfea_num;
            CorX=[startCol,startCol+1]; CorX=CorX.*cellsize-cellsize/2+LeftX; Sfea_ChiMap(Sfea_num).X=[CorX,NaN];
            CorY=[startRow,startRow+1]; CorY=UpY-CorY.*cellsize+cellsize/2; Sfea_ChiMap(Sfea_num).Y=[CorY,NaN];
            Sfea_ChiMap(Sfea_num).ChiValue=ChiMap(startRow+1,startCol+1);
            Sfea_num=Sfea_num+1;

        end
        if startRow<m && FlowDir(startRow+1,startCol)==2^6 && FlowAcc(startRow+1,startCol)>Ac
            % �õ� �·� ����õ㣬���ڻ�ˮ������ޣ���ջ
            stepPoint=[stepPoint;[startRow+1,startCol,FlowAcc(startRow+1,startCol),startRow,startCol]];
            ChiMap(startRow+1,startCol)=ChiMap(startRow,startCol)+(1/(FlowAcc(startRow,startCol).*cellsize^2))^mVSn*cellsize;
            
            Sfea_ChiMap(Sfea_num).Geometry='Line';
            Sfea_ChiMap(Sfea_num).ID=Sfea_num;
            CorX=[startCol,startCol];   CorX=CorX.*cellsize-cellsize/2+LeftX; Sfea_ChiMap(Sfea_num).X=[CorX,NaN];
            CorY=[startRow,startRow+1]; CorY=UpY-CorY.*cellsize+cellsize/2; Sfea_ChiMap(Sfea_num).Y=[CorY,NaN];
            Sfea_ChiMap(Sfea_num).ChiValue=ChiMap(startRow+1,startCol);
            Sfea_num=Sfea_num+1;

        end
        if startRow<m && startCol>1 && FlowDir(startRow+1,startCol-1)==2^7 && FlowAcc(startRow+1,startCol-1)>Ac
            % �õ� ���·� ����õ㣬���ڻ�ˮ������ޣ���ջ
            stepPoint=[stepPoint;[startRow+1,startCol-1,FlowAcc(startRow+1,startCol-1),startRow,startCol]];
            ChiMap(startRow+1,startCol-1)=ChiMap(startRow,startCol)+(1/(FlowAcc(startRow,startCol).*cellsize^2))^mVSn*cellsize*sqrt(2);
            
            Sfea_ChiMap(Sfea_num).Geometry='Line';
            Sfea_ChiMap(Sfea_num).ID=Sfea_num;
            CorX=[startCol,startCol-1]; CorX=CorX.*cellsize-cellsize/2+LeftX; Sfea_ChiMap(Sfea_num).X=[CorX,NaN];
            CorY=[startRow,startRow+1]; CorY=UpY-CorY.*cellsize+cellsize/2; Sfea_ChiMap(Sfea_num).Y=[CorY,NaN];
            Sfea_ChiMap(Sfea_num).ChiValue=ChiMap(startRow+1,startCol-1);
            Sfea_num=Sfea_num+1;

        end
        if startCol>1 && FlowDir(startRow,startCol-1)==2^0 && FlowAcc(startRow,startCol-1)>Ac
            % �õ� �� ����õ㣬���ڻ�ˮ������ޣ���ջ
            stepPoint=[stepPoint;[startRow,startCol-1,FlowAcc(startRow,startCol-1),startRow,startCol]];
            ChiMap(startRow,startCol-1)=ChiMap(startRow,startCol)+(1/(FlowAcc(startRow,startCol).*cellsize^2))^mVSn*cellsize;
            
            Sfea_ChiMap(Sfea_num).Geometry='Line';
            Sfea_ChiMap(Sfea_num).ID=Sfea_num;
            CorX=[startCol,startCol-1]; CorX=CorX.*cellsize-cellsize/2+LeftX; Sfea_ChiMap(Sfea_num).X=[CorX,NaN];
            CorY=[startRow,startRow];   CorY=UpY-CorY.*cellsize+cellsize/2; Sfea_ChiMap(Sfea_num).Y=[CorY,NaN];
            Sfea_ChiMap(Sfea_num).ChiValue=ChiMap(startRow,startCol-1);
            Sfea_num=Sfea_num+1;

        end
        if startRow>1 && startCol>1 && FlowDir(startRow-1,startCol-1)==2^1 && FlowAcc(startRow-1,startCol-1)>Ac
            % �õ� ���Ϸ� ����õ㣬���ڻ�ˮ������ޣ���ջ
            stepPoint=[stepPoint;[startRow-1,startCol-1,FlowAcc(startRow-1,startCol-1),startRow,startCol]];
            ChiMap(startRow-1,startCol-1)=ChiMap(startRow,startCol)+(1/(FlowAcc(startRow,startCol).*cellsize^2))^mVSn*cellsize*sqrt(2);
            
            Sfea_ChiMap(Sfea_num).Geometry='Line';
            Sfea_ChiMap(Sfea_num).ID=Sfea_num;
            CorX=[startCol,startCol-1]; CorX=CorX.*cellsize-cellsize/2+LeftX; Sfea_ChiMap(Sfea_num).X=[CorX,NaN];
            CorY=[startRow,startRow-1]; CorY=UpY-CorY.*cellsize+cellsize/2; Sfea_ChiMap(Sfea_num).Y=[CorY,NaN];
            Sfea_ChiMap(Sfea_num).ChiValue=ChiMap(startRow-1,startCol-1);
            Sfea_num=Sfea_num+1;

        end
        if startRow>1 && FlowDir(startRow-1,startCol)==2^2 && FlowAcc(startRow-1,startCol)>Ac
            % �õ� �Ϸ� ����õ㣬���ڻ�ˮ������ޣ���ջ
            stepPoint=[stepPoint;[startRow-1,startCol,FlowAcc(startRow-1,startCol),startRow,startCol]];
            ChiMap(startRow-1,startCol)=ChiMap(startRow,startCol)+(1/(FlowAcc(startRow,startCol).*cellsize^2))^mVSn*cellsize;
            
            Sfea_ChiMap(Sfea_num).Geometry='Line';
            Sfea_ChiMap(Sfea_num).ID=Sfea_num;
            CorX=[startCol,startCol];   CorX=CorX.*cellsize-cellsize/2+LeftX; Sfea_ChiMap(Sfea_num).X=[CorX,NaN];
            CorY=[startRow,startRow-1]; CorY=UpY-CorY.*cellsize+cellsize/2; Sfea_ChiMap(Sfea_num).Y=[CorY,NaN];
            Sfea_ChiMap(Sfea_num).ChiValue=ChiMap(startRow-1,startCol);
            Sfea_num=Sfea_num+1;

        end
        if startRow>1 && startCol<n && FlowDir(startRow-1,startCol+1)==2^3 && FlowAcc(startRow-1,startCol+1)>Ac
            % �õ� ���Ϸ� ����õ㣬���ڻ�ˮ������ޣ���ջ
            stepPoint=[stepPoint;[startRow-1,startCol+1,FlowAcc(startRow-1,startCol+1),startRow,startCol]];
            ChiMap(startRow-1,startCol+1)=ChiMap(startRow,startCol)+(1/(FlowAcc(startRow,startCol).*cellsize^2))^mVSn*cellsize*sqrt(2);
            
            Sfea_ChiMap(Sfea_num).Geometry='Line';
            Sfea_ChiMap(Sfea_num).ID=Sfea_num;
            CorX=[startCol,startCol+1]; CorX=CorX.*cellsize-cellsize/2+LeftX; Sfea_ChiMap(Sfea_num).X=[CorX,NaN];
            CorY=[startRow,startRow-1]; CorY=UpY-CorY.*cellsize+cellsize/2; Sfea_ChiMap(Sfea_num).Y=[CorY,NaN];
            Sfea_ChiMap(Sfea_num).ChiValue=ChiMap(startRow-1,startCol+1);
            Sfea_num=Sfea_num+1;

        end
        % ÿһ��������ɺ�ԭ������ʼ�� ��ջ
        stepPoint(1,:)=[]; % ��һ����Ϊ��һ����������
        empty_stepPoint=~isempty(stepPoint);
    end
end


str2='_ChiMap';str3='.shp';
str=strcat(str1,str2,str3);
outpath=str;
shapewrite(Sfea_ChiMap,str);

str2='_ChiMap';str3='.txt';
str=strcat(str1,str2,str3);
fid=fopen(str,'w');
for j=1:m
    for k=1:n-1
        fprintf(fid,'%f\t',ChiMap(j,k));
    end
    fprintf(fid,'%f\n',ChiMap(j,n));
end
fclose(fid);

end