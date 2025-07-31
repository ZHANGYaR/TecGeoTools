function [outpath,ChiMap]=Auto_ChiMap(FlowAcc,FlowDir,LeftX,DownY,cellsize,startX,startY,Ac,mVSn,str1)
% 使用方法:
%   [outpath,ChiMap]=Auto_ChiMap(FlowAcc,FlowDir,LeftX,DownY,cellsize,startX,startY,Ac,mVSn,str1);
%
% 描述:
%   该函数用于计算 chi - map，并生成矢量文件（参考 willitte 2014 Science）。
%   它通过给定的汇水面积矩阵、流向矩阵以及相关的地理坐标信息，从指定的起始点开始搜索，
%   计算每个符合条件点的 chi 值，并将结果保存为矢量文件和文本文件。
%
% 必要输入:
%   FlowAcc - 汇水面积矩阵，包含了每个栅格单元的汇水面积信息。
%   FlowDir - 流向矩阵，指示了每个栅格单元的水流方向。
%   LeftX - 栅格图像起始左下点的 X 坐标。
%   DownY - 栅格图像起始左下点的 Y 坐标。
%   cellsize - 栅格图像的步长，即每个栅格单元的边长。
%   startX - 搜索起始点的 X 坐标数组。
%   startY - 搜索起始点的 Y 坐标数组。
%   Ac - 汇水面积阈值，单位为平方米。
%   mVSn - 用于计算 chi 值的参数。
%   str1 - 输出文件的名称前缀。
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数编写者: 王一舟- 更新日期: 2025年2月4日
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


format long;
[m,n]=size(FlowAcc);
ChiMap=zeros(m,n);
% ChiMap=zeros(m,n);
UpY=DownY+cellsize*m;
Ac=floor(Ac/(cellsize^2)); % Ac位汇水面积临界值

len_startX=length(startX);
stepPoint=[];

Sfea_ChiMap=struct(); %生成polyline的结构体
Sfea_num=1; % chimap矢量文件的累积增加

for ilen=1:len_startX
    % 起始搜索点的行、列号
    startRow=floor((UpY-startY(ilen))/cellsize)+1; 
    startCol=floor((startX(ilen)-LeftX)/cellsize)+1; 
    ChiMap(startRow,startCol)=0;

    stepPoint=[stepPoint;[startRow,startCol,FlowAcc(startRow,startCol),0,0]];
    empty_stepPoint=1;
    
    % 搜索点行、列号；该点的汇水面积；该点汇入区域的行、列号（0,0表示无汇入）
    while empty_stepPoint
        startRow=stepPoint(1,1);
        startCol=stepPoint(1,2);
        
        if startCol<n && FlowDir(startRow,startCol+1)==2^4 && FlowAcc(startRow,startCol+1)>Ac
            % 该点 右侧 汇入该点，大于汇水面积下限，入栈
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
             % 该点 右下方 汇入该点，大于汇水面积下限，入栈
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
            % 该点 下方 汇入该点，大于汇水面积下限，入栈
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
            % 该点 左下方 汇入该点，大于汇水面积下限，入栈
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
            % 该点 左方 汇入该点，大于汇水面积下限，入栈
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
            % 该点 左上方 汇入该点，大于汇水面积下限，入栈
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
            % 该点 上方 汇入该点，大于汇水面积下限，入栈
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
            % 该点 右上方 汇入该点，大于汇水面积下限，入栈
            stepPoint=[stepPoint;[startRow-1,startCol+1,FlowAcc(startRow-1,startCol+1),startRow,startCol]];
            ChiMap(startRow-1,startCol+1)=ChiMap(startRow,startCol)+(1/(FlowAcc(startRow,startCol).*cellsize^2))^mVSn*cellsize*sqrt(2);
            
            Sfea_ChiMap(Sfea_num).Geometry='Line';
            Sfea_ChiMap(Sfea_num).ID=Sfea_num;
            CorX=[startCol,startCol+1]; CorX=CorX.*cellsize-cellsize/2+LeftX; Sfea_ChiMap(Sfea_num).X=[CorX,NaN];
            CorY=[startRow,startRow-1]; CorY=UpY-CorY.*cellsize+cellsize/2; Sfea_ChiMap(Sfea_num).Y=[CorY,NaN];
            Sfea_ChiMap(Sfea_num).ChiValue=ChiMap(startRow-1,startCol+1);
            Sfea_num=Sfea_num+1;

        end
        % 每一轮搜索完成后，原来的起始点 出栈
        stepPoint(1,:)=[]; % 第一个点为上一步的搜索点
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