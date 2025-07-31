function  outpath=Auto_ksMap(Elev,FlowDir,ChiMap,cellsize,headX,headY,LeftX,DownY,LenWin,str1)
%
% 使用方法:
%   outpath=Auto_ksMap(Elev,FlowDir,ChiMap,cellsize,headX,headY,LeftX,DownY,LenWin,str1);
%
% 描述:
%   该函数用于直接生成矢量与栅格形式的 ksMap（河道陡峭度指数图），并按照高程间隔进行计算。函数从指定的河流起始点（headX, headY）开始搜索，根据流向矩阵（FlowDir）和 chi-map（ChiMap）等信息，计算每个符合条件点的局部 ksn 值，并将结果保存为矢量文件和文本文件。
%
% 必要输入:
%   Elev - 高程矩阵，包含了每个栅格单元的高程信息。
%   FlowDir - 流向矩阵，指示了每个栅格单元的水流方向，取值范围为 1（East）、2（SE）、4（South）、8（SW）、16（West）、32（NW）、64（North）、128（NE）。
%   ChiMap - chi-map 矩阵，存储了每个栅格单元的 chi 值。
%   cellsize - 栅格图像的步长，即每个栅格单元的边长。
%   headX - 计算 pntKs 的起始点（河流起点）的 X 坐标数组。
%   headY - 计算 pntKs 的起始点（河流起点）的 Y 坐标数组。
%   LeftX - 栅格图像起始左下点的 X 坐标。
%   DownY - 栅格图像起始左下点的 Y 坐标。
%   LenWin - 搜索的距离限制，当搜索的河段长度达到该值时，开始计算局部 ksn 值。
%   str1 - 输出文件的名称前缀。
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数编写者: 王一舟 - 更新日期: 2025年2月4日
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
[m,n]=size(Elev);lenheadX=length(headX); % RiverPro=cell(lenheadX,1);
ksMap=zeros(m,n); % ks分布数据
ksValue=[];
% RasterToASCII工具中得到的txt的左下角坐标，但行列号开始的是左上
upY=DownY+m*cellsize;
max_ksn=1000; % 当超过最大值，则全部赋这个值
% MulNum=zeros(lenheadX,3);

Sfea_ksnMap=struct(); Sfea_num=1;
%%
for j=1:lenheadX
%% midMtrix存储搜索到的栅格，信息：行号，列号，高程，chi；"kk"记录midMtrix行数；midMtrix_len记录该河段长度；kk_sfea记录河段中，实际纳入矢量显示的行数
    midMtrix=[]; kk=0; midMtrix_len=0; kk_sfea=0; 
    % 从midRow,midCol开始寻找 
    midRow=floor((abs(upY-headY(j)))/cellsize)+1; midCol=floor((abs(headX(j)-LeftX))/cellsize)+1; 
    % 起始的搜索点越界，直接退出本次循环
    if midRow<=1 || midRow>=m || midCol<=1 || midCol>=n;    continue;    end
    
    % midRow,midCol虽然不越界，但这些点有可能ChiMap为0 （点在河道上，但面积有限，未被计算ChiMap）
    while (midRow>1 && midRow<m && midCol>1 && midCol<n && ChiMap(midRow,midCol)==0)  % 找到ChiMap的起始点，即ChiMap~=0，开始记录
        switch FlowDir(midRow,midCol)
            case 1;   midRow=midRow;   midCol=midCol+1; % East
            case 2;   midRow=midRow+1; midCol=midCol+1; % SE
            case 4;   midRow=midRow+1; midCol=midCol;   % South
            case 8;   midRow=midRow+1; midCol=midCol-1; % SW
            case 16;  midRow=midRow;   midCol=midCol-1; % West
            case 32;  midRow=midRow-1; midCol=midCol-1; % NW
            case 64;  midRow=midRow-1; midCol=midCol;   % North
            case 128; midRow=midRow-1; midCol=midCol+1; % NE
        end
    end
    
    % 从开始搜索ChiMap非0点，但一直到越界了，都没找到，直接退出本次循环
    if midRow<=1 || midRow>=m || midCol<=1 || midCol>=n;    continue;    end
    
    % 初始搜索点
    midMtrix=[midMtrix;midRow,midCol,Elev(midRow,midCol),ChiMap(midRow,midCol)]; kk=kk+1; kk_sfea=kk_sfea+1;

%%   
    while (midRow>1 && midRow<m && midCol>1 && midCol<n && ChiMap(midRow,midCol)>0)  % 该点是否值得搜索：必须不越矩阵的界，且有ChiMap值 
        % while中的判断，可以保证，即便在行列号追加之后，还能不越界；因为，河流出山口的ChiMap是0值，所以这里需要的判别条件是ChiMap>0
        
%         midMtrix=[midMtrix;midRow,midCol,Elev(midRow,midCol),ChiMap(midRow,midCol)]; kk=kk+1; % kk_sfea=kk_sfea+1;
%         if ksMap(midRow,midCol)==0;            kk_sfea=kk_sfea+1;        end
        
        if FlowDir(midRow,midCol)==1
            midRow=midRow;   midCol=midCol+1;    % East
            midMtrix=[midMtrix;midRow,midCol,Elev(midRow,midCol),ChiMap(midRow,midCol)]; kk=kk+1; midMtrix_len=midMtrix_len+cellsize; % 不管ksnmap是否被计算，都纳入计算局部ksn的河段
            if ksMap(midRow,midCol)==0;  kk_sfea=kk_sfea+1; end % 但是，只有当ksnmap未被计算，才纳入矢量显示
        elseif FlowDir(midRow,midCol)==2
            midRow=midRow+1; midCol=midCol+1;    % SE
            midMtrix=[midMtrix;midRow,midCol,Elev(midRow,midCol),ChiMap(midRow,midCol)]; kk=kk+1; midMtrix_len=midMtrix_len+cellsize*sqrt(2);  
            if ksMap(midRow,midCol)==0;  kk_sfea=kk_sfea+1; end
        elseif FlowDir(midRow,midCol)==4
            midRow=midRow+1; midCol=midCol;      % South
            midMtrix=[midMtrix;midRow,midCol,Elev(midRow,midCol),ChiMap(midRow,midCol)]; kk=kk+1; midMtrix_len=midMtrix_len+cellsize;
            if ksMap(midRow,midCol)==0;  kk_sfea=kk_sfea+1; end
        elseif FlowDir(midRow,midCol)==8
            midRow=midRow+1; midCol=midCol-1;    % SW
            midMtrix=[midMtrix;midRow,midCol,Elev(midRow,midCol),ChiMap(midRow,midCol)]; kk=kk+1; midMtrix_len=midMtrix_len+cellsize*sqrt(2);
            if ksMap(midRow,midCol)==0;  kk_sfea=kk_sfea+1; end
        elseif FlowDir(midRow,midCol)==16
            midRow=midRow;   midCol=midCol-1;    % West
            midMtrix=[midMtrix;midRow,midCol,Elev(midRow,midCol),ChiMap(midRow,midCol)]; kk=kk+1; midMtrix_len=midMtrix_len+cellsize;  
            if ksMap(midRow,midCol)==0;  kk_sfea=kk_sfea+1; end
        elseif FlowDir(midRow,midCol)==32
            midRow=midRow-1; midCol=midCol-1;    % NW
            midMtrix=[midMtrix;midRow,midCol,Elev(midRow,midCol),ChiMap(midRow,midCol)]; kk=kk+1; midMtrix_len=midMtrix_len+cellsize*sqrt(2); 
            if ksMap(midRow,midCol)==0;  kk_sfea=kk_sfea+1; end
        elseif FlowDir(midRow,midCol)==64
            midRow=midRow-1; midCol=midCol;      % North
            midMtrix=[midMtrix;midRow,midCol,Elev(midRow,midCol),ChiMap(midRow,midCol)]; kk=kk+1; midMtrix_len=midMtrix_len+cellsize;      
            if ksMap(midRow,midCol)==0;  kk_sfea=kk_sfea+1; end
        elseif FlowDir(midRow,midCol)==128
            midRow=midRow-1; midCol=midCol+1;    % NE
            midMtrix=[midMtrix;midRow,midCol,Elev(midRow,midCol),ChiMap(midRow,midCol)]; kk=kk+1; midMtrix_len=midMtrix_len+cellsize*sqrt(2); 
            if ksMap(midRow,midCol)==0;  kk_sfea=kk_sfea+1; end
        end
        
        %%
        if midMtrix_len>=LenWin % 搜索点未过界，就已经超过距离限制
            % 计算河段的局部ksn值；midMtrix存储搜索到的栅格，信息：行号，列号，高程，chi；
            [Ks,Ks_error]=regress(midMtrix(:,3),[ones(kk,1),midMtrix(:,4)]); 
            if abs(Ks(2))>max_ksn; ksValue_1=max_ksn; else; ksValue_1=abs(Ks(2)); end
            ksValue=[ksValue;[ksValue_1,abs(Ks_error(2,2)-Ks(2))]];
            
            % 生成栅格文件 
            for i_kk=1:kk_sfea; ksMap(midMtrix(i_kk,1),midMtrix(i_kk,2))=ksValue_1;end
            
            % 生成矢量文件  midMtrix 行号，列号，高程，chi   
            Sfea_ksnMap(Sfea_num).Geometry='Line';
            Sfea_ksnMap(Sfea_num).ID=Sfea_num;
            if kk>kk_sfea % 在距离限制内，栅格被计算的ksn值
                CorX=midMtrix(1:(kk_sfea+1),2); CorX=CorX.*cellsize-cellsize/2+LeftX; Sfea_ksnMap(Sfea_num).X=[CorX;NaN];%如果出错，就去掉"+1"
                CorY=midMtrix(1:(kk_sfea+1),1); CorY=upY-CorY.*cellsize+cellsize/2;   Sfea_ksnMap(Sfea_num).Y=[CorY;NaN];
                Sfea_ksnMap(Sfea_num).ksn=ksValue_1;
                Sfea_num=Sfea_num+1;
                midMtrix=[]; kk=0; midMtrix_len=0; kk_sfea=0;
                break; % 当在距离限制内，下段栅格已经被计算ksn，一般出现在河流汇合处，退出由这个初始点开始while语句的河道搜索
            else % 在距离限制内，栅格未被计算的ksn值
                CorX=midMtrix(1:(kk_sfea),2); CorX=CorX.*cellsize-cellsize/2+LeftX; Sfea_ksnMap(Sfea_num).X=[CorX;NaN]; 
                CorY=midMtrix(1:(kk_sfea),1); CorY=upY-CorY.*cellsize+cellsize/2;   Sfea_ksnMap(Sfea_num).Y=[CorY;NaN];
                Sfea_ksnMap(Sfea_num).ksn=ksValue_1;
                Sfea_num=Sfea_num+1;
                midMtrix=[]; kk=0; midMtrix_len=0; kk_sfea=0;
            end
%             midMtrix=[midMtrix;midRow,midCol,Elev(midRow,midCol),ChiMap(midRow,midCol)]; kk=kk+1; % midMtrix_len=midMtrix_len+cellsize*sqrt(2); 
            
        else % 搜索点没有超过距离限制，进行是否越界判别
            if midRow<=1 || midRow>=m || midCol<=1 || midCol>=n %|| ChiMap(midRow,midCol)<=0
            ChiMap_text=ChiMap(midRow,midCol);
            if ChiMap_text==0  % 搜索点还没有超过距离限制，就已经过界
            % 生成栅格文件 
                len_ksValue=length(ksValue(:,1)); ksValue_1=ksValue(len_ksValue,1); 
                for i_kk=1:kk_sfea;  ksMap(midMtrix(i_kk,1),midMtrix(i_kk,2))=ksValue_1;  end
            
            % 生成矢量文件  midMtrix 行号，列号，高程，chi   
                Sfea_ksnMap(Sfea_num).Geometry='Line';
                Sfea_ksnMap(Sfea_num).ID=Sfea_num;
                CorX=midMtrix(1:(kk_sfea),2); CorX=CorX.*cellsize-cellsize/2+LeftX; Sfea_ksnMap(Sfea_num).X=[CorX;NaN];
                CorY=midMtrix(1:(kk_sfea),1); CorY=upY-CorY.*cellsize+cellsize/2;   Sfea_ksnMap(Sfea_num).Y=[CorY;NaN];
                Sfea_ksnMap(Sfea_num).ksn=ksValue_1;
                Sfea_num=Sfea_num+1;
                midMtrix=[]; kk=0; midMtrix_len=0; kk_sfea=0;
                break;
            end
            end
        end

%         if (ksMap(midRow,midCol)~=0 || Elev(midRow,midCol)<0) 
%             break; % 该点是否值得搜索：ksMap（空值=0）已经被计算 或者 越高程界；那么下一个不必再寻找，直接退出while循环，寻找下一个起始点   
%         end 
%         if ChiMap(midRow,midCol)>0 % 判断搜索出的点，是否满足要求（ChiMap非空，也是针对起始点的问题而言的）
%             midMtrix=[midMtrix;midRow,midCol,Elev(midRow,midCol),ChiMap(midRow,midCol)];kk=kk+1;
%         end
    end
end

str2='_KsnMap';str3='.shp';
str=strcat(str1,str2,str3);
outpath=str;
shapewrite(Sfea_ksnMap,str);

str2='_KsnMap';str3='.txt';
str=strcat(str1,str2,str3);
fid=fopen(str,'w');
for i=1:m
    for j=1:n-1; fprintf(fid,'%f ',ksMap(i,j));end
    fprintf(fid,'%f\n',ksMap(i,n));
end
fclose(fid);
end