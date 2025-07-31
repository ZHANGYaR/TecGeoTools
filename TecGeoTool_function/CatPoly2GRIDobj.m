function [OUT,look_table]=CatPoly2GRIDobj(DEM,poly_shape,field,varargin)
	%
% 用法：
%   [GRIDobj,look_up_table]=CatPoly2GRIDobj(DEM,poly_shape,field);
%
% 描述：
%   本函数将分类矢量数据（如地质图数字化成果）转换为与DEM匹配的GRIDobj格式，
%   适用于'ProcessRiverBasins'处理流程
%
% 必需输入：
%   DEM - 目标数字高程模型，用于定义输出栅格的空间参考
%   poly_shape - 分类矢量文件路径或名称（shapefile格式）
%   field - 矢量文件中包含分类数据的字段名称
%
% 可选输入：
%   table_in - 自定义查找表，需为包含两列的表格：
%       第一列为整型数值（n x 1数组），第二列为分类名称（n x 1单元格数组）
%       分类名称必须与矢量文件'field'字段中的条目完全匹配，否则输出栅格值全为0
%
% 输出：
%   OUT - 与DEM空间一致的GRIDobj，数值对应look_table中的分类编码
%   look_table - 查找表（n x 2表格），包含'Numbers'和'Categories'两列，
%       用于数字编码与原始分类的对应关系
%
% 示例：
%   [GEO,geo_table]=CatPoly2GRIDobj(DEM,'地质图.shp','岩性类型');
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 函数作者：Adam M. Forte - 更新日期：2019年7月5日 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% 解析输入参数
	p = inputParser;
	p.FunctionName = 'CatPoly2GRIDobj';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'poly_shape',@(x) ischar(x));
	addRequired(p,'field',@(x) ischar(x));

	addParameter(p,'table_in',[],@(x) istable(x));

	parse(p,DEM,poly_shape,field,varargin{:});
	DEM=p.Results.DEM;
	poly_shape=p.Results.poly_shape;
	field=p.Results.field;

	table_in=p.Results.table_in;

	% 验证自定义查找表格式
	if ~isempty(table_in)
		if size(table_in,2)~=2
			error('提供的"table_in"必须包含两列');
		elseif table_in.(1)(1)~=0
			error('"table_in"第一列的首个元素必须为0');
		elseif ~iscell(table_in.(2))
			error('"table_in"的第二列必须为单元格数组');
		elseif ~ischar(table_in.(2){1})
			error('"table_in"第二列的元素必须为字符类型');
		end
	end
			

	% 读取矢量文件并转换为表格
	PS=shaperead(poly_shape);
	TS=struct2table(PS);

	% 提取目标字段数据
	Foi=TS.(field);

	if isempty(table_in)
		% 生成唯一类别列表及查找表
		Categories=unique(Foi);
		Numbers=[1:numel(Categories)]; Numbers=Numbers';
		% 添加'未定义'类别处理可能存在的读取误差
		Numbers=vertcat(0,Numbers);
		Categories=vertcat('未定义',Categories);
		look_table=table(Numbers,Categories);
	else 
		% 加载用户提供的查找表
		Categories=table_in.(2);
		Numbers=table_in.(1);
		look_table=table(Numbers,Categories);
	end

	% 将分类名称替换为数字编码
	for ii=1:numel(PS)
		Eoi=PS(ii,1).(field);
		ix=find(strcmp(Categories,Eoi));
		PS(ii,1).replace_number=Numbers(ix);
	end

	% 创建缓冲栅格以避免边缘0值问题
	DEMp=GRIDobj(DEM,'logical');
	DEMp.Z(:,:)=true;
	DEMp=pad(DEMp,1,false);
	[OUT]=polygon2GRIDobj(DEMp,PS,'replace_number');
	OUT=crop(OUT,DEMp);

	% 清理查找表中未使用的分类
	pres=unique(OUT.Z(:));
	ix=ismember(look_table.Numbers,pres);
	look_table=look_table(ix,:);
end