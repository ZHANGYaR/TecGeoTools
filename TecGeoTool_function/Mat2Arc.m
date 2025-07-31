function Mat2Arc(mat_file,file_prefix,varargin)
%
	% 用法:
	%	Mat2Arc(mat_file,file_prefix);
	%
	% 描述:
	% 	此函数将包含在mat文件中的所有有效topotoolbox文件转换为兼容ArcGIS的输出。
	%	具体来说，将任何GRIDobj转换为ASCII文件，将任何STREAMobj转换为shapefile，
	%	将任何FLOWobj转换为ArcGIS流向栅格并保存为ASCII文件，将任何有效的地图结构转换为shapefile。
	%
	% 必须输入参数:
	%	mat_file - mat文件的名称或路径
	%	file_prefix - 要添加到所有输出文件前缀的字符
	%
	% 可选输入参数:
	%	raster_type ['ascii'] - 指定栅格导出的格式，合法输入为'tif'或'ascii'
	%
	% 示例:
	%	包含以下内容的mat文件 'Data.mat'：
	%		DEM - GRIDobj
	%		FD - FLOWobj
	%		S - STREAMobj
	%		MS_KSN - 包含ksn数据的Mapstructure
	%		
	%	Mat2Arc('Data.mat','Example'); 会生成:
	%		Example_DEM.txt - DEM的ASCII文件
	%		Example_FD.txt - Arc流向栅格的ASCII文件
	%		Example_S.shp - STREAMobj的河流shapefile
	%		Example_MS_KSN.shp - 包含ksn的河流数据shapefile
	%
	%	Mat2Arc('Data.mat','Example','raster_type','tif'); 会生成:
	%		Example_DEM.tif - DEM的GeoTiff文件
	%		Example_FD.tif - Arc流向栅格的GeoTiff文件
	%		Example_S.shp - STREAMobj的河流shapefile
	%		Example_MS_KSN.shp - 包含ksn的河流数据shapefile
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% 函数由Adam M. Forte编写 - 更新日期: 06/18/18 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'Mat2Arc';
	addRequired(p,'mat_file',@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.mat'))));
	addRequired(p,'file_prefix',@(x) ischar(x));

	addParameter(p,'raster_type','ascii',@(x) ischar(validatestring(x,{'ascii','tif'})));
	addParameter(p,'out_dir',[],@(x) isdir(x));

	parse(p,mat_file,file_prefix,varargin{:});
	mat_file=p.Results.mat_file;
	file_prefix=p.Results.file_prefix;

	raster_type=p.Results.raster_type;
	out_dir=p.Results.out_dir;

	if isempty(out_dir)
		out_dir=pwd;
	end

	file_prefix=[out_dir filesep file_prefix];

	d=whos('-file',mat_file);

	num_vars=numel(d);

	for ii=1:num_vars

		classOI=d(ii,1).class;
		varNM=d(ii,1).name;

		if strcmp(classOI,'GRIDobj')
			varOI=load(mat_file,varNM);
			switch raster_type
			case 'ascii'
				out_name=[file_prefix '_' varNM '.txt'];
				GRIDobj2ascii(varOI.(varNM),out_name);
			case 'tif'
				out_name=[file_prefix '_' varNM '.tif'];
				GRIDobj2geotiff(varOI.(varNM),out_name);
			end
		elseif strcmp(classOI,'STREAMobj')
			varOI=load(mat_file,varNM);
			out_name=[file_prefix '_' varNM '.shp'];
			MS=STREAMobj2mapstruct(varOI.(varNM));
			shapewrite(MS,out_name);
		elseif strcmp(classOI,'FLOWobj')
			varOI=load(mat_file,varNM);
			G=FLOWobj2GRIDobj(varOI.(varNM));
			switch raster_type
			case 'ascii'
				out_name=[file_prefix '_' varNM '.txt'];
				GRIDobj2ascii(G,out_name);
			case 'tif'
				out_name=[file_prefix '_' varNM '.tif'];
				GRIDobj2geotiff(G,out_name);
			end
		elseif strcmp(classOI,'struct')
			varOI=load(mat_file,varNM);
			st=varOI.(varNM);
			fn=fieldnames(st);
			if any(strcmp(fn,'Geometry')) && any(strcmp(fn,'X')) && any(strcmp(fn,'Y'))
				out_name=[file_prefix '_' varNM '.shp'];
				shapewrite(st,out_name);
			end
		end
	end
end
		