function CheckDependencies()
	%
	% 用法:
	%	CheckTAKDependencies;
	%
	% 描述:
	% 	本函数用于检查运行所需的MATLAB工具箱
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% 函数作者：Adam M. Forte - 更新日期：2018年6月18日 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	p=ver; % 获取已安装工具箱信息

	warn_flag=false; % 初始化警告标志

	% 检查TopoToolbox安装情况
	ix=find(strcmp(cellstr(char(p.Name)),'TopoToolbox'));
	if isempty(ix)
		warning('致命错误：未检测到TopoToolbox安装或未添加至路径，请从 https://github.com/wschwanghart/topotoolbox/ 下载');
	end

	% 验证TopoToolbox版本
	if str2num(p(ix).Version)<2.3
		warning('TopoToolbox版本低于2.3，请下载最新版本：https://github.com/wschwanghart/topotoolbox/');
	end

	% 检查图像处理工具箱
	ix=find(strcmp(cellstr(char(p.Name)),'Image Processing Toolbox'));
	if isempty(ix)
		warning('致命错误：未检测到图像处理工具箱授权，TopoToolbox无法正常运行');
		warn_flag=true;
	end

	% 检查地图工具箱
	ix=find(strcmp(cellstr(char(p.Name)),'Mapping Toolbox'));
	if isempty(ix)
		warning('致命错误：未检测到地图工具箱授权，TopoToolbox无法正常运行');
		warn_flag=true;
	end

	% 检查曲线拟合工具箱
	ix=find(strcmp(cellstr(char(p.Name)),'Curve Fitting Toolbox'));
	if isempty(ix)
		warning('警告：未检测到曲线拟合工具箱授权，部分功能受限');
		warn_flag=true;
	end	

	% 检查优化工具箱
	ix=find(strcmp(cellstr(char(p.Name)),'Optimization Toolbox'));
	if isempty(ix)
		warning('警告：未检测到优化工具箱授权，部分功能受限');
		warn_flag=true;
	end	

	% 检查统计与机器学习工具箱
	ix=find(strcmp(cellstr(char(p.Name)),'Statistics and Machine Learning Toolbox'));
	if isempty(ix)
		warning('警告：未检测到统计与机器学习工具箱授权，部分功能受限');
		warn_flag=true;
	end	

	% 检查并行计算工具箱
	ix=find(strcmp(cellstr(char(p.Name)),'Parallel Computing Toolbox'));
	if isempty(ix)
		warning('警告：未检测到并行计算工具箱授权，部分功能受限');
		warn_flag=true;
	end

	% 检查MATLAB版本
	ix=find(strcmp(cellstr(char(p.Name)),'MATLAB'));
	if str2num(p(ix).Version)<9.4
		warning('提示：工具箱使用MATLAB 2021a编写，旧版可能兼容异常');
	end	

	% 汇总警告信息
	if warn_flag
		warning('系统缺少部分必备工具箱，建议使用预编译版本保证功能完整性');
	end

end