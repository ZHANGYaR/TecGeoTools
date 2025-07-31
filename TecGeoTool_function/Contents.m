% TecGeoTools V1.0 
%
%
% 如在出版物中使用或修改对应功能，请引用论文:
% Forte, A. M., & Whipple, K. X. (2019). Short communication: 
% The Topographic Analysis Kit (TAK) for TopoToolbox. Earth Surface Dynamics,
% 7(1), 87-95. doi:10.5194/esurf-7-87-2019
% 
% Schwanghart, W., & Scherler, D. (2014). Short Communication: TopoToolbox 2 – 
% MATLAB-based software for topographic analysis and modeling in Earth surface 
% sciences. Earth Surface Dynamics, 2(1), 1-7. doi:10.5194/esurf-2-1-2014
% 
% 王一舟, 郑德文, & 张会平. (2022). 河流高程剖面分析的方法与程序实现——基于
% Matlab平台编写的开源函数集RiverProAnalysis. 中国科学:地球科学, 
% 52(10), 22.doi:10.1360/SSTe-2021-0261

% 所需MATLAB工具箱列表：
%
% Image Processing Toolbox 图像处理工具箱
% Mapping Toolbox 地图工具箱
% Statistics and Machine Learning Toolbox 统计与机器学习工具箱
% Curve Fitting Toolbox 曲线拟合工具箱
% Optimization Toolbox 优化工具箱
% Parallel Computing Toolbox 并行计算工具箱
%
% 建议使用MATLAB 2021a或更高版本以获得完整功能支持。
%
% 主要功能列表： 
% Auto_ChiMap            自动生成ksn图
% Auto_ksMap             自动生成chi图
% AutoKsnProfiler        自动河道剖面分析
% Basin2Raster           流域转栅格
% Basin2Shape            流域转矢量
% BasinPicker            交互式流域选择
% BasinStatsPlots        流域统计绘图
% BestFitSmallCircle     最佳拟合小圆
% CatPoly2GRIDobj        面要素转栅格对象
% CheckTAKDepedencies    检查TAK依赖项
% ClassifyKnicks         河道拐点分类
% CompileBasinStats      流域统计编译
% ConditionDEM           地形数据预处理
% CrossSwath             交叉地形剖面
% DippingBedFinder       地层倾斜识别
% EroGrid                侵蚀速率网格生成
% FindBasinKnicks        流域拐点识别
% FindCentroid           流域质心定位
% FindThreshold          阈值查找
% HackRelationship       河网分级关系分析
% InspectJunction        节点检查
% JunctionAngle          河道交汇角度计算
% JunctionLinks          节点连接分析
% KsnChiBatch            批量χ-ksn分析
% ksncolor               着色带
% KsnProfiler            河道陡峭指数剖面生成
% linear_uplift_uniformChi 线性模拟隆升历史
% LinearUp_MisFit        寻找最佳模拟chi值
% MakeCombinedSwath      合并剖面带
% MakeSerialSwath        创建连续剖面
% MakeStreams            河网生成
% MakeTopoSwath          生成地形剖面带
% Mat2Arc                MATLAB数据转ArcGIS格式
% PlotIndividualBasins   单流域可视化
% ProcessRiverBasins     流域处理流程
% ProjectedIncision      投影切割分析
% ProjectGPSOntoSwath    GPS数据投影至剖面
% ProjectOntoSwath       通用数据剖面投影
% ProjectSmallCircleOntoSwath 小圆投影至剖面
% RemoveFlats            平坦区域处理
% SegmentPicker          交互式河段选择
% SegmentPlotter         河段可视化
% SegmentProjector       河段投影分析
% SubDivideBigBasins     细分大流域