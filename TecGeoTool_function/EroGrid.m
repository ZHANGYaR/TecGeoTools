function [ERO,varargout]=EroGrid(DEM,KSN,rel_type,varargin)
% 用法：
%   [ERO] = EroGrid(DEM, KSN, rel_type);
%   [ERO] = EroGrid(DEM, KSN, rel_type, 'name', value, ...);
%   [ERO, ERO_P, ERO_M] = EroGrid(DEM, KSN, rel_type, 'name', value, ...);
%
% 描述：
%   该函数基于归一化河道陡度（ksn）与侵蚀率（E）之间的经验关系生成侵蚀速率图，
%   关系形式为 ksn = C * E ^ phi，或使用Lague等(2005)描述的复杂随机阈值模型。
%   该关系通常源自集水区平均宇宙成因侵蚀率与流域平均归一化河道陡度的相关性。
%   可为不同区域（通过可选输入'VAL'定义）建立不同ksn-E关系，并考虑ksn插值不确定性
%   和/或拟合参数不确定性（如幂律模型的C和phi），计算侵蚀速率的最小/最大估计值。
%
% 必需输入：
%   DEM - DEM GRID对象
%   KSN - 连续ksn GRID对象（如KsnChiBatch输出的'ksngrid'）或ksn地图结构体
%         （如KsnChiBatch输出的'ksn'）。若KSN GRID对象由其他方法生成，
%         需确保与DEM维度/像素大小一致（validatealignment(DEM, KSN)返回true）
%   rel_type - ksn与E的关系类型，可选：
%         'power'       - 简单幂律关系 ksn = C * E ^ phi
%         'stochastic_threshold' - Lague等(2005)和DiBiase与Whipple(2011)描述的
%                                   随机阈值模型（DiBiase & Whipple 2011方程10）
%
% 幂律模型可选参数：
%   C []       - 幂律系数（ksn = C * E^phi）。可为单值或数组，数组时表示不同区域关系
%   phi []     - 幂律指数（ksn = C * E^phi）。格式同C
%   
% 随机阈值模型可选参数：
%   k_e [1e-12]       - 侵蚀效率常数（m^3/(N·yr)）
%   tau_crit [45]     - 临界剪切应力（Pa）
%   Rb [1]            - 平均径流量（mm/day）
%   k [0.5]           - 气候变异参数（反Gamma分布参数）
%   k_w [15]          - 河道宽度-流量关系的幅度因子（m^(1-omega_a)）
%   f [0.08313]       - 达西-韦斯巴赫摩擦因子（无量纲）
%   omega_a [0.55]    - 河道宽度下游尺度指数（无量纲）
%   omega_s [0.25]    - 局部流量-宽度尺度指数（无量纲）
%   alpha_val [2/3]   - 流量的摩擦指数（无量纲）
%   beta_val [2/3]    - 坡度的摩擦指数（无量纲）
%   a [1.5]           - 剪切应力指数（无量纲）
%
% 其他可选参数：
%   radius [5000]     - 若KSN为地图结构体，生成空间平均ksn网格的半径（m）
%   KSNstd            - ksn标准差GRID对象（如KsnChiBatch输出的'ksngrid'），
%                       用于量化ksn平滑不确定性。若为true且KSN为地图结构体，
%                       将计算KSNstd并纳入误差分析
%   C_std []          - 幂律系数C的标准差（需与C维度一致）
%   phi_std []        - 幂律指数phi的标准差（需与phi维度一致）
%   VAL []            - 定义分区的GRID对象（如降水图区分不同ksn-E关系）
%   edges []          - VAL的分界值数组，条目数应为C/phi条目数+1
%   resample_method ['nearest'] - VAL与DEM尺寸不匹配时的重采样方法
%                       可选'nearest'/'bilinear'/'bicubic'
%   plot_result [false] - 是否绘制结果图
%
% 注意：
%   'power'关系类型下，KSN参数不限于河道陡度，但'stochastic_threshold'必须使用陡度。
%   随机阈值模型涉及复杂数值解，不支持参数不确定性分析。
%   'power'模型的侵蚀速率单位由拟合数据决定，'stochastic_threshold'输出单位为m/Myr
%
% 示例：
%   % 单一ksn-E关系
%   [ERO] = EroGrid(DEM, KSN, 'power', 'C', 100, 'phi', 0.5);
%
%   % 基于降水量的分区关系（PRECIP为降水量GRID对象，单位m/yr）
%   [ERO] = EroGrid(DEM, KSN, 'power', 'C', [100 316 1000], 'phi', [0.5 0.5 0.5],...
%                   'VAL', PRECIP, 'edges', [0 1 2 4]);
%
%   % 计算含ksn标准差的不确定性范围
%   [ERO, ERO_P, ERO_M] = EroGrid(DEM, KSN, 'power', 'C', 100, 'phi', 0.5, 'KSNstd', KSNstd);
%
%   % 计算含参数不确定性的范围
%   [ERO, ERO_P, ERO_M] = EroGrid(DEM, KSN, 'power', 'C', 100, 'phi', 0.5,...
%                   'C_std', 5, 'phi_std', 0.05);
%
% 输出：
%   ERO    - 侵蚀速率图（单位：mm/yr或m/Myr）
%   ERO_P  - 侵蚀速率上限估计
%   ERO_M  - 侵蚀速率下限估计
%
% 示例调用：
%   % 计算并绘制侵蚀速率图
%   [ERO] = EroGrid(DEM, KSN, 'power', 'C', 100, 'phi', 0.5, 'plot_result', true);
%

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% 函数作者：Adam M. Forte - 最近更新：2019/10/30        %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% 解析输入参数
	p = inputParser;         
	p.FunctionName = 'EroGrid';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'KSN',@(x) isa(x,'GRIDobj') | isstruct(x));
	addRequired(p,'rel_type',@(x) ischar(validatestring(x,{'power','stochastic_threshold'})));

	addParameter(p,'C',[],@(x) isnumeric(x));
	addParameter(p,'phi',[],@(x) isnumeric(x));
	addParameter(p,'k_e',1e-12,@(x) isnumeric(x));
	addParameter(p,'tau_crit',45,@(x) isnumeric(x));
	addParameter(p,'Rb',1,@(x) isnumeric(x));
	addParameter(p,'k',0.5,@(x) isnumeric(x));
	addParameter(p,'k_w',15,@(x) isnumeric(x) && isscalar(x));
	addParameter(p,'f',0.08313,@(x) isnumeric(x) && isscalar(x));	
	addParameter(p,'omega_a',0.55,@(x) isnumeric(x) && isscalar(x));
	addParameter(p,'omega_s',0.25,@(x) isnumeric(x) && isscalar(x));
	addParameter(p,'alpha_val',2/3,@(x) isnumeric(x) && isscalar(x));
	addParameter(p,'beta_val',2/3,@(x) isnumeric(x) && isscalar(x));
	addParameter(p,'a',1.5,@(x) isnumeric(x) && isscalar(x));
	addParameter(p,'radius',5000,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'KSNstd',false,@(x) isa(x,'GRIDobj') | islogical(x));
	addParameter(p,'C_std',[],@(x) isnumeric(x) || isempty(x));
	addParameter(p,'phi_std',[],@(x) isnumeric(x)|| isempty(x));
	addParameter(p,'VAL',[],@(x) isa(x,'GRIDobj') || isempty(x));
	addParameter(p,'edges',[],@(x) isnumeric(x) || isempty(x));
	addParameter(p,'resample_method','nearest',@(x) ischar(validatestring(x,{'nearest','bilinear','bicubic'})));
	addParameter(p,'error_type','std',@(x) ischar(validatestring(x,{'std','std_error'})));
	addParameter(p,'plot_result',false,@(x) islogical(x) && isscalar(x));

	parse(p,DEM,KSN,rel_type,varargin{:});
	DEM=p.Results.DEM;
	KSN=p.Results.KSN;
	rel_type=p.Results.rel_type;

	C=p.Results.C;
	phi=p.Results.phi;
	k_e=p.Results.k_e;
	tau_crit=p.Results.tau_crit;
	Rb=p.Results.Rb;
	k=p.Results.k;
	k_w=p.Results.k_w;
	f=p.Results.f;
	omega_a=p.Results.omega_a;
	omega_s=p.Results.omega_s;
	alpha_val=p.Results.alpha_val;
	beta_val=p.Results.beta_val;
	a=p.Results.a;
	radius=p.Results.radius;
	ksn_std=p.Results.KSNstd;
	C_std=p.Results.C_std;
	phi_std=p.Results.phi_std;
	VAL=p.Results.VAL;
	edges=p.Results.edges;
	resample_method=p.Results.resample_method;
	error_type=p.Results.error_type;
	plot_result=p.Results.plot_result;

	switch rel_type
	case 'power'

		% 检查必要输入
		if isempty(C) | isempty(phi)
			if isdeployed
				errordlg('关系类型为"power"，必须提供"C"和"phi"参数')
			end
			error('关系类型为"power"，必须提供"C"和"phi"参数');
		end

		% 输入维度校验
		if numel(C) ~= numel(phi)
			if isdeployed
				errordlg('"C"参数的数量必须与"phi"参数的数量一致')
			end
			error('"C"参数的数量必须与"phi"参数的数量一致')
		end

		if numel(C_std) ~= numel(phi_std)
			if isdeployed
				errordlg('"C_std"和"phi_std"的条目数必须相同')
			end
			error('"C_std"和"phi_std"的条目数必须相同');
		end

		if ~isempty(C_std) & numel(C_std) ~= numel(C)
			if isdeployed
				errordlg('"C_std"和"C"的条目数必须相同')
			end
			error('"C_std"和"C"的条目数必须相同');
		end

		if ~isempty(phi_std) & numel(phi_std) ~= numel(phi)
			if isdeployed
				errordlg('"phi_std"和"phi"的条目数必须相同')
			end
			error('"phi_std"和"phi"的条目数必须相同');
		end

		% 分区参数校验
		if ~isempty(edges)
			if isempty(VAL)
				if isdeployed
					errordlg('使用"edges"时必须提供对应的"VAL"参数')
				end
				error('使用"edges"时必须提供对应的"VAL"参数')
			elseif numel(edges) ~= numel(C)+1
				if isdeployed
					errordlg('"edges"条目数与系数数量不匹配')
				end
				error('"edges"条目数与系数数量不匹配')
			elseif min(edges)> min(VAL.Z(:),[],'omitnan')
				if isdeployed
					warndlg('"edges"最小值超过"VAL"最小值，结果可能出现空值区域')
				end
				warning('"edges"最小值超过"VAL"最小值，结果可能出现空值区域')
			elseif max(edges)< max(VAL.Z(:),[],'omitnan')
				if isdeployed
					warndlg('"edges"最大值小于"VAL"最大值，结果可能出现空值区域')
				end
				warning('"edges"最大值小于"VAL"最大值，结果可能出现空值区域')			
			end
		end

		% 校验VAL与DEM对齐
		if ~isempty(VAL)
			if ~validatealignment(VAL,DEM)
				disp(['使用' resample_method '方法将VAL重采样至DEM的分辨率和尺寸']);
				VAL=resample(VAL,DEM,resample_method);
			end
		end

		% 处理KSN输入格式
		if isstruct(KSN)
			disp('生成连续KSN栅格');
			[KSN,KSNstd]=KsnAvg(DEM,KSN,radius,error_type);
			disp('连续KSN栅格生成完成');
		end

		% 确定不确定性计算标志
		if isa(ksn_std,'GRIDobj')
			ksn_std_flag=true;
			KSNstd=ksn_std;
		elseif islogical(ksn_std) && ksn_std && isa(KSNstd,'GRIDobj')
			ksn_std_flag=true;
		else
			ksn_std_flag=false;
		end

		if ~isempty(C_std)
			fit_unc_flag=true;
		else
			fit_unc_flag=false;
		end

		if fit_unc_flag | ksn_std_flag
			std_flag=true;
		else
			std_flag=false;
		end

		% 开始计算
		if isempty(edges)
			% 单一关系模型
			n=1/phi;
			K=C^(-n);

			ERO=K.*KSN.^n;

			% 处理不确定性
			if ksn_std_flag && ~fit_unc_flag
				ERO_P=K.*(KSN+KSNstd).^n;
				ERO_M=K.*(KSN-KSNstd).^n;
			elseif ~ksn_std_flag && fit_unc_flag
				n_m=1/(phi+phi_std);
				K_m=(C+C_std)^-n_m;
				n_p=1/(phi-phi_std);
				K_p=(C-C_std)^-n_p;

				ERO_P=(K_p).*KSN.^n_p;
				ERO_M=(K_m).*KSN.^n_m;
			elseif ksn_std_flag && fit_unc_flag
				n_m=1/(phi+phi_std);
				K_m=(C+C_std)^-n_m;
				n_p=1/(phi-phi_std);
				K_p=(C-C_std)^-n_p;

				ERO_P=(K_p).*(KSN+KSNstd).^n_p;
				ERO_M=(K_m).*(KSN-KSNstd).^n_m;			
			end
		else
			% 分区计算
			ERO=GRIDobj(DEM);

			if ksn_std_flag
				ERO_P=GRIDobj(DEM);
				ERO_M=GRIDobj(DEM);
			end		

			n=1./phi;
			K=C.^(-n);

			num_bins=numel(edges)-1;
			for ii=1:num_bins
				IDX = VAL>=edges(ii) & VAL<edges(ii+1);
				ERO.Z(IDX.Z)=K(ii).*KSN.Z(IDX.Z).^n(ii);

				% 处理各分区不确定性
				if ksn_std_flag && ~fit_unc_flag
					ERO_P.Z(IDX.Z)=K(ii).*(KSN.Z(IDX.Z)+KSNstd.Z(IDX.Z)).^n(ii);
					ERO_M.Z(IDX.Z)=K(ii).*(KSN.Z(IDX.Z)-KSNstd.Z(IDX.Z)).^n(ii);
				elseif ~ksn_std_flag && fit_unc_flag
					n_m(ii)=1/(phi(ii)+phi_std(ii));
					K_m(ii)=(C(ii)+C_stdi(ii))^-n_m;
					n_p(ii)=1/(phi(ii)-phi_std(ii));
					K_p(ii)=(C(ii)-C_std(ii))^-n_p;

					ERO_P.Z(IDX.Z)=K_p(ii).*KSN.Z(IDX.Z).^n_p(ii);
					ERO_M.Z(IDX.Z)=K_m(ii).*KSN.Z(IDX.Z).^n_m(ii);
				elseif ksn_std_flag && fit_unc_flag 
					n_m(ii)=1/(phi(ii)+phi_std(ii));
					K_m(ii)=(C(ii)+C_std(ii))^-n_m;
					n_p(ii)=1/(phi(ii)-phi_std(ii));
					K_p(ii)=(C(ii)-C_std(ii))^-n_p;

					ERO_P.Z(IDX.Z)=K_p(ii).*(KSN.Z(IDX.Z)+KSNstd.Z(IDX.Z)).^n_p(ii);
					ERO_M.Z(IDX.Z)=K_m(ii).*(KSN.Z(IDX.Z)-KSNstd.Z(IDX.Z)).^n_m(ii);
				end
			end
		end

		% 处理无效值
		IDX=GRIDobj(DEM,'logical');
		IDX.Z(isnan(DEM.Z))=true;

		ERO.Z(IDX.Z)=NaN;

		% 清除虚数值
		ERO.Z(imag(ERO.Z)~=0)=NaN;

		if std_flag
			ERO_P.Z(IDX.Z)=NaN;
			ERO_M.Z(IDX.Z)=NaN;

			ERO_P.Z(imag(ERO_P.Z)~=0)=NaN;
			ERO_M.Z(imag(ERO_M.Z)~=0)=NaN;

			varargout{1}=ERO_P;
			varargout{2}=ERO_M;
		end

		% 绘图
		if plot_result
			f1=figure(1);
			clf 
			set(f1,'unit','normalized','position',[0.1 0.1 0.8 0.8]);

			ksn_min=min(KSN.Z(:),[],'omitnan');
			ksn_max=max(KSN.Z(:),[],'omitnan');
			ksn_vec=linspace(ksn_min,ksn_max,100);

			sbplt1=subplot(3,3,[1:6]);
			hold on 
			imageschs(DEM,ERO,'colorbarlabel','侵蚀速率');
			disableDefaultInteractivity(sbplt1);
			hold off

			sbplt2=subplot(3,3,[7:9]);
			hold on 
			if fit_unc_flag & isempty(edges)
				E_vec=K.*ksn_vec.^n;
				plot(E_vec,ksn_vec,'-k','LineWidth',2);
				E_vec_m=K_m.*ksn_vec.^n_m;
				E_vec_p=K_p.*ksn_vec.^n_p;
				plot(E_vec_m,ksn_vec,':k','LineWidth',1);
				plot(E_vec_p,ksn_vec,':k','LineWidth',1);
			elseif fit_unc_flag & ~isempty(edges)
				for ii=1:num_bins
					E_vec=K(ii).*ksn_vec.^n(ii);
					plt(ii)=plot(E_vec,ksn_vec,'-','LineWidth',2);
					E_vec_m=K_m(ii).*ksn_vec.^n_m(ii);
					E_vec_p=K_p(ii).*ksn_vec.^n_p(ii);
					plot(E_vec_m,ksn_vec,':','LineWidth',1);
					plot(E_vec_p,ksn_vec,':','LineWidth',1);
					leg{ii}=['分区' num2str(ii)];
				end
				legend(plt,leg,'location','best');
			elseif ~isempty(edges) & ~fit_unc_flag
				for ii=1:num_bins
					E_vec=K(ii).*ksn_vec.^n(ii);
					plt(ii)=plot(E_vec,ksn_vec,'-','LineWidth',2);
					leg{ii}=['分区' num2str(ii)];
				end
				legend(plt,leg,'location','best');				
			else
				E_vec=K.*ksn_vec.^n;
				plot(E_vec,ksn_vec,'-k','LineWidth',2);
			end
			xlabel('侵蚀速率');
			ylabel('K_{sn}');
			disableDefaultInteractivity(sbplt1);
			hold off
		end

	case 'stochastic_threshold'

		% 参数校验
		if numel(k_e) ~= numel(tau_crit) | numel(k_e) ~= numel(k) | numel(k_e) ~= numel(Rb)
			if isdeployed
				errordlg('"k_e"、"tau_crit"、"k"和"Rb"参数的数量必须相同')
			end
			error('"k_e"、"tau_crit"、"k"和"Rb"参数的数量必须相同');
		end

		% 分区参数校验
		if ~isempty(edges)
			if isempty(VAL)
				if isdeployed
					errordlg('使用"edges"时必须提供对应的"VAL"参数')
				end
				error('使用"edges"时必须提供对应的"VAL"参数')
			elseif numel(edges) ~= numel(k_e)+1
				if isdeployed
					errordlg('"edges"条目数与参数数量不匹配')
				end
				error('"edges"条目数与参数数量不匹配')
			elseif min(edges)> min(VAL.Z(:),[],'omitnan')
				if isdeployed
					warndlg('"edges"最小值超过"VAL"最小值，结果可能出现空值区域')
				end
				warning('"edges"最小值超过"VAL"最小值，结果可能出现空值区域')
			elseif max(edges)< max(VAL.Z(:),[],'omitnan')
				if isdeployed
					warndlg('"edges"最大值小于"VAL"最大值，结果可能出现空值区域')
				end
				warning('"edges"最大值小于"VAL"最大值，结果可能出现空值区域')			
			end
		end

		% 校验VAL与DEM对齐
		if ~isempty(VAL)
			if ~validatealignment(VAL,DEM)
				disp(['使用' resample_method '方法将VAL重采样至DEM的分辨率和尺寸']);
				VAL=resample(VAL,DEM,resample_method);
			end
		end

		% 处理KSN输入格式
		if isstruct(KSN)
			disp('生成连续KSN栅格');
			[KSN,KSNstd]=KsnAvg(DEM,KSN,radius,error_type);
			disp('连续KSN栅格生成完成');
		end

		% 确定ksn标准差标志
		if isa(ksn_std,'GRIDobj')
			ksn_std_flag=true;
			KSNstd=ksn_std;
		elseif islogical(ksn_std) && ksn_std && isa(KSNstd,'GRIDobj')
			ksn_std_flag=true;
		else
			ksn_std_flag=false;
		end

		% 开始计算
		if isempty(edges)
			% 全局计算
			min_ksn=min(KSN.Z(:),[],'omitnan');
			max_ksn=max(KSN.Z(:),[],'omitnan');
			ERO=GRIDobj(DEM);

			% 数值求解
			[E,Ks]=stoch_thresh(min_ksn,max_ksn,k_e,tau_crit,Rb,k,k_w,f,omega_a,omega_s,alpha_val,beta_val,a);

			% 分配侵蚀速率
			ix=discretize(KSN.Z,Ks);
			w1=waitbar(0,'正在生成侵蚀速率栅格...');
			for ii=1:numel(Ks)-1
				E_val=mean([E(ii) E(ii+1)]);
				idx=ix==ii;
				ERO.Z(idx)=E_val;
				waitbar(ii/(numel(Ks)-1));
			end
			close(w1);

			% 处理不确定性
			if ksn_std_flag

				KSN_MAX=KSN+KSNstd;
				max_ksn_min=min(KSN_MAX.Z(:),[],'omitnan');
				max_ksn_max=max(KSN_MAX.Z(:),[],'omitnan');

				KSN_MIN=KSN-KSNstd;
				min_ksn_min=min(KSN_MIN.Z(:),[],'omitnan');
				min_ksn_max=max(KSN_MIN.Z(:),[],'omitnan');

				ERO_P=GRIDobj(DEM);
				ERO_M=GRIDobj(DEM);

				% 计算上下限
				[EMax,KsMax]=stoch_thresh(max_ksn_min,max_ksn_max,k_e,tau_crit,Rb,k,k_w,f,omega_a,omega_s,alpha_val,beta_val,a);
				[EMin,KsMin]=stoch_thresh(min_ksn_min,min_ksn_max,k_e,tau_crit,Rb,k,k_w,f,omega_a,omega_s,alpha_val,beta_val,a);

				% 分配上下限值
				ix_max=discretize(KSN_MAX.Z,KsMax);
				ix_min=discretize(KSN_MIN.Z,KsMin);

				w1=waitbar(0,'正在生成侵蚀速率上下限栅格...');
				for ii=1:numel(Ks)-1
					E_val_max=mean([EMax(ii) EMax(ii+1)]);
					idx_max=ix_max==ii;
					ERO_P.Z(idx_max)=E_val_max;

					E_val_min=mean([EMin(ii) EMin(ii+1)]);
					idx_min=ix_min==ii;
					ERO_M.Z(idx_min)=E_val_min;	
					waitbar(ii/(numel(Ks)-1));				
				end
				close(w1);
			end

		else
			% 分区计算
			ERO=GRIDobj(DEM);

			if ksn_std_flag
				ERO_P=GRIDobj(DEM);
				ERO_M=GRIDobj(DEM);
			end		

			num_bins=numel(edges)-1;
			for kk=1:num_bins
				IDX = VAL>=edges(kk) & VAL<edges(kk+1);

				min_ksn=min(KSN.Z(IDX.Z),[],'omitnan');
				max_ksn=max(KSN.Z(IDX.Z),[],'omitnan');	

				% 数值求解
				[E,Ks]=stoch_thresh(min_ksn,max_ksn,k_e(kk),tau_crit(kk),Rb(kk),k(kk),k_w,f,omega_a,omega_s,alpha_val,beta_val,a);

				% 分配侵蚀速率
				KSN_TEMP=KSN;
				KSN_TEMP.Z(~IDX.Z)=NaN;
				ix=discretize(KSN_TEMP.Z,Ks);
				w1=waitbar(0,['正在处理分区' num2str(kk) '...']);
				for ii=1:numel(Ks)-1
					E_val=mean([E(ii) E(ii+1)]);
					idx=ix==ii;
					ERO.Z(idx)=E_val;
					waitbar(ii/(numel(Ks)-1));	
				end	

				% 存储绘图数据
				if plot_result
					Eout{kk,1}=E;
					Ksout{kk,1}=Ks;	
				end					
				close(w1);	

				% 处理分区不确定性
				if ksn_std_flag
					KSN_MAX=KSN+KSNstd;
					max_ksn_min=min(KSN_MAX.Z(IDX.Z),[],'omitnan');
					max_ksn_max=max(KSN_MAX.Z(IDX.Z),[],'omitnan');

					KSN_MIN=KSN-KSNstd;
					min_ksn_min=min(KSN_MIN.Z(IDX.Z),[],'omitnan');
					min_ksn_max=max(KSN_MIN.Z(IDX.Z),[],'omitnan');

					% 计算上下限
					[EMax,KsMax]=stoch_thresh(max_ksn_min,max_ksn_max,k_e(kk),tau_crit(kk),Rb(kk),k(kk),k_w,f,omega_a,omega_s,alpha_val,beta_val,a);
					[EMin,KsMin]=stoch_thresh(min_ksn_min,min_ksn_max,k_e(kk),tau_crit(kk),Rb(kk),k(kk),k_w,f,omega_a,omega_s,alpha_val,beta_val,a);

					% 分配上下限值
					KSN_MAX.Z(~IDX.Z)=NaN;
					KSN_MIN.Z(~IDX.Z)=NaN;
					ix_max=discretize(KSN_MAX.Z,KsMax);
					ix_min=discretize(KSN_MIN.Z,KsMin);

					w1=waitbar(0,['处理分区' num2str(kk) '的上下限...']);
					for ii=1:numel(Ks)-1
						E_val_max=mean([EMax(ii) EMax(ii+1)]);
						idx_max=ix_max==ii;
						ERO_P.Z(idx_max)=E_val_max;

						E_val_min=mean([EMin(ii) EMin(ii+1)]);
						idx_min=ix_min==ii;
						ERO_M.Z(idx_min)=E_val_min;	
						waitbar(ii/(numel(Ks)-1));				
					end
					close(w1);
				end
			end
		end

		% 处理无效值
		IDX=GRIDobj(DEM,'logical');
		IDX.Z(isnan(DEM.Z))=true;

		ERO.Z(IDX.Z)=NaN;

		% 清除虚数值
		ERO.Z(imag(ERO.Z)~=0)=NaN;

		if ksn_std_flag
			ERO_P.Z(IDX.Z)=NaN;
			ERO_M.Z(IDX.Z)=NaN;

			ERO_P.Z(imag(ERO_P.Z)~=0)=NaN;
			ERO_M.Z(imag(ERO_M.Z)~=0)=NaN;

			varargout{1}=ERO_P;
			varargout{2}=ERO_M;
		end

		% 绘图
		if plot_result
			f1=figure(1);
			clf 
			set(f1,'unit','normalized','position',[0.1 0.1 0.8 0.8]);

			sbplt1=subplot(3,3,[1:6]);
			hold on 
			imageschs(DEM,ERO,'colorbarlabel','侵蚀速率 [m/Myr]');
			disableDefaultInteractivity(sbplt1);
			hold off

			if isempty(edges)
				sbplt2=subplot(3,3,[7:9]);
				hold on 
				plot(E,Ks,'-k','LineWidth',2);
				xlabel('侵蚀速率 [m/Myr]');
				ylabel('K_{sn}');
				disableDefaultInteractivity(sbplt1);
				hold off
			else
				sbplt2=subplot(3,3,[7:9]);
				hold on 
				for ii=1:num_bins
					plot(Eout{ii},Ksout{ii},'-','LineWidth',2);
					leg{ii}=['分区' num2str(ii)];
				end
				legend(leg,'location','best');
				xlabel('侵蚀速率 [m/Myr]');
				ylabel('K_{sn}');
				disableDefaultInteractivity(sbplt1);
				hold off				
			end
		end

	end

end

function [E,Ks]=stoch_thresh(ksn_min,ksn_max,k_e,tau_crit,Rb,k,k_w,f,omega_a,omega_s,alpha_val,beta_val,a)
	% Lague等(2005)方程16的数值积分实现，代码基于Roman DiBiase的原始版本改编

	% 将径流量Rb转换为k_q参数
	k_q=Rb/(24 * 60 * 60 * 10 * 100);

	% 派生参数
	k_t = 0.5 * 1000*(9.81^(2/3))*(f^(1/3));      % 设为1000，遵循Tucker 2004设定
	y = a*alpha_val*(1-omega_s);                % gamma指数
	m = a*alpha_val*(1-omega_a);                % 侵蚀定律中的m参数
	n = a*beta_val;                             % 侵蚀定律中的n参数
	psi_crit = k_e*tau_crit^a;                 % 侵蚀定律中的临界剪切应力项
	K = k_e*(k_t^a)*(k_w^(-a*alpha_val));       % 侵蚀效率系数

	% 生成Ksn范围
	Ks=linspace(floor(ksn_min),ceil(ksn_max),1000);
	E = zeros(size(Ks));        
	Q_starc = zeros(size(Ks));

	% 设置积分参数
	q_min = 0.00368*k;           % 最小流量（对应频率>1e-8）
	q_max = 1000000*exp(-k);     % 最大流量（对应频率<1e-8）	

	% 对Ks范围进行数值积分
	w1=waitbar(0,'正在数值积分...');
	for ii = 1:length(Ks)
	    % 计算每个Ks对应的临界流量（方程27）
	    Q_starc(ii) = ((K./psi_crit).*(Ks(ii).^(n))*(k_q^m)).^(-1./y);
	    if Q_starc(ii) < q_min
	        Q_starc(ii) = q_min;        % 防止下溢
	    elseif Q_starc(ii) > q_max
	        Q_starc(ii) = q_max - 1;    % 防止上溢
	    end
	    
	    Er =    @erosion_law;           % 侵蚀定律函数（方程13）
	    PDF =   @inv_gamma;             % 流量概率分布函数（方程3）
	    % 构造被积函数：侵蚀功率 × 概率密度
	    ErPow = @(ks,q,k,kq,kw,ke,tc,f,omega_a,omega_s,alpha_val,beta_val,a) Er(ks,q,kq,kw,ke,tc,f,omega_a,omega_s,alpha_val,beta_val,a).*PDF(q,k);       

	    % 执行数值积分
	    E(ii) = integral(@(x) ErPow(Ks(ii),x,k,k_q,k_w,k_e,tau_crit,f,omega_a,omega_s,alpha_val,beta_val,a),Q_starc(ii),q_max);
	    waitbar(ii/length(Ks));
	end
	close(w1);
end

function [I] = erosion_law(Ks,Q_star,k_q,k_w,k_e,tau_crit,f,omega_a,omega_s,alpha_val,beta_val,a)
%EROSION_LAW 计算"日"侵蚀率，实现Lague 2005方程13或Tucker 2004方程10

sec_per_yr = 31556926;     % 年秒数
Ma =  1000000;             % 百万年数

% 派生参数
k_t = 0.5 * 1000*(9.81^(2/3))*(f^(1/3));      % 遵循Tucker 2004设定
y = a*alpha_val*(1-omega_s);                % gamma指数
m = a*alpha_val*(1-omega_a);                % 侵蚀定律中的m参数
n = a*beta_val;                             % 侵蚀定律中的n参数
psi_crit = k_e*tau_crit^a;                 % 临界剪切应力项
K = k_e*(k_q^m)*(k_t^a)*(k_w^(-a*alpha_val)); % 侵蚀效率系数

% 主计算过程
I = sec_per_yr*Ma*(K*(Ks.^n).*(Q_star.^y)-psi_crit); % 侵蚀速率，单位m/Myr

end


function [f] = inv_gamma(Q_star,k)
%INV_GAMMA 计算逆Gamma概率分布函数，用于Lague等(2005)方程3

% 最近检查者：Roman DiBiase 2019/10/23

f=exp(-k./Q_star).*((k^(k+1))*Q_star.^(-(2+k)))/(gamma(k+1));

end

function [KSNGrid,KSNstdGrid] = KsnAvg(DEM,ksn_ms,radius,er_type)
% 生成空间平均的KSN网格及其标准差/标准误网格

	% 计算半径（像素单位）
	radiuspx = ceil(radius/DEM.cellsize);
	SE = strel('disk',radiuspx,0);  % 创建圆形结构元素

	% 记录当前NaN区域的掩膜
	MASK=isnan(DEM.Z);

	% 创建沿河道赋值的网格
	KSNGrid=GRIDobj(DEM);
	KSNGrid.Z(:,:)=NaN;
	% 将河道节点的KSN值填入网格
	for ii=1:numel(ksn_ms)
		ix=coord2ind(DEM,ksn_ms(ii).X,ksn_ms(ii).Y);
		ix(isnan(ix))=[];
		KSNGrid.Z(ix)=ksn_ms(ii).ksn;
	end

	% 基于半径的局部均值计算
	ISNAN=isnan(KSNGrid.Z);
    [~,L] = bwdist(~ISNAN,'e');    % 计算最近非NaN索引
    ksng = KSNGrid.Z(L);           % 用最近邻填充NaN
    FLT   = fspecial('disk',radiuspx); % 创建圆形滤波器
    ksng   = imfilter(ksng,FLT,'symmetric','same','conv'); % 执行均值滤波

    % 计算局部标准差
    nhood   = getnhood(SE);
    ksnstd   = stdfilt(ksng,nhood); 

    % 处理误差类型
    switch er_type
    case 'std_error'  % 计算标准误
    	II=~MASK; II=single(II);
    	avg_num=imfilter(II,FLT,'symmetric','same','conv'); % 计算有效像素数
    	num_nhood_pix=sum(SE.Neighborhood(:)); % 邻域总像素数
    	num_pix=avg_num.*num_nhood_pix;        % 实际参与计算的像素数
    	ksnstder=ksnstd./sqrt(num_pix);        % 计算标准误
    	ksnstder(MASK)=NaN;                   % 恢复原始NaN区域
    end

    % 恢复原始NaN区域
    ksng(MASK)=NaN;
    ksnstd(MASK)=NaN;

    % 输出结果
    KSNGrid.Z=ksng;

    switch er_type
    case 'std'        % 标准差网格
	    KSNstdGrid=GRIDobj(DEM);
	    KSNstdGrid.Z=ksnstd;
	case 'std_error'  % 标准误网格
	    KSNstdGrid=GRIDobj(DEM);
	    KSNstdGrid.Z=ksnstder;	
    end	
end