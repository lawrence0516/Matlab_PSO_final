%opensource
clear all; clc;
warning('off');

% --- 啟動平行運算 ---
if isempty(gcp('nocreate'))
    parpool; 
end

% --- 實驗設定 ---
D = 30;             
Xmin = -100; Xmax = 100;
fes_max = 10000 * D; 
pop_size = 60;    
Max_Gen = ceil(fes_max / pop_size);
fhd = str2func('cec17_func');
jingdu = 0;
func_list = [1]; 

scenarios = {
    0.2, '--g', '$c_{i,1} = 0.2$', 0.25 ;  
    0.5, '-.b', '$c_{i,1} = 0.5$', 0.25 ;  
    0.8, ':k',  '$c_{i,1} = 0.8$', 0.25;  
    [],  '-r',  '$c_{i,1}\ in\ MFCPSO$',0.25  ; 
};

results_curves = cell(1, 4); 
fprintf('開始運算\n');
tic; 

parfor s = 1:4
    % 取得參數
    fixed_c1_val = scenarios{s, 1};
    % label_name = scenarios{s, 3};
    fixed_p_val  = scenarios{s, 4}; 

    [~, curve] = Gaussian_p_freedom(jingdu, func_list(1), fhd, D, pop_size, Max_Gen, fes_max, Xmin, Xmax, fixed_c1_val, fixed_p_val, func_list(1));
    
    % 數據截斷處理
    curve = max(curve, 1e-8); % 避免 log(0) 或過小
    
    results_curves{s} = curve;
end

elapsed_time = toc;

% 1. 設定畫布大小 (稍微加大一點，讓字體不擁擠)
figure('Name', 'Fig 4 Reproduction', 'Color', 'w', 'Position', [150, 150, 700, 550]); 
hold on;

% 2. 畫線 (Line Plotting)
for s = 1:4
    curve_to_plot = results_curves{s}; 
    
    style = scenarios{s, 2};
    
    % X 軸轉換
    X = (1:length(curve_to_plot)) * pop_size;
    
    % 畫圖 (線寬加粗到 2.0，這樣縮小看才清楚
    semilogy(X, curve_to_plot, style, 'LineWidth', 1.5); 
end

% 3. 座標軸美化 (Axes Customization)
xlim([0, 300000]); 
ylim([1e0, 1e10]);

% X 軸刻度：只顯示 0, 0.5, 1... 簡潔有力
xticks([0, 0.5e5, 1e5, 1.5e5, 2e5, 2.5e5, 3e5]);
xticklabels({'0', '0.5', '1', '1.5', '2', '2.5', '3'});


set(gca, 'YScale', 'log');
ylim([1e-2, 1e6]); 
yticks([1e-2 1e0 1e2 1e4 1e6]);
% 使用 LaTeX 格式顯示指數
set(gca, 'YTickLabel', {'$10^{-2}$','$10^{0}$','$10^{2}$','$10^{4}$','$10^{6}$'}, 'TickLabelInterpreter', 'latex');

set(gca, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 14, ...
    'LineWidth', 1.2, ...           % 座標軸框線加粗
    'TickLabelInterpreter', 'latex', ... % 刻度標籤使用 LaTeX
    'YScale', 'log', ...
    'Box', 'on');

% X軸與Y軸標籤
xlabel('FEs', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 16); 
ylabel('Fitness', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 16);

text(1.02, -0.06, '$\times 10^5$', ...
    'Units', 'normalized', ...
    'Interpreter', 'latex', ...
    'FontSize', 14, ...
    'FontName', 'Times New Roman');

% 6. 圖例設定 (Legend)
lgd = legend(scenarios(:,3), 'Location', 'northeast');
set(lgd, ...
    'Box', 'off', ...      
    'Interpreter', 'latex', ...     
    'FontName', 'Times New Roman', ...
    'FontSize', 14);                

hold off;

fprintf('繪圖完成！請檢查字體是否滿意。\n');
