clear all; clc;
warning('off');

if isempty(gcp('nocreate'))
    parpool; 
end

D = 30;             
Xmin = -100; Xmax = 100;
fes_max = 10000 * D; 
pop_size = 60;    
Max_Gen = ceil(fes_max / pop_size);
fhd = str2func('cec17_func');
jingdu = 0;
func_list = [1];


scenarios = {
    0.2, '--g', '$c_{i,1} = 0.2$', 0.5 ;  
    0.5, '-.b', '$c_{i,1} = 0.5$', 0.5 ; 
    0.8, ':k',  '$c_{i,1} = 0.8$', 0.5 ; 
    [],  '-r',  '$c_{i,1}\ in\ MFCPSO$', 0.5  ;  
};

results_curves = cell(1, 4); 
fprintf('開始運算\n');
tic; 

parfor s = 1:4
    fixed_c1_val = scenarios{s, 1};
    fixed_p_val  = scenarios{s, 4}; 
    
    [~, curve] = MFCPSO_func_PaperLogic(jingdu, func_list(1), fhd, D, pop_size, Max_Gen, fes_max, Xmin, Xmax, fixed_c1_val, fixed_p_val, func_list(1));

    curve = max(curve, 1e-8); % 避免 log(0) 或過小
    
    results_curves{s} = curve;
end

elapsed_time = toc;
fprintf('=== 計算完成！總耗時：%.2f 秒 ===\n', elapsed_time);

fprintf('正在繪製圖表...\n');

figure('Name', 'Fig 4 Reproduction', 'Color', 'w', 'Position', [150, 150, 700, 550]); 
hold on;

for s = 1:4
    curve_to_plot = results_curves{s}; 
    
    style = scenarios{s, 2};
    
    X = (1:length(curve_to_plot)) * pop_size;
    
    % 畫圖 (線寬加粗到 2.0，這樣縮小看才清楚)
    semilogy(X, curve_to_plot, style, 'LineWidth', 2.0); 
end

% 3. 座標軸美化 (Axes Customization)
xlim([0, 300000]); 
ylim([1e-2, 1e6]);

xticks([0, 0.5e5, 1e5, 1.5e5, 2e5, 2.5e5, 3e5]);
xticklabels({'0', '0.5', '1', '1.5', '2', '2.5', '3'});
yticks([1e-2 1e0 1e2 1e4 1e6]);
yticklabels({'$10^{-2}$', '$10^{0}$', '$10^{2}$', '$10^{4}$','$10^{6}$'}); % 強制用 LaTeX 數學格式

set(gca, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 14, ...
    'LineWidth', 1.2, ...          
    'TickLabelInterpreter', 'latex', ... % 刻度標籤使用 LaTeX
    'YScale', 'log', ...
    'Box', 'on');

% X軸與Y軸標籤
xlabel('FEs', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 16); 
ylabel('Fitness', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 16);

% 手動加入右下角的 x10^5 (位置微調過)
text(1.02, -0.06, '$\times 10^5$', ...
    'Units', 'normalized', ...
    'Interpreter', 'latex', ...
    'FontSize', 14, ...
    'FontName', 'Times New Roman');

% 6. 圖例設定 (Legend)
lgd = legend(scenarios(:,3), 'Location', 'northeast');
set(lgd, ...
    'Box', 'off', ...                % 去除圖例邊框
    'Interpreter', 'latex', ...      % 圖例文字使用 LaTeX (這樣 c_{i,1} 下標才會漂亮)
    'FontName', 'Times New Roman', ...
    'FontSize', 14);                 % 圖例字體加大

hold off;

fprintf('繪圖完成！請檢查字體是否滿意。\n');
