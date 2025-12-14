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
func_list = [1]; % 注意：論文的 F1 應該對應 func_num=1 (除非你的 cec17_func 索引不同)

% 定義情境： {Fixed_c1, LineStyle, Label, Fixed_p}
% 【關鍵修改】新增第四欄：Fixed_p
scenarios = {
    0.2, '--g', '$c_{i,1} = 0.2$', 0.25 ;  % 虛線組：
    0.5, '-.b', '$c_{i,1} = 0.5$', 0.25 ;  % 虛線組：
    0.8, ':k',  '$c_{i,1} = 0.8$', 0.25 ;  % 虛線組：
    [],  '-r',  '$c_{i,1}\ in\ MFCPSO$', 0.25  ;  % 紅線組：Fuzzy p (給空值)
};

results_curves = cell(1, 4); 
fprintf('開始運算\n');
tic; 

% ★★★ 第一階段：平行運算 ★★★
parfor s = 1:4
    % 取得參數
    fixed_c1_val = scenarios{s, 1};
    % label_name = scenarios{s, 3};
    fixed_p_val  = scenarios{s, 4}; % 【關鍵修改】讀取 Fixed_p
    
    % ★★★ 關鍵修改：呼叫時傳入 Fixed_p ★★★
    % 請確認你的 MFCPSO_func_pure 輸入參數順序是否為：
    % (..., Fixed_c1, Fixed_p, func_id)
    [~, curve] = MFCPSO_func_PaperLogic(jingdu, func_list(1), fhd, D, pop_size, Max_Gen, fes_max, Xmin, Xmax, fixed_c1_val, fixed_p_val, func_list(1));
    
    % 數據截斷處理
    curve = max(curve, 1e-8); % 避免 log(0) 或過小
    
    results_curves{s} = curve;
end

elapsed_time = toc;
fprintf('=== 計算完成！總耗時：%.2f 秒 ===\n', elapsed_time);

% ... (前面的運算部分不用動) ...

% ★★★ 第二階段：統一繪圖 (Visual Optimization) ★★★
fprintf('正在繪製圖表...\n');

% 1. 設定畫布大小 (稍微加大一點，讓字體不擁擠)
figure('Name', 'Fig 4 Reproduction', 'Color', 'w', 'Position', [150, 150, 700, 550]); 
hold on;

% 2. 畫線 (Line Plotting)
for s = 1:4
    % 這裡要對應你上面儲存的變數名稱
    % 如果你是跑 51 run 取中位數，變數可能是 representative_curves
    % 如果你是跑 1 run，變數可能是 results_curves
    % 這裡預設用你程式碼最後寫的 results_curves，若不同請自行修改
    curve_to_plot = results_curves{s}; 
    
    style = scenarios{s, 2};
    
    % X 軸轉換
    X = (1:length(curve_to_plot)) * pop_size;
    
    % 畫圖 (線寬加粗到 2.0，這樣縮小看才清楚)
    semilogy(X, curve_to_plot, style, 'LineWidth', 2.0); 
end

% 3. 座標軸美化 (Axes Customization)
xlim([0, 300000]); 
ylim([1e-2, 1e6]);

% X 軸刻度：只顯示 0, 0.5, 1... 簡潔有力
xticks([0, 0.5e5, 1e5, 1.5e5, 2e5, 2.5e5, 3e5]);
xticklabels({'0', '0.5', '1', '1.5', '2', '2.5', '3'});

% Y 軸刻度：設定為 10^0, 10^5, 10^10
yticks([1e-2 1e0 1e2 1e4 1e6]);
yticklabels({'$10^{-2}$', '$10^{0}$', '$10^{2}$', '$10^{4}$','$10^{6}$'}); % 強制用 LaTeX 數學格式

% 4. 全局字體設定 (Global Font Settings)
% 這行最重要！把所有座標軸文字改成 Times New Roman，大小設為 14
set(gca, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 14, ...
    'LineWidth', 1.2, ...           % 座標軸框線加粗
    'TickLabelInterpreter', 'latex', ... % 刻度標籤使用 LaTeX
    'YScale', 'log', ...
    'Box', 'on');

% 5. 標籤與特殊文字 (Labels & Text)
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