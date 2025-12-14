%opensource
clear all; clc;
warning('off');

% --- 2. 實驗參數設定 (黃金組合) ---
Runs = 51;          
D = 30;             
Xmin = -100; Xmax = 100;
fes_max = 10000 * D; 

pop_size = 60;   
Max_Gen = ceil(fes_max / pop_size);

fhd = str2func('cec17_func');
jingdu = 0;
func_num = 1; 

% --- 3. 定義情境 (補回 Fixed_p) ---
% 格式: {Fixed_c1, LineStyle, Label, Color, Fixed_p}
% 【關鍵修正 B】虛線組 p=0.5 (產生肩膀/浮起效果), 紅線 p=[] (Fuzzy)
scenarios = {
    0.2, '--', 'c_{i,1} = 0.2',     [0 1 0], 0.5;
    0.5, '-.', 'c_{i,1} = 0.5',     [0 0 1], 0.5;
    0.8, ':',  'c_{i,1} = 0.8',     [0 0 0], 0.5;
    [],  '-',  'c_{i,1} in MFCPSO', [1 0 0], 0.5
};

representative_curves = cell(1, 4); 

fprintf('=== 開始執行 51 次實驗 (取中位數) ===\n');

for s = 1:4
    % 讀取參數
    fixed_c1_val = scenarios{s, 1};
    label        = scenarios{s, 3};
    fixed_p_val  = scenarios{s, 5}; % 讀取 Fixed_p
    
    fprintf('正在運行情境: %s (p=%s)...\n', label, num2str(fixed_p_val));
    
    all_runs_matrix = zeros(Runs, Max_Gen);
    
    % --- 平行運算 51 次 ---
    parfor r = 1:Runs
        % 【關鍵修正 C】呼叫時傳入 fixed_p_val
        [~, curve] = Gaussian_rule(jingdu, func_num, fhd, D, pop_size, Max_Gen, fes_max, Xmin, Xmax, fixed_c1_val, fixed_p_val, func_num);
        
        % 確保長度一致 (補齊或裁切)
        len = length(curve);
        if len < Max_Gen
            curve = [curve, repmat(curve(end), 1, Max_Gen - len)];
        elseif len > Max_Gen
            curve = curve(1:Max_Gen);
        end
        
        all_runs_matrix(r, :) = curve;
    end
    
    % --- 計算中位數 (Median Run) ---
    % 邏輯：根據"最終 Fitness"排序，選中間那一次的完整曲線
    % 這比單純計算每一代的平均值更好，因為能保留曲線的"特徵形狀" (如階梯狀)
    final_values = all_runs_matrix(:, end);
    [~, sorted_idx] = sort(final_values);
    median_run_idx = sorted_idx(ceil(Runs/2)); 
    
    representative_curves{s} = all_runs_matrix(median_run_idx, :);
    
    fprintf('  >> 完成。中位數最終值: %e\n', final_values(median_run_idx));
end

% --- 4. 繪圖 ---
fprintf('正在繪圖...\n');
figure('Name', 'Median Convergence (51 Runs)', 'Color', 'w', 'Position', [100, 100, 700, 550]);
hold on;

for s = 1:4
    curve_data = representative_curves{s}; 
    style = scenarios{s, 2};
    color = scenarios{s, 4};
    
    % X 軸轉換: FEs
    X = (1:length(curve_data)) * pop_size;
    
    % 避免 log(0)
    curve_to_plot = max(curve_data, 1e-10); 
    
    semilogy(X, curve_to_plot, style, 'Color', color, 'LineWidth', 1.5);
end

% --- 5. 視覺美化 (跟論文一致) ---
xlim([0, 300000]); 
xticks([0, 0.5e5, 1e5, 1.5e5, 2e5, 2.5e5, 3e5]);
xticklabels({'0', '0.5', '1', '1.5', '2', '2.5', '3'});

% 手動加入 x10^5
text(1.02, -0.08, '\times 10^5', 'Units', 'normalized', 'FontSize', 11, 'FontName', 'Times New Roman');

xlabel('FEs', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Fitness', 'FontName', 'Times New Roman', 'FontSize', 12);

set(gca, 'YScale', 'log');
ylim([1e-2, 1e6]); 
yticks([1e-2 1e0 1e2 1e4 1e6]);
% 使用 LaTeX 格式顯示指數
set(gca, 'YTickLabel', {'$10^{-2}$','$10^{0}$','$10^{2}$','$10^{4}$','$10^{6}$'}, 'TickLabelInterpreter', 'latex');

grid off;
box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.0);

% Legend
legend(scenarios(:,3), 'Location', 'northeast', 'Box', 'off', 'FontName', 'Times New Roman');

hold off;
fprintf('全部完成！\n');