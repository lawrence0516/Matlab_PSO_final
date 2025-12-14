clear all; clc;
warning('off');

% --- 啟動平行運算 (強制開啟 16 核心火力全開) ---
% 1. 檢查當前是否有已開啟的 Pool
current_pool = gcp('nocreate'); 

% 2. 判斷邏輯：
%    - 如果沒開 -> 開啟 16 個
%    - 如果有開，但數量不是 16 個 -> 關掉重開 16 個
%    - 如果已經開了 16 個 -> 繼續用
if isempty(current_pool)
    parpool('local', 16); 
elseif current_pool.NumWorkers ~= 16
    delete(current_pool); % 刪除舊的 (例如只開了 4 個的)
    parpool('local', 16); % 重新建立 16 個 workers
end

fprintf('平行運算已就緒，目前 Worker 數量: %d\n', 16);

% --- 實驗設定 ---
Runs = 51;          % 設定跑 51 次
D = 30;             
Xmin = -100; Xmax = 100;
fes_max = 10000 * D; 
pop_size = 60;    
Max_Gen = ceil(fes_max / pop_size);
fhd = str2func('cec17_func');
jingdu = 0;
func_list = [2]; 

% --- 定義情境 Scenarios ---
% 格式：{Fixed_c1, LineStyle, Label, Color, Fixed_p}
% 我把顏色獨立出來 (第四欄)，這樣畫圖控制更精準
scenarios = {
    0.2, '--', '$c_{i,1} = 0.2$',     [0 1 0], 0.5;   % 綠色虛線
    0.5, '-.', '$c_{i,1} = 0.5$',     [0 0 1], 0.5;   % 藍色點劃線
    0.8, ':',  '$c_{i,1} = 0.8$',     [0 0 0], 0.5;   % 黑色點線
    [],  '-',  '$c_{i,1}\ in\ MFCPSO$', [1 0 0], []   % 紅色實線
};

% 用來儲存最後要畫的那 4 條「中位數曲線」
representative_curves = cell(1, 4); 

fprintf('開始執行 51 次實驗取中位數 (尋找第 26 名)...\n');
tic; 

% ★★★ 第一階段：針對 4 種情境，各跑 51 次 ★★★
for s = 1:4
    % 讀取參數
    fixed_c1_val = scenarios{s, 1};
    label_name   = scenarios{s, 3};
    fixed_p_val  = scenarios{s, 5}; 
    
    fprintf('正在執行情境 %d/4 : %s ...\n', s, label_name);
    
    % 暫存該情境的 51 條曲線
    all_runs_matrix = zeros(Runs, Max_Gen);
    
    % --- 平行運算 51 次 (把最耗時的部分平行化) ---
    parfor r = 1:Runs
        % 呼叫函式
        [~, curve] = MFCPSO_func_plot(jingdu, func_list(1), fhd, D, pop_size, Max_Gen, fes_max, Xmin, Xmax, fixed_c1_val, fixed_p_val, func_list(1));
        
        % 數據處理：避免 log(0)
        curve = max(curve, 1e-10);
        
        % 確保長度一致 (補齊或裁切，防止維度不對)
        len = length(curve);
        if len < Max_Gen
            % 如果提前結束，後面補最後一個值
            curve = [curve, repmat(curve(end), 1, Max_Gen - len)];
        elseif len > Max_Gen
            curve = curve(1:Max_Gen);
        end
        
        all_runs_matrix(r, :) = curve;
    end
    
    % --- 計算中位數 (Median Run) ---
    % 1. 抓出每一次的「最終收斂值」(最後一代)
    final_values = all_runs_matrix(:, end);
    
    % 2. 排序
    [~, sorted_idx] = sort(final_values);
    
    % 3. 找出中間名次的索引 (第26名)
    median_run_idx = sorted_idx(ceil(Runs/2)); 
    
    % 4. 儲存這條冠軍曲線
    representative_curves{s} = all_runs_matrix(median_run_idx, :);
    
    fprintf('  >> 完成！中位數結果為: %e\n', final_values(median_run_idx));
end

elapsed_time = toc;
fprintf('=== 全部計算完成！總耗時：%.2f 秒 ===\n', elapsed_time);


% ★★★ 第二階段：統一繪圖 (Visual Optimization) ★★★
fprintf('正在繪製圖表 (中位數版)...\n');

% 1. 設定畫布大小
figure('Name', 'Fig 4 Reproduction (Median)', 'Color', 'w', 'Position', [150, 150, 700, 550]); 
hold on;

% 2. 畫線
for s = 1:4
    curve_to_plot = representative_curves{s}; 
    style = scenarios{s, 2};
    line_color = scenarios{s, 4}; % 讀取 RGB 顏色
    
    % X 軸轉換
    X = (1:length(curve_to_plot)) * pop_size;
    
    % 畫圖 (使用指定的顏色和樣式，線寬 1.0)
    semilogy(X, curve_to_plot, style, 'Color', line_color, 'LineWidth', 2.0); 
end

% 3. 座標軸美化
xlim([0, 300000]); 
ylim([1e-2, 1e6]);

% X 軸刻度
xticks([0, 0.5e5, 1e5, 1.5e5, 2e5, 2.5e5, 3e5]);
xticklabels({'0', '0.5', '1', '1.5', '2', '2.5', '3'});

% Y 軸刻度：設定為 10^0, 10^5, 10^10
yticks([1e-2 1e0 1e2 1e4 1e6]);
yticklabels({'$10^{-2}$', '$10^{0}$', '$10^{2}$', '$10^{4}$','$10^{6}$'}); % 強制用 LaTeX 數學格式

% 4. 全局字體設定
set(gca, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 14, ...
    'LineWidth', 1.0, ...           
    'TickLabelInterpreter', 'latex', ... 
    'YScale', 'log', ...
    'Box', 'on');

% 5. 標籤與特殊文字
xlabel('FEs', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 16); 
ylabel('Fitness', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 16);

% 右下角的 x10^5
text(1.02, -0.06, '$\times 10^5$', ...
    'Units', 'normalized', ...
    'Interpreter', 'latex', ...
    'FontSize', 14, ...
    'FontName', 'Times New Roman');
% (B) ★新增：標註實驗次數 (Run 51 times)★
% 位置設在左下角 (x=0.05, y=0.08)，這裡通常是空白的
text(0.05, 0.08, '\textbf{(run 51 times)}', ... 
    'Units', 'normalized', ...       % 使用相對座標 (0~1)，不怕座標軸範圍改變
    'Interpreter', 'latex', ...      % 使用 LaTeX 渲染
    'FontName', 'Times New Roman', ...
    'FontSize', 15, ...              % 字體稍微大一點，清楚顯眼
    'FontWeight', 'bold', ...        % 加粗強調
    'Color', 'k');                   % 黑色文字

% 6. 圖例設定
lgd = legend(scenarios(:,3), 'Location', 'northeast');
set(lgd, ...
    'Box', 'off', ...                
    'Interpreter', 'latex', ...      
    'FontName', 'Times New Roman', ...
    'FontSize', 14);                 

hold off;
fprintf('繪圖完成！這是 51 次實驗的中位數結果。\n');