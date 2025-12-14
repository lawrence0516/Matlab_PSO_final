clear all; % 修正了原本的拼寫錯誤
clc;
warning('off');

% --- 1. 平行運算設定 (16 Workers) ---
desired_workers = 16; 
current_pool = gcp('nocreate'); 
if isempty(current_pool)
    parpool('local', desired_workers);
elseif current_pool.NumWorkers ~= desired_workers
    delete(current_pool);
    parpool('local', desired_workers);
end

p = gcp; 
if p.NumWorkers == desired_workers
    fprintf('Successfully Connected to parallel pool with %d workers.\n', p.NumWorkers);
else
    fprintf('Warning: Connected to %d workers (Target was %d).\n', p.NumWorkers, desired_workers);
end

% --- 2. 實驗參數設定 ---
D = 100; % 設定維度 (你可以改回 10, 30, 50)
Xmin = -100;
Xmax = 100;
fes_max = 10000 * D; 
runtimes = 51; % 跑 51 次
fhd = str2func('cec17_func');
fbias = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,...
         1100,1200,1300,1400,1500, 1600,1700,1800,1900,2000,...
         2100,2200,2300,2400,2500, 2600,2700,2800,2900,3000];
jingdu = 0;
funset = [1:30]; % 跑 F1 到 F30

% --- 【修改 1】設定輸出檔名 (加上 _PaperLogic 以區分) ---
filename_data = sprintf('CEC2017_Result_D%d_PaperLogic.xlsx', D);       
filename_summary = sprintf('CEC2017_Summary_D%d_PaperLogic.xlsx', D); 

% 初始化表格
All_Data_Table = []; 
Summary_Table = [];

for fun = 1:length(funset)
    func_num = funset(fun);
    fprintf('[%dD] Running CEC2017 Function F%d (Paper Logic Version)...\n', D, func_num);
    
    temp_fbest = zeros(runtimes, 1);
    
    % --- 平行運算迴圈 ---
    parfor runs = 1:runtimes
        % 粒子數設定 (依照論文 Table III: 2*D)
        current_pop_size = 2 * D;   
        current_iter_max = ceil(fes_max / current_pop_size);
    
        % --- 【修改 2】呼叫 "PaperLogic" 版本的函式 ---
        % 這裡使用的是你修改過 (pos + 無雜訊) 的程式碼
        [gbest, gbestval, FES] = MFCPSO_func_PaperLogic(jingdu, func_num, fhd, D, current_pop_size, current_iter_max, fes_max, Xmin, Xmax, func_num);
        
        temp_fbest(runs) = gbestval;
        
        % fprintf('Run %d: Best Fitness = %1.4e\n', runs, gbestval);
    end
    % -------------------------------
    
    % --- 統計計算 ---
    f_mean = mean(temp_fbest);
    f_std = std(temp_fbest);
    
    fprintf('\nFunction F%d Result (PaperLogic):\nMean = %1.4e\nStd  = %1.4e\n', func_num, f_mean, f_std);
    fprintf('--------------------------------------------------\n');
    
    % --- 儲存詳細數據 ---
    Func_Col = repmat(func_num, runtimes, 1);      
    Run_Col = (1:runtimes)';
    
    Best_Str = compose('%1.4e', temp_fbest);
    Mean_Str = compose('%1.4e', repmat(f_mean, runtimes, 1));
    Std_Str  = compose('%1.4e', repmat(f_std, runtimes, 1));
    
    Temp_Table = table(Func_Col, Run_Col, Best_Str, Mean_Str, Std_Str, ...
        'VariableNames', {'Function', 'Run', 'BestFitness', 'Mean', 'SD'});
    
    if isempty(All_Data_Table)
        All_Data_Table = Temp_Table;
    else
        All_Data_Table = [All_Data_Table; Temp_Table];
    end
    
    % --- 儲存總結數據 ---
    S_Func = func_num;
    S_Mean = compose('%1.4e', f_mean);
    S_Std  = compose('%1.4e', f_std);
    
    Temp_Summary = table(S_Func, S_Mean, S_Std, ...
        'VariableNames', {'Function', 'Mean', 'SD'});
        
    if isempty(Summary_Table)
        Summary_Table = Temp_Summary;
    else
        Summary_Table = [Summary_Table; Temp_Summary];
    end
    
    % --- 即時寫入 Excel ---
    writetable(All_Data_Table, filename_data);
    writetable(Summary_Table, filename_summary);
end

fprintf('All calculations complete!\n');
fprintf('1. Detailed Data saved to: %s\n', filename_data);
fprintf('2. Summary Report saved to: %s\n', filename_summary);