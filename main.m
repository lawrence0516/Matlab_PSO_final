clear all
clc
warning('off')

% 16 workers---
desired_workers = 16; 
current_pool = gcp('nocreate'); 
if isempty(current_pool)
    parpool('local', desired_workers);
elseif current_pool.NumWorkers ~= desired_workers
    delete(current_pool);
    parpool('local', desired_workers);
end
% --------------------------------------------
% ¡¾New¡¿Confirm and Print Message
p = gcp; % Get current pool object
if p.NumWorkers == desired_workers
    fprintf('Successfully Connected to parallel pool with %d workers.\n', p.NumWorkers);
else
    fprintf('Warning: Connected to %d workers (Target was %d).\n', p.NumWorkers, desired_workers);
end
% --------------------------------------------
% Dimension Setting 
D = 100; 
Xmin = -100;
Xmax = 100;
fes_max = 10000 * D; 
runtimes = 51; % Total runs 51 times 

fhd = str2func('cec17_func');
fbias = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,...
         1100,1200,1300,1400,1500, 1600,1700,1800,1900,2000,...
         2100,2200,2300,2400,2500, 2600,2700,2800,2900,3000];
jingdu = 0;
funset = [1:30];

% Define Output Filenames
filename_data = sprintf('CEC2017_original_Result_D%d.xlsx', D);       % Detailed Data (All Runs)
filename_summary = sprintf('CEC2017_original_Summary_D%d.xlsx', D); % Summary Report (Mean & Std)

% Initialize Tables
All_Data_Table = []; 
Summary_Table = [];

for fun = 1:length(funset)
    func_num = funset(fun);
    fprintf('[%dD] Running CEC2017 Function F%d ...\n', D, func_num);
    
    % Temporary storage for parfor
    temp_fbest = zeros(runtimes, 1);
    
    % --- Parallel Computing Loop ---
    parfor runs = 1:runtimes
        % Population Size = 2 * Dimension (Paper Setting)
        current_pop_size = 2 * D;   
        
        current_iter_max = ceil(fes_max / current_pop_size);
    
        % Execute Algorithm
        [gbest, gbestval, FES] = MFCPSO_func_PaperLogic(jingdu, func_num, fhd, D, current_pop_size, current_iter_max, fes_max, Xmin, Xmax, func_num);
        
        % Record Result
        temp_fbest(runs) = gbestval;
        
        % Output to Command Window (English)
        fprintf('Run %d: Best Fitness = %1.4e\n', runs, gbestval);
    end
    % -------------------------------
    
    % --- Statistics ---
    f_mean = mean(temp_fbest);
    f_std = std(temp_fbest);
    
    fprintf('\nFunction F%d Result:\nMean = %1.4e\nStd  = %1.4e\n', func_num, f_mean, f_std);
    fprintf('--------------------------------------------------\n');
    
    % --- 1. Prepare Detailed Data Table ---
    Func_Col = repmat(func_num, runtimes, 1);      
    Run_Col = (1:runtimes)';
    
    % Convert numbers to Scientific Notation Strings
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
    
    % --- 2. Prepare Summary Table (F1-F30) ---
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

    % --- Real-time Save ---
    writetable(All_Data_Table, filename_data);
    writetable(Summary_Table, filename_summary);
end

fprintf('All calculations complete!\n');
fprintf('1. Detailed Data saved to: %s\n', filename_data);
fprintf('2. Summary Report saved to: %s\n', filename_summary);