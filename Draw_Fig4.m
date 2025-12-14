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
func_list = [23];

scenarios = {
    0.2, '--g', '$c_{i,1} = 0.2$', 0.5 ;  % 虛線組：
    0.5, '-.b', '$c_{i,1} = 0.5$', 0.5 ;  % 虛線組：
    0.8, ':k',  '$c_{i,1} = 0.8$', 0.5 ;  % 虛線組：
    [],  '-r',  '$c_{i,1}\ in\ MFCPSO$', 0.5  ;  % 紅線組：Fuzzy p (給空值)
};

results_curves = cell(1, 4); 
fprintf('開始運算\n');
tic; 

parfor s = 1:4
    fixed_c1_val = scenarios{s, 1};
    % label_name = scenarios{s, 3};
    fixed_p_val  = scenarios{s, 4};
    
    [~, curve] = MFCPSO_func_plot(jingdu, func_list(1), fhd, D, pop_size, Max_Gen, fes_max, Xmin, Xmax, fixed_c1_val, fixed_p_val, func_list(1));
    
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
    
    semilogy(X, curve_to_plot, style, 'LineWidth', 1.5); 
end

xlim([0, 300000]); 
ylim([1e0, 1e10]);

xticks([0, 0.5e5, 1e5, 1.5e5, 2e5, 2.5e5, 3e5]);
xticklabels({'0', '0.5', '1', '1.5', '2', '2.5', '3'});


set(gca, 'YScale', 'log');
ylim([300, 1000]); 
yticks([300 400 500 600 700 800 900 1000]);
set(gca, 'YTickLabel', {'300','400','500','600','700','800','900','1000'}, 'TickLabelInterpreter', 'latex');

set(gca, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 14, ...
    'LineWidth', 1.2, ...       
    'TickLabelInterpreter', 'latex', ... 
    'YScale', 'log', ...
    'Box', 'on');

xlabel('FEs', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 16); 
ylabel('Fitness', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 16);

text(1.02, -0.06, '$\times 10^5$', ...
    'Units', 'normalized', ...
    'Interpreter', 'latex', ...
    'FontSize', 14, ...
    'FontName', 'Times New Roman');

lgd = legend(scenarios(:,3), 'Location', 'northeast');
set(lgd, ...
    'Box', 'off', ...               
    'Interpreter', 'latex', ...    
    'FontName', 'Times New Roman', ...
    'FontSize', 14);              

hold off;

fprintf('繪圖完成！請檢查字體是否滿意。\n');
