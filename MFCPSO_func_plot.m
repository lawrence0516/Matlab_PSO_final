function [gbestval, Convergence_Curve] = MFCPSO_func_plot(jingdu,func_num,fhd,Dimension,Particle_Number,Max_Gen,Max_FES,VRmin,VRmax, Fixed_c1, Fixed_p, func_id)
    
    % --- 參數設定 ---
    fbias = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,...
             1100,1200,1300,1400,1500, 1600,1700,1800,1900,2000,...
             2100,2200,2300,2400,2500, 2600,2700,2800,2900,3000];
    
    % Fuzzy 規則庫
    p_KB= [ 7, 6, 5, 4, 4, 4, 4; 7, 6, 5, 4, 4, 4, 3; 6, 5, 4, 4, 4, 4, 3; 6, 5, 4, 4, 4, 3, 2; 5, 4, 4, 4, 4, 3, 2; 5, 4, 4, 4, 3, 2, 1; 4, 4, 4, 4, 3, 2, 1 ];
    c1_KB= [1, 2, 3, 4, 4, 4, 4; 1, 2, 3, 4, 4, 4, 5; 2, 3, 4, 4, 4, 4, 5; 2, 3, 4, 4, 4, 5, 6; 3, 4, 4, 4, 4, 5, 6; 3, 4, 4, 4, 5, 6, 7; 4, 4, 4, 4, 5, 6, 7 ];
    
    % 初始化亂數
    try
        rand('state',sum(100*clock));
    catch
        rng('shuffle');
    end
    
    ps = Particle_Number;
    D = Dimension;
    
    if length(VRmin) == 1
        VRmin = repmat(VRmin,1,D);
        VRmax = repmat(VRmax,1,D);
    end
    lu = [VRmin; VRmax];
    mv = 0.2.*(VRmax-VRmin);
    Vmin = repmat(-mv,ps,1);
    Vmax = -Vmin;
    
    % 初始化
    pos = repmat(lu(1, :), ps, 1) + rand(ps, D) .* (repmat(lu(2, :) - lu(1, :), ps, 1));
    vel = Vmin+2.*Vmax.*rand(ps,D);
    
    fit = (feval(fhd,pos',func_id)-fbias(func_num))';
    fitcount = ps;
    
    pbest = pos;
    pbestval = fit;
    [gbestval, IX] = min(pbestval);
    gbest = pbest(IX,:);    
    
    % 收斂曲線紀錄
    Convergence_Curve = zeros(1, Max_Gen);
    Convergence_Curve(1) = gbestval; 
    
    % Fuzzy 初始化
    avgFit = mean(fit); avgNovel = 0;
    stdFit=avgFit-(1:Max_FES).*(avgFit./Max_FES);
    stdNovel=avgNovel-(1:Max_FES).*(avgNovel./Max_FES); 
    eFit_Last = 0; eNovel_Last = 0;
    
    % c1, c2 初始化
    c1 = 0.5 * ones(ps,1);
    c2 = 0.5 * ones(ps,1);
    
    c1_DeltaRange = [0.0, 1.0];
    p_ratio_DeltaRange = [0.05, 0.5];
    
    % Novelty 初始化 (呼叫乾淨版 Computing_Novelty)
    num_neighbors = ceil(ps/20);
    % 這裡傳入 0 (false) 或任何值都沒差，因為下面的函數已經被淨化，會完全忽略這個參數
    [novel, pos] = Computing_Novelty(ps, D, pbest, gbest, num_neighbors, lu, 0); 
    stdNovel=mean(novel)-(1:Max_FES).*(mean(novel)./Max_FES);
    
    iwt = 0.90-(1:Max_Gen).*(0.5./Max_Gen);
    gen = 1;
    
    % 主迴圈 
    while fitcount < Max_FES && gen < Max_Gen
        gen = gen + 1;
        
        gap = 1;
        maxFit = max(fit); minFit = min(fit);
        
        if mod(gen,gap)==0 && (maxFit <= stdFit(fitcount) || minFit >= stdFit(fitcount))   
            avgFit = mean(fit);  
            stdFit(fitcount:Max_FES) = avgFit-(fitcount:Max_FES).*(avgFit./Max_FES);
        end
        eFit = stdFit(fitcount) - fit;
        ecFit = eFit - eFit_Last;
        eFit_Range = [min(eFit), max(eFit)];
        ecFit_Range = [min(ecFit), max(ecFit)];
        eFit_Last = eFit; 
        
        maxNovel = max(novel); minNovel = min(novel);
        if mod(gen,gap)==0 && (maxNovel <= stdNovel(fitcount) || minNovel >= stdNovel(fitcount))   
            avgNovel = mean(novel);    
            stdNovel(fitcount:Max_FES) = avgNovel-(fitcount:Max_FES).*(avgNovel./Max_FES);      
        end
        eNovel = novel - stdNovel(fitcount);
        ecNovel = eNovel - eNovel_Last;
        eNovel_Range = [min(eNovel), max(eNovel)];
        ecNovel_Range = [min(ecNovel), max(ecNovel)];
        eNovel_Last = eNovel; 
        
        number_MF = 6;
        
        for i=1:ps
            % A. 計算 p
            p_fuzzy_val = Compute_c1c2(eFit(i), eFit_Range, eNovel(i), eNovel_Range, number_MF, p_KB, p_ratio_DeltaRange);
            
            if ~isempty(Fixed_p)
                p_ratio_fit(1,i) = Fixed_p;
            else
                p_ratio_fit(1,i) = p_fuzzy_val;
            end
            % B. 計算 c1
            c1_fuzzy_val = Compute_c1c2(eFit(i), eFit_Range, eNovel(i), eNovel_Range, number_MF, c1_KB, c1_DeltaRange);
            if ~isempty(Fixed_c1)
                c1(i,1) = Fixed_c1;
            else
                c1(i,1) = c1_fuzzy_val;
            end
            
            % C. 計算 c2
            c2(i,1) = 1.0 - c1(i,1);
        end
        
        % 計算 Maverick 比例
        if ~isempty(Fixed_p)
            % 若為固定參數，使用直觀邏輯 (1-p)
            p_ratio_novelty = 1.0 - p_ratio_fit; 
        else
            % 若為 Fuzzy，維持作者公式 (0.45-p)
            p_ratio_novelty = p_ratio_DeltaRange(2) - p_ratio_DeltaRange(1) - p_ratio_fit;
        end
        
        pbest_fit = pbest;
        pbest_novel = pbest;
        
        [~,indFit]=sort(pbestval,'ascend');
        [~, indNovel] = sort(novel, 'descend');
        
      for i=1:ps
            idx1 = ceil(ps * p_ratio_fit(1,i) * rand); 
            idx1 = max(1, min(ps, idx1));
            
            % 菁英 (Elite) = pbest (這是 PSO 標準)
            pbest_fit(i,:) = pbest(indFit(idx1), :); 
            
            % 計算 Maverick 索引
            pNP_novelty = p_ratio_novelty(1,i) * rand; 
            idx2 = ceil(ps * pNP_novelty); 
            idx2 = max(1, min(ps, idx2)); 
            
            % ★關鍵修改：Maverick 使用 pos (當前位置)★
            % 這符合論文 "Position" 的描述，而不是 pbest
            % 這會導致粒子跟隨一個不穩定的目標，造成震盪
            pbest_novel(i,:) = pos(indNovel(idx2), :); 
        end
        
        aa = repmat(c1,1,D).*rand(ps,D).*(pbest_fit-pos) + repmat(c2,1,D).*rand(ps,D).*(pbest_novel-pos);
        vel=iwt(gen)*vel+aa;
        vel=(vel>Vmax).*Vmax+(vel<=Vmax).*vel;
        vel=(vel<Vmin).*Vmin+(vel>=Vmin).*vel;
        pos=pos+vel;
        pos=(pos>VRmax).*VRmax+(pos<=VRmax).*pos; 
        pos=(pos<VRmin).*VRmin+(pos>=VRmin).*pos;
        
        fit = (feval(fhd,pos',func_id)-fbias(func_num))';
        fitcount = fitcount + ps;
        
        improved=(pbestval>fit);
        temp=repmat(improved,1,D);
        pbest=temp.*pos+(1-temp).*pbest;
        pbestval=improved.*fit+(1-improved).*pbestval;
        [gbestval, gbestid] = min(pbestval);
        gbest = pbest(gbestid,:);
        
        % 呼叫乾淨版 Novelty (無雜訊)
        [novel, pos] = Computing_Novelty(ps, D, pos, gbest, 2, lu, 0); 
        
        Convergence_Curve(gen) = gbestval;
        
    end
    
    Convergence_Curve = Convergence_Curve(1:gen);
end 

% === 完全淨化版 Computing_Novelty (No Noise, No Update) ===
% 即使上層傳入 enable_noise，這裡也完全忽略它，確保純淨
function [Novelty,pos] = Computing_Novelty(NP, D, pos, gbest, num_neighbors, lu, ~)
     for i = 1:NP
        % 1. 計算距離
        distance = zeros(1, NP);
        for j = 1:NP
           distance(1,j) = norm(pos(i,:)-pos(j,:));
        end        
        
        % 2. 正常過濾：只刪除完全重疊(距離為0)的
        distance(distance == 0) = []; 
        
        % 3. 純粹計算 Novelty，不做任何雜訊注入
        % 移除了 normrnd，也移除了 is_stuck 判斷
        if ~isempty(distance)
            [Novel_val, ~] = sort(distance,'ascend');
            Novelty(i,1) = mean(Novel_val(1:min(num_neighbors,length(Novel_val))));
        else
            Novelty(i,1) = 0;
        end
     end
end

% --- 以下 Fuzzy 函數保持不變 ---
function ef = Compute_MF(e, Range, number_MF)
    f = (Range(2)-Range(1)) / number_MF;   
    for i=1:number_MF+1
        r(i) = Range(1) + f*(i-1);  
    end 
    ef(1)=FTraL(e,r(1),r(3));        
    ef(7)=FTraR(e,r(5),r(7));        
    ef(2)=FTri(e,r(1),r(2),r(4));
    ef(3)=FTri(e,r(1),r(3),r(5));
    ef(4)=FTri(e,r(2),r(4),r(6));
    ef(5)=FTri(e,r(3),r(5),r(7));
    ef(6)=FTri(e,r(4),r(6),r(7));
end
function c_Exact = DeFuzzy(MaxX, MaxY, form, KB, OutRange)
    f = (OutRange(2)-OutRange(1)) / 6;    
    for i=1:7
        e(i) = OutRange(1) + f*(i-1);   
    end 
    c_fuzzy = form(MaxX,MaxY);
    tmp = KB(MaxX,MaxY);
    if tmp == 1
        c_Exact = uFTraL(c_fuzzy,e(1),e(3));
    elseif tmp == 2 
		c_Exact = uFTri(c_fuzzy,e(1),e(2),e(4));
    elseif tmp == 3 
		c_Exact = uFTri(c_fuzzy,e(1),e(3),e(5));
    elseif tmp == 4 
		c_Exact = uFTri(c_fuzzy,e(2),e(4),e(6));
    elseif tmp == 5 
		c_Exact = uFTri(c_fuzzy,e(3),e(5),e(7));
    elseif tmp == 6 
		c_Exact = uFTri(c_fuzzy,e(4),e(6),e(7));
    elseif tmp == 7 
		c_Exact = uFTraR(c_fuzzy,e(5),e(7));
    end
end
function c = Compute_c1c2(e, e_Range, ec, ec_Range, number_MF, c_KB, c_Range)
    ef = Compute_MF(e,e_Range,number_MF);
    ecf = Compute_MF(ec,ec_Range,number_MF);
    for j=1:7	
        for k=1:7		
            form(j,k)=fAND(ef(1,j),ecf(1,k));
        end
    end
    MaxX=1; MaxY=1;
    for j=1:7	
        for k=1:7	
           if  form(MaxX,MaxY)<form(j,k)
                MaxX=j;
                MaxY=k;
           end
        end
    end
    c = DeFuzzy(MaxX, MaxY, form, c_KB, c_Range);   
end
function result = FTri(  x, a, b, c)  
	if x<=a
		result = 0;
    elseif a<x && x<=b
		result = (x-a)/(b-a);
    elseif((b<x)&&(x<=c))
		result = (c-x)/(c-b);
    elseif(x>c)
		result = 0;
	else
		result = 0;       
	end
end
function result = FTraL(  x, a, b)  
	if x<=a
		result = 1;
    elseif a<x && x<=b
		result = (b-x)/(b-a);
    elseif x>b
		result = 0;    
	else
		result = 0;       
	end
end
function result = FTraR(  x, a, b)  
	if x<=a
		result = 0;
    elseif a<x && x<=b
		result = (x-a)/(b-a);
    elseif x>=b
		result = 1;    
	else
		result = 1;       
	end
end
function result = uFTri( x, a, b, c)  
	z = (b-a)*x+a;
	y = c-(c-b)*x;
	result = (y+z)/2; 
end
function result = uFTraL(  x, a, b)  	
	result = b-(b-a)*x; 
end
function result = uFTraR(  x, a, b)  	
	result = (b-a)*x + a; 
end
function result = fAND(a, b)  
    if a<b
        result = a;
    else
        result = b;
    end
end
function result = fOR(a, b)  
    if a<b
        result = b;
    else
        result = a;
    end
end