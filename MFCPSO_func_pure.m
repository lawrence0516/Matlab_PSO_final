function [gbestval, Convergence_Curve] = MFCPSO_func_pure(jingdu,func_num,fhd,Dimension,Particle_Number,Max_Gen,Max_FES,VRmin,VRmax, Fixed_c1, Fixed_p, func_id)
    % 保留核心演算法與 Fuzzy 控制
    

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
    
    % Novelty 初始化
    num_neighbors = ceil(ps/20);
    [novel, pos] = Computing_Novelty(ps, D, pbest, gbest, num_neighbors, lu, 1);
    stdNovel=mean(novel)-(1:Max_FES).*(mean(novel)./Max_FES);
    
    iwt = 0.90-(1:Max_Gen).*(0.5./Max_Gen);
    gen = 1;
    
    % 主迴圈 
    while fitcount < Max_FES && gen < Max_Gen
        gen = gen + 1;
        
        %% === 核心邏輯：Fuzzy 計算 ===
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
            p_ratio_fit(1,i) = Compute_c1c2(eFit(i), eFit_Range, eNovel(i), eNovel_Range, number_MF, p_KB, p_ratio_DeltaRange);
            % B. 計算 c1
            c1_fuzzy_val = Compute_c1c2(eFit(i), eFit_Range, eNovel(i), eNovel_Range, number_MF, c1_KB, c1_DeltaRange);
            % C. 決定 c1
            if ~isempty(Fixed_c1)
                c1(i,1) = Fixed_c1;
            else
                c1(i,1) = c1_fuzzy_val;
            end
            % D. 計算 c2
            c2(i,1) = 1.0 - c1(i,1);
        end
        
        p_ratio_novelty = 1 - p_ratio_fit; 
        
        pbest_fit = pbest;
        pbest_novel = pbest;
        
        %% === 粒子更新 ===
        [~,indFit]=sort(pbestval,'ascend');
        [~, indNovel] = sort(novel, 'descend');
        
        for i=1:ps
            idx1 = ceil(ps * p_ratio_fit(1,i) * rand); 
            idx1 = max(1, min(ps, idx1));
            pbest_fit(i,:) = pbest(indFit(idx1), :);
            
            idx2 = ceil(ps * p_ratio_novelty(1,i) * rand);
            idx2 = max(1, min(ps, idx2));
            pbest_novel(i,:) = pbest(indNovel(idx2), :);
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
        
        [novel, pos] = Computing_Novelty(ps, D, pos, gbest, 2, lu, 1);
        
        %% === 紀錄收斂曲線 ===
        Convergence_Curve(gen) = gbestval;
        
    end
    
    Convergence_Curve = Convergence_Curve(1:gen);

end 

% === 以下是輔助函數 ===
function [Novelty,pos] = Computing_Novelty(NP, D, pos, gbest, num_neighbors, lu, type)
     for i = 1:NP
        for j = 1:NP
           distance(1,j) = norm(pos(i,:)-pos(j,:));
        end        
        distance(distance==0)=[];
        if ~isempty(distance)
            [Novel_val, ~] = sort(distance,'ascend');
            Novelty(i,1) = mean(Novel_val(1:min(num_neighbors,length(Novel_val))));
        else
            Novelty(i,1) = 0;
        end
     end
end

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

function result = FZ(  x, a, b)  
	if x<=a
		result = 1;
    elseif a<x && x<=(a+b)/2
		result = 1-2*((x-a)/(b-a))*((x-a)/(b-a));
    elseif(((a+b)/2<x)&&(x<=b))
		result = 2*((x-b)/(b-a))*((x-b)/(b-a));
    elseif(x>b)
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

function result = uFZ(  x, a, b)  
	if x>=0.5
		result = sqrt((1-x)/2)*(b-a) + a;
    elseif x<0.5
        result = b - sqrt((x)/2)*(b-a);
	end
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