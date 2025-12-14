function [gbest,gbestval,fitcount,suc,suc_fes]= Gaussian_p_freedom(jingdu,func_num,fhd,Dimension,Particle_Number,Max_Gen,Max_FES,VRmin,VRmax,varargin)
%%
% 將 FuzzyPSO_func 中的結構體編程替換為常規方法。
%%
fbias=[100, 200, 300, 400, 500,...
       600, 700, 800, 900, 1000,...
       1100,1200,1300,1400,1500,...
       1600,1700,1800,1900,2000,...
       2100,2200,2300,2400,2500,...
       2600,2700,2800,2900,3000 ];

% 定義符號常量
%NB = sym('1');
%NM = sym('2');
%NS = sym('3');
% ZO = sym('4');
% PS = sym('5');
% PM = sym('6');
% PB = sym('7');

% 模糊規則，或知識庫
p_KB= [ 7, 6, 5, 4, 4, 4, 4;
	    7, 6, 5, 4, 4, 4, 3;
	    6, 5, 4, 4, 4, 4, 3;
	    6, 5, 4, 4, 4, 3, 2;
	    5, 4, 4, 4, 4, 3, 2;
	    5, 4, 4, 4, 3, 2, 1;
	    4, 4, 4, 4, 3, 2, 1 ];
    
c1_KB= [1, 2, 3, 4, 4, 4, 4;
	    1, 2, 3, 4, 4, 4, 5;
	    2, 3, 4, 4, 4, 4, 5;
	    2, 3, 4, 4, 4, 5, 6;
	    3, 4, 4, 4, 4, 5, 6;
	    3, 4, 4, 4, 5, 6, 7;
	    4, 4, 4, 4, 5, 6, 7 ];
    
c2_KB=[ 7, 7, 7, 7, 6, 5, 4;
	    7, 7, 7, 6, 5, 4, 4;
	    7, 7, 6, 5, 4, 4, 4;
	    7, 6, 5, 4, 4, 4, 3;
	    6, 5, 4, 4, 4, 3, 3;
	    5, 4, 4, 4, 3, 3, 2;
	    4, 4, 4, 3, 3, 2, 1 ];  % 第一版裡面使用的
    
c3_KB=[ 5, 4, 4, 4, 3, 2, 1;
        5, 5, 4, 4, 4, 3, 2;
        6, 5, 5, 4, 4, 4, 3;
        6, 6, 5, 5, 4, 4, 4;
        7, 6, 6, 5, 5, 4, 4;
        7, 7, 6, 6, 5, 5, 4;
        7, 7, 7, 6, 6, 5, 5 ];
    
    

rand('state',sum(100*clock));
me = Max_Gen;
ps = Particle_Number;
D = Dimension;

cc = [2.0 2.0];   %acceleration constants
iwt = 0.90-(1:me).*(0.5./me);

recorded = 0;  % 達到精度時紀錄相關信息
suc = 0;
suc_fes = 0;

if length(VRmin) == 1
    VRmin = repmat(VRmin,1,D);
    VRmax = repmat(VRmax,1,D);
end
lu = [VRmin; VRmax];

mv = 0.2.*(VRmax-VRmin);
Vmin = -mv;
Vmax = -Vmin;

fitcount = 0;
gen = 0;

%% Initialize the main population and velocity

% Initialize position and velocity
pos = repmat(lu(1, :), ps, 1) + rand(ps, D) .* (repmat(lu(2, :) - lu(1, :), ps, 1));
mv = 0.10*(lu(2,:) - lu(1,:));
Vmin = repmat(-mv,ps,1);
Vmax = -Vmin;
vel = Vmin+2.*Vmax.*rand(ps,D);

% Evaluate the population
fit = (feval(fhd,pos',varargin{:})-fbias(func_num))';
fitcount = fitcount + ps;
%gen = gen + 1;

% Initialize the pbest and the pbestval
pbest = pos;
pbestval = fit;

% Initialize the gbest and the gbestval
[val,IX] = min(pbestval);
gbestval = val;
gbest = pbest(IX,:);    
gbestrep=repmat(gbest,ps,1);%update the gbest

% Compute the average fitness and novelty of the population
avgFit = mean(fit);
maxFit = max(fit);
minFit = min(fit);

num_neighbors = ceil(ps/20);
Novelty_type = 1;
    % 1. 個體的 pbest 與最接近的 num_neighbors 個個體 pbest 的距離來表示 Novelty
    % 2. 個體的 pbest 與 gbest 的距離來表示 Novelty  
    % 3. 個體的 pbest 與 pbest 的中心點的距離來表示 Novelty
    % 4. 個體的 pbest 與 gbest 及 pbest 的中心點的距離來表示 Novelty  
[novel, pos] = Computing_Novelty(ps, D, pbest, gbest, num_neighbors, lu, Novelty_type); 
avgNovel = mean(novel);
maxNovel = max(novel);
minNovel = min(novel);

%% ============ 确定模糊控制的標準量 ================
change_type = 1;  

stdFit=avgFit-(1:Max_FES).*(avgFit./Max_FES);    % 利用 Fuzzy 進行控制時，fitness 在每代中的標準值 
stdNovel=avgNovel-(1:Max_FES).*(avgNovel./Max_FES);  % 利用 Fuzzy 進行控制時，novelty 在每代中的標準值

eFit_Last = stdFit(fitcount) - fit;
eNovel_Last = novel - stdNovel(fitcount);

%% Initial c1 and c2. Each particle has its own c1 and c2. 

c1 = ones(ps,1);   % c1 的標準量
c2 = ones(ps,1);   % c1 的標準量
%======== 這裡的 switch 是用作測試 ==============%
c1_DeltaRange = [0.0, 1.0];  % c1 的控制量的波動範圍   type_range = 1;
c2_DeltaRange = c1_DeltaRange; 
p_ratio_DeltaRange = [0.05, 0.5];  % c1 的控制量的波動範圍

%% ================== 以下為畫圖用 ====================%%   
% % ------------- 適應值變化 -------------------%
% old = 1; new = fitcount;
% fitness(old:new) = gbestval;
% maxfitness(old:new) = maxFit;
% minfitness(old:new) = minFit;
% maxnovelty(old:new) = maxNovel;
% minnovelty(old:new) = minNovel;
% old = new;   

%% ================== 畫圖結束 ====================%% 
p_val = 0; gbestval_initial = gbestval;   pos_initial = pos; vel_initial=vel;  pbestval_initial = pbestval;

 fitcount = ps;    gen = 0;  
 p_val = 4; % 0: p_val=0.2;  1: p_val=0.5;  3:p_val=0.8;  4:p_val controlled by Fuzzy
 old = 1; new = fitcount;
 pos = pos_initial; vel=vel_initial;  pbestval = pbestval_initial;
 gbestval = gbestval_initial; fitness(old:new) = gbestval;
while fitcount < Max_FES && gen < Max_Gen
    gen = gen + 1; 

    %% =================== PSO: Update position and velocity ================
    if gen == 1
       p_ratio_fit = 0.5*ones(1,ps);  % 利用族群適應值及多樣性的平均值來確定此參數          
    else
       for i=1:ps
           p_ratio_fit(1,i) = Compute_c1c2(eFit(i), eFit_Range, eNovel(i), eNovel_Range, number_MF, p_KB, p_ratio_DeltaRange);
       end
    end
    % 修改後的代碼 (加入 max 保護機制)
% 邏輯：原本的總額度是 0.45，如果你設 0.9，這裡會變負的。
% 我們用 max(0, ...) 強制讓它最小為 0，這樣就不會出錯。
p_ratio_novelty = max(0, p_ratio_DeltaRange(2) - p_ratio_DeltaRange(1) - p_ratio_fit);  

    %% =================== Find many solutions with higher fitness ==========  
   [~,indFit]=sort(pbestval,'ascend');   % 對適應值進行排序 
   [~, indNovel] = sort(novel, 'descend');   % 對新穎性進行排序 
   for i=1:ps
        pNP_fit(1,i) = p_ratio_fit(1,i)*rand; 
        index_fit(1,i) = ceil(ps * pNP_fit(1,i)); % select from [1, 2, 3, ..., pNP]
        index_fit(1,i) = max(1, index_fit(1,i)); % to avoid the problem that rand = 0 and thus ceil(rand) = 0        
        IX1 = indFit(index_fit(1,i));
        pbest_fit(i,:) = pbest(IX1, :);
        
        pNP_novelty(1,i) = p_ratio_novelty(1,i)*rand;
        index_novelty(1,i) = ceil(ps * pNP_novelty(1,i)); % select from [1, 2, 3, ..., pNP]
        index_novelty(1,i) = max(1, index_novelty(1,i)); % to avoid the problem that rand = 0 and thus ceil(rand) = 0          
        IX2 = indNovel(index_novelty(1,i));
        pbest_novel(i,:) = pbest(IX2, :);
   end
    %% =================== Find many solutions with higher novelty ==========
    % 利用模糊 c1,c2 的 fitness-novelty-based 策略來更新vel
    aa = repmat(c1,1,D).*rand(ps,D).*(pbest_fit-pos) + repmat(c2,1,D).*rand(ps,D).*(pbest_novel-pos);
    
    vel=iwt(gen)*vel+aa;
    vel=(vel>Vmax).*Vmax+(vel<=Vmax).*vel;
    vel=(vel<Vmin).*Vmin+(vel>=Vmin).*vel;
    pos=pos+vel;
    pos=(pos>VRmax).*VRmax+(pos<=VRmax).*pos; 
    pos=(pos<VRmin).*VRmin+(pos>=VRmin).*pos;
    fit = (feval(fhd,pos',varargin{:})-fbias(func_num))';
    fitcount = fitcount + ps;
     % 更新 pbest 及 pbestval
    improved=(pbestval>fit);     % individuals are improved or not
    temp=repmat(improved,1,D);
    pbest=temp.*pos+(1-temp).*pbest;
    pbestval=improved.*fit+(1-improved).*pbestval;      % update the pbest
    
    % 更新 gbest 及 gbestval
    [gbestval,gbestid] = min(pbestval);  % 注意方括號不能少
    gbest = pbest(gbestid,:); %initialize the gbest and the gbest's fitness value
    gbestrep=repmat(gbest,ps,1);%update the gbest
    
    num_neighbors = 2;     
    [novel, pos] = Computing_Novelty(ps, D, pos, gbest, num_neighbors, lu, Novelty_type); 

    %% ======================= Fuzzy 控制 ===========================%
   
    % ================= 計算各模糊控制的參數 ===================
    % 如果當前族群的適應值都處於標準控制量的一側，則改變標準控制量
    gap = 1; % 測試 gap = 10, 100, 500, 1000, 2000
    maxFit = max(fit); minFit = min(fit);
    if mod(gen,gap)==0 && (maxFit <= stdFit(fitcount) || minFit >= stdFit(fitcount))   
        avgFit = mean(fit);  
        stdFit(fitcount:Max_FES) = avgFit-(fitcount:Max_FES).*(avgFit./Max_FES);    % 利用 Fuzzy 進行控制時，fitness 在每代中的標準值       
    end
    eFit = stdFit(fitcount) - fit;  % 個體適應值與當前標準適應值的偏差
    ecFit = eFit - eFit_Last;       % 個體適應值與當前標準適應值的偏差的變化
    eFit_Range = [min(eFit), max(eFit)];
    ecFit_Range = [min(ecFit), max(ecFit)];
        
    % 如果當前除群最的多樣性都處於標準控制量的一側，則改變標準控制量
    maxNovel = max(novel); minNovel = min(novel);
    if mod(gen,gap)==0 && (maxNovel <= stdNovel(fitcount) || minNovel >= stdNovel(fitcount))   
        avgNovel = mean(novel);    
        stdNovel(fitcount:Max_FES) = avgNovel-(fitcount:Max_FES).*(avgNovel./Max_FES);  % 利用 Fuzzy 進行控制時，novelty 在每代中的標準值      
    end
    eNovel = novel - stdNovel(fitcount);  % 個體新穎度與當前標準新穎度的偏差
    ecNovel = eNovel - eNovel_Last;    % 個體新穎度與當前標準新穎度的偏差的變化
    eNovel_Range = [min(eNovel), max(eNovel)];
    ecNovel_Range = [min(ecNovel), max(ecNovel)];
        
    number_MF = 6;
    
    % ------------------ 模糊調整 c1,c2 -------------------- %
    for i=1:ps  
        c1(i,1) = Compute_c1c2(eFit(i), eFit_Range, eNovel(i), eNovel_Range, number_MF, c1_KB, c1_DeltaRange);
        c2(i,1) = c1_DeltaRange(1) + c1_DeltaRange(2) - c1(i,1);  %
    end

%% ================== 以下為畫圖用 ====================%%   
% % ------------- 適應值變化 -------------------%
%     new = fitcount;
%     fitness(old:new) = gbestval;
%     maxfitness(old:new) = maxFit;
%     minfitness(old:new) = minFit;
%     maxnovelty(old:new) = maxNovel;
%     minnovelty(old:new) = minNovel;
%     old = new;   

%% ================== 畫圖结束 ====================%% 

   %%%%%%%%%%%%%%%
     if gbestval <= jingdu && recorded == 0
         recorded = 1;
         suc = 1;
         suc_fes = fitcount;
     end
     %%%%%%%%%%%%%%%

    if fitcount>=Max_FES
        break;
    end
    
end



end


%========================= UpdatePop ==========================================
function [pos, fit, Novelty, pbest, pbestval, gbestval, gbest, vel] = UpdatePop(gbest, pbestfit, pbestnovel, w, c1, c2, pos, vel, pbest, pbestval, N, D, fhd, lu, Vmax, Vmin, func_num,fbias)
    
    c1 = repmat(c1, 1, D);
    c2 = repmat(c2, 1, D);
 
    aa = c1.*rand(N,D).*(pbestfit-pos) + c2.*rand(N,D).*(pbestnovel-pos);
    vel = w*vel+aa;
    vel = (vel>Vmax).*Vmax + (vel<=Vmax).*vel;
    vel = (vel<Vmin).*Vmin + (vel>=Vmin).*vel;
    pos = pos+vel;  
    if rand>0.5
        pos = (pos>lu(2,:)).*lu(2,:) + (pos<=lu(2,:)).*pos; 
        pos = (pos<lu(1,:)).*lu(1,:) + (pos>=lu(1,:)).*pos;
    else
        pos = ((pos>=lu(1,:))&(pos<=lu(2,:))).*pos ...
            + (pos<lu(1,:)).*(lu(1,:)+0.2.*(lu(2,:)-lu(1,:)).*rand(N,D)) + (pos>lu(2,:)).*(lu(2,:)-0.2.*(lu(2,:)-lu(1,:)).*rand(N,D));
    end
    
    fit(:,1) = feval(fhd,(pos)',func_num)-fbias(func_num);  

    % 更新 pbest 及 pbestval
    improved=(pbestval>fit);     % individuals are improved or not
    temp=repmat(improved,1,D);
    pbest=temp.*pos+(1-temp).*pbest;
    pbestval=improved.*fit+(1-improved).*pbestval;      % update the pbest
    
    % 更新 gbest 及 gbestval
    [gbestval,gbestid] = min(pbestval);  % 注意方括號不能少
    gbest = pbest(gbestid,:); %initialize the gbest and the gbest's fitness value

    num_neighbors = 2;  
    Novelty_type = 1;
%     pop(i).novel = novel;
    [Novelty,pos] = Computing_Novelty(N, D, pos, gbest, num_neighbors, lu, Novelty_type); 

end

%========================= Computing_Novelty ==========================================
function [Novelty,pos] = Computing_Novelty(NP, D, pos, gbest, num_neighbors, lu, type)
    switch type
    case 1
    % 個體的 pbest 與最鄰近的 num_neighbors 個個體 pbest 的距離來表示 Novelty
         for i = 1:NP
            for j = 1:NP
               distance(1,j) = norm(pos(i,:)-pos(j,:));
            end        

            distance(distance==0)=[]; % 刪掉 0 元素，即不考慮個體與其自身的距離
            if ~isempty(distance)
                [Novel_val, Novel_IX] = sort(distance,'ascend');

                if length(Novel_IX) == 1   % 當族群已經收斂時，對 pos 随機擾動
                    Novelty(i,1) = mean(Novel_val(Novel_IX(1)));
                    sigma = (lu(2)-lu(1))/1000;
                    pos(i,:) = pos(i,:) + normrnd(0,sigma,1,D);
                
                else
                    Novelty(i,1) = mean(Novel_val(Novel_IX(1:min(num_neighbors,length(Novel_IX)))));
                end               
            else
                Novelty(i,1) = 0;
            end
         end

    case 2
         % 個體的 pbest 與gbest 的距離來表示 Novelty     
         for i = 1:NP        
            distance(i) = norm(pos(i,:)-gbest(1,:));              
         end

         if sum(distance)==0  % 族群已收斂，對 pbest 進行擾動        
                sigma = (lu(2)-lu(1))/200;
                pos(i,:) = pos(i,:) + normrnd(0,sigma,1,D); 
                for i = 1:NP        
                    distance(i) = norm(pos(i,:)-gbest(1,:));              
                end
         end
         Novelty = distance';
            
     case 3
        % 個體的 pbest 與 pbest 的中心點的距離來表示 Novelty
         cbest = mean(pos);
         for i = 1:NP        
            distance(i) = norm(pos(i,:)-cbest(1,:));              
         end

         if sum(distance)==0  % 族群已收斂，對 pbest 進行擾動        
                sigma = (lu(2)-lu(1))/200;
                pos(i,:) = pos(i,:) + normrnd(0,sigma,1,D); 
                for i = 1:NP        
                    distance(i) = norm(pos(i,:)-gbest(1,:));              
                end
         end
         Novelty = distance'; 

      case 4
          % 個體的 pbest 與 gbest 及 pbest 的中心點的距離來表示 Novelty
          cbest = mean(pos);
          for i = 1:NP        
            distance_1(i) = norm(pos(i,:)-gbest(1,:));
            distance_2(i) = norm(pos(i,:)-cbest(1,:));
            distance(i) = distance_1(i)+distance_1(i);
          end

          if sum(distance)==0  % 除群已收斂，對 pbest 進行擾動        
                sigma = (lu(2)-lu(1))/200;
                pos(i,:) = pos(i,:) + normrnd(0,sigma,1,D); 
                for i = 1:NP        
                    distance(i) = norm(pos(i,:)-gbest(1,:));              
                end
          end
         Novelty = distance';

    end
            
end

%% ======================= 計算隸屬度 ===================================
function ef = Compute_MF(e, Range, number_MF)
    % 定義符號常量
    %     NB = sym('1');
    %     NM = sym('2');
    %     NS = sym('3');
    %     ZO = sym('4');
    %     PS = sym('5');
    %     PM = sym('6');
    %     PB = sym('7');
    % NB = 1; NM = 2; NS = 3; ZO = 4; PS = 5; PM = 6; PB = 7;
    f = (Range(2)-Range(1)) / number_MF;  % 輸入變量值的論域   
    for i=1:number_MF+1
        r(i) = Range(1) + f*(i-1);  % 
    end 
        
    ef(1)=FTraL(e,r(1),r(3));        
    ef(7)=FTraR(e,r(5),r(7));        
    ef(2)=FTri(e,r(1),r(2),r(4));
    ef(3)=FTri(e,r(1),r(3),r(5));
    ef(4)=FTri(e,r(2),r(4),r(6));
    ef(5)=FTri(e,r(3),r(5),r(7));
    ef(6)=FTri(e,r(4),r(6),r(7));
   
end

%% ======================= 解模糊 ===================================
function c_Exact = DeFuzzy(MaxX, MaxY, form, KB, OutRange)
     % 定義符號常量
    %     NB = sym('1');
    %     NM = sym('2');
    %     NS = sym('3');
    %     ZO = sym('4');
    %     PS = sym('5');
    %     PM = sym('6');
    %     PB = sym('7');
    % NB = 1; NM = 2; NS = 3; ZO = 4; PS = 5; PM = 6; PB = 7;
    
    f = (OutRange(2)-OutRange(1)) / 6;  % 輸入變量值的論域   
    for i=1:7
        e(i) = OutRange(1) + f*(i-1);  % 
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


%% ======================= 計算 c1 和 c2 的值 ==========================
function c = Compute_c1c2(e, e_Range, ec, ec_Range, number_MF, c_KB, c_Range)
                        
  % 計算 eFit 和 ecFit 隸屬度
    ef = Compute_MF(e,e_Range,number_MF);
    ecf = Compute_MF(ec,ec_Range,number_MF);

    % 計算隸屬度表，确定 e 和 ec 相關聯後表格各項隸屬度的值
    % NB = 1; NM = 2; NS = 3; ZO = 4; PS = 5; PM = 6; PB = 7;
    for j=1:7	
        for k=1:7		
            form(j,k)=fAND(ef(1,j),ecf(1,k));
        end
    end

    % 取出具有最大隸屬度的那一項
    MaxX=1; MaxY=1;
    for j=1:7	
        for k=1:7	
           if  form(MaxX,MaxY)<form(j,k)
                MaxX=j;
                MaxY=k;
           end
        end
    end

    % 進行模糊推理，並去模糊
    c = DeFuzzy(MaxX, MaxY, form, c_KB, c_Range);   
end


%% =================== 確定隸屬度 =========================
function result = FTri(  x, a, b, c)  % FuzzyTriangle - 三角形
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


function result = FZ(  x, a, b)  % FuzzyZ - Z形
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


function result = FTraL(  x, a, b)  % FuzzyTrapezoidLeft - 梯形左
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


function result = FTraR(  x, a, b)  % FuzzyTrapezoidLeft - 梯形右
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


%% ====================== 解模糊 =========================
function result = uFTri( x, a, b, c)  % FuzzyTriangle - 三角形
	z = (b-a)*x+a;
	y = c-(c-b)*x;
	result = (y+z)/2; 
end

function result = uFZ(  x, a, b)  % FuzzyZ - Z形
	if x>=0.5
		result = sqrt((1-x)/2)*(b-a) + a;
    elseif x<0.5
        result = b - sqrt((x)/2)*(b-a);
	end
end


function result = uFTraL(  x, a, b)  % FuzzyTrapezoidLeft - 梯形左	
	result = b-(b-a)*x; 
end


function result = uFTraR(  x, a, b)  % FuzzyTrapezoidLeft - 梯形右	
	result = (b-a)*x + a; 
end


%% ================== 求交、並集 =======================
function result = fAND(a, b)  % 交集
    if a<b
        result = a;
    else
        result = b;
    end
end

function result = fOR(a, b)  % 並集
    if a<b
        result = b;
    else
        result = a;
    end
end