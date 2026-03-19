%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
%  Rift-Drift Field Optimizer (RDFO) source codes 
%  
%  Developed in:	MATLAB (R2024a)
%  
%  Original paper:	Rift-Drift Field Optimizer: A 
%                   geophysics-inspired metaheuristic for global 
%                   and engineering optimization
%   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bestx, bestf, ConvergenceCurve] = RDFO(fun, nvars, lb, ub, options)

%% Initialization

% --- Population and Iteration settings  ---
if isfield(options, 'PopulationSize')
    N = options.PopulationSize;
elseif isfield(options, 'PopSize')
    N = options.PopSize;
else
    error('RDFO:MissingPopulationSize', 'options.PopulationSize or options.PopSize is required.');
end

if isfield(options, 'MaxIterations')
    T = options.MaxIterations;
elseif isfield(options, 'MaxFEs')
    T = max(1, ceil(options.MaxFEs / N));
else
    error('RDFO:MissingBudget', 'options.MaxFEs or options.MaxIterations is required.');
end

% --- Core RDFO Parameters ---
Pe = 0.95;               % Elite ratio for Latent Stress Field 
eta = 0.8;               % Learning rate for the Latent Stress Field 
Cmax = 265;              % Max stress tolerance for Fission Reinjection 

% Dynamic step sizes
alphaD_max = 0.9;        % Maximum drift step size 
alphaR_max = 0.5;        % Maximum rift step size 


ConvergenceCurve = zeros(1, T);

% Initialize population positions
[x, lb, ub] = Initialization_rdfo(N, nvars, ub, lb);
f = inf(N, 1);

H = zeros(1, nvars);           % Latent Stress Field (H), initialized to zero

bestf = inf;
bestx = zeros(1, nvars);
gBest_prev_score = inf;
Cstress = 0;                   % Stress counter for stagnation (C_stress)

%% Main Loop
for t = 1 : T
    
    % Evaluation & Boundary enforcement
    for i = 1:N
        x(i, :) = max(x(i, :), lb);
        x(i, :) = min(x(i, :), ub);
        f(i) = fun(x(i, :));
    end
    
    [current_best_score, best_idx] = min(f);

    if current_best_score < bestf
        bestf = current_best_score;
        bestx = x(best_idx, :);
    end

    % ---------------------------------------------------------
    % Stress-Triggered Fission Reinjection 
    % ---------------------------------------------------------
    if abs(bestf - gBest_prev_score) < 1e-10
        Cstress = Cstress + 1;
    else
        Cstress = 0;
        gBest_prev_score = bestf;
    end
    
    if Cstress >= Cmax
        % --- Fission Reset Triggered ---
        pop_mean = mean(x, 1);
        Vesc = normalize_vector_rdfo(pop_mean - bestx); % Escape vector (V_esc)
        
        % Select underperforming agents (S_reset) for reset
        [~, sorted_indices] = sort(f);
        Sreset_indices = sorted_indices(round(N/2):end);
        
        Rfission = (ub - lb) * (1 - t/T); % Fission range diminishes over time
        
        for i = 1:length(Sreset_indices)
            idx = Sreset_indices(i);
            % Project along escape vector
            projection_point = x(idx,:) + norm(x(idx,:) - bestx) * Vesc;
            % Add scaled random noise
            x(idx, :) = projection_point + Rfission .* (rand(1, nvars) - 0.5);
        end
        
        H = zeros(1, nvars); % Purge stale directional memory
        Cstress = 0;         % Reset stress counter
        
        % Re-evaluate after fission reinjection
        for i = 1:N
            x(i, :) = max(x(i, :), lb);
            x(i, :) = min(x(i, :), ub);
            f(i) = fun(x(i, :));
        end
    end
    
    % ---------------------------------------------------------
    % Rank-Adaptive Complementary Phase Gates & Update
    % ---------------------------------------------------------
    [sorted_fitness, sorted_indices] = sort(f);
    
    % Cache old states (X_old, f_old) for field assimilation
    X_old = x;
    f_old = f;
    
    % Dual-Phase attenuation schedule (Eq. 4-5)
    alphaD = alphaD_max * (1 - (t-1)/T);
    alphaR = alphaR_max * (1 - (t-1)/T);

    next_x = x;
    successful_moves = zeros(N, nvars);
    successful_gains = zeros(N, 1);
    
    for i = 1:N
        original_index = sorted_indices(i);
        pi_rank = (i-1) / (N-1); % Normalized rank p_i in [0, 1]
        
        % Phase Gates (Eq. 7-8)
        GD = cos(pi_rank * pi / 2);
        GR = sin(pi_rank * pi / 2);
        
        % Rift Shear Vector computation (Exploration)
        r1 = randi([1, N]); while r1 == original_index, r1 = randi([1, N]); end
        r2 = randi([1, N]); while r2 == original_index || r2 == r1, r2 = randi([1, N]); end
        
        VR = x(r1, :) - x(r2, :); % Rift Shear Vector
        
        % Core Position Update (Eq. 10)
        step = alphaD * GD * H + alphaR * GR * VR;
        
        xcand = x(original_index, :) + step;

        % Greedy Selection (Eq. 11)
        xcand = max(xcand, lb);
        xcand = min(xcand, ub);
        fcand = fun(xcand);

        % Record successful improvements for field assimilation
        if fcand < f(original_index)
            next_x(original_index, :) = xcand;
            f(original_index) = fcand;
            successful_moves(original_index, :) = xcand - X_old(original_index, :);
            successful_gains(original_index) = f_old(original_index) - fcand;
        end
    end
    x = next_x;
    
    % ---------------------------------------------------------
    % Latent Stress Field Evolution (Gain-weighted assimilation)
    % ---------------------------------------------------------
    [~, sorted_indices_old] = sort(f_old);
    num_elites = floor(N * Pe);
    elite_indices = sorted_indices_old(1:num_elites); % Elite subset E(t)
    
    H_raw = zeros(length(elite_indices), nvars);
    weights = zeros(length(elite_indices), 1);
    
    valid_elite_moves = 0;
    for i = 1:length(elite_indices)
        idx = elite_indices(i);
        if successful_gains(idx) > 0
            valid_elite_moves = valid_elite_moves + 1;
            H_raw(valid_elite_moves, :) = successful_gains(idx) * successful_moves(idx, :);
            weights(valid_elite_moves) = successful_gains(idx);
        end
    end
    
    if sum(weights) > 1e-10
        H_consensus = sum(H_raw, 1) / sum(weights);
        % Smooth update of the Latent Stress Field H (Eq. 13)
        H = (1 - eta) * H + eta * normalize_vector_rdfo(H_consensus);
    end
    
    % Record convergence curve
    ConvergenceCurve(t) = bestf;

end

end

% -------------------------------------------------------------------------
% Helper function for Initialization
% -------------------------------------------------------------------------
function[Positions, lb_extended, ub_extended] = Initialization_rdfo(SearchAgents_no, dim, ub, lb)
    Boundary_no = size(ub, 2);
    if Boundary_no == 1
        Positions = rand(SearchAgents_no, dim) .* (ub - lb) + lb;
        lb_extended = repmat(lb, 1, dim);
        ub_extended = repmat(ub, 1, dim);
    else
        Positions = zeros(SearchAgents_no, dim);
        for i = 1:dim
            ub_i = ub(i);
            lb_i = lb(i);
            Positions(:, i) = rand(SearchAgents_no, 1) .* (ub_i - lb_i) + lb_i;
        end
        lb_extended = lb;
        ub_extended = ub;
    end
end

% -------------------------------------------------------------------------
% Helper function to normalize a vector
% -------------------------------------------------------------------------
function v_norm = normalize_vector_rdfo(v)
    norm_val = norm(v);
    if norm_val > 1e-10
        v_norm = v / norm_val;
    else
        v_norm = v;
    end
end