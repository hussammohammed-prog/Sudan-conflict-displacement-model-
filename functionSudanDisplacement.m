function functionSudanDisplacement_GlobalSensitivity
%% ================================================================
% SUDAN DISPLACEMENT DYNAMICAL MODEL
% Standalone MATLAB Version - NO EXCEL FILE REQUIRED
% Khartoum - Port Sudan
%
% Save this file exactly as:
% functionSudanDisplacement_GlobalSensitivity.m
%
% Run by typing:
% functionSudanDisplacement_GlobalSensitivity
%% ================================================================

clear;
clc;
close all;

fprintf('\n============================================================\n');
fprintf('       SUDAN DISPLACEMENT DYNAMICAL MODEL\n');
fprintf('       Khartoum - Port Sudan\n');
fprintf('============================================================\n\n');

%% 1. OBSERVED DISPLACEMENT DATA
Date = datetime({ ...
    '2023-05-01'
    '2023-06-01'
    '2023-07-01'
    '2023-08-01'
    '2023-09-01'
    '2023-10-01'
    '2023-11-01'
    '2023-12-01'
    '2024-01-01'
    '2024-02-01'
    '2024-03-01'
    '2024-04-01'},'InputFormat','yyyy-MM-dd');

D_K_obs = [19585;22975;34750;39503;41538;60100;64030;44769;NaN;NaN;54355;69057];
D_R_obs = [17631;36835;96328;107942;112609;124470;129785;182238;239027;249555;247012;247874];
D_P_obs = [NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;231554;230701];

%% 2. CONFLICT DATA
Conflict_Start = datetime({ ...
    '2023-04-15'
    '2023-05-20'
    '2023-06-17'
    '2023-07-15'
    '2023-08-05'
    '2023-10-28'
    '2023-11-25'
    '2024-01-06'
    '2024-02-10'},'InputFormat','yyyy-MM-dd');

Conflict_End = datetime({ ...
    '2023-05-19'
    '2023-06-16'
    '2023-07-14'
    '2023-08-04'
    '2023-09-01'
    '2023-11-24'
    '2024-01-05'
    '2024-02-09'
    '2024-03-08'},'InputFormat','yyyy-MM-dd');

Conflict_Events = [215;182;200;260;300;240;440;440;157];

Conflict_Type = { ...
    'Observed'
    'Derived'
    'LowerBound'
    'LowerBound'
    'LowerBound'
    'LowerBound'
    'LowerBound'
    'LowerBound'
    'LowerBound'};

%% 3. MONTHLY CONFLICT INTENSITY
nMonths = length(Date);
Conflict_Raw = zeros(nMonths,1);

for i = 1:nMonths
    monthStart = Date(i);
    monthEnd = dateshift(Date(i),'end','month');

    for j = 1:length(Conflict_Events)
        overlapStart = max(monthStart,Conflict_Start(j));
        overlapEnd   = min(monthEnd,Conflict_End(j));

        if overlapStart <= overlapEnd
            totalDays = days(Conflict_End(j)-Conflict_Start(j))+1;
            overlapDays = days(overlapEnd-overlapStart)+1;
            Conflict_Raw(i) = Conflict_Raw(i) + ...
                Conflict_Events(j)*(overlapDays/totalDays);
        end
    end
end

if max(Conflict_Raw) > 0
    W_K = Conflict_Raw./max(Conflict_Raw);
else
    W_K = zeros(size(Conflict_Raw));
end

% Proxy for shelling severity because the supplied dataset has no
% independent shelling-severity measurement.
S_K = sqrt(W_K);

%% 4. TIME AND INITIAL CONDITIONS
t = (0:nMonths-1)';

idxK0 = find(~isnan(D_K_obs),1,'first');
D_K0 = D_K_obs(idxK0);

% No Port Sudan observation exists before March 2024.
% A small positive numerical seed is used.
D_P0 = 1;

%% 5. INITIAL PARAMETER VALUES
% p = [lambda_KP lambda_PK mu_K mu_P alpha beta K_P q_K q_P]

p0 = [0.005 0.002 0.001 0.001 0.50 0.20 400000 5000 20000];

%% 6. PARAMETER BOUNDS
lb = [0 0 0 0 0 0 max(D_P_obs(~isnan(D_P_obs)))+1 0 0];
ub = [1 1 1 1 10 10 2000000 200000 200000];

%% 7. OBSERVATION WEIGHTS
wK = ones(nMonths,1);
wP = zeros(nMonths,1);
wP(~isnan(D_P_obs)) = 2;

%% 8. PARAMETER ESTIMATION
fprintf('Estimating model parameters...\n');
fprintf('Please wait...\n\n');

% IMPORTANT:
% We use a bounded transformation so the optimization is numerically
% stable and all parameters remain inside their physical bounds.

z0 = parameterTransformForward(p0,lb,ub);

objectiveZ = @(z) modelObjective( ...
    parameterTransformBackward(z,lb,ub), ...
    t,D_K_obs,D_P_obs,W_K,S_K,D_K0,D_P0,wK,wP);

options = optimset( ...
    'Display','off', ...
    'MaxIter',5000, ...
    'MaxFunEvals',20000, ...
    'TolX',1e-8, ...
    'TolFun',1e-8);

[z_est,fval,exitflag,output] = fminsearch(objectiveZ,z0,options);
p_est = parameterTransformBackward(z_est,lb,ub);

%% 9. ESTIMATED PARAMETERS
lambda_KP = p_est(1);
lambda_PK = p_est(2);
mu_K      = p_est(3);
mu_P      = p_est(4);
alpha     = p_est(5);
beta      = p_est(6);
K_P       = p_est(7);
q_K       = p_est(8);
q_P       = p_est(9);

fprintf('\n============================================================\n');
fprintf('             ESTIMATED PARAMETERS\n');
fprintf('============================================================\n');
fprintf('lambda_KP = %.10f\n',lambda_KP);
fprintf('lambda_PK = %.10f\n',lambda_PK);
fprintf('mu_K      = %.10f\n',mu_K);
fprintf('mu_P      = %.10f\n',mu_P);
fprintf('alpha     = %.10f\n',alpha);
fprintf('beta      = %.10f\n',beta);
fprintf('K_P       = %.4f\n',K_P);
fprintf('q_K       = %.4f\n',q_K);
fprintf('q_P       = %.4f\n',q_P);
fprintf('============================================================\n\n');

%% 10. MODEL SIMULATION
[t_model,Y_model] = simulateModel(p_est,t,D_K0,D_P0,W_K,S_K);

D_K_model = Y_model(:,1);
D_P_model = Y_model(:,2);

%% 11. MODEL PERFORMANCE
validK = ~isnan(D_K_obs);
validP = ~isnan(D_P_obs);

RMSE_K = sqrt(mean((D_K_model(validK)-D_K_obs(validK)).^2));
MAE_K  = mean(abs(D_K_model(validK)-D_K_obs(validK)));

if any(validP)
    RMSE_P = sqrt(mean((D_P_model(validP)-D_P_obs(validP)).^2));
    MAE_P  = mean(abs(D_P_model(validP)-D_P_obs(validP)));
else
    RMSE_P = NaN;
    MAE_P  = NaN;
end

SSresK = sum((D_K_obs(validK)-D_K_model(validK)).^2);
SStotK = sum((D_K_obs(validK)-mean(D_K_obs(validK))).^2);
R2_K = 1-SSresK/SStotK;

if sum(validP) >= 2
    SSresP = sum((D_P_obs(validP)-D_P_model(validP)).^2);
    SStotP = sum((D_P_obs(validP)-mean(D_P_obs(validP))).^2);
    R2_P = 1-SSresP/SStotP;
else
    R2_P = NaN;
end

fprintf('============================================================\n');
fprintf('             MODEL PERFORMANCE\n');
fprintf('============================================================\n');
fprintf('Khartoum RMSE  = %.4f\n',RMSE_K);
fprintf('Khartoum MAE   = %.4f\n',MAE_K);
fprintf('Khartoum R^2   = %.4f\n',R2_K);
fprintf('Port Sudan RMSE = %.4f\n',RMSE_P);
fprintf('Port Sudan MAE  = %.4f\n',MAE_P);
fprintf('Port Sudan R^2  = %.4f\n',R2_P);
fprintf('Objective       = %.8f\n',fval);
fprintf('============================================================\n\n');

%% 12. RESULTS TABLE
Results = table(Date,D_K_obs,D_K_model,D_R_obs,D_P_obs,D_P_model,W_K,S_K, ...
    'VariableNames',{'Date','Khartoum_Observed','Khartoum_Model', ...
    'RedSea_Observed','PortSudan_Observed','PortSudan_Model', ...
    'Conflict_Intensity','Shelling_Proxy'});

disp(Results);

%% 13. KHARTOUM PLOT
figure('Name','Khartoum Displacement');
plot(Date,D_K_obs,'o','LineWidth',1.5);
hold on;
plot(Date,D_K_model,'-','LineWidth',2);
xlabel('Date');
ylabel('Displaced Population');
title('Khartoum: Observed vs Model');
legend('Observed','Model','Location','best');
grid on;

%% 14. PORT SUDAN PLOT
figure('Name','Port Sudan Displacement');
plot(Date,D_P_obs,'o','LineWidth',1.5);
hold on;
plot(Date,D_P_model,'-','LineWidth',2);
xlabel('Date');
ylabel('Displaced Population');
title('Port Sudan: Observed vs Model');
legend('Observed','Model','Location','best');
grid on;

%% 15. CONFLICT INTENSITY PLOT
figure('Name','Conflict Intensity');
bar(Date,W_K);
xlabel('Date');
ylabel('Normalized Conflict Intensity');
title('Khartoum Conflict Intensity Index');
grid on;

%% 16. DISPLACEMENT PRESSURE INDEX
P_K = alpha.*W_K + beta.*S_K;

if max(P_K)>0
    DPI = P_K./max(P_K);
else
    DPI = zeros(size(P_K));
end

figure('Name','Displacement Pressure Index');
plot(Date,DPI,'-','LineWidth',2);
hold on;
yline(0.80,'--','Critical Threshold');
xlabel('Date');
ylabel('DPI');
title('Displacement Pressure Index');
grid on;

CriticalThreshold = 0.80;
critical_idx = find(DPI >= CriticalThreshold);

fprintf('============================================================\n');
fprintf('        DISPLACEMENT PRESSURE INDEX\n');
fprintf('============================================================\n');

for i = 1:length(Date)
    fprintf('%s   DPI = %.4f\n',datestr(Date(i),'yyyy-mm'),DPI(i));
end

fprintf('\nCritical threshold = %.2f\n',CriticalThreshold);

if isempty(critical_idx)
    fprintf('No period reached the critical threshold.\n');
else
    fprintf('Critical periods:\n');
    for i = 1:length(critical_idx)
        fprintf('%s   DPI = %.4f\n', ...
            datestr(Date(critical_idx(i)),'yyyy-mm'), ...
            DPI(critical_idx(i)));
    end
end

fprintf('============================================================\n\n');

%% 17. EQUILIBRIUM ANALYSIS
% Equilibria are checked against the frozen autonomous system obtained by
% fixing the conflict intensity at W_eq = max(W_K). Because the calibrated
% model may contain a parameter very close to zero, a numerical fsolve
% result is NOT accepted unless its residual is sufficiently small.
W_eq = max(W_K);
S_eq = sqrt(W_eq);

fprintf('============================================================\n');
fprintf('              EQUILIBRIUM ANALYSIS\n');
fprintf('============================================================\n');
fprintf('Frozen conflict level: W_eq = %.6f, S_eq = %.6f\n',W_eq,S_eq);

F_equilibrium = @(Y) displacementRHS(Y,p_est,W_eq,S_eq);

initial_guesses = [ ...
    D_K0  D_P0
    0     0
    50000 50000
    100000 100000
    200000 200000
    500000 300000];

Equilibria = [];
EquilibriumResiduals = [];
StabilityLabels = cell(0,1);
Eigenvalues_All = cell(0,1);

if exist('fsolve','file') == 2
    for i = 1:size(initial_guesses,1)
        try
            optionsF = optimoptions('fsolve', ...
                'Display','off', ...
                'FunctionTolerance',1e-10, ...
                'StepTolerance',1e-10);

            Yeq = fsolve(F_equilibrium,initial_guesses(i,:)',optionsF);
            residual = norm(F_equilibrium(Yeq),2);
            scaleResidual = max(1,norm([q_K;q_P],2));
            relativeResidual = residual/scaleResidual;

            % Accept only genuine numerical roots.
            if all(isfinite(Yeq)) && all(Yeq >= 0) && relativeResidual <= 1e-6
                if isempty(Equilibria)
                    Equilibria = Yeq(:)';
                    EquilibriumResiduals = residual;
                else
                    difference = vecnorm(Equilibria-Yeq(:)',2,2);
                    if all(difference > 1)
                        Equilibria = [Equilibria;Yeq(:)'];
                        EquilibriumResiduals = [EquilibriumResiduals;residual];
                    end
                end
            end
        catch
        end
    end
else
    fprintf('fsolve is not available. Numerical equilibrium search skipped.\n');
end

if isempty(Equilibria)
    fprintf('\nNo verified numerical equilibrium was detected for the frozen system.\n');
    fprintf(['This can occur when the exogenous source q_P remains positive while ' ...
             'mu_P is near zero.\n']);
    fprintf(['Therefore, no asymptotic-stability conclusion is assigned from ' ...
             'the failed/near-zero eigenvalue calculation.\n']);
else
    for i = 1:size(Equilibria,1)
        DK_eq = Equilibria(i,1);
        DP_eq = Equilibria(i,2);

        fprintf('\nVerified Equilibrium %d:\n',i);
        fprintf('D_K* = %.8f\n',DK_eq);
        fprintf('D_P* = %.8f\n',DP_eq);
        fprintf('Residual norm = %.6e\n',EquilibriumResiduals(i));

        J = numericalJacobian(F_equilibrium,[DK_eq;DP_eq]);
        eigenvalues = eig(J);

        fprintf('Jacobian:\n');
        disp(J);
        fprintf('Eigenvalues:\n');
        disp(eigenvalues);

        stabilityTol = 1e-8;
        if all(real(eigenvalues) < -stabilityTol)
            stabilityLabel = 'ASYMPTOTICALLY STABLE';
        elseif any(real(eigenvalues) > stabilityTol)
            stabilityLabel = 'UNSTABLE';
        else
            stabilityLabel = 'NON-HYPERBOLIC / INCONCLUSIVE';
        end

        fprintf('Stability: %s\n',stabilityLabel);
        StabilityLabels{i,1} = stabilityLabel;
        Eigenvalues_All{i,1} = eigenvalues;
    end
end

fprintf('============================================================\n\n');

%% 18. LOCAL SENSITIVITY ANALYSIS (OAT: +/-10%)
% One-at-a-time local sensitivity analysis.
% The analysis measures the effect of each parameter on the
% time-integrated displacement in Khartoum, Port Sudan, and both cities.

fprintf('\n============================================================\n');
fprintf('             LOCAL SENSITIVITY ANALYSIS\n');
fprintf('============================================================\n');

parameterNames = { ...
    'lambda_KP', 'lambda_PK', 'mu_K', 'mu_P', ...
    'alpha', 'beta', 'K_P', 'q_K', 'q_P'};

nParameters = length(p_est);
perturbation = 0.10;

% Baseline integrated outputs.
Baseline_K = trapz(t_model,D_K_model);
Baseline_P = trapz(t_model,D_P_model);
Baseline_Total = Baseline_K + Baseline_P;

% Use the same physical bounds employed during calibration.
% This protects the sensitivity runs from leaving the admissible region.
if exist('lb','var') && exist('ub','var')
    sensitivity_lb = lb(:);
    sensitivity_ub = ub(:);
else
    sensitivity_lb = zeros(nParameters,1);
    sensitivity_ub = inf(nParameters,1);
end

Sensitivity_K = NaN(nParameters,1);
Sensitivity_P = NaN(nParameters,1);
Sensitivity_Total = NaN(nParameters,1);
AbsoluteSensitivity_K = NaN(nParameters,1);
AbsoluteSensitivity_P = NaN(nParameters,1);
AbsoluteSensitivity_Total = NaN(nParameters,1);

for i = 1:nParameters
    p_plus = p_est;
    p_minus = p_est;

    % Symmetric perturbation, clipped to the calibration bounds.
    p_plus(i) = min(sensitivity_ub(i), ...
        max(sensitivity_lb(i),p_est(i)*(1 + perturbation)));
    p_minus(i) = min(sensitivity_ub(i), ...
        max(sensitivity_lb(i),p_est(i)*(1 - perturbation)));

    [t_plus,Y_plus] = simulateModel(p_plus,t,D_K0,D_P0,W_K,S_K);
    [t_minus,Y_minus] = simulateModel(p_minus,t,D_K0,D_P0,W_K,S_K);

    K_plus = trapz(t_plus,Y_plus(:,1));
    P_plus = trapz(t_plus,Y_plus(:,2));
    Total_plus = K_plus + P_plus;

    K_minus = trapz(t_minus,Y_minus(:,1));
    P_minus = trapz(t_minus,Y_minus(:,2));
    Total_minus = K_minus + P_minus;

    delta_p = p_plus(i) - p_minus(i);

    if abs(delta_p) > eps
        % Absolute finite-difference sensitivities.
        AbsoluteSensitivity_K(i) = (K_plus-K_minus)/delta_p;
        AbsoluteSensitivity_P(i) = (P_plus-P_minus)/delta_p;
        AbsoluteSensitivity_Total(i) = (Total_plus-Total_minus)/delta_p;

        % Dimensionless elasticities. These are left as NaN when the
        % estimated parameter is approximately zero or the baseline is zero.
        if abs(p_est(i)) > 1e-12 && abs(Baseline_K) > eps
            Sensitivity_K(i) = AbsoluteSensitivity_K(i)*p_est(i)/Baseline_K;
        end
        if abs(p_est(i)) > 1e-12 && abs(Baseline_P) > eps
            Sensitivity_P(i) = AbsoluteSensitivity_P(i)*p_est(i)/Baseline_P;
        end
        if abs(p_est(i)) > 1e-12 && abs(Baseline_Total) > eps
            Sensitivity_Total(i) = AbsoluteSensitivity_Total(i)*p_est(i)/Baseline_Total;
        end
    end
end

% Calculate an overall dimensionless sensitivity while ignoring NaN values.
Overall_Sensitivity = NaN(nParameters,1);
for i = 1:nParameters
    values = [abs(Sensitivity_K(i)), ...
              abs(Sensitivity_P(i)), ...
              abs(Sensitivity_Total(i))];
    values = values(isfinite(values));
    if ~isempty(values)
        Overall_Sensitivity(i) = mean(values);
    end
end

Sensitivity_Results = table( ...
    parameterNames',p_est(:),Sensitivity_K, ...
    Sensitivity_P,Sensitivity_Total, ...
    AbsoluteSensitivity_K,AbsoluteSensitivity_P, ...
    AbsoluteSensitivity_Total,Overall_Sensitivity, ...
    'VariableNames',{'Parameter','Estimated_Value', ...
    'Sensitivity_Khartoum','Sensitivity_Port_Sudan', ...
    'Sensitivity_Total','Absolute_Khartoum', ...
    'Absolute_Port_Sudan','Absolute_Total', ...
    'Overall_Sensitivity'});

% Sort a copy of the sensitivity table without changing the bar-plot order.
SortValues = Overall_Sensitivity;
SortValues(~isfinite(SortValues)) = -Inf;
[~,sortIndex] = sort(SortValues,'descend');
Sensitivity_Results_Ordered = Sensitivity_Results(sortIndex,:);

fprintf('\nSensitivity results ordered by absolute overall sensitivity:\n');
disp(Sensitivity_Results_Ordered);

figure('Name','Parameter Sensitivity');
bar(Overall_Sensitivity);
set(gca,'XTick',1:nParameters);
set(gca,'XTickLabel',parameterNames);
if exist('xtickangle','file') == 2
    xtickangle(45);
end
ylabel('Mean Absolute Elasticity');
title('Local Parameter Sensitivity (+/-10%)');
grid on;

save('Sudan_Displacement_Sensitivity.mat', ...
    'Sensitivity_Results','Sensitivity_Results_Ordered', ...
    'Sensitivity_K','Sensitivity_P','Sensitivity_Total', ...
    'AbsoluteSensitivity_K','AbsoluteSensitivity_P', ...
    'AbsoluteSensitivity_Total','Overall_Sensitivity', ...
    'Baseline_K','Baseline_P','Baseline_Total','perturbation');

fprintf('Sensitivity results saved to Sudan_Displacement_Sensitivity.mat\n');
fprintf('============================================================\n\n');

%% 19. GLOBAL SENSITIVITY ANALYSIS - LATIN HYPERCUBE SAMPLING
% Global sensitivity is performed over the full admissible parameter
% ranges rather than only +/-10%% around the calibrated point.
% The primary global measure is Partial Rank Correlation (PRCC).
% A manual LHS generator is used so that the analysis does not require
% the Statistics and Machine Learning Toolbox.

fprintf('\n============================================================\n');
fprintf('       GLOBAL SENSITIVITY ANALYSIS (LHS + PRCC)\n');
fprintf('============================================================\n');

% Number of Latin Hypercube samples.
% 1000 is suitable for the final study. For a quick mobile test,
% this can be reduced to 300-500.
nLHS = 1000;

rng(20260922,'twister');

nPar = length(p_est);
U = zeros(nLHS,nPar);

% Latin Hypercube construction: one randomized point in each stratum.
for j = 1:nPar
    U(:,j) = ((0:nLHS-1)' + rand(nLHS,1))/nLHS;
    U(:,j) = U(randperm(nLHS),j);
end

LHS_Parameters = zeros(nLHS,nPar);
for j = 1:nPar
    LHS_Parameters(:,j) = sensitivity_lb(j) + ...
        U(:,j).*(sensitivity_ub(j)-sensitivity_lb(j));
end

% Force the calibrated parameter vector into the sampling design by
% replacing the first row. This gives a direct baseline reference.
LHS_Parameters(1,:) = p_est(:)';

% Outputs used in the global analysis.
LHS_Final_K = NaN(nLHS,1);
LHS_Final_P = NaN(nLHS,1);
LHS_Final_Total = NaN(nLHS,1);
LHS_Cumulative_K = NaN(nLHS,1);
LHS_Cumulative_P = NaN(nLHS,1);
LHS_Cumulative_Total = NaN(nLHS,1);
LHS_Objective = NaN(nLHS,1);

fprintf('Running %d global parameter samples...\n',nLHS);

for s = 1:nLHS
    p_s = LHS_Parameters(s,:)';
    try
        [t_s,Y_s] = simulateModel(p_s,t,D_K0,D_P0,W_K,S_K);
        DK_s = Y_s(:,1);
        DP_s = Y_s(:,2);

        if all(isfinite(DK_s)) && all(isfinite(DP_s))
            LHS_Final_K(s) = DK_s(end);
            LHS_Final_P(s) = DP_s(end);
            LHS_Final_Total(s) = DK_s(end) + DP_s(end);
            LHS_Cumulative_K(s) = trapz(t_s,DK_s);
            LHS_Cumulative_P(s) = trapz(t_s,DP_s);
            LHS_Cumulative_Total(s) = LHS_Cumulative_K(s) + LHS_Cumulative_P(s);
            LHS_Objective(s) = modelObjective(p_s,t,D_K_obs,D_P_obs,W_K,S_K,D_K0,D_P0,wK,wP);
        end
    catch
        % Failed trajectories remain NaN and are excluded from PRCC.
    end

    if mod(s,max(1,round(nLHS/10))) == 0
        fprintf('  %d%% completed\n',round(100*s/nLHS));
    end
end

validGlobal = isfinite(LHS_Final_Total) & ...
              isfinite(LHS_Final_K) & isfinite(LHS_Final_P);

nValidGlobal = sum(validGlobal);
fprintf('Valid global simulations: %d of %d\n',nValidGlobal,nLHS);

% Rank-correlation / PRCC analysis.
Global_Spearman_Final_K = NaN(nPar,1);
Global_Spearman_Final_P = NaN(nPar,1);
Global_Spearman_Final_Total = NaN(nPar,1);
Global_PRCC_Final_K = NaN(nPar,1);
Global_PRCC_Final_P = NaN(nPar,1);
Global_PRCC_Final_Total = NaN(nPar,1);
Global_PRCC_Cumulative_Total = NaN(nPar,1);

if nValidGlobal >= max(20,nPar+5)
    X = LHS_Parameters(validGlobal,:);
    YK = LHS_Final_K(validGlobal);
    YP = LHS_Final_P(validGlobal);
    YT = LHS_Final_Total(validGlobal);
    YCT = LHS_Cumulative_Total(validGlobal);

    XR = zeros(size(X));
    for j = 1:nPar
        XR(:,j) = rankData(X(:,j));
    end

    YKR = rankData(YK);
    YPR = rankData(YP);
    YTR = rankData(YT);
    YCTR = rankData(YCT);

    for j = 1:nPar
        Global_Spearman_Final_K(j) = simpleCorrelation(XR(:,j),YKR);
        Global_Spearman_Final_P(j) = simpleCorrelation(XR(:,j),YPR);
        Global_Spearman_Final_Total(j) = simpleCorrelation(XR(:,j),YTR);

        Global_PRCC_Final_K(j) = partialRankCorrelation(XR,YKR,j);
        Global_PRCC_Final_P(j) = partialRankCorrelation(XR,YPR,j);
        Global_PRCC_Final_Total(j) = partialRankCorrelation(XR,YTR,j);
        Global_PRCC_Cumulative_Total(j) = partialRankCorrelation(XR,YCTR,j);
    end
end

Global_Overall_PRCC = NaN(nPar,1);
for j = 1:nPar
    v = [abs(Global_PRCC_Final_K(j)), ...
         abs(Global_PRCC_Final_P(j)), ...
         abs(Global_PRCC_Final_Total(j)), ...
         abs(Global_PRCC_Cumulative_Total(j))];
    v = v(isfinite(v));
    if ~isempty(v)
        Global_Overall_PRCC(j) = mean(v);
    end
end

Global_Sensitivity_Results = table( ...
    parameterNames', ...
    Global_Spearman_Final_K, ...
    Global_Spearman_Final_P, ...
    Global_Spearman_Final_Total, ...
    Global_PRCC_Final_K, ...
    Global_PRCC_Final_P, ...
    Global_PRCC_Final_Total, ...
    Global_PRCC_Cumulative_Total, ...
    Global_Overall_PRCC, ...
    'VariableNames',{'Parameter','Spearman_Final_Khartoum', ...
    'Spearman_Final_Port_Sudan','Spearman_Final_Total', ...
    'PRCC_Final_Khartoum','PRCC_Final_Port_Sudan', ...
    'PRCC_Final_Total','PRCC_Cumulative_Total', ...
    'Overall_PRCC'});

SortGlobal = Global_Overall_PRCC;
SortGlobal(~isfinite(SortGlobal)) = -Inf;
[~,globalSortIndex] = sort(SortGlobal,'descend');
Global_Sensitivity_Results_Ordered = Global_Sensitivity_Results(globalSortIndex,:);

fprintf('\nGlobal sensitivity results ordered by absolute overall PRCC:\n');
disp(Global_Sensitivity_Results_Ordered);

% Direct comparison between the local elasticity analysis and the global
% PRCC analysis. The parameter order is kept identical in both columns.
Sensitivity_Comparison = table( ...
    parameterNames',Overall_Sensitivity,Global_Overall_PRCC, ...
    'VariableNames',{'Parameter','Local_Overall_Elasticity','Global_Overall_PRCC'});

fprintf('\nLOCAL vs GLOBAL SENSITIVITY:\n');
disp(Sensitivity_Comparison);

figure('Name','Global PRCC Sensitivity');
bar(Global_Overall_PRCC);
set(gca,'XTick',1:nPar);
set(gca,'XTickLabel',parameterNames);
if exist('xtickangle','file') == 2
    xtickangle(45);
end
ylabel('Mean Absolute PRCC');
title('Global Parameter Sensitivity - Latin Hypercube Sampling');
grid on;

save('Sudan_Displacement_Global_Sensitivity.mat', ...
    'nLHS','LHS_Parameters','validGlobal','nValidGlobal', ...
    'LHS_Final_K','LHS_Final_P','LHS_Final_Total', ...
    'LHS_Cumulative_K','LHS_Cumulative_P','LHS_Cumulative_Total', ...
    'LHS_Objective','Global_Sensitivity_Results', ...
    'Global_Sensitivity_Results_Ordered', ...
    'Global_Spearman_Final_K','Global_Spearman_Final_P', ...
    'Global_Spearman_Final_Total','Global_PRCC_Final_K', ...
    'Global_PRCC_Final_P','Global_PRCC_Final_Total', ...
    'Global_PRCC_Cumulative_Total','Global_Overall_PRCC', ...
    'Sensitivity_Comparison');

fprintf('Global sensitivity results saved to Sudan_Displacement_Global_Sensitivity.mat\n');
fprintf('============================================================\n\n');

%% 20. SAVE RESULTS
save('Sudan_Displacement_Model_Results.mat', ...
    'Date','D_K_obs','D_R_obs','D_P_obs', ...
    'Conflict_Start','Conflict_End','Conflict_Events','Conflict_Type', ...
    'Conflict_Raw','W_K','S_K', ...
    'D_K_model','D_P_model','p_est', ...
    'lambda_KP','lambda_PK','mu_K','mu_P','alpha','beta','K_P','q_K','q_P', ...
    'Results','RMSE_K','MAE_K','RMSE_P','MAE_P','R2_K','R2_P', ...
    'DPI','CriticalThreshold','Equilibria', ...
    'StabilityLabels','Eigenvalues_All','EquilibriumResiduals', ...
    'Global_Sensitivity_Results','Global_Sensitivity_Results_Ordered', ...
    'Global_Overall_PRCC','Sensitivity_Comparison','nLHS','nValidGlobal', ...
    'exitflag','output');

fprintf('Results saved to Sudan_Displacement_Model_Results.mat\n');
fprintf('\nMODEL RUN COMPLETED.\n');

end


%% ================================================================
% LOCAL / SUBFUNCTIONS
%% ================================================================

function J = numericalJacobian(F,Y)
n = length(Y);
J = zeros(n,n);

for i = 1:n
    h = 1e-5*max(1,abs(Y(i)));
    Y_plus = Y;
    Y_minus = Y;
    Y_plus(i) = Y_plus(i)+h;
    Y_minus(i) = Y_minus(i)-h;
    J(:,i) = (F(Y_plus)-F(Y_minus))/(2*h);
end
end


function cost = modelObjective(p,t,D_K_obs,D_P_obs,W_K,S_K,D_K0,D_P0,wK,wP)

try
    if any(~isfinite(p)) || any(p<0)
        cost = 1e30;
        return;
    end

    [~,Y] = simulateModel(p,t,D_K0,D_P0,W_K,S_K);

    DK = Y(:,1);
    DP = Y(:,2);

    if any(~isfinite(DK)) || any(~isfinite(DP))
        cost = 1e30;
        return;
    end

    validK = ~isnan(D_K_obs);
    errK = D_K_obs(validK)-DK(validK);
    costK = sum(wK(validK).*errK.^2);

    validP = ~isnan(D_P_obs);

    if any(validP)
        errP = D_P_obs(validP)-DP(validP);
        costP = sum(wP(validP).*errP.^2);
    else
        costP = 0;
    end

    scaleK = max(D_K_obs(validK));

    if any(validP)
        scaleP = max(D_P_obs(validP));
    else
        scaleP = 1;
    end

    cost = costK/(scaleK^2) + costP/(scaleP^2);

catch
    cost = 1e30;
end
end


function [t_out,Y] = simulateModel(p,t,D_K0,D_P0,W_K,S_K)

Y0 = [D_K0;D_P0];

odefun = @(tt,Y) displacementRHSInterpolated(tt,Y,p,t,W_K,S_K);

options = odeset( ...
    'RelTol',1e-7, ...
    'AbsTol',1e-8, ...
    'NonNegative',[1 2]);

[t_out,Y] = ode45(odefun,t,Y0,options);
end


function dYdt = displacementRHSInterpolated(tt,Y,p,t,W_K,S_K)

W = interp1(t,W_K,tt,'linear','extrap');
S = interp1(t,S_K,tt,'linear','extrap');

W = max(0,min(1,W));
S = max(0,min(1,S));

dYdt = displacementRHS(Y,p,W,S);
end


function dYdt = displacementRHS(Y,p,W_K,S_K)

lambda_KP = p(1);
lambda_PK = p(2);
mu_K      = p(3);
mu_P      = p(4);
alpha     = p(5);
beta      = p(6);
K_P       = p(7);
q_K       = p(8);
q_P       = p(9);

D_K = max(0,Y(1));
D_P = max(0,Y(2));

P_K = alpha*W_K + beta*S_K;

A_eff = max(0,1-D_P/K_P);

F_KP = lambda_KP*D_K*P_K*A_eff;

F_PK = lambda_PK*D_P*(1-W_K);

E_K = mu_K*D_K*P_K;
E_P = mu_P*D_P;

% Exogenous displacement sources (persons per month).
B_K = q_K;
B_P = q_P;

dD_K = B_K - F_KP + F_PK - E_K;
dD_P = B_P + F_KP - F_PK - E_P;

dYdt = [dD_K;dD_P];
end



function r = rankData(x)
% Rank data with average ranks for ties. No Statistics Toolbox required.
x = x(:);
[sorted,order] = sort(x);
r = zeros(size(x));
i = 1;
n = length(x);
while i <= n
    j = i;
    while j < n && sorted(j+1) == sorted(i)
        j = j + 1;
    end
    avgRank = (i+j)/2;
    r(order(i:j)) = avgRank;
    i = j + 1;
end
end


function r = simpleCorrelation(x,y)
x = x(:); y = y(:);
mx = mean(x); my = mean(y);
dx = x-mx; dy = y-my;
den = sqrt(sum(dx.^2)*sum(dy.^2));
if den <= eps
    r = NaN;
else
    r = sum(dx.*dy)/den;
end
end


function prcc = partialRankCorrelation(XR,yR,j)
% Partial rank correlation of parameter j with output y, controlling
% for all remaining parameters.
[n,m] = size(XR);
if m <= 1
    prcc = simpleCorrelation(XR(:,j),yR);
    return;
end

others = setdiff(1:m,j);
A = [ones(n,1),XR(:,others)];

try
    betaY = A\yR;
    betaX = A\XR(:,j);
    resY = yR - A*betaY;
    resX = XR(:,j) - A*betaX;
    prcc = simpleCorrelation(resX,resY);
catch
    prcc = NaN;
end
end

function z = parameterTransformForward(p,lb,ub)

q = (p-lb)./(ub-lb);
q = min(max(q,1e-8),1-1e-8);
z = log(q./(1-q));
end


function p = parameterTransformBackward(z,lb,ub)

q = 1./(1+exp(-z));
p = lb + q.*(ub-lb);
end
