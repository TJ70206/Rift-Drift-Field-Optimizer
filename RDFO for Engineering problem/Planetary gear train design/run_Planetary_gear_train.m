clc
clear
close all

%% Planetary gear train experiment
baseDir = fileparts(mfilename('fullpath'));
addpath(genpath(baseDir));

selectedProblemName = 'Planetary_gear_train_design';
problemNumber = 1;

%% Experiment settings
rng('default');
rng(42);

options.PopulationSize = 30;
options.MaxFEs = 1 * 10^5;
curveLength = getCurveLength(options);
runs = 30;
baseSeed = 42;

resultsRoot = fullfile(baseDir, 'results');
if ~exist(resultsRoot, 'dir')
    mkdir(resultsRoot);
end

%% Run RDFO experiment
[dim, lb, ub, Vio, GloMin, Obj, EvalObj] = ProbInfo(problemNumber);
fobj = @(x) localPenaltyObjective(x, Vio, Obj);
[lbVec, ubVec] = normalizeBounds(lb, ub, dim);

useScientificFormat = contains(lower(selectedProblemName), 'gear_train');

problemDir = fullfile(resultsRoot, selectedProblemName);
if ~exist(problemDir, 'dir')
    mkdir(problemDir);
end

bestScores = NaN(runs, 1);
runtimes = NaN(runs, 1);
allCurves = NaN(runs, curveLength);
bestPositions = NaN(runs, dim);
successMask = false(runs, 1);
runSeeds = baseSeed + (1:runs);

fprintf('\nRunning %s (D=%d)\n', strrep(selectedProblemName, '_', ' '), dim);

for r = 1:runs
    try
        rng(runSeeds(r), 'twister');
        [bestScore, bestPos, cgCurve, runtime] = runSingleAlgorithm(options, lbVec, ubVec, fobj);
        bestScores(r) = bestScore;
        runtimes(r) = runtime;
        allCurves(r, :) = reshape(cgCurve, 1, []);
        bestPositions(r, :) = reshape(bestPos, 1, []);
        successMask(r) = isfinite(bestScore);
        fprintf('  Run %02d/%02d: %s | time = %.2fs\n', r, runs, formatScoreForDisplay(bestScore, useScientificFormat), runtime);
    catch ME
        fprintf('  Run %02d/%02d failed: %s\n', r, runs, ME.message);
    end
end

% Evaluate run results and save one submission file.
[objectiveValues, penaltyValues, maxViolations, feasibleFlags] = evaluateRunResults(bestPositions, successMask, Vio, Obj, EvalObj);
resultSummary = summarizeRuns(bestScores, runtimes, objectiveValues, penaltyValues, maxViolations, feasibleFlags, successMask);
bestRunDetail = buildRunDetail(resultSummary.BestRun, resultSummary.BestScore, bestPositions, Vio, Obj, EvalObj);

fprintf('  best = %s | mean time = %.2fs\n', formatScoreForDisplay(resultSummary.Best, useScientificFormat), resultSummary.AverageTime);

save(fullfile(problemDir, 'result.mat'), 'problemNumber', 'selectedProblemName', 'dim', 'lbVec', 'ubVec', 'GloMin', ...
    'options', 'curveLength', 'runs', 'baseSeed', 'runSeeds', 'bestScores', 'runtimes', 'allCurves', ...
    'bestPositions', 'successMask', 'objectiveValues', 'penaltyValues', 'maxViolations', 'feasibleFlags', ...
    'resultSummary', 'bestRunDetail');

%% Local functions
% Run RDFO once and normalize the saved convergence curve.
function [bestScore, bestPos, cgCurve, runtime] = runSingleAlgorithm(options, lb, ub, fobj)
curveLength = getCurveLength(options);
dim = numel(lb);
tic
[bestPos, bestScore, rawCurve] = RDFO(fobj, dim, lb, ub, options);
runtime = toc;

cgCurve = normalizeCurve(rawCurve, curveLength);
if isnan(bestScore) && any(isfinite(cgCurve))
    bestScore = cgCurve(end);
end
end

function curve = normalizeCurve(rawCurve, maxIter)
if isempty(rawCurve)
    curve = NaN(1, maxIter);
    return
end
rawCurve = reshape(rawCurve, 1, []);
finiteMask = isfinite(rawCurve);
if ~any(finiteMask)
    curve = NaN(1, maxIter);
    return
end
lastFiniteIdx = find(finiteMask, 1, 'last');
rawCurve = rawCurve(1:lastFiniteIdx);
if isscalar(rawCurve)
    curve = repmat(rawCurve, 1, maxIter);
else
    sampleIdx = round(linspace(1, numel(rawCurve), maxIter));
    sampleIdx(sampleIdx < 1) = 1;
    curve = rawCurve(sampleIdx);
end
firstFiniteIdx = find(isfinite(curve), 1, 'first');
if isempty(firstFiniteIdx)
    curve = NaN(1, maxIter);
    return
end
if firstFiniteIdx > 1
    curve(1:firstFiniteIdx - 1) = curve(firstFiniteIdx);
end
for k = 2:maxIter
    if ~isfinite(curve(k))
        curve(k) = curve(k - 1);
    end
end
curve = cummin(curve);
end

function curveLength = getCurveLength(options)
curveLength = max(1, ceil(options.MaxFEs / options.PopulationSize));
end

function [lbVec, ubVec] = normalizeBounds(lb, ub, dim)
if isscalar(lb)
    lbVec = repmat(lb, 1, dim);
else
    lbVec = reshape(lb, 1, []);
end
if isscalar(ub)
    ubVec = repmat(ub, 1, dim);
else
    ubVec = reshape(ub, 1, []);
end
end

function value = localPenaltyObjective(x, Vio, Obj)
[value, ~] = CostFunction(x, Vio, Obj);
end

% Recover objective values and constraint violations for all successful runs.
function [objectiveValues, penaltyValues, maxViolations, feasibleFlags] = evaluateRunResults(bestPositions, successMask, Vio, Obj, EvalObj)
numRuns = size(bestPositions, 1);
objectiveValues = NaN(numRuns, 1);
penaltyValues = NaN(numRuns, 1);
maxViolations = NaN(numRuns, 1);
feasibleFlags = false(numRuns, 1);
for r = 1:numRuns
    if successMask(r)
        bestPos = reshape(bestPositions(r, :), 1, []);
        if all(isfinite(bestPos))
            [objectiveValue, gValues, hValues] = EvalObj(bestPos);
            [~, evalData] = CostFunction(bestPos, Vio, Obj);
            gValues = normalizeConstraintVector(gValues);
            hValues = normalizeConstraintVector(hValues);
            allViolations = [max(0, gValues), abs(hValues)];
            if isempty(allViolations)
                maxViolation = 0;
            else
                maxViolation = max(allViolations);
            end
            objectiveValues(r) = objectiveValue;
            penaltyValues(r) = evalData.v;
            maxViolations(r) = maxViolation;
            feasibleFlags(r) = maxViolation <= 1e-8;
        end
    end
end
end

function resultSummary = summarizeRuns(bestScores, runtimes, objectiveValues, penaltyValues, maxViolations, feasibleFlags, successMask)
resultSummary = struct('Best', NaN, 'Worst', NaN, 'Mean', NaN, 'Std', NaN, 'AverageTime', NaN, ...
    'FeasibleRate', NaN, 'BestRun', NaN, 'BestScore', NaN, 'BestObjective', NaN, 'BestPenalty', NaN, ...
    'BestMaxViolation', NaN, 'BestIsFeasible', false);

validRuns = find(successMask);
if isempty(validRuns)
    return
end

validTimes = runtimes(validRuns);
validTimes = validTimes(isfinite(validTimes));
if ~isempty(validTimes)
    resultSummary.AverageTime = mean(validTimes);
end
resultSummary.FeasibleRate = mean(feasibleFlags(validRuns));

if any(feasibleFlags(validRuns))
    feasibleRuns = validRuns(feasibleFlags(validRuns));
    feasibleObjectives = objectiveValues(feasibleRuns);
    [resultSummary.Best, bestLocalIdx] = min(feasibleObjectives);
    resultSummary.BestRun = feasibleRuns(bestLocalIdx);
    resultSummary.Worst = max(feasibleObjectives);
    resultSummary.Mean = mean(feasibleObjectives);
    resultSummary.Std = std(feasibleObjectives);
else
    fallbackBasis = [maxViolations(validRuns), bestScores(validRuns)];
    [~, sortIdx] = sortrows(fallbackBasis, [1 2]);
    resultSummary.BestRun = validRuns(sortIdx(1));
    validScores = bestScores(validRuns);
    resultSummary.Best = bestScores(resultSummary.BestRun);
    resultSummary.Worst = max(validScores);
    resultSummary.Mean = mean(validScores);
    resultSummary.Std = std(validScores);
end

resultSummary.BestScore = bestScores(resultSummary.BestRun);
resultSummary.BestObjective = objectiveValues(resultSummary.BestRun);
resultSummary.BestPenalty = penaltyValues(resultSummary.BestRun);
resultSummary.BestMaxViolation = maxViolations(resultSummary.BestRun);
resultSummary.BestIsFeasible = feasibleFlags(resultSummary.BestRun);
end

% Build detailed information for the best run of the current problem.
function bestRunDetail = buildRunDetail(bestRunIndex, bestScore, bestPositions, Vio, Obj, EvalObj)
bestRunDetail = struct('BestRun', NaN, 'BestScore', NaN, 'BestPosition', [], 'Variables', struct(), ...
    'ObjectiveValue', NaN, 'PenaltyValue', NaN, 'MaxViolation', NaN, 'IsFeasible', false, ...
    'ConstraintValues', [], 'Constraints', struct(), 'EqualityValues', [], 'Equalities', struct());
if ~isfinite(bestRunIndex)
    return
end

bestPos = reshape(bestPositions(bestRunIndex, :), 1, []);
if ~all(isfinite(bestPos))
    return
end

[objectiveValue, gValues, hValues] = EvalObj(bestPos);
[~, evalData] = CostFunction(bestPos, Vio, Obj);
gValues = normalizeConstraintVector(gValues);
hValues = normalizeConstraintVector(hValues);
allViolations = [max(0, gValues), abs(hValues)];
if isempty(allViolations)
    maxViolation = 0;
else
    maxViolation = max(allViolations);
end

bestRunDetail.BestRun = bestRunIndex;
bestRunDetail.BestScore = bestScore;
bestRunDetail.BestPosition = bestPos;
bestRunDetail.Variables = vectorToNamedStruct(bestPos, 'x');
bestRunDetail.ObjectiveValue = objectiveValue;
bestRunDetail.PenaltyValue = evalData.v;
bestRunDetail.MaxViolation = maxViolation;
bestRunDetail.IsFeasible = maxViolation <= 1e-8;
bestRunDetail.ConstraintValues = gValues;
bestRunDetail.Constraints = vectorToNamedStruct(gValues, 'g');
bestRunDetail.EqualityValues = hValues;
bestRunDetail.Equalities = vectorToNamedStruct(hValues, 'h');
end

function variableStruct = vectorToNamedStruct(values, prefix)
if nargin < 2
    prefix = 'x';
end
variableStruct = struct();
values = reshape(values, 1, []);
for idx = 1:numel(values)
    variableStruct.(sprintf('%s%d', prefix, idx)) = values(idx);
end
end

function values = normalizeConstraintVector(values)
if isempty(values)
    values = [];
    return
end
values = reshape(values, 1, []);
if isscalar(values) && abs(values) <= eps
    values = [];
end
end

function textValue = formatScoreForDisplay(value, useScientificFormat)
if ~isfinite(value)
    textValue = 'NaN';
elseif useScientificFormat
    textValue = sprintf('%.4E', value);
else
    textValue = sprintf('%.4f', value);
end
end
