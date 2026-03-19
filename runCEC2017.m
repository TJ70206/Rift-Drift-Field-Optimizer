clc;
clear;
close all;

baseDir = fileparts(mfilename('fullpath'));
oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));
cd(baseDir);
addpath(genpath(baseDir));

% Manual setting
functionId = 1;
dim = 30;
runs = 30;
popSize = 50;
maxIter = 3000;
baseSeed = 42;

if ~isscalar(functionId) || functionId < 1 || functionId > 30 || functionId ~= fix(functionId)
    error('functionId must be an integer between 1 and 30.');
end

% Prepare benchmark executable
prepareCec17Mex(baseDir);

options.PopulationSize = popSize;
options.MaxIterations = maxIter;
resultsRoot = fullfile(baseDir, 'results');
if ~exist(resultsRoot, 'dir')
    mkdir(resultsRoot);
end

[lb, ub, dimUsed, fobj] = Get_Functions_cec2017(functionId, dim);
lb = reshape(lb, 1, []);
ub = reshape(ub, 1, []);

funcDir = fullfile(resultsRoot, sprintf('F%d_D%d', functionId, dimUsed));
if ~exist(funcDir, 'dir')
    mkdir(funcDir);
end

bestScores = NaN(runs, 1);
runtimes = NaN(runs, 1);
allCurves = NaN(runs, maxIter);
bestPositions = NaN(runs, dimUsed);
runSeeds = baseSeed + functionId * 1000 + (1:runs);

% Run RDFO
fprintf('\nRunning CEC2017 F%d (D=%d)\n', functionId, dimUsed);
for r = 1:runs
    rng(runSeeds(r), 'twister');
    try
        tic;
        [bestPos, bestScore, rawCurve] = RDFO(fobj, dimUsed, lb, ub, options);
        runtime = toc;
        curve = normalizeCurve(rawCurve, maxIter);
        if ~isfinite(bestScore) && any(isfinite(curve))
            bestScore = curve(find(isfinite(curve), 1, 'last'));
        end
        bestScores(r) = bestScore;
        runtimes(r) = runtime;
        allCurves(r, :) = curve;
        bestPositions(r, :) = reshape(bestPos, 1, []);
        fprintf('  Run %02d/%02d: %.2E | time = %.2fs\n', r, runs, bestScore, runtime);
    catch ME
        fprintf('  Run %02d/%02d failed: %s\n', r, runs, ME.message);
    end
end

resultSummary = buildSummary(bestScores, runtimes);

% Save results
save(fullfile(funcDir, 'result.mat'), ...
    'functionId', 'dimUsed', 'lb', 'ub', 'options', 'runs', 'runSeeds', ...
    'bestScores', 'runtimes', 'allCurves', 'bestPositions', 'resultSummary');

fprintf('  best = %.2E | mean time = %.2fs\n', ...
    resultSummary.Best, resultSummary.MeanTime);

function curve = normalizeCurve(rawCurve, targetLength)
curve = NaN(1, targetLength);
if isempty(rawCurve)
    return;
end
rawCurve = reshape(rawCurve, 1, []);
finiteMask = isfinite(rawCurve);
if ~any(finiteMask)
    return;
end
rawCurve = rawCurve(1:find(finiteMask, 1, 'last'));
if isscalar(rawCurve)
    curve(:) = rawCurve;
else
    sampleIdx = round(linspace(1, numel(rawCurve), targetLength));
    sampleIdx(sampleIdx < 1) = 1;
    curve = rawCurve(sampleIdx);
end
firstFinite = find(isfinite(curve), 1, 'first');
if isempty(firstFinite)
    curve = NaN(1, targetLength);
    return;
end
if firstFinite > 1
    curve(1:firstFinite - 1) = curve(firstFinite);
end
for k = 2:targetLength
    if ~isfinite(curve(k))
        curve(k) = curve(k - 1);
    end
end
curve = cummin(curve);
end

function resultSummary = buildSummary(bestScores, runtimes)
validScoreMask = isfinite(bestScores);
validTimeMask = isfinite(runtimes);
validScores = bestScores(validScoreMask);
validTimes = runtimes(validTimeMask);
resultSummary.Best = safeStat(validScores, @min);
resultSummary.Worst = safeStat(validScores, @max);
resultSummary.Mean = safeStat(validScores, @mean);
resultSummary.Std = safeStat(validScores, @std);
resultSummary.MeanTime = safeStat(validTimes, @mean);
resultSummary.SuccessRuns = sum(validScoreMask);
end

function value = safeStat(data, funcHandle)
if isempty(data)
    value = NaN;
else
    value = funcHandle(data);
end
end

function prepareCec17Mex(baseDir)
cecMexFile = fullfile(baseDir, ['cec17_func.', mexext]);
if exist(cecMexFile, 'file') == 2
    return;
end

workspaceRoot = fileparts(fileparts(fileparts(baseDir)));
candidateFiles = dir(fullfile(workspaceRoot, '**', ['cec17_func.', mexext]));
candidateFiles = candidateFiles(~contains(string({candidateFiles.folder}), string(baseDir)));
if ~isempty(candidateFiles)
    addpath(candidateFiles(1).folder, '-begin');
    return;
end

cppFile = fullfile(baseDir, 'cec17_func.cpp');
if exist(cppFile, 'file') ~= 2
    error('cec17_func executable was not found, and cec17_func.cpp is also missing.');
end

compilerConfig = mex.getCompilerConfigurations('C++', 'Installed');
if isempty(compilerConfig)
    error(['No supported C/C++ compiler was detected, and no prebuilt ', ['cec17_func.', mexext], ...
        ' was found. Please copy an existing compiled cec17_func file into this folder.']);
end

mex('-output', 'cec17_func', cppFile);
end
