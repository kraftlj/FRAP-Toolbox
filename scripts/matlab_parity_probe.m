function matlab_parity_probe()
%MATLAB_PARITY_PROBE Collect MATLAB-only parity evidence for FRAP-Toolbox.
%
% Run from the repository root with:
%   addpath(fullfile(pwd, 'scripts'));
%   matlab_parity_probe;
%
% The probe writes only to scratch/matlab-parity-output/.

repoRoot = resolveRepoRoot();
addpath(repoRoot);

outputDir = fullfile(repoRoot, 'scratch', 'matlab-parity-output');
ensureDir(outputDir);

diaryPath = fullfile(outputDir, 'matlab_parity_probe_console.txt');
if exist(diaryPath, 'file') == 2
    delete(diaryPath);
end
diary(diaryPath);
cleanup = onCleanup(@() diary('off')); %#ok<NASGU>

fprintf('FRAP-Toolbox MATLAB parity probe\n');
fprintf('Repository root: %s\n', repoRoot);
fprintf('Output directory: %s\n\n', outputDir);

sectionStatus = {};

try
    writeEnvironmentFacts(outputDir, repoRoot);
    sectionStatus{end + 1} = 'environment: ok'; %#ok<AGROW>
catch err
    sectionStatus{end + 1} = ['environment: failed - ', err.message]; %#ok<AGROW>
    reportSectionError('environment', err);
end

try
    writeDiffusionRoiMasks(outputDir);
    sectionStatus{end + 1} = 'diffusion ROI masks: ok'; %#ok<AGROW>
catch err
    sectionStatus{end + 1} = ['diffusion ROI masks: failed - ', err.message]; %#ok<AGROW>
    reportSectionError('diffusion ROI masks', err);
end

try
    runDiffusionOptimizerProbes(outputDir, repoRoot);
    sectionStatus{end + 1} = 'diffusion optimizer traces: ok'; %#ok<AGROW>
catch err
    sectionStatus{end + 1} = ['diffusion optimizer traces: failed - ', err.message]; %#ok<AGROW>
    reportSectionError('diffusion optimizer traces', err);
end

try
    copied = copyOptionalReactionExports(outputDir, repoRoot);
    sectionStatus{end + 1} = sprintf('optional reaction exports copied: %d', copied); %#ok<AGROW>
catch err
    sectionStatus{end + 1} = ['optional reaction exports: failed - ', err.message]; %#ok<AGROW>
    reportSectionError('optional reaction exports', err);
end

writeProbeReadme(outputDir, sectionStatus);
fprintf('\nProbe complete. Review %s\n', fullfile(outputDir, 'README.txt'));
end

function repoRoot = resolveRepoRoot()
scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    scriptDir = pwd;
else
    scriptDir = fileparts(scriptPath);
end

candidates = {fileparts(scriptDir), scriptDir, pwd};
repoRoot = '';
for idx = 1:numel(candidates)
    candidate = candidates{idx};
    if exist(fullfile(candidate, 'DiffusionModel_2.m'), 'file') == 2 ...
            && exist(fullfile(candidate, 'ROIinitialization_Diffusion.m'), 'file') == 2
        repoRoot = candidate;
        break;
    end
end

if isempty(repoRoot)
    error('Could not locate the FRAP-Toolbox repository root.');
end
end

function ensureDir(pathToCreate)
if exist(pathToCreate, 'dir') ~= 7
    mkdir(pathToCreate);
end
end

function reportSectionError(sectionName, err)
fprintf('\nERROR in %s:\n%s\n', sectionName, err.message);
for idx = 1:numel(err.stack)
    fprintf('  at %s line %d\n', err.stack(idx).name, err.stack(idx).line);
end
fprintf('\n');
end

function writeEnvironmentFacts(outputDir, repoRoot)
symbols = {'lsqcurvefit', 'lsqnonlin', 'poly2mask', 'roipoly', 'bfopen'};
versionText = version;
versionDetails = ver; %#ok<NASGU>
defaultLsqcurvefitOptions = optimset('lsqcurvefit'); %#ok<NASGU>
defaultLsqnonlinOptions = optimset('lsqnonlin'); %#ok<NASGU>
whichOutput = cell(numel(symbols), 2);

textPath = fullfile(outputDir, 'environment.txt');
fid = fopen(textPath, 'w');
if fid < 0
    error('Could not write %s', textPath);
end

fprintf(fid, 'MATLAB version: %s\n', versionText);
fprintf(fid, 'Computer: %s\n', computer);
fprintf(fid, 'Repository root: %s\n\n', repoRoot);

fprintf(fid, 'Installed toolbox versions:\n');
for idx = 1:numel(versionDetails)
    release = '';
    if isfield(versionDetails, 'Release')
        release = versionDetails(idx).Release;
    end
    fprintf(fid, '- %s | Version %s | %s | %s\n', ...
        versionDetails(idx).Name, versionDetails(idx).Version, release, versionDetails(idx).Date);
end

fprintf(fid, '\nwhich -all output:\n');
for idx = 1:numel(symbols)
    symbol = symbols{idx};
    try
        text = evalc(['which -all ', symbol]);
    catch err
        text = ['ERROR: ', err.message];
    end
    whichOutput{idx, 1} = symbol;
    whichOutput{idx, 2} = text;
    fprintf(fid, '\n--- %s ---\n%s\n', symbol, text);
end
fclose(fid);

save(fullfile(outputDir, 'environment.mat'), ...
    'versionText', 'versionDetails', 'whichOutput', ...
    'defaultLsqcurvefitOptions', 'defaultLsqnonlinOptions');

writeStructText(fullfile(outputDir, 'default_lsqcurvefit_options.txt'), defaultLsqcurvefitOptions);
writeStructText(fullfile(outputDir, 'default_lsqnonlin_options.txt'), defaultLsqnonlinOptions);
end

function writeStructText(path, value)
fid = fopen(path, 'w');
if fid < 0
    error('Could not write %s', path);
end
fields = fieldnames(value);
for idx = 1:numel(fields)
    field = fields{idx};
    text = formatValue(value.(field));
    fprintf(fid, '%s: %s\n', field, text);
end
fclose(fid);
end

function text = formatValue(value)
if isempty(value)
    text = '<empty>';
elseif isnumeric(value)
    if isscalar(value)
        text = sprintf('%.17g', value);
    else
        sizeText = sprintf('%dx', size(value));
        text = ['numeric[', sizeText(1:end - 1), ']'];
    end
elseif ischar(value)
    text = value;
elseif isa(value, 'function_handle')
    text = func2str(value);
elseif iscell(value)
    text = sprintf('cell[%d]', numel(value));
else
    text = class(value);
end
end

function writeDiffusionRoiMasks(outputDir)
x0 = 256;
y0 = 23;
R0 = 9;
imageShape = [512, 512]; %#ok<NASGU>
t = 0:pi/20:2*pi; %#ok<NASGU>
xi = R0 * cos(t) + x0; %#ok<NASGU>
yi = R0 * sin(t) + y0; %#ok<NASGU>
bleachroimask = poly2mask(xi, yi, imageShape(1), imageShape(2)); %#ok<NASGU>
adjacentXi = R0 * cos(t) + x0 + R0 * 2.5; %#ok<NASGU>
adjacentYi = yi; %#ok<NASGU>
adjacentroimask = poly2mask(adjacentXi, adjacentYi, imageShape(1), imageShape(2)); %#ok<NASGU>

save(fullfile(outputDir, 'roi_masks.mat'), ...
    'x0', 'y0', 'R0', 'imageShape', 't', 'xi', 'yi', ...
    'bleachroimask', 'adjacentXi', 'adjacentYi', 'adjacentroimask');

writeMaskCoordinates(fullfile(outputDir, 'diffusion_bleach_mask_coordinates.csv'), bleachroimask);
writeMaskCoordinates(fullfile(outputDir, 'diffusion_adjacent_mask_coordinates.csv'), adjacentroimask);
writeRowIntervals(fullfile(outputDir, 'diffusion_bleach_mask_row_intervals.csv'), bleachroimask);
writeRowIntervals(fullfile(outputDir, 'diffusion_adjacent_mask_row_intervals.csv'), adjacentroimask);

summaryPath = fullfile(outputDir, 'roi_mask_summary.txt');
fid = fopen(summaryPath, 'w');
if fid < 0
    error('Could not write %s', summaryPath);
end
fprintf(fid, 'Diffusion ROI masks from ROIinitialization_Diffusion.m\n');
fprintf(fid, 'x0=%g, y0=%g, R0=%g, image shape=%dx%d\n', x0, y0, R0, imageShape(1), imageShape(2));
fprintf(fid, 'bleach nnz=%d\n', nnz(bleachroimask));
fprintf(fid, 'adjacent nnz=%d\n', nnz(adjacentroimask));
fprintf(fid, '\nBleach row intervals:\n');
writeRowIntervalsToOpenFile(fid, bleachroimask);
fprintf(fid, '\nAdjacent row intervals:\n');
writeRowIntervalsToOpenFile(fid, adjacentroimask);
fclose(fid);
end

function writeMaskCoordinates(path, mask)
[rows, cols] = find(mask);
fid = fopen(path, 'w');
if fid < 0
    error('Could not write %s', path);
end
fprintf(fid, 'row,col\n');
for idx = 1:numel(rows)
    fprintf(fid, '%d,%d\n', rows(idx), cols(idx));
end
fclose(fid);
end

function writeRowIntervals(path, mask)
fid = fopen(path, 'w');
if fid < 0
    error('Could not write %s', path);
end
fprintf(fid, 'row,start_col,end_col,count\n');
writeRowIntervalsToOpenFile(fid, mask);
fclose(fid);
end

function writeRowIntervalsToOpenFile(fid, mask)
for row = 1:size(mask, 1)
    cols = find(mask(row, :));
    if ~isempty(cols)
        fprintf(fid, '%d,%d,%d,%d\n', row, cols(1), cols(end), numel(cols));
    end
end
end

function runDiffusionOptimizerProbes(outputDir, repoRoot)
diffusionDir = fullfile(repoRoot, 'test-data', 'Diffusion');
paramsPath = fullfile(diffusionDir, 'Venus_Cytoplasm_Diffusion_Fit_Parameters.txt');
frapPath = fullfile(diffusionDir, 'Venus_Cytoplasm_Diffusion_FRAP_datasets.txt');

if exist(paramsPath, 'file') ~= 2
    error('Missing parameter table: %s', paramsPath);
end
if exist(frapPath, 'file') ~= 2
    error('Missing FRAP export: %s', frapPath);
end

params = readParameterTable(paramsPath);
frapVectors = readVectorExport(frapPath);
targets = {'Venus_Cytoplasm_3.lsm', 'Venus_Cytoplasm_5.lsm', 'Venus_Cytoplasm_9.lsm'};

summaryPath = fullfile(outputDir, 'diffusion_optimizer_summary.csv');
fid = fopen(summaryPath, 'w');
if fid < 0
    error('Could not write %s', summaryPath);
end
fprintf(fid, ['file,guide_D,guide_MF,guide_SS,computed_guide_SS,normal_D,normal_MF,', ...
    'normal_SS,normal_exitflag,traced_D,traced_MF,traced_SS,traced_exitflag,trace_rows\n']);

results = struct([]);
for idx = 1:numel(targets)
    result = runOneDiffusionProbe(outputDir, params, frapVectors, targets{idx});
    results = appendStruct(results, result); %#ok<AGROW>
    fprintf(fid, '%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%.17g,%.17g,%.17g,%d,%d\n', ...
        result.file, result.guide.D, result.guide.MF, result.guide.SS, result.guide.computedSS, ...
        result.normal.params(1), result.normal.params(2), result.normal.resnorm, result.normal.exitflag, ...
        result.traced.params(1), result.traced.params(2), result.traced.resnorm, result.traced.exitflag, ...
        result.traceRowCount);
end
fclose(fid);

save(fullfile(outputDir, 'diffusion_optimizer_results.mat'), 'results');
end

function out = appendStruct(in, value)
if isempty(in)
    out = value;
else
    out = in;
    out(end + 1) = value;
end
end

function result = runOneDiffusionProbe(outputDir, params, frapVectors, fileName)
postBleachImage = 21;
voxelSizeX = 0.10986328426900892;
nominalRadius = 9 * voxelSizeX;

guide.k = getTableValue(params, fileName, 'k');
guide.re = getTableValue(params, fileName, 're');
guide.D = getTableValue(params, fileName, 'D');
guide.MF = getTableValue(params, fileName, 'MF');
guide.SS = getTableValue(params, fileName, 'SS');

time = getVector(frapVectors, fileName, 'Time');
correctedFrap = getVector(frapVectors, fileName, 'Corrected FRAP');
fitTime = time(postBleachImage:end) - time(postBleachImage);
f = correctedFrap(postBleachImage:end);
weights = sum(f);
weightedTarget = f ./ (fitTime + weights);
guide.computedSS = diffusionWeightedSse(guide.D, guide.MF, fitTime, weightedTarget, ...
    guide.k, guide.re, nominalRadius, f(1), weights);

options = optimset('lsqcurvefit');
options.Display = 'off';
model = @(p, t) diffusionWeightedModel(t, p, guide.k, guide.re, nominalRadius, f(1), weights);

[normalParams, normalResnorm, normalResidual, normalExitflag, normalOutput] = ...
    lsqcurvefit(model, [10, 1], fitTime, weightedTarget, [0, 0], [Inf, 2], options);

traceRows = struct('state', {}, 'x', {}, 'optimValues', {});
traceOptions = optimset(options, 'OutputFcn', @traceOutputFcn);
[tracedParams, tracedResnorm, tracedResidual, tracedExitflag, tracedOutput] = ...
    lsqcurvefit(model, [10, 1], fitTime, weightedTarget, [0, 0], [Inf, 2], traceOptions);

safeName = regexprep(fileName, '[^A-Za-z0-9_]', '_');
traceCsvPath = fullfile(outputDir, [safeName, '_lsqcurvefit_trace.csv']);
writeTraceCsv(traceCsvPath, traceRows);

matPath = fullfile(outputDir, [safeName, '_lsqcurvefit_probe.mat']);
normal.params = normalParams;
normal.resnorm = normalResnorm;
normal.residual = normalResidual;
normal.exitflag = normalExitflag;
normal.output = normalOutput;
traced.params = tracedParams;
traced.resnorm = tracedResnorm;
traced.residual = tracedResidual;
traced.exitflag = tracedExitflag;
traced.output = tracedOutput;
save(matPath, ...
    'fileName', 'postBleachImage', 'voxelSizeX', 'nominalRadius', ...
    'guide', 'time', 'correctedFrap', 'fitTime', 'f', 'weights', ...
    'weightedTarget', 'options', 'normal', 'traced', 'traceRows');

result.file = fileName;
result.guide = guide;
result.normal = normal;
result.traced = traced;
result.traceRowCount = numel(traceRows);
result.traceCsvPath = traceCsvPath;
result.matPath = matPath;

    function stop = traceOutputFcn(x, optimValues, state)
        stop = false;
        rowIdx = numel(traceRows) + 1;
        traceRows(rowIdx).state = state;
        traceRows(rowIdx).x = x(:)';
        traceRows(rowIdx).optimValues = optimValues;
    end
end

function y = diffusionWeightedModel(t, p, k, re, rn, f0, weights)
D = p(1);
MF = p(2);
fit = kangFrapProbe(t, re, rn, D, k) .* MF + (1 - MF) .* f0;
y = fit ./ (t + weights);
end

function sse = diffusionWeightedSse(D, MF, t, target, k, re, rn, f0, weights)
prediction = diffusionWeightedModel(t, [D, MF], k, re, rn, f0, weights);
residual = prediction - target;
sse = sum(residual .^ 2);
end

function y = kangFrapProbe(t, re, rn, D, K)
t = t(:);
m = (0:10)';
a = ((-K) .^ m .* re .^ 2) ./ factorial(m);
b = re .^ 2 + m * (8 .* D .* t' + rn .^ 2);
c = (a * ones(1, numel(t))) ./ b;
y = sum(c, 1)';
end

function writeTraceCsv(path, traceRows)
fid = fopen(path, 'w');
if fid < 0
    error('Could not write %s', path);
end
fprintf(fid, 'row,state,D,MF,iteration,funccount,resnorm,firstorderopt,stepsize,procedure\n');
for idx = 1:numel(traceRows)
    x = traceRows(idx).x;
    if numel(x) < 2
        x = [NaN, NaN];
    end
    optimValues = traceRows(idx).optimValues;
    procedure = safeOptimText(optimValues, 'procedure');
    procedure = strrep(procedure, ',', ';');
    fprintf(fid, '%d,%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%s\n', ...
        idx, traceRows(idx).state, x(1), x(2), ...
        safeOptimNumber(optimValues, 'iteration'), ...
        safeOptimNumber(optimValues, 'funccount'), ...
        safeOptimNumber(optimValues, 'resnorm'), ...
        safeOptimNumber(optimValues, 'firstorderopt'), ...
        safeOptimNumber(optimValues, 'stepsize'), ...
        procedure);
end
fclose(fid);
end

function value = safeOptimNumber(optimValues, field)
value = NaN;
if isstruct(optimValues) && isfield(optimValues, field)
    raw = optimValues.(field);
    if isnumeric(raw) && isscalar(raw)
        value = raw;
    end
end
end

function text = safeOptimText(optimValues, field)
text = '';
if isstruct(optimValues) && isfield(optimValues, field)
    raw = optimValues.(field);
    if ischar(raw)
        text = raw;
    elseif isnumeric(raw) && isscalar(raw)
        text = sprintf('%.17g', raw);
    else
        text = class(raw);
    end
end
end

function table = readParameterTable(path)
lines = readTextLines(path);
if isempty(lines)
    error('Parameter table is empty: %s', path);
end
headers = splitTabLine(lines{1});
headers = cleanTextCells(headers);
names = {};
values = [];
for idx = 2:numel(lines)
    cells = splitTabLine(lines{idx});
    if isempty(cells) || isempty(strtrim(cells{1}))
        continue;
    end
    row = NaN(1, numel(headers) - 1);
    for col = 2:numel(headers)
        if col <= numel(cells)
            text = strtrim(cells{col});
            if ~isempty(text)
                row(col - 1) = str2double(text);
            end
        end
    end
    names{end + 1} = cleanText(cells{1}); %#ok<AGROW>
    values(end + 1, :) = row; %#ok<AGROW>
end
table.headers = headers(2:end);
table.names = names;
table.values = values;
end

function vectors = readVectorExport(path)
lines = readTextLines(path);
vectors = struct('name', {}, 'label', {}, 'values', {});
for idx = 1:numel(lines)
    cells = splitTabLine(lines{idx});
    if numel(cells) < 3
        continue;
    end
    values = [];
    for col = 3:numel(cells)
        text = strtrim(cells{col});
        if ~isempty(text)
            values(end + 1, 1) = str2double(text); %#ok<AGROW>
        end
    end
    vectorIdx = numel(vectors) + 1;
    vectors(vectorIdx).name = cleanText(cells{1});
    vectors(vectorIdx).label = cleanText(cells{2});
    vectors(vectorIdx).values = values;
end
end

function lines = readTextLines(path)
fid = fopen(path, 'r');
if fid < 0
    error('Could not read %s', path);
end
lines = {};
line = fgetl(fid);
while ischar(line)
    lines{end + 1} = line; %#ok<AGROW>
    line = fgetl(fid);
end
fclose(fid);
end

function cells = splitTabLine(line)
cells = regexp(line, '\t', 'split');
end

function cells = cleanTextCells(cells)
for idx = 1:numel(cells)
    cells{idx} = cleanText(cells{idx});
end
end

function text = cleanText(text)
text = strtrim(text);
if ~isempty(text) && double(text(1)) == 65279
    text = text(2:end);
    text = strtrim(text);
end
end

function value = getTableValue(table, rowName, headerName)
rowIdx = findTextIndex(table.names, rowName);
colIdx = findTextIndex(table.headers, headerName);
if isempty(rowIdx)
    error('Missing row %s in parameter table.', rowName);
end
if isempty(colIdx)
    error('Missing column %s in parameter table.', headerName);
end
value = table.values(rowIdx, colIdx);
end

function values = getVector(vectors, rowName, label)
values = [];
for idx = 1:numel(vectors)
    if strcmpi(cleanText(vectors(idx).name), cleanText(rowName)) ...
            && strcmpi(cleanText(vectors(idx).label), cleanText(label))
        values = vectors(idx).values;
        return;
    end
end
error('Missing vector row: %s / %s', rowName, label);
end

function idx = findTextIndex(values, target)
idx = [];
target = cleanText(target);
for valueIdx = 1:numel(values)
    if strcmpi(cleanText(values{valueIdx}), target)
        idx = valueIdx;
        return;
    end
end
end

function copied = copyOptionalReactionExports(outputDir, repoRoot)
copied = 0;
targetDir = fullfile(outputDir, 'optional-reaction-exports');
ensureDir(targetDir);

reactionDirs = { ...
    fullfile(repoRoot, 'test-data', 'Reaction 1'), ...
    fullfile(repoRoot, 'test-data', 'Reaction 2') ...
};
patterns = {'*_Reaction_FRAP_datasets.txt', '*_Reaction2_FRAP_datasets.txt'};

for dirIdx = 1:numel(reactionDirs)
    reactionDir = reactionDirs{dirIdx};
    if exist(reactionDir, 'dir') ~= 7
        continue;
    end
    for patternIdx = 1:numel(patterns)
        matches = dir(fullfile(reactionDir, patterns{patternIdx}));
        for matchIdx = 1:numel(matches)
            source = fullfile(reactionDir, matches(matchIdx).name);
            destination = fullfile(targetDir, matches(matchIdx).name);
            copyfile(source, destination);
            copied = copied + 1;
        end
    end
end
end

function writeProbeReadme(outputDir, sectionStatus)
readmePath = fullfile(outputDir, 'README.txt');
fid = fopen(readmePath, 'w');
if fid < 0
    error('Could not write %s', readmePath);
end
fprintf(fid, 'FRAP-Toolbox MATLAB parity probe output\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now));
fprintf(fid, 'Section status:\n');
for idx = 1:numel(sectionStatus)
    fprintf(fid, '- %s\n', sectionStatus{idx});
end
fprintf(fid, '\nExpected key files:\n');
fprintf(fid, '- environment.txt\n');
fprintf(fid, '- default_lsqcurvefit_options.txt\n');
fprintf(fid, '- roi_masks.mat\n');
fprintf(fid, '- roi_mask_summary.txt\n');
fprintf(fid, '- diffusion_optimizer_summary.csv\n');
fprintf(fid, '- Venus_Cytoplasm_3_lsm_lsqcurvefit_probe.mat and trace CSV\n');
fprintf(fid, '- Venus_Cytoplasm_5_lsm_lsqcurvefit_probe.mat and trace CSV\n');
fprintf(fid, '- Venus_Cytoplasm_9_lsm_lsqcurvefit_probe.mat and trace CSV\n');
fprintf(fid, '\nReturn this entire directory to the FRAP-Toolbox modernization owner.\n');
fclose(fid);
end
