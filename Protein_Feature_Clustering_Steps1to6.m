%% Protein Feature Extraction and Clustering (Steps 1-6)
% Input: ProteinData.xlsx
% Required columns: ProteinID, Sequence, PfamDomains
% Output: ProteinClusteringResults.mat, SelectedFeatures.xlsx
%
% Notes:
% - Uses standard amino acids and overlapping 3-mers (8000 possible features).
% - Standard Laplacian Score favors LOWER scores; this script selects the
%   92 features with the smallest finite scores.
% - Exact k-nearest-neighbor search can be expensive for very large datasets.

clear; clc; close all;

%% STEP 1: Input data integration
proteinFile = 'ProteinData.xlsx';
proteins = readtable(proteinFile, 'VariableNamingRule', 'preserve');

requiredColumns = {'ProteinID', 'Sequence', 'PfamDomains'};
if ~all(ismember(requiredColumns, proteins.Properties.VariableNames))
    error('ProteinData.xlsx must contain columns: ProteinID, Sequence, PfamDomains.');
end

proteins.ProteinID = string(proteins.ProteinID);
proteins.Sequence = upper(string(proteins.Sequence));
proteins.PfamDomains = string(proteins.PfamDomains);

fprintf('Imported %d protein records.\\n', height(proteins));

%% STEP 2: Data cleaning and redundancy removal
validRows = ~ismissing(proteins.ProteinID) & ...
    strlength(strtrim(proteins.ProteinID)) > 0 & ...
    ~ismissing(proteins.Sequence) & ...
    strlength(strtrim(proteins.Sequence)) > 0;

proteins = proteins(validRows, :);
proteins.Sequence = regexprep(proteins.Sequence, '\\s+', '');

% Keep sequences composed only of the 20 standard amino acids.
validAA = cellfun(@(s) all(ismember(s, 'ACDEFGHIKLMNPQRSTVWY')), ...
    cellstr(proteins.Sequence));
proteins = proteins(validAA, :);

% Remove duplicate sequences, then duplicate identifiers.
[~, idx] = unique(proteins.Sequence, 'stable');
proteins = proteins(idx, :);
[~, idx] = unique(proteins.ProteinID, 'stable');
proteins = proteins(idx, :);

fprintf('%d proteins remain after cleaning.\\n', height(proteins));

%% STEP 3: Feature extraction using overlapping 3-mers
sequences = cellstr(proteins.Sequence);
alphabet = 'ACDEFGHIKLMNPQRSTVWY';
k = 3;
numFeatures = numel(alphabet)^k;
numProteins = numel(sequences);

% Generate all 8000 possible 3-mers in a fixed order.
allKmers = strings(numFeatures, 1);
idx = 1;
for a = 1:numel(alphabet)
    for b = 1:numel(alphabet)
        for c = 1:numel(alphabet)
            allKmers(idx) = string([alphabet(a), alphabet(b), alphabet(c)]);
            idx = idx + 1;
        end
    end
end

% Store normalized 3-mer frequencies as single precision to reduce memory.
X = zeros(numProteins, numFeatures, 'single');
aaMap = zeros(1, 256);
for a = 1:numel(alphabet)
    aaMap(double(alphabet(a)) + 1) = a;
end

for i = 1:numProteins
    seq = sequences{i};
    if numel(seq) < k
        continue;
    end

    counts = zeros(1, numFeatures, 'single');
    for j = 1:(numel(seq) - k + 1)
        chars = double(seq(j:j+k-1)) + 1;
        a = aaMap(chars(1));
        b = aaMap(chars(2));
        c = aaMap(chars(3));

        % Skip any unexpected characters safely.
        if a == 0 || b == 0 || c == 0
            continue;
        end

        featureIdx = (a-1)*numel(alphabet)^2 + ...
                     (b-1)*numel(alphabet) + c;
        counts(featureIdx) = counts(featureIdx) + 1;
    end

    validWindows = numel(seq) - k + 1;
    X(i, :) = counts / max(validWindows, 1);

    if mod(i, 1000) == 0
        fprintf('Processed 3-mers for %d/%d proteins.\\n', i, numProteins);
    end
end

fprintf('Feature matrix: %d proteins x %d 3-mer features.\\n', ...
    size(X,1), size(X,2));

%% STEP 4: Laplacian Score using a sparse k-nearest-neighbor graph
% NOTE: Exact knnsearch on hundreds of thousands of proteins is costly.
% Test on a smaller subset first, or replace with an approximate NN method.
numNeighbors = 10;

% knnsearch requires Statistics and Machine Learning Toolbox.
% Convert to double for compatibility with common MATLAB versions.
Xdouble = double(X);
neighborIdx = knnsearch(Xdouble, Xdouble, 'K', numNeighbors + 1);

n = size(Xdouble, 1);
edgeCount = n * numNeighbors;
rowIdx = zeros(edgeCount, 1);
colIdx = zeros(edgeCount, 1);
weights = zeros(edgeCount, 1);
pos = 1;

for i = 1:n
    neighbors = neighborIdx(i, 2:end);
    for j = 1:numNeighbors
        nb = neighbors(j);
        delta = Xdouble(i,:) - Xdouble(nb,:);
        distSquared = sum(delta .* delta);
        rowIdx(pos) = i;
        colIdx(pos) = nb;
        weights(pos) = exp(-distSquared);
        pos = pos + 1;
    end
    if mod(i, 1000) == 0
        fprintf('Built neighbor edges for %d/%d proteins.\\n', i, n);
    end
end

W = sparse(rowIdx, colIdx, weights, n, n);
W = max(W, W');  % undirected graph
degree = full(sum(W, 2));
degree(degree == 0) = eps;

weightedMean = (degree' * Xdouble) / sum(degree);
Xcentered = Xdouble - weightedMean;

laplacianScore = inf(size(Xdouble, 2), 1);
for j = 1:size(Xdouble, 2)
    f = Xcentered(:, j);
    numerator = 0.5 * sum(sum(W .* ((f - f').^2)));
    denominator = sum(degree .* (f.^2));

    if denominator > eps
        laplacianScore(j) = numerator / denominator;
    end
end

clear Xdouble Xcentered W neighborIdx rowIdx colIdx weights;

%% STEP 5: Feature selection
numSelectedFeatures = 92;
[sortedScores, sortedIdx] = sort(laplacianScore, 'ascend');
validIdx = sortedIdx(isfinite(sortedScores));

if numel(validIdx) < numSelectedFeatures
    error('Fewer than %d features have finite Laplacian Scores.', ...
        numSelectedFeatures);
end

selectedIdx = validIdx(1:numSelectedFeatures);
Xselected = X(:, selectedIdx);
selectedKmers = allKmers(selectedIdx);

FeatureResults = table(selectedKmers, laplacianScore(selectedIdx), ...
    'VariableNames', {'Kmer', 'LaplacianScore'});
writetable(FeatureResults, 'SelectedFeatures.xlsx');

fprintf('Selected %d features with the lowest Laplacian Scores.\\n', ...
    numSelectedFeatures);

%% STEP 6: K-means clustering
numClusters = 8;
rng(42);

% kmeans may require Statistics and Machine Learning Toolbox.
[clusterLabels, clusterCenters] = kmeans(double(Xselected), ...
    numClusters, 'Distance', 'sqeuclidean', ...
    'Replicates', 5, 'MaxIter', 300);

proteins.Cluster = clusterLabels;

clusterCounts = accumarray(clusterLabels, 1, [numClusters, 1]);
disp(table((1:numClusters)', clusterCounts, ...
    'VariableNames', {'Cluster', 'ProteinCount'}));

% Save all outputs required by the second script.
save('ProteinClusteringResults.mat', ...
    'proteins', 'Xselected', 'selectedIdx', 'selectedKmers', ...
    'laplacianScore', 'clusterLabels', 'clusterCenters', ...
    'numClusters', '-v7.3');

fprintf('Steps 1-6 complete. Saved ProteinClusteringResults.mat.\\n');
