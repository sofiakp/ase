function bigWig2mat(infile, outfile, chrLenFile, chunkLen, append, tempDir, binPath) 
%bigWig2mat Converts bigWig file to mat signal file in the format used by extractSignal.
%bigWig2mat(infile, outfile, chrLenFile, chunkLen, append, tempDir, binPath)
%Input arguments:
%   infile: input bigWig file
%   outfile: output mat file
%   chrLenFile: file with chromosome lengths. First column is chr name, second is length in bp.
%   chunkLen: chunk length for output signal
%   append: append to existing files (T) or overwrite (F, default).
%   tempDir: directory for temporary files (default: directory of output file).   
%   binPath: If bigWigSummary is not in your path, then provide the full
%   path to bigWigSummary here.

[outDir, outBase, ~] = fileparts(outfile);

if nargin < 7
    binPath = 'bigWigSummary';
    if nargin < 6
        tempDir = outDir;
        if nargin < 5
            append = false;
        end
    end
end

tempFile = fullfile(tempDir, strcat(outBase, '_TMP', num2str(randi(10000, 1))));

% Read chromosome lengths
f = fopen(chrLenFile, 'r');
chrData = textscan(f, '%s%d', 'HeaderLines', 0);
chrNames = chrData{1};
chrLen = chrData{2};
fclose(f);

if ~exist(outfile, 'file')
    append = false;
end

if ~append
    % Set parameters
    globalParameters.oChunkLen = chunkLen;
    globalParameters.program = 'bigWig2mat';    
    [specificParameters.chrNames, idx] = sort(chrNames');
    specificParameters.chrLen = double(chrLen(idx));
    save(outfile, 'globalParameters', 'specificParameters');
end

% Get existing file contents in case we want to append to file
contents = whos('-file', outfile);

outPref = 'signal';

for c = 1:length(chrLen)
    numChunks = ceil(double(chrLen(c)) / chunkLen); % Make sure this is not an integer division
    numVals = 0;
    
    disp(['Starting ', chrNames{c}, ', chunks ', num2str(numChunks)]);
    for i = 1:numChunks
        if mod(i, 10) == 0
            disp(['Finished chunk ', num2str(i)]);
        end
        outVar = genvarname([outPref, '_', chrNames{c}, '_', num2str(i)]);
        if append && ismember(outVar, {contents.name})
            continue;
        end
        
        % Start and end coordinates in BED format
        chunkStart = (i - 1) * chunkLen;
        chunkEnd = min(i * chunkLen, chrLen(c));
        trueLen = chunkEnd - chunkStart;
        
        % bigWigSummary will output the mean value in each position, which is the same as the
        % value at that position. The output file will have n/a for missing
        % values. Replace these with NaN.
        command = [binPath, ' ', infile, ' ', chrNames{c}, ' ', num2str(chunkStart), ...
            ' ', num2str(chunkEnd), ' ', num2str(trueLen), ' | sed ''s/n\/a/NaN/g'' > ', tempFile];
        try
            system(command);
            signal = single(dlmread(tempFile));
        catch exc
            % bigWigSummary might throw an exception if for example there
            % are no data in the region specified. Replace the values with
            % NaNs.
            signal = single(nan(trueLen, 1));
        end
        signal = reshape(signal, [trueLen, 1]); % Make sure it's a column vector
        assert(length(signal) == trueLen);
        numVals = numVals + length(signal);
        eval([outVar '= signal;']);
        save(outfile, outVar, '-append');
    end
    disp(['Lines read: ', num2str(numVals)]);
end
delete(tempFile);
end
