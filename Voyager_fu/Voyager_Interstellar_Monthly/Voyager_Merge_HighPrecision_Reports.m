function report = Voyager_Merge_HighPrecision_Reports(outputFolder)
%Voyager_Merge_HighPrecision_Reports Merge parallel segment manifests.
%   Chunk manifests are combined, sorted by spacecraft and UTC start time,
%   and the PNG files are renamed to give each spacecraft one continuous
%   segment number sequence. No science values are modified.

if nargin < 1 || isempty(outputFolder)
    outputFolder = fullfile(fileparts(mfilename('fullpath')), ...
        'HighestPrecision_Segments');
end
outputFolder = char(outputFolder);
files = dir(fullfile(outputFolder, ...
    'HighestPrecision_Segments_manifest_*.csv'));
if isempty(files)
    error('VoyagerMerge:NoChunkManifests', ...
        'No chunk manifests were found in %s.', outputFolder);
end

parts = cell(numel(files), 1);
for ii = 1:numel(files)
    parts{ii} = readtable(fullfile(files(ii).folder, files(ii).name), ...
        'TextType', 'string');
end
report = vertcat(parts{:});
report = sortrows(report, {'Spacecraft', 'StartUTC'});

for spacecraft = unique(report.Spacecraft).'
    rows = find(report.Spacecraft == spacecraft);
    for jj = 1:numel(rows)
        row = rows(jj);
        source = char(report.FigureFile(row));
        [folder, name, extension] = fileparts(source);
        newName = regexprep(name, 'segment_\d+_', ...
            sprintf('segment_%04d_', jj), 'once');
        destination = fullfile(folder, [newName, extension]);
        if ~strcmpi(source, destination)
            if isfile(destination)
                error('VoyagerMerge:DestinationExists', ...
                    'Rename destination already exists: %s', destination);
            end
            if ~isfile(source)
                error('VoyagerMerge:SourceMissing', ...
                    'Segment PNG is missing: %s', source);
            end
            [ok, message] = movefile(source, destination);
            if ~ok, error('VoyagerMerge:RenameFailed', '%s', message); end
        end
        report.Segment(row) = jj;
        report.FigureFile(row) = string(destination);
    end
end

writetable(report, fullfile(outputFolder, ...
    'HighestPrecision_Segments_manifest.csv'));
save(fullfile(outputFolder, 'HighestPrecision_Segments_report.mat'), ...
    'report', '-v7.3');
fprintf('Merged %d high-precision segments.\n', height(report));
end
