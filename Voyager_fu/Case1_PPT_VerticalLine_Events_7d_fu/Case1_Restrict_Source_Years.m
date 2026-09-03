function files = Case1_Restrict_Source_Years(files,catalog,contextDays)
%Case1_Restrict_Source_Years Read only annual files needed by event windows.
files = string(files(:));
years = [];
for ii = 1:height(catalog)
    a = catalog.StartUTC(ii)-days(contextDays);
    b = catalog.EndUTCExclusive(ii)+days(contextDays)-seconds(1);
    years = [years,year(a):year(b)]; %#ok<AGROW>
end
keep = false(size(files));
for ii = 1:numel(files)
    token = regexp(char(files(ii)),'_(\d{4})0101_v','tokens','once');
    assert(~isempty(token),'Unexpected annual source filename.');
    keep(ii) = ismember(str2double(token{1}),years);
end
files = files(keep);
end
