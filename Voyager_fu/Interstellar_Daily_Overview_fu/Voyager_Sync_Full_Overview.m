function audit=Voyager_Sync_Full_Overview(dataRoot,spacecraft)
% Inventory official NASA directories and download missing original CDFs.
base=sprintf('https://cdaweb.gsfc.nasa.gov/pub/data/voyager/voyager%d/',spacecraft);
root=fullfile(dataRoot,sprintf('voyager%d',spacecraft));
audit=struct('VerifiedOnline',false,'CheckedUTC',datetime('now','TimeZone','UTC'), ...
    'Listings',struct('URL',{},'HTML',{}),'Files',table);
cohoURL=[base,'coho1hr_magplasma/'];
html=fetchText(cohoURL);
audit.Listings(end+1)=struct('URL',cohoURL,'HTML',html);
t=regexp(html,'href="((?:19|20)\d{2})/"','tokens');
years=unique(cellfun(@(x)str2double(x{1}),t)); years=years(years>=1977);
assert(~isempty(years),'Official COHO directory has no recognized years.');
for yy=years
    url=sprintf('%s%d/',cohoURL,yy); html=fetchText(url);
    audit.Listings(end+1)=struct('URL',url,'HTML',html);
    names=cdfNames(html,sprintf('voyager%d_coho1hr_merged_mag_plasma_',spacecraft));
    for k=1:numel(names)
        name=char(names(k)); token=regexp(name,'_(\d{4})(\d{2})\d{2}_v','tokens','once');
        folder=fullfile(root,'coho','1hr','l2','merged_mag_plasma',token{1},token{2});
        audit.Files=[audit.Files;ensureFile([url,name],folder,name)]; %#ok<AGROW>
    end
end
products={'lev-1-rates','lev-2-daily-avg'};
for pp=1:2
    url=[base,'particle/lecp/final-cdf/',products{pp},'/']; html=fetchText(url);
    audit.Listings(end+1)=struct('URL',url,'HTML',html);
    names=cdfNames(html,[sprintf('voyager-%d_lecp_',spacecraft),products{pp},'_']);
    for k=1:numel(names)
        name=char(names(k)); token=regexp(name,'_(\d{4})0101_v','tokens','once');
        if str2double(token{1})<1977, continue, end
        if pp==1, folder=fullfile(root,'lecp','native','l1','sectored_rates',token{1});
        else, folder=fullfile(root,'lecp','1d','l2','sectored_flux',token{1}); end
        audit.Files=[audit.Files;ensureFile([url,name],folder,name)]; %#ok<AGROW>
    end
end
audit.VerifiedOnline=true;
out = fullfile(dataRoot,'source_verification',sprintf('V%d_full_overview',spacecraft));
if ~isfolder(out), mkdir(out); end
save(fullfile(out,'official_archive_inventory.mat'),'audit');
end

function names=cdfNames(html,prefix)
t=regexp(html,['href="(',regexptranslate('escape',prefix),'[^"/]+\.cdf)"'],'tokens');
names=sort(unique(string(cellfun(@(x)x{1},t,'UniformOutput',false))));
assert(~isempty(names),'No recognized original CDF files in official listing.');
end

function html=fetchText(url)
[status,html]=system(sprintf('curl.exe -6 --noproxy "*" --fail --silent --show-error --retry 3 --retry-all-errors --connect-timeout 15 --max-time 45 "%s"',url));
assert(status==0,'Official archive listing failed: %s\n%s',url,html);
end

function row=ensureFile(url,folder,name)
if ~isfolder(folder), mkdir(folder); end
file=fullfile(folder,name); existed=isfile(file);
if ~existed
    fprintf('Downloading original CDF: %s\n',url);
    partial=[file,'.download'];
    [status,message]=system(sprintf(['curl.exe -6 --noproxy "*" --fail --silent --show-error ', ...
        '--retry 3 --connect-timeout 15 --max-time 180 --output "%s" "%s"'],partial,url));
    assert(status==0,'Download failed; partial retained: %s\n%s',partial,message);
    obj=dataobj(partial); assert(~isempty(obj.Variables),'Downloaded file is not a readable CDF.');
    movefile(partial,file);
end
info=dir(file);
row=table(string(url),string(file),info.bytes,existed, ...
    'VariableNames',{'URL','SourceFile','Bytes','ExistedBefore'});
end



