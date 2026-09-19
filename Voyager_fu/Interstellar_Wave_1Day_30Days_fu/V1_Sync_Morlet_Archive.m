function audit = V1_Sync_Morlet_Archive(dataRoot,startUTC)
% Verify the official COHO/VIM inventories and the last valid magnetic record.
% Source listings/download metadata are archived beside original CDFs.
base='https://cdaweb.gsfc.nasa.gov/pub/data/voyager/voyager1/';
checked=datetime('now','TimeZone','UTC');
stamp=char(datetime(checked,'Format','yyyyMMdd_HHmmss'));
metadataRoot=fullfile(dataRoot,'source_verification','V1_daily_morlet_latest',stamp);
if ~isfolder(metadataRoot), mkdir(metadataRoot); end
audit=struct('VerifiedOnline',false,'CheckedUTC',checked,'MetadataRoot',metadataRoot, ...
    'Dataset','VOYAGER1_COHO1HR_MERGED_MAG_PLASMA','Listings',table,'Files',table);
url=[base 'coho1hr_magplasma/'];
[html,row]=listing(url,fullfile(metadataRoot,'coho_root.html'));
audit.Listings=[audit.Listings;row];
tokens=regexp(html,'href="(20\d{2})/"','tokens');
years=unique(cellfun(@(x)str2double(x{1}),tokens));
years=years(years>=year(startUTC));
assert(~isempty(years),'Official archive years were not found.');
rows=table;
for yy=years
    yearURL=sprintf('%s%d/',url,yy);
    [html,row]=listing(yearURL,fullfile(metadataRoot,sprintf('coho_%d.html',yy)));
    audit.Listings=[audit.Listings;row]; %#ok<AGROW>
    tokens=regexp(html,'href="(voyager1_coho1hr_merged_mag_plasma_(\d{8})_v(\d+)\.cdf)"','tokens');
    assert(~isempty(tokens),'No recognized CDFs in official year %d.',yy);
    for k=1:numel(tokens)
        z=tokens{k};
        parts=regexp(z{1},'_(\d{8})_v(\d+)\.cdf$','tokens','once');
        date=datetime(parts{1},'InputFormat','yyyyMMdd','TimeZone','UTC');
        if date+calmonths(1)<=startUTC, continue; end
        assert(year(date)==yy,'CDF filename year differs from listing.');
        name=string(z{1}); version=str2double(parts{2});
        folder=fullfile(dataRoot,'voyager1','coho','1hr','l2','merged_mag_plasma', ...
            sprintf('%04d',year(date)),sprintf('%02d',month(date)));
        one=table(date,version,name,string([yearURL char(name)]),string(fullfile(folder,name)), ...
            'VariableNames',{'MonthUTC','Version','Name','URL','SourceFile'});
        rows=[rows;one]; %#ok<AGROW>
    end
end
rows=sortrows(rows,{'MonthUTC','Version'});
[~,last]=unique(rows.MonthUTC,'last');
rows=rows(sort(last),:);
existed=false(height(rows),1); bytes=zeros(height(rows),1);
for k=1:height(rows)
    [existed(k),bytes(k)]=ensureCDF(rows.URL(k),rows.SourceFile(k));
end
rows.ExistedBefore=existed; rows.Bytes=bytes;
audit.Files=rows;
months=(dateshift(startUTC,'start','month'):calmonths(1):max(rows.MonthUTC)).';
audit.VerifiedMissingMonths=string(datestr(setdiff(months,rows.MonthUTC),'yyyymmdd'));
audit.VerifiedMissingMonths=audit.VerifiedMissingMonths(:);

% A filename endpoint can contain fill; inspect the actual scalar payload.
found=false;
for k=height(rows):-1:1
    q=Voyager_Read_CDF_Product(rows.SourceFile(k),'coho');
    if ~isfield(q,'ABS_B') || ~isfield(q,'Epoch'), continue; end
    valid=isfinite(q.ABS_B(:)) & ~isnat(q.Epoch(:));
    if ~any(valid), continue; end
    [lastUTC,j]=max(q.Epoch(valid));
    values=q.ABS_B(valid);
    audit.LastValidUTC=lastUTC;
    audit.LastValidB_nT=values(j);
    audit.EndpointSourceFile=rows.SourceFile(k);
    audit.EndpointSHA256=Case1_File_SHA256(rows.SourceFile(k));
    found=true; break
end
assert(found,'No valid official scalar magnetic record.');
audit.StopExclusiveUTC=dateshift(audit.LastValidUTC,'start','day')+days(1);

% Independently check the latest original reviewed 48-second MAG product.
vimURL=[base 'magnetic_fields_cdaweb/vim_48secmag/'];
[html,row]=listing(vimURL,fullfile(metadataRoot,'vim_root.html'));
audit.Listings=[audit.Listings;row];
tokens=regexp(html,'href="(voyager1_48s_mag-vim_(\d{8})_v(\d+)\.cdf)"','tokens');
assert(~isempty(tokens),'No recognized reviewed VIM products.');
names=sort(unique(string(cellfun(@(x)x{1},tokens,'UniformOutput',false))));
name=names(end); token=regexp(char(name),'_(\d{4})\d{4}_v','tokens','once');
file=fullfile(dataRoot,'voyager1','mag','48s','reviewed_vim',token{1},name);
ensureCDF([vimURL char(name)],file);
q=Voyager_Read_CDF_Product(file,'mag48s');
valid=isfinite(q.F1(:)) & ~isnat(q.Epoch(:));
assert(any(valid),'Latest reviewed VIM file has no valid scalar magnetic data.');
audit.VIMLastValidUTC=max(q.Epoch(valid));
audit.VIMSourceFile=string(file);
audit.VIMSHA256=Case1_File_SHA256(file);
assert(dateshift(audit.VIMLastValidUTC,'start','day')<= ...
    dateshift(audit.LastValidUTC,'start','day'), ...
    'Reviewed 48-second MAG extends beyond COHO; inspect product choice before extending.');
audit.VerifiedOnline=true;
save(fullfile(metadataRoot,'official_archive_inventory.mat'),'audit');
writetable(rows,fullfile(metadataRoot,'official_cdf_inventory.csv'));
fprintf('Last valid hourly |B|: %s UTC\n',string(audit.LastValidUTC));
fprintf('Latest reviewed 48-second |B|: %s UTC\n',string(audit.VIMLastValidUTC));
fprintf('Officially absent monthly CDFs: %s\n',strjoin(audit.VerifiedMissingMonths,', '));
end

function [html,row]=listing(url,file)
[status,message]=system(sprintf(['curl.exe --fail --silent --show-error --location ', ...
    '--retry 2 --connect-timeout 15 --max-time 60 "%s" --output "%s"'],url,file));
assert(status==0,'Official listing failed: %s\n%s',url,message);
html=fileread(file);
row=table(string(url),string(file),Case1_File_SHA256(file), ...
    'VariableNames',{'URL','SavedHTML','SHA256'});
end

function [existed,bytes]=ensureCDF(url,file)
file=char(file); url=char(url); existed=isfile(file);
if ~existed
    folder=fileparts(file); if ~isfolder(folder), mkdir(folder); end
    partial=[file '.download'];
    fprintf('Downloading official CDF: %s\n',url);
    [status,message]=system(sprintf(['curl.exe --fail --silent --show-error --location ', ...
        '--retry 2 --connect-timeout 15 --max-time 180 "%s" --output "%s"'],url,partial));
    assert(status==0,'CDF download failed: %s\n%s',url,message);
    check=dataobj(partial);
    assert(~isempty(check.Variables),'Downloaded file is not a readable CDF.');
    movefile(partial,file);
end
info=dir(file); bytes=info.bytes;
end


