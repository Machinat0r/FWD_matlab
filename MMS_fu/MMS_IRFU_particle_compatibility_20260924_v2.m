%% IRFU现有接口兼容修正：不新增科学处理function
% FEEPS未启用探头的CDF值为-2147483648，与其FILLVAL=-1e31不一致。
% 仅将该精确占位值设为NaN，继续使用get_data原有的探头平均。
Source=which('mms.get_data');s=fileread(Source);
if ~contains(s,'inactive FEEPS eye pad')
    Backup=fullfile(RecordDir,'get_data_before_feeps_pad_fix.m');
    if ~isfile(Backup),copyfile(Source,Backup);end
    s=strrep(s,'energies = energies + energy_top_tmp.data;','energies = cat(3,energies,energy_top_tmp.data);');
    s=strrep(s,'energies = energies + energy_bottom_tmp.data;','energies = cat(3,energies,energy_bottom_tmp.data);');
    s=strrep(s,'energies = energies / 2 / nSensors;',['energies(energies==double(intmin(''int32''))) = NaN; % inactive FEEPS eye pad' newline '          energies = mean(energies,3,''omitnan'');']);
    s=strrep(s,'top = mms.db_get_ts(dsetName, [dsetPref ''_top_'' suf], Tint);',['top = mms.db_get_ts(dsetName, [dsetPref ''_top_'' suf], Tint);' newline '          top.data(top.data==double(intmin(''int32''))) = NaN; % inactive FEEPS eye pad']);
    s=strrep(s,'bot = mms.db_get_ts(dsetName,[dsetPref ''_bottom_'' suf],Tint);',['bot = mms.db_get_ts(dsetName,[dsetPref ''_bottom_'' suf],Tint);' newline '          bot.data(bot.data==double(intmin(''int32''))) = NaN; % inactive FEEPS eye pad']);
    fid=fopen(Source,'w','n','UTF-8');fprintf(fid,'%s',s);fclose(fid);
    clear mms.get_data
end
