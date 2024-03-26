function Identification(Path)
FilenamePath = [Path '.m'];
fid = fopen(FilenamePath);
chars = char(fread(fid)');
expr = '-+.+Wending Fu.+-+';
loc = regexpi(chars,expr, 'once');
if isempty(loc)
    error('DO NOT CHANGE THE SIGNATURE!!! (*￣︿￣)')
end
end