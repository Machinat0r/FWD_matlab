function value = Case1_File_SHA256(fileName)
%Case1_File_SHA256 Hash source bytes for the reproducibility audit.
%   No equivalent SHA-256 helper was found in the local IRFU package.
fid = fopen(fileName, 'rb');
if fid < 0, error('VoyagerCase1:SourceOpen', '%s', fileName); end
cleanup = onCleanup(@() fclose(fid));
digest = java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid)
    bytes = fread(fid, 1024*1024, '*uint8');
    digest.update(typecast(bytes, 'int8'));
end
value = string(lower(reshape(dec2hex(typecast(digest.digest(), ...
    'uint8'), 2).', 1, [])));
end
