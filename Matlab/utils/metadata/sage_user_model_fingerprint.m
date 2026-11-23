function value=sage_user_model_fingerprint(folder)
% Identify numerical source and metadata, excluding platform-specific MEX
% files.
names={'user_model.cpp','user_model.hpp','user_model_prepare.cpp', ...
    'user_model_prepare.hpp','crr_user_model_mex.cpp','delta_n.hpp', ...
        'user_model_info.mat'};
md=java.security.MessageDigest.getInstance('SHA-256');
for k=1:numel(names)
    file=fullfile(folder,names{k});assert(isfile(file), ...
        'SAGE:UserModelIdentity','Missing model definition: %s',file);
    fid=fopen(file,'rb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
    bytes=fread(fid,Inf,'*uint8');clear cleanup
    md.update(unicode2native([names{k} char(0)],'UTF-8'));md.update(bytes);
end
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
