function value = sage_relocate_paths(value)
% Translate installation paths in memory; never rewrite saved scientific
% files.
matlabRoot=fileparts(fileparts(fileparts(mfilename('fullpath'))));
sageRoot=fileparts(matlabRoot);
if isstruct(value)
 for k=1:numel(value)
  fields=fieldnames(value(k));
  for j=1:numel(fields), value(k).(fields{j})=sage_relocate_paths(value( ...
      k).(fields{j})); end
 end
elseif iscell(value)
 for k=1:numel(value), value{k}=sage_relocate_paths(value{k}); end
elseif isstring(value)
 for k=1:numel(value), value(k)=string(sage_relocate_paths(char(value( ...
     k)))); end
elseif ischar(value) && ...
    isrow(value)
 normalized=strrep(value,'\','/');prefix=strrep(sageRoot,'\','/');
 if strcmpi(normalized,prefix), value=matlabRoot;return,end
 if startsWith(normalized,[prefix '/'],'IgnoreCase',true)
  relative=extractAfter(normalized,strlength(prefix)+1);
  parts=split(string(relative),'/');
  if ~ismember(lower(parts(1)),["matlab","python"])
   value=fullfile(matlabRoot,char(relative));
  end
 end
end
end
