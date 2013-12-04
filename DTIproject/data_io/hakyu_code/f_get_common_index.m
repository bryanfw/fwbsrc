
function index0_v = f_get_common_index(loadopts,attr_m,fnames,fname_in_attr_index)
%[f_get_common_index] finds common index in attr_m based on all loadopts fields.
% USAGE:
%   index0_v = f_get_common_index(loadopts,attr_m,fnames,fname_in_attr_index)
%
%
% Last modified
% 2011.04.16.

for ind=1:length(fname_in_attr_index)
    v = attr_m(:,fname_in_attr_index(ind)+1);
    index_v = find(v==loadopts.(sprintf('%s',fnames{ind})));
    if ind==1
        index0_v=index_v;
    end
    index1_v = intersect(index0_v,index_v);
    index0_v = index1_v;
    clear  index1_v  index_v  v
end

return