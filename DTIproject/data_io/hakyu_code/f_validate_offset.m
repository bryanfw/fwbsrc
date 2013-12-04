
function tmp_attr_m = f_validate_offset(tmp_attr_m,tmp_count,a,idx_first_wrong_offset,count)
%[f_validate_offset]
% Check if offset value is correct.
% For reading a large LIST file, offset values are read wrong. This
% problem has been solved not using 'single' precision attr matrix.
% With 'single' precision, offset values are sometimes read wrong.
%
% 2011.10.31.
%
% Ha-Kyu

siz = a(19);
offset = a(20);
if tmp_count>1
    siz_prev = tmp_attr_m(tmp_count-1,20);
    offset_prev = tmp_attr_m(tmp_count-1,21);
    if offset ~= (offset_prev+siz_prev)
        tmp_attr_m(tmp_count,21) = offset_prev+siz_prev;
        if isempty(idx_first_wrong_offset)
            idx_first_wrong_offset=count;
            fprintf('  WRONG first OFFSET FOUND\n')
            fprintf('  siz[%d],offset[%d],siz_prev[%d],offset_prev[%d],offset corrected[%d]\n',...
                siz,offset,siz_prev,offset_prev,(offset_prev+siz_prev))
        end
    end
else
    error('tmp_count must be > 1')
end