function [idx] = strfindcell(cellarray,string,idxtype)

if nargin == 2
    idxtype = 'logical';
end

switch idxtype
    case 'find'
        idx = find(not(cellfun('isempty',strfind(cellarray,string))));
    case 'logical'
        idx = not(cellfun('isempty',strfind(cellarray,string)));
end

end          