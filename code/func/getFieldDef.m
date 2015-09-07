function out = getFieldDef(S, field, default)
%getFieldDef - returns either S.<field> or default
%out = getFieldDef(S, field, default)
%
% If string 'field' is a field of structure S, it returns this value
% otherwise, it returns default in out.
%
% Input:
%   S       - any structure
%   field   - string with possible field of S
%   default - default to assign to 'out', if 'field' doesn't exist
% Output:
%   out     - any variable output
%
% Created  By: Alex Harner (1),     alexh@atr.jp      06/07/05
% Modified By: Alex Harner (1),     alexh@atr.jp      06/07/12
% Modified By: Satoshi MURATA (1),  satoshi-m@atr.jp  08/08/21
% (1) ATR Intl. Computational Neuroscience Labs, Decoding Group


%% Check and get pars:
if exist('default','var')==0 || isempty(default)
    default = [];
end

if exist('S','var')==0 || isempty(S) || exist('field','var')==0 || isempty(field)
    out = default;
    return;
end


%% Return value:
out = default;
if isfield(S,field)
    out = S.(field);
end
