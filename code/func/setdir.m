function [dpath,id]=setdir(dpath)
% setdir - output dpath and make directory if the directory did not exist.
% function [dpath,id]=setdir(dpath)
% 
% 
% 
% 
% 
%  Created by Tomoyasu Horikawa horikawa-t@atr.jp 2011/01/21
% 
% 
id=exist(dpath,'dir');
if ~id
    mkdir(dpath)
end

