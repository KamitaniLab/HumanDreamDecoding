function [vec,r,c]=getTri(mat,lu)
% getTri -- extract upper triangle componet with diagonal from matrix
%
% [Input]
%   -mat: square matrix
%   -lu: 'lower' or 'upper', default: get lower components
%
% [Output]
%   -vec:vector
%   -r:row index
%   -c:column index
%
%
%
% [note]
%   each component was extracted in the ? direction, not -> direction
%   if mat is symmetry, use squareform(vec) to recover original matrix
%   (except diag. comp.) from vec
% 
% 
% [e.g.,]
% mat=[1,2,3,4;5,6,7,8;9,10,11,12;13,14,15,16]
% [v,r,c]=getTri(mat)
% squareform(v)
%
% Created By Tomoyasu Horikawa horikawa-t@atr.jp 2009/10/27
%
r=zeros(size(mat,2)*(size(mat,2)-1)/2,1);
c=zeros(size(mat,2)*(size(mat,2)-1)/2,1);
if iscell(mat(1))
vec=cell(size(mat,2)*(size(mat,2)-1)/2,1);
counter=0;
if ~exist('lu','var')
    lu='l';
end
switch lu(1)
    case 'l'
        for i=1:(size(mat,2)-1)
            for j=(i+1):size(mat,2)
                counter=counter+1;
                vec{counter}=mat{j,i};
                r(counter)=j;
                c(counter)=i;
            end
        end
    case 'u'
        for i=1:(size(mat,2)-1)
            for j=(i+1):size(mat,2)
                counter=counter+1;
                vec{counter}=mat{i,j};
                r(counter)=j;
                c(counter)=i;
            end
        end
    otherwise
        error('Invalid 2nd inputs.')
end
else
vec=zeros(size(mat,2)*(size(mat,2)-1)/2,1);
counter=0;
if ~exist('lu','var')
    lu='l';
end
switch lu(1)
    case 'l'
        for i=1:(size(mat,2)-1)
            for j=(i+1):size(mat,2)
                counter=counter+1;
                vec(counter)=mat(j,i);
                r(counter)=j;
                c(counter)=i;
            end
        end
    case 'u'
        for i=1:(size(mat,2)-1)
            for j=(i+1):size(mat,2)
                counter=counter+1;
                vec(counter)=mat(i,j);
                r(counter)=j;
                c(counter)=i;
            end
        end
    otherwise
        error('Invalid 2nd inputs.')
end
end
end
