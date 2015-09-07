function writeHDF5FromStruct(dataPath, mStruct, rowMajor )
% 
% function writeHDF5FromStruct(dataPath, mStruct, compress)
% Input:
%   dataPath      - a path to the output HDF file.
%                   If file exists, data is overwritten.
%   mStruct       - < 1x1 struct > 
%
%   rowMajor      - If 1 (default), change order of storing
%                   multidimensional arrays.
%                   This is for compatibility with data format stored using 
%                   row-major based languages such as python, C/C++, ...
%
%
%   Note: If this script works on Matlab version >= 7.12 (2011a) and data size > 10MB, 
%         numerical data will be compressed. 
%
% ex.) 
%      data.data1=1;data.group1.data2={'a','b'};
%      writeHDF5FromStruct( 'sampleHDF5.h5', data );
%  
%      datasets  '/data1' and '/group1/data2' are created
%      and saved to 'sampleHDF5.h5'
%    
%
if ~exist('dataPath', 'var')||isempty(dataPath)
    error('Specify output HDF file');    
end
if ~exist('mStruct', 'var')
    error('Specify an input Matlab struct');
end
if ~exist('rowMajor', 'var') || isempty(rowMajor)
    rowMajor = 1;
end

% If this script works on Matlab version >= 7.12 (2011a),
% large data is compressed.
compress = 0;
matlabInfo = ver('MATLAB');
[ verIdx1, remain ] = strtok( matlabInfo.Version, '.');
if str2num(verIdx1) >= 8
    compress = 1;
else
    verIdx2 = strtok( remain, '.') ;
    if str2num(verIdx1) == 7 && str2num(verIdx2) >= 12
        compress = 1;
    end
end

if exist(dataPath, 'file')
    delete(dataPath);
end

groupHierachy = [];
if ~isstruct(mStruct)
    evalString = [ inputname(2) ' = mStruct ;' ];
    eval(evalString);
    eval([ 'writeHDF5Hierachy(dataPath,' inputname(2) ', groupHierachy, compress, rowMajor) ;' ] );
else
    writeHDF5Hierachy(dataPath, mStruct, groupHierachy, compress, rowMajor);
end

end
%%
function writeHDF5Hierachy(dataPath, mStruct, groupHierachy, compress, rowMajor)

if ~isstruct(mStruct)
    data = mStruct ;
    clear mStruct
    mStruct.(inputname(2)) = data;
end
fields = fieldnames(mStruct);
nField = length(fields);
for fieldIdx = 1:nField
    data = mStruct.(fields{fieldIdx});
    gh   = [ groupHierachy, '/', fields{fieldIdx} ];
    if isstruct(data)        

        createGroup(dataPath, gh);
        writeHDF5Hierachy(dataPath, data, gh, compress, rowMajor);                
    else
        % check if data has empty element []
        if ~isempty(data)
            if rowMajor
                dataSize = size(data);
                nDim     = sum(dataSize >= 1);
                if ~ischar(data)
                    data = permute(data, [ nDim:-1:1 ]);
                end
            end           
            % check if data is cell array of cells
            if iscell(data)
                if ~isscalar(data{1}) && ~ischar(data{1})
                    if rowMajor
                        dataSize  = size(data);
                        nElements = numel(data);
                        
                        dataTemp = reshape(data, 1, nElements);
                        for idx = 1:nElements
                            nDim = sum(dataSize >= 1) ;
                            dataTemp{idx} = permute(dataTemp{idx}, [ nDim:-1:1 ]) ;
                        end
                        data     = reshape(dataTemp, dataSize);
                        
                    end
                end
            end
        end
        
        dataInfo = whos('data');
        % compress if data size > 10000000 bytes (10MB);
        if dataInfo.bytes > 10000000 && isnumeric(data) && compress 
            chunkSize = floor(size(data)/10) + 1 ;
            h5create(dataPath, gh, size(data), 'Deflate',9,'ChunkSize',chunkSize);
            h5write(dataPath, gh, data);
        elseif numel(data) == 1 && isnumeric(data)
            writeScaler(dataPath, gh, data) ;
        elseif isnumeric(data)
            writeNumericData(dataPath, gh, data, rowMajor);
        elseif ischar(data)
            %if data is 1D string
            write1DString(dataPath, gh, data) ;
        else
            if exist(dataPath, 'file')
                hdf5write(dataPath, gh , data, 'WriteMode', 'append');
            else
                hdf5write(dataPath, gh , data);
            end
        end
    
    end
        
end

end

%% sub functions
function writeScaler(dataPath, gh, data)
opt1  = 'H5P_DEFAULT';

type_id  = H5T.copy('H5T_NATIVE_DOUBLE');
space_id = H5S.create('H5S_SCALAR');

if exist(dataPath, 'file')
    fid     = H5F.open(dataPath, 'H5F_ACC_RDWR', opt1);
else
    fcpl = H5P.create('H5P_FILE_CREATE');
    fapl = H5P.create('H5P_FILE_ACCESS');
    fid = H5F.create( dataPath,'H5F_ACC_TRUNC',fcpl,fapl);
end
dset_id = H5D.create(fid, gh, type_id, space_id, opt1);
H5D.write(dset_id,'H5ML_DEFAULT','H5S_ALL','H5S_ALL', opt1, data);

H5S.close(space_id);
H5T.close(type_id);
H5D.close(dset_id);
H5F.close(fid);

end

function writeNumericData(dataPath, gh, data, rowMajor)
opt1  = 'H5P_DEFAULT';
dataSize = size(data);
if rowMajor
    dataSize = fliplr(dataSize);
end
type_id  = H5T.copy('H5T_NATIVE_DOUBLE');
%space_id = H5S.create('H5S_SCALAR');
space_id = H5S.create_simple( length(dataSize), dataSize, dataSize);

if exist(dataPath, 'file')
    fid     = H5F.open(dataPath, 'H5F_ACC_RDWR', opt1);
else
    fcpl = H5P.create('H5P_FILE_CREATE');
    fapl = H5P.create('H5P_FILE_ACCESS');
    fid = H5F.create( dataPath,'H5F_ACC_TRUNC',fcpl,fapl);
end
dset_id = H5D.create(fid, gh, type_id, space_id, opt1);
H5D.write(dset_id,'H5ML_DEFAULT','H5S_ALL','H5S_ALL', opt1, data);

H5S.close(space_id);
H5T.close(type_id);
H5D.close(dset_id);
H5F.close(fid);

end


function write1DString(dataPath, gh, data)

stringLength = length(data);
opt1    = 'H5P_DEFAULT';

type_id  = H5T.copy('H5T_C_S1');
H5T.set_size(type_id, stringLength);

space_id = H5S.create('H5S_SCALAR');

if exist(dataPath, 'file')
    fid     = H5F.open(dataPath,'H5F_ACC_RDWR', opt1);
else
    fcpl = H5P.create('H5P_FILE_CREATE');
    fapl = H5P.create('H5P_FILE_ACCESS');
    fid = H5F.create( dataPath,'H5F_ACC_TRUNC',fcpl,fapl);
end
dset_id = H5D.create(fid, gh, type_id, space_id, opt1);
H5D.write(dset_id,'H5ML_DEFAULT','H5S_ALL','H5S_ALL', opt1, data);

H5S.close(space_id);
H5T.close(type_id);
H5D.close(dset_id);
H5F.close(fid);

end

function createGroup( dataPath, gh )
opt1 = 'H5P_DEFAULT';
if exist(dataPath, 'file')
    fid     = H5F.open(dataPath,'H5F_ACC_RDWR',opt1);
else
    fcpl = H5P.create('H5P_FILE_CREATE');
    fapl = H5P.create('H5P_FILE_ACCESS');
    fid = H5F.create( dataPath,'H5F_ACC_TRUNC',fcpl,fapl);
end

group_id  = H5G.create(fid, gh, opt1,opt1,opt1);
H5G.close(group_id);

end















