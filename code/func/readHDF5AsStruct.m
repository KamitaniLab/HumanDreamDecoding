function [ data, info ] = readHDF5AsStruct(hdf5Path, rowMajor)
%
% Read an HDF5 file and convert it to a Matlab struct.
%
% ex.) [ data, info ] = readHDF5AsStruct('sampleData.hdf5')
%
% Input:
%   hdf5Path   - path to the HDF5 file.
%   rowMajor   - if 1 (default), data is assumed to be saved in row-major
%                 order ( C/C++ or python style), so dimension of matrix is permutated. 
%                 if 0, dimension of matrix does'nt change.
%
% Output:
%   data       - < 1x1 struct > data converted from the HDF5 file
%   info       - < 1x1 struct > Information about the HDF5 file
%                ( output of hdf5info(hdf5Path) ).
%
%%
if ~exist( 'hdf5Path', 'var') || isempty(hdf5Path)
    help(mfilename)
    error('Specify a path to an HDF5 file.')
end
if ~exist('rowMajor', 'var')|| isempty(rowMajor)
    rowMajor = 1;
end

newVersion = 0 ; % 1 if MATLAB version >= ver. 2011a
matlabInfo = ver('MATLAB');
[ verIdx1, remain ] = strtok( matlabInfo.Version, '.');
if str2num(verIdx1) >= 8
    newVersion = 1;
else
    verIdx2 = strtok( remain, '.') ;
    if str2num(verIdx1) == 7 && str2num(verIdx2) >= 12
        newVersion = 1;
    end
end


% get info
info   = hdf5info(hdf5Path);
groupHierarchy = info.GroupHierarchy ;

% load data
data = [];
data = readDataSet(data, groupHierarchy, newVersion, rowMajor );

end


%% readDataSet
%  read DataSet and Groups in a hierarchical way
function data = readDataSet(data, groupHierarchy, newVersion, rowMajor)

if ~isempty(groupHierarchy.Groups)
    nGroup = length(groupHierarchy.Groups);
    for groupIdx = 1:nGroup
        data = readDataSet(data, groupHierarchy.Groups(groupIdx), newVersion, rowMajor);
    end
end

% load dataset
if ~isempty(groupHierarchy.Datasets)
    nDatasets = length(groupHierarchy.Datasets);
    for dataIdx = 1:nDatasets
        dataSetName = groupHierarchy.Datasets(dataIdx).Name;
        % replace '.' to '_' if exist
        dataSetName( ismember(dataSetName, '.') ) = '_';
        % replace '/' to '.'
        dataSetName( ismember(dataSetName, '/') ) = '.';
        
        dataSet = hdf5read(groupHierarchy.Datasets(dataIdx)) ;

        % check if the type of the dataSet is a cell or string array or not.
        if strcmp( groupHierarchy.Datasets(dataIdx).Datatype.Class, 'H5T_ARRAY') || strcmp(groupHierarchy.Datasets(dataIdx).Datatype.Class, 'H5T_STRING')
            
            nElements = numel(dataSet);
            cellDataTemp = cell(1,nElements);
            dataSetTemp  = reshape(dataSet, 1, nElements);
            for idx = 1:nElements
                if rowMajor
                    if ~isscalar(dataSetTemp(idx).Data) && ~ischar(dataSetTemp(idx).Data)
                        nDim     = sum(size(dataSetTemp(idx).Data) >= 1);
                        cellDataTemp{idx} = permute(dataSetTemp(idx).Data, [ nDim:-1:1 ] );
                    else
                        cellDataTemp{idx} = dataSetTemp(idx).Data;
                    end
                else
                    cellDataTemp{idx} = dataSetTemp(idx).Data;
                end
            end
            cellData = reshape(cellDataTemp, size(dataSet));
            
            % write
            evalString = [ 'data', dataSetName, ' = ', 'cellData ;' ] ;
            eval(evalString);
            
            if numel(cellData) == 1
                evalString = [ 'data', dataSetName, ' = ', 'cell2mat(data', dataSetName,  ') ;'    ] ;
                eval(evalString);
                
            elseif rowMajor
                evalString = [ 'dataSize = size(data', dataSetName, ') ;' ];
                eval(evalString);
                dataDim = length(dataSize);
                if dataDim > 2 || ( dataSize(1) > 1 && dataSize(2) > 1 );
                    evalString = [ 'data', dataSetName, ' = permute(data', dataSetName, ',[ ', num2str(dataDim:-1:1), '] );'] ;
                    eval(evalString);
                end
            end
                        
        else 
            evalString = [ 'data', dataSetName, ' = ',  'dataSet ;' ];
            eval(evalString);
            
            if rowMajor
                evalString = [ 'dataSize = size(data', dataSetName, ') ;' ];
                eval(evalString);
                dataDim = length(dataSize) ;                
                if dataDim > 2 || ( dataSize(1) > 1 && dataSize(2) > 1 ) || isnumeric(dataSet)
                    evalString = [ 'data', dataSetName, ' = permute(data', dataSetName, ',[ ', num2str(dataDim:-1:1),  '] );' ];
                    eval(evalString);
                end
            end
            
        end
    end
    
end

% load attribute
if ~isempty(groupHierarchy.Attributes)
    nAttributes = length(groupHierarchy.Attributes);
    for attrIdx = 1:nAttributes
        attrLocation = groupHierarchy.Attributes(attrIdx).Location;
        attrName     = groupHierarchy.Attributes(attrIdx).Shortname;
        % replace '/' to '.'
        attrLocation( ismember(attrLocation, '/') ) = '.';
        if length(attrLocation) > 1
            attrLocation = [ attrLocation, '.Attributes.', attrName ];
        else
            attrLocation = [ attrLocation, 'Attributes.', attrName ];
        end
        attrLocation( ismember(attrLocation, ' ') ) = '_';
        
        if isnumeric(groupHierarchy.Attributes(attrIdx).Value)
            attrData = groupHierarchy.Attributes(attrIdx).Value;
        else
            dataSize = length( groupHierarchy.Attributes(attrIdx).Value );
            attrData = cell(dataSize, 1);
            for dataSizeIdx = 1:dataSize
                attrData{dataSizeIdx} = groupHierarchy.Attributes(attrIdx).Value(dataSizeIdx).Data;
            end
        end
        
        evalString = [ 'data', attrLocation, ' = ', 'attrData ;' ];
        eval(evalString);
    end
    
end

end







