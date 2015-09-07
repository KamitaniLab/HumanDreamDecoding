% This make.m is used under Windows

% add -largeArrayDims on 64-bit machines

if ispc
    mex -largeArrayDims -O -c svm.cpp
    mex -largeArrayDims -O -c svm_model_matlab.c
    mex -largeArrayDims -O svmtrain.c svm.obj svm_model_matlab.obj
    mex -largeArrayDims -O svmpredict.c svm.obj svm_model_matlab.obj
    mex -largeArrayDims -O libsvmread.c
    mex -largeArrayDims -O libsvmwrite.c
else
    mex -largeArrayDims -O -c svm.cpp
    mex -largeArrayDims -O -c svm_model_matlab.c
    mex -largeArrayDims -O svmtrain.c svm.o svm_model_matlab.o
    mex -largeArrayDims -O svmpredict.c svm.o svm_model_matlab.o
    mex -largeArrayDims -O libsvmread.c
    mex -largeArrayDims -O libsvmwrite.c
end
