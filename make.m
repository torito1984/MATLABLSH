% Compiles the LSH function to be used by 
% Matlab for search of nearest neighbors

status = rmdir('./lib', 's');
mkdir ./lib

out = fullfile(pwd, '/lib');

cd ./sources/

% The definition of the type is only needed in mac (in other OS could be ommited)

mex('-outdir', out, '-c', '-I.', '-DREAL_DOUBLE', '-Dchar16_t=uint16_T', 'BucketHashing.cpp', 'Geometry.cpp', 'GlobalVars.cpp')
mex('-outdir', out, '-c', '-I.', '-DREAL_DOUBLE', '-Dchar16_t=uint16_T', 'LocalitySensitiveHashing.cpp', 'NearNeighbors.cpp')
mex('-outdir', out, '-c', '-I.', '-DREAL_DOUBLE', '-Dchar16_t=uint16_T', 'Random.cpp', 'SelfTuning.cpp', 'Util.cpp')

mex('-outdir', out,  '-Dchar16_t=uint16_T', '-I.',  '-DREAL_DOUBLE', 'lshfind.cpp', '../lib/BucketHashing.o', ...
'../lib/Geometry.o', '../lib/GlobalVars.o', '../lib/LocalitySensitiveHashing.o', '../lib/NearNeighbors.o', '../lib/Random.o',...
'../lib/SelfTuning.o', '../lib/Util.o'); 


cd ..