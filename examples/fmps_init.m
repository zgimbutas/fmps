DataFolderName = ['..' filesep 'data'];
InfoFileName   = 'info.txt';

% folder with utilites (functions, mex-files ...)
UtilFolderName = ['..' filesep 'fmps'];
addpath(UtilFolderName); % add folder to MatLab path

% triangulations
UtilFolderName = [pwd filesep 'aqctri'];
addpath(UtilFolderName); % add folder to MatLab path

% iterative methods: gmres, conjugate gradients, etc.
UtilFolderName = [pwd filesep 'iterate'];
addpath(UtilFolderName); % add folder to MatLab path

