%% Installation
% To install |BiMat|, copy the downloaded zip file to a directory of interest and unzip it. Next, you
% will need to add |BiMat| to the |Matlab| path either temporally or permanently:
%
% *Temporal path:* Add the |BiMat| directory (and sub-directories) to the
% \matlab path. The user will need to type the next lines everytime he start a new |Matlab| session:
%Replace next line with your appropiate path'
bimat_user_path = 'mypackages/bimat';
g=genpath(bimat_user_path);
addpath(g);	 
%%
% *Permanent path:* Alternatively, the user can add |BiMat| to the |Matlab| search path permanently. 
% Instructions about how to do that can be found on: 
% <http://www.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html Matlab Search Path>
%