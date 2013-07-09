% Add the code to the matlab path
g = genpath('.');
addpath(g);

% Create an object of the class BiWeb. By just creating it,
% nestedness will be calculated;
w = Reading.CREATE_FROM_MATRIX_WITH_LABELS('rodents.web'); %You can create your network directly from a text file
%matrix = [1 1 0 1 0 1 1;1 0 1 0 1 0 0;1 0 0 1 0 1 1;1 1 1 0 0 1 0;1 1 1 1 1 1 1];
%w = Bipartite(matrix); %You can send a matlab matrix directly to the constructor

%Accesing the values of nestedness
nodf = w.nestedness.nodf;
nodf_up = w.nestedness.nodf_rows;
nodf_low = w.nestedness.nodf_cols;

fprintf('\nValues of Nestedness using NODF algorithm:\n');
fprintf('\tNODF:                 \t%f\n', nodf);
fprintf('\tNODF in rows:         \t%f\n', nodf_up);
fprintf('\tNODF in columns:      \t%f\n', nodf_low);

%Accesing the values of modularity using LP&BRIM algorithm
w.modules.Detect();
N = w.modules.N;
Qb = w.modules.Qb;
Qr = w.modules.Qr;

fprintf('\nValues of Modularity using LP&BRIM algorithm:\n');
fprintf('\tNumber of modules:    \t%f\n', N);
fprintf('\tBipartite modularity: \t%f\n', Qb);
fprintf('\tRealized modularity:  \t%f\n', Qr);

%Create statistics for modularity and nestedness;
n_trials = 10; %Increase this value for a best accurate result
w.statistics.DoNulls(@NullModels.NULL_1,n_trials); %Using bernoulli null model
w.statistics.Nestedness();
w.statistics.Modularity();
w.printer.NetworkLevel();
w.statistics.DoNulls(@NullModels.NULL_2_,n_trials); %Using Average in columns and rows
w.statistics.Nestedness();
w.statistics.Modularity();
%Two more null models exist
w.printer.NetworkLevel(); %Add a line with results to networklevel.txt file. The headers in that file mean:
% Basic properties
%-- Name = Name (by default the name of the input data).
%-- Conn = Connectance (Number of interactions/Size).
%-- Size = n_rows * n_cols
%-- n_rows = Number of rows in the matrix
%-- n_cols = Number of columns in the matrix)
% Nestedness Results
%-- nodf = NODF nestedness result in the entire matrix (0 - 100)
%-- nodf_r = NODF nestedness result in the rows of the matrix (0 - 100)
%-- nodf_c = NODF nestedness result in the columns of the matrix (0 100)
%I have never used the next four results, but you can consult on
%http://nicolasmouquet.free.fr/publications/Poisot_et_al_2012_MEE.pdf what
%they mean
%-- aspe = Specificity
%-- arr = Resource range
%-- resp = ? 
%-- inc = ?
% Modularity results.
%-- mN = Number of modules of the modularity algorithm (if the module detections was not used a - appears).
%-- mQb = Standard modularity result
%-- mQr = Fraction of interactions that are inside modules.
% Statistical results - I suggest in focus only in tReps, null_name and all
% the z-scores (n_z, mQb_z).
%-- n_sim, n_pval, nicLow, nicUp, n_z = Statistical measures of nestedness for the entire matrix
%-- nco_sim nco_pval ncoiclow, ncoicup nco_z = Statistical measures of nestedness for rows in the matrix
%-- nro_sim nro_pval nroiclow enroicup nro_z  = Statistical measures of nestedness for columns in the matrix
%-- mQr_sim mQr_val mQrIclow mQrIcup mQr_z = Statistical measures for modularity Qr (mQr)
%-- mQb_sim mQb_val mQbIclow mQbIcup mQb_z = Statistical measures for modularity Qb (mQr)
%-- tReps Number of random networks that was used for calculating the statistical measures
%-- null_name Null model that was used for creating the random networks.


%Plot network as matrices
w.plotter.use_type_interaction = 0; %Change to 1 only if your matrix has values 0,1,2
figure(1);
w.plotter.PlotMatrix();
figure(2);
w.plotter.PlotNestedMatrix();
figure(3);
w.plotter.PlotModularMatrix();

%Plot network as beads (graph)
figure(4);
w.plotter.PlotBMatrix();
figure(5);
w.plotter.PlotBNestedMatrix();
figure(6);
w.plotter.PlotBModularMatrix();
