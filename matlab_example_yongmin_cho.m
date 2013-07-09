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
n_trials = 100; %Increase this value for a best accurate result
%Do the test analysis using a Bernoulli Random matrix as null model
w.statistics.DoCompleteAnalysis(n_trials, @NullModels.NULL_1);
%Do the test analsysi using a Degree random matrix as null model
w.statistics.DoCompleteAnalysis(n_trials, @NullModels.NULL_2);

%Plot network as matrices
w.plotter.use_module_format = 0;
w.plotter.use_isocline = 0;
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
