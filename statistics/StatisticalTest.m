% StatisticalTest Statistical analysis class for a bipartite object in terms of
% modularity and nestedness
%
% StatisticalTest Properties:
%     bipweb - A bipartite network object in which the analysis will be done
%     nulls - A set of random matrices used as null model
%     tempvals - Results for the NTC value
%     eigvals - Results for the espectral Radius algorithm value
%     nestvals - Results for the NODF algorithm value
%     nestvals_rows - Results for the NODF algorithm in rows
%     nestvals_cols - Results for the NODF algorithm in columns
%     qb_vals - Results for the standard modularity value
%     qr_vals - Results for the ratio of interactions inside modules
%     nest_row_contrib - Row NODF contributions
%     nest_col_contrib - Column Nodf Contributions
%     model - Null Model that will be used for creating the random networks nulls
%     replicates - Number of random networks used for the null model
%     modul_done - The analysis in modularity was performed
%     nest_done - The analysis in NODF was performed
%     temp_done - The analysis in NTC was performed
%     eig_done - The analysis en expectral radius was performed
%     nest_contrib_done - The analysis on nestedness contribution was performed
%     print_output - Flag to print output       
%
% StatisticalTest Methods:
%     StatisticalTest - Main Constructor
%     DoNulls - Create random matrices for the statistical analysis
%     DoCompleteAnalysis - Perform the entire modularity and nestedness analysis
%     Nestedness - Perform the NODF Statistical Analysis
%     Temperature - Perform the NTC Statistical Analysis
%     MaxEigenvalue - Perform the spectral radius Statistical Analysis
%     Modularity - Perform the Modularity Statistical Analysis
%     NestednessContributions - Perform the nestedness contribution Statistical Analysis
%     GET_DEV_MODUL - Perform a Modularity Statistical Analysis
%     GET_DEV_NEST - Perform a NODF Statistical Analysis
%     GET_DEV_TEMP - Perform a NTC Statistical Analysis
%     GET_DEV_EIG - Perform a spectral radius Statistical Analysis
%     GET_NEST_CONTRIBUTIONS - Get NODF nestedness contributions of a
%
% See also:
%     BipartiteModularity, NODF, NestednessBINMATNEST
classdef StatisticalTest < handle

    properties(GetAccess = 'public', SetAccess = 'public')
        bipweb       = {};     % A bipartite network object in which the analysis will be done
        nulls        = {};     % A set of random matrices used as null model
        tempvals      = [];    % Results for the NTC value
        eigvals       = [];    % Results for the espectral Radius algorithm value
        nestvals      = [];    % Results for the NODF algorithm value
        nestvals_rows = [];    % Results for the NODF algorithm in rows
        nestvals_cols = [];    % Results for the NODF algorithm in columns
        qb_vals       = [];    % Results for the standard modularity value
        qr_vals       = [];    % Results for the ratio of interactions inside modules
        nest_row_contrib = []; % Row NODF contributions
        nest_col_contrib = []; % Column Nodf Contributions
        model        = {};     % Null Model that will be used for creating the random networks nulls
        replicates   = Options.REPLICATES;    % Number of random networks used for the null model
        modul_done   = 0;      % The analysis in modularity was performed
        nest_done    = 0;      % The analysis in NODF was performed
        temp_done    = 0;      % The analysis in NTC was performed
        eig_done     = 0;      % The analysis en expectral radius was performed
        nest_contrib_done = 0; % The analysis on nestedness contribution was performed
        print_output = 1;      % Flag to print output
    end
    
    methods

        function obj = StatisticalTest(webbip)
        % StatisticalTest - Main Constructor
        %   obj = StatisticalTest(webbip) Create a StatisticalTest object that makes
        %   reference to the Bipartite object webbip
        
            obj.bipweb = webbip;
        end
        
        function obj = DoNulls(obj,nullmodel,replic)
        % DoNulls - Create random matrices for the statistical analysis
        %   obj = DoNulls(obj) Create Options.REPLICATES random matrices using the
        %   EQUIPROBABLE null model.
        %
        %   obj = DoNulls(obj,nullmodel) Create Options.REPLICATES random matrices using the
        %   indicated Null Model
        %
        %   obj = DoNulls(obj,nullmodel,replic) Create replic random
        %   matrices using the null model indicated in the variable
        %   nullmodel
        
            %obj.modul_done = 0;
            %obj.nest_done = 0;
            
            if(nargin == 1)
                obj.model = Options.DEFAULT_NULL_MODEL;
                obj.replicates = Options.REPLICATES;
            elseif(nargin == 2)
                obj.model = nullmodel;
                obj.replicates = Options.REPLICATES;
            else
                obj.model = nullmodel;
                obj.replicates = replic;
            end
            
            obj.nulls = NullModels.NULL_MODEL(obj.bipweb.matrix,obj.model,obj.replicates);
            
        end
        
        function obj = DoCompleteAnalysis(obj, replic, nullmodel)
        % DoCompleteAnalysis - Perform the entire modularity and nestedness analysis
        %   obj = DoCompleteAnalysis(obj) Perform the entire analysis for
        %   nestedness and modularity using the EQUIPROBABLE null model and Options.REPLICATES random matrices.
        %
        %   obj = DoCompleteAnalysis(obj,replic) Perform the entire analysis for
        %   nestedness and modularity using the EQUIPROBABLE null model and a total of replic random matrices.
        %
        %   obj = DoCompleteAnalysis(obj,replic,nullmodel) Perform the entire analysis for
        %   nestedness and modularity using the the specified Null model
        %   and a total of replic random matrices
        
            
            if(nargin == 1)
                nullmodel = Options.DEFAULT_NULL_MODEL;
                replic = Options.REPLICATES;
            elseif(nargin == 2)
                nullmodel = Options.DEFAULT_NULL_MODEL;
            end
            
            obj.DoNulls(nullmodel,replic);
            obj.Nestedness();
            obj.Modularity();
            obj.Temperature();
            %obj.MaxEigenvalue();
            
            if(obj.print_output == 1)
                fprintf('Null Model: %s\n', func2str(obj.model));
                fprintf('Trials: %i\n', replic);
                fprintf('Modularity\n');
                fprintf('\tQb: %f\n', obj.qb_vals.value);
                fprintf('\tz-score: %f\n', obj.qb_vals.zscore);
                fprintf('\tpercentage: %f\n', obj.qb_vals.percent);
                fprintf('Nestedness\n');
                fprintf('\tNodf: %f\n', obj.nestvals.value);
                fprintf('\tz-score: %f\n', obj.nestvals.zscore);
                fprintf('\tpercentage: %f\n', obj.nestvals.percent);
                fprintf('Temperature\n');
                fprintf('\tNTC: %f\n', obj.tempvals.value);
                fprintf('\tz-score: %f\n', obj.tempvals.zscore);
                fprintf('\tpercentage: %f\n', obj.tempvals.percent);
                %fprintf('Eigenvalue Nestedness\n');
                %fprintf('\tMax Eigenvalue: %f\n', obj.eigvals.maxe);
                %fprintf('\tz-score: %f\n', obj.eigvals.zscore);
                %fprintf('\tpercentage: %f\n', obj.eigvals.percent);
            end
            
        end
        
        function obj = Nestedness(obj)
        % Nestedness - Perform the NODF Statistical Analysis
        %   obj = Nestedness(obj) Perform the entire NODF analsysis. Be
        %   sure to create the random matrices before calling this
        %   function. Otherwise only Options.REPLICATES equiprobable random matrices will
        %   be used for the analysis
        
            if(isempty(obj.nulls))
                obj.DoNulls();
            end
            [obj.nestvals obj.nestvals_rows obj.nestvals_cols] = StatisticalTest.GET_DEV_NEST(obj.bipweb,obj.nulls);
            
            obj.nestvals.replicates = obj.replicates;
            obj.nestvals.model = obj.model;
            
            obj.nest_done = 1;
            
        end
        
        function obj = Temperature(obj)
        % Temperature - Perform the NTC Statistical Analysis
        %   obj = Temperature(obj) Perform the NTC Statistical analsysis. Be
        %   sure to create the random matrices before calling this
        %   function. Otherwise only Options.REPLICATES equiprobable random matrices will
        %   be used for the analysis
        
            if(isempty(obj.nulls))
                obj.DoNulls();
            end
            [obj.tempvals] = StatisticalTest.GET_DEV_TEMP(obj.bipweb,obj.nulls);
            
            obj.tempvals.replicates = obj.replicates;
            obj.tempvals.model = obj.model;
            
            obj.temp_done = 1;
        end
        
        function obj = MaxEigenvalue(obj)
        % MaxEigenvalue - Perform the spectral radius Statistical Analysis
        %   obj = MaxEigenvalue(obj) Perform the NTC Statistical analsysis. Be
        %   sure to create the random matrices before calling this
        %   function. Otherwise only Options.REPLICATES equiprobable random matrices will
        %   be used for the analysis
            
            if(isempty(obj.nulls))
                obj.DoNulls();
            end
            
            [obj.eigvals] = StatisticalTest.GET_DEV_EIG(obj.bipweb,obj.nulls);
            
            obj.eig_done = 1;
        end
        
        function obj = Modularity(obj)
        % Modularity - Perform the Modularity Statistical Analysis
        %   obj = Temperature(obj) Perform the NTC Statistical analsysis. Be
        %   sure to create the random matrices before calling this
        %   function. Otherwise only Options.REPLICATES equiprobable random matrices will
        %   be used for the analysis
        
            if(isempty(obj.nulls))
                obj.DoNulls();
            end
            
            %Calculate the modularity of the Bipartite object
            if(obj.bipweb.modules.done == 0)
                obj.bipweb.modules.Detect(100);
            end
            
            [obj.qb_vals obj.qr_vals] = StatisticalTest.GET_DEV_MODUL(obj.bipweb, obj.nulls);
            obj.qb_vals.algorithm = class(obj.bipweb.modules);
            obj.qb_vals.replicates = obj.replicates;
            obj.qb_vals.model = obj.model;
            obj.modul_done = 1;
        end
        
        function obj = NestednessContributions(obj)
        % NestednessContributions - Perform the nestedness contribution Statistical Analysis
        %   obj = NestednessContributions(obj) Perform the nestedness contribution Statistical analsysis. Be
        %   sure to create the random matrices before calling this
        %   function. Otherwise only Options.REPLICATES equiprobable random matrices will
        %   be used for the analysis
        
            if(isempty(obj.nulls))
                obj.DoNulls();
            end
            [obj.nest_row_contrib obj.nest_col_contrib] = StatisticalTest.GET_NEST_CONTRIBUTIONS(obj.bipweb.matrix,obj.nulls);
            
            obj.nest_row_contrib.replicates = obj.replicates;
            obj.nest_row_contrib.model = obj.model;
            
            obj.nest_col_contrib.replicates = obj.replicates;
            obj.nest_col_contrib.model = obj.model;
            
            obj.nest_contrib_done = 1;
        end
        
    end
    
    methods(Static)
        
        function [out_b out_r] = GET_DEV_MODUL(webbip,rmatrices)
        % GET_DEV_MODUL - Perform a Modularity Statistical Analysis
        %   [out_b out_r] =  GET_DEV_MODUL(webbip,rmatrices) Perform t-test and
        %   z-test in the modularity value of the bipartite object webbip using
        %   the ser of random matrices rmatrices. Return an
        %   structure for both Qb (standard modularity definition) and Qr
        %   (ratio of inside module vs total interactions). These
        %   structures contain the next values:
        %      value   - The value in the empirical matrix
        %      p       - p-value of the performed t-test
        %      ci      - Confidence interval of the performet t-test
        %      percent - The percent of random networks for which  the
        %                empirical value is bigger than the value of the random
        %                networks
        %      z       - z-score of the empirical value
        
            wQr = webbip.modules.Qr;
            wQb = webbip.modules.Qb;
            n = length(rmatrices);
            
            Qb_random = zeros(n,1);
            Qr_random = zeros(n,1);
            
            modul_class = str2func(class(webbip.modules));
            n_trials = webbip.modules.trials;
            parfor i = 1:n
                modularity = feval(modul_class,rmatrices{i});
                modularity.trials = n_trials;
                modularity.Detect(10);
                Qb_random(i) = modularity.Qb;
                Qr_random(i) = modularity.Qr; 
                %fprintf('Trial %i:\n', i);
            end
            
            out_b = StatisticalTest.GET_STATISTICAL_VALUES(wQb,Qb_random);
            out_r = StatisticalTest.GET_STATISTICAL_VALUES(wQr,Qr_random);
            
        end
        
        function [out out_row out_col] = GET_DEV_NEST(webbip,rmatrices) 
        % GET_DEV_NEST - Perform a NODF Statistical Analysis
        %   [out out_row out_col] =  GET_DEV_NEST(webbip,rmatrices) Perform t-test and
        %   z-test in the NODF value of the bipartite object webbip using
        %   the ser of random matrices rmatrices. Return a
        %   structure for nodf values in the entire matrix(out), rows (out_row)
        %   and columns (out_col) with the next elements:
        %      value   - The value in the empirical matrix
        %      p       - p-value of the performed t-test
        %      ci      - Confidence interval of the performet t-test
        %      percent - The percent of random networks for which  the
        %                empirical value is bigger than the value of the random
        %                networks
        %      z       - z-score of the empirical value    
            n = length(rmatrices);
            expect = zeros(n,1);
            expect_row = zeros(n,1);
            expect_col = zeros(n,1);
            
            parfor i = 1:n
                Nodf = NODF(rmatrices{i});
                expect(i) = Nodf.nodf;
                expect_row(i) = Nodf.nodf_rows;
                expect_col(i) = Nodf.nodf_cols;
            end
            
            out = StatisticalTest.GET_STATISTICAL_VALUES(webbip.nodf.nodf,expect);
            out_row = StatisticalTest.GET_STATISTICAL_VALUES(webbip.nodf.nodf_rows,expect_row);
            out_col = StatisticalTest.GET_STATISTICAL_VALUES(webbip.nodf.nodf_cols,expect_col);
            
        end
        
        function [out] = GET_DEV_TEMP(webbip,rmatrices)
        % GET_DEV_TEMP - Perform a NTC Statistical Analysis
        %   [out out_row out_col] =  GET_DEV_TEMP(webbip,rmatrices) Perform t-test and
        %   z-test in the NTC value of the bipartite object webbip using
        %   the ser of random matrices rmatrices. Return a
        %   structure for NTC statistical values in the entire matrix(out)
        %   with the next elements:
        %      value   - The value in the empirical matrix
        %      p       - p-value of the performed t-test
        %      ci      - Confidence interval of the performet t-test
        %      percent - The percent of random networks for which  the
        %                empirical value is bigger than the value of the random
        %                networks
        %      z       - z-score of the empirical value                 
            n = length(rmatrices);
            nestedness = NestednessBINMATNEST(webbip.webmatrix>0);
            nestedness.CalculateNestedness();
            ntc = nestedness.N;
            
            expect = zeros(n,1);
            
            for i = 1:n
                nestedness.SetMatrix(rmatrices{i});
                nestedness.DoGeometry = 0;
                nestedness.CalculateNestedness();
                expect(i) = nestedness.N;
                %fprintf('Trial %i:\n', i);
            end
            
            out = StatisticalTest.GET_STATISTICAL_VALUES(ntc,expect);
        end
        
        function [out] = GET_DEV_EIG(webbip,rmatrices)
        % GET_DEV_EIG - Perform a spectral radius Statistical Analysis
        %   [out out_row out_col] =  GET_DEV_EIG(webbip,rmatrices) Perform t-test and
        %   z-test in the NTC value of the bipartite object webbip using
        %   the ser of random matrices rmatrices. Return a
        %   structure for NTC statistical values in the entire matrix(out)
        %   with the next elements:
        %      value   - The value in the empirical matrix
        %      p       - p-value of the performed t-test
        %      ci      - Confidence interval of the performet t-test
        %      percent - The percent of random networks for which  the
        %                empirical value is bigger than the value of the random
        %                networks
        %      z       - z-score of the empirical value      
            n = length(rmatrices);
            expect = zeros(n,1);
            
            maxe = MatrixNull.GetBiggestEigenvalue(webbip.matrix);
            
            for i = 1:n
                expect(i) = MatrixNull.GetBiggestEigenvalue(rmatrices{i});
            end
            
            out = StatisticalTest.GET_STATISTICAL_VALUES(maxe,expect);
            
        end
        
        function [c_rows c_cols] = GET_NEST_CONTRIBUTIONS(matrix, rmatrices)
        % GET_NEST_CONTRIBUTIONS - Get NODF nestedness contributions of a
        % matrix
        %   [c_rows c_cols] = GET_NEST_CONTRIBUTIONS(matrix, rmatrices) Get
        %   the NODF contributions for both rows and columns. The
        %   contribution is calculated as the z-score of each individually
        %   randomly permuted row and columns
            nodf = NODF(matrix);
            N = nodf.nodf;
            [n_rows n_cols] = size(matrix);
            nn = length(rmatrices);
            row_contrib = zeros(nn,n_rows);
            col_contrib = zeros(nn,n_cols);
            %orig_matrix = matrix;
            for i = 1:n_rows
                temp_row = matrix(i,:);
                for j = 1:nn
                    matrix(i,:) = rmatrices{j}(i,:);
                    nodf = NODF(matrix,0);
                    row_contrib(j,i) = nodf.nodf;
                    if(isnan(nodf.nodf))
                        continue;
                    end
                end
                matrix(i,:) = temp_row;
            end
            
            for i = 1:n_cols
                temp_col = matrix(:,i);
                for j = 1:nn
                    matrix(:,i) = rmatrices{j}(:,i);
                    nodf = NODF(matrix,0);
                    col_contrib(j,i) = nodf.nodf;
                end
                matrix(:,i) = temp_col;
            end
            std_rows = std(row_contrib);
            std_cols = std(col_contrib);
            mean_rows = mean(row_contrib);
            mean_cols = mean(col_contrib);
            
            c_rows = (N - mean_rows)./std_rows;
            c_cols = (N - mean_cols)./std_cols;
            
        end
        
        function out = GET_STATISTICAL_VALUES(real_value, random_values)
            [~, p ci] = ttest(random_values, real_value);
            n = length(random_values);
            me = mean(random_values);
            zscore = (real_value - mean(random_values))/std(random_values);
            stad = std(random_values);
            percent = 100.0*sum(real_value>random_values)/n;
            
            out.value = real_value; out.mean = me; out.std = stad;
            out.p = p; out.ci = ci; out.zscore = zscore; out.percent = percent;
            out.random_values = sort(random_values);
            
        end
        
    end

end
