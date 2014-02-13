% StatisticalTest Statistical analysis class for a bipartite object in terms of
% modularity and nestedness
%
% StatisticalTest Properties:
%     bipweb - A bipartite network object in which the analysis will be done
%     nulls - A set of random matrices used as null model
%     tempvals - Results for the NTC value
%     eigvals - Results for the espectral Radius algorithm value
%     nodf_vals - Results for the NODF algorithm value
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
%
% See also:
%     BipartiteModularity, NODF, NestednessNTC
classdef StatisticalTest < handle

    properties(GetAccess = 'public', SetAccess = 'public')
        bipweb       = {};     % A bipartite network object in which the analysis will be done
        nulls        = {};     % A set of random matrices used as null model
        tempvals      = [];    % Results for the NTC value
        eigvals       = [];    % Results for the espectral Radius algorithm value
        nodf_vals      = [];    % Results for the NODF algorithm value
        nestvals_rows = [];    % Results for the NODF algorithm in rows
        nestvals_cols = [];    % Results for the NODF algorithm in columns
        qb_vals       = [];    % Results for the standard modularity value
        qr_vals       = [];    % Results for the ratio of interactions inside modules
        nest_row_contrib = []; % Row NODF contributions
        nest_col_contrib = []; % Column Nodf Contributions
        model        = Options.DEFAULT_NULL_MODEL; % Null Model that will be used for creating the random networks nulls
        replicates   = Options.REPLICATES;    % Number of random networks used for the null model
        modul_done   = 0;      % The analysis in modularity was performed
        nest_done    = 0;      % The analysis in NODF was performed
        temp_done    = 0;      % The analysis in NTC was performed
        eig_done     = 0;      % The analysis en expectral radius was performed
        nest_contrib_done = 0; % The analysis on nestedness contribution was performed
        print_output = 1;      % Flag to print output
        print_status = 1;      % Print status of the statistical test (print the random matrix that is tested)
    end
    
    methods

        function obj = StatisticalTest(webbip)
        % StatisticalTest - Main Constructor
        %   obj = StatisticalTest(webbip) Create a StatisticalTest object that makes
        %   reference to the Bipartite object webbip
        
            obj.bipweb = webbip;
        end
        
        function obj = DoNulls(obj,replic,nullmodel)
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
                obj.model = Options.DEFAULT_NULL_MODEL;
                obj.replicates = replic;
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
            
            fprintf('Creating %i null random matrices...\n', replic);
            obj.DoNulls(replic,nullmodel);
            fprintf('Performing NODF statistical analysis...\n');
            obj.TestNODF();
            fprintf('Performing Modularity statistical analysis...\n');
            obj.TestModularity();
            fprintf('Performing NTC statistical analysis...\n');
            obj.TestNTC();
            %obj.MaxEigenvalue();
            
        end
        
        function obj = TestNODF(obj)
        % TestNODF - Perform the NODF Statistical Analysis
        %
        %   obj = Nestedness(obj) Perform the entire NODF analsysis. Be
        %   sure to create the random matrices before calling this
        %   function. Otherwise only Options.REPLICATES equiprobable random matrices will
        %   be used for the analysis
        
            if(isempty(obj.nulls))
                obj.DoNulls();
            end
            
            nodf = obj.bipweb.nodf.nodf;
            nodf_rows = obj.bipweb.nodf.nodf_rows;
            nodf_cols = obj.bipweb.nodf.nodf_cols;
            n = obj.replicates;
            
            expect = zeros(n,1);
            expect_row = zeros(n,1);
            expect_col = zeros(n,1);
            
            print_stat = floor(linspace(1,n+1,11));
            jp = 1;
            
            for i = 1:n
                if(obj.print_status && print_stat(jp)==i)
                    fprintf('Evaluating NODF in random matrices: %i - %i\n', print_stat(jp),print_stat(jp+1)-1);
                    jp = jp+1;
                end
                Nodf = NODF(obj.nulls{i});
                expect(i) = Nodf.nodf;
                expect_row(i) = Nodf.nodf_rows;
                expect_col(i) = Nodf.nodf_cols;
            end
            
            obj.nodf_vals = StatisticalTest.GET_STATISTICAL_VALUES(nodf,expect);
            obj.nestvals_rows = StatisticalTest.GET_STATISTICAL_VALUES(nodf_rows,expect_row);
            obj.nestvals_cols = StatisticalTest.GET_STATISTICAL_VALUES(nodf_cols,expect_col);
            
            obj.nodf_vals.replicates = obj.replicates;
            obj.nodf_vals.model = obj.model;
            
            obj.nest_done = 1;
            
        end
        
        function obj = TestNTC(obj)
        % TestNTC - Perform the NTC Statistical Analysis
        %
        %   obj = Temperature(obj) Perform the NTC Statistical analsysis. Be
        %   sure to create the random matrices before calling this
        %   function. Otherwise only Options.REPLICATES equiprobable random matrices will
        %   be used for the analysis
        
            if(isempty(obj.nulls))
                obj.DoNulls();
            end
            
            n = length(obj.nulls);
            
            nestedness = obj.bipweb.ntc;
            if(nestedness.done == 0)
                nestedness.Detect();
                ntc = obj.bipweb.ntc.N;
            else
                ntc = obj.bipweb.ntc.N;
            end
            
            expect = zeros(n,1);
            
            print_stat = floor(linspace(1,n+1,11));
            jp = 1;
            
            for i = 1:n
                if(obj.print_status && print_stat(jp)==i)
                    fprintf('Evaluating NTC in random matrices: %i - %i\n', print_stat(jp),print_stat(jp+1)-1);
                    jp = jp+1;
                end
                nestedness.SetMatrix(obj.nulls{i});
                nestedness.do_geometry = 0;
                nestedness.Detect();
                expect(i) = nestedness.N;
                %fprintf('Trial %i:\n', i);
            end
                       
            obj.tempvals = StatisticalTest.GET_STATISTICAL_VALUES(ntc,expect);
            
            obj.tempvals.replicates = obj.replicates;
            obj.tempvals.model = obj.model;
            
            obj.temp_done = 1;
        end
        
        function obj = MaxEigenvalue(obj)
        % MaxEigenvalue - Perform the spectral radius Statistical Analysis
        %
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
        
        function obj = TestModularity(obj)
        % TestModularity - Perform the Modularity Statistical Analysis
        %
        %   obj = TestModularity(obj) Perform the Modularity Statistical analsysis. Be
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
            
            wQr = obj.bipweb.modules.Qr;
            wQb = obj.bipweb.modules.Qb;
            n = obj.replicates;
            
            Qb_random = zeros(n,1);
            Qr_random = zeros(n,1);
            
            if(isprop(obj.bipweb.modules,'DoKernighanLinTunning'))
                exist_kernig = 1;
                do_kernig = obj.bipweb.modules.DoKernighanLinTunning;
            else
                exist_kernig = 0;
                do_kernig = 0;
            end
            
            modul_class = str2func(class(obj.bipweb.modules));
            n_trials = obj.bipweb.modules.trials;
            
            print_stat = floor(linspace(1,n+1,11));
            jp = 1;
            
            for i = 1:n
                if(obj.print_status && print_stat(jp)==i)
                    fprintf('Evaluating Modularity in random matrices: %i - %i\n', print_stat(jp),print_stat(jp+1)-1);
                    jp = jp+1;
                end
                modularity = feval(modul_class,obj.nulls{i});
                modularity.trials = n_trials;
                if(exist_kernig == 1)
                    modularity.DoKernighanLinTunning = do_kernig;
                end
                modularity.Detect(10);
                Qb_random(i) = modularity.Qb;
                Qr_random(i) = modularity.Qr; 
            end
            
            obj.qb_vals = StatisticalTest.GET_STATISTICAL_VALUES(wQb,Qb_random);
            obj.qr_vals = StatisticalTest.GET_STATISTICAL_VALUES(wQr,Qr_random);
            
            obj.qb_vals.algorithm = class(obj.bipweb.modules);
            obj.qb_vals.replicates = obj.replicates;
            obj.qb_vals.model = obj.model;
            obj.modul_done = 1;
        end
        
        function obj = NestednessContributions(obj)
        % NOT TESTED!!!!!!
        % NestednessContributions - Perform the nestedness contribution Statistical Analysis
        %
        %   obj = NestednessContributions(obj) Perform the nestedness contribution Statistical analsysis. Be
        %   sure to create the random matrices before calling this
        %   function. Otherwise default number of replicates and null model
        %   will be used (see main/Options.m)
        
            if(isempty(obj.nulls))
                obj.DoNulls();
            end
            [obj.nest_row_contrib,obj.nest_col_contrib] = StatisticalTest.GET_NEST_CONTRIBUTIONS(obj.bipweb.matrix,obj.nulls);
            
            obj.nest_row_contrib.replicates = obj.replicates;
            obj.nest_row_contrib.model = obj.model;
            
            obj.nest_col_contrib.replicates = obj.replicates;
            obj.nest_col_contrib.model = obj.model;
            
            obj.nest_contrib_done = 1;
        end
        
        function obj = Print(obj,filename)
        
            str = '';
            if(obj.modul_done == 1)
                str = 'Modularity\n';
                str = [str, '\t Used algorithm:\t', sprintf('%30s',obj.qb_vals.algorithm), '\n'];
                str = [str, '\t Null model:    \t', sprintf('%30s',func2str(obj.qb_vals.model)), '\n'];
                str = [str, '\t Replicates:    \t', sprintf('%30i',obj.qb_vals.replicates), '\n'];
                str = [str, '\t Qb value:      \t', sprintf('%30.4f',obj.qb_vals.value), '\n'];
                str = [str, '\t     mean:      \t', sprintf('%30.4f',obj.qb_vals.mean), '\n'];
                str = [str, '\t     std:       \t', sprintf('%30.4f',obj.qb_vals.std), '\n'];
                str = [str, '\t     z-score:   \t', sprintf('%30.4f',obj.qb_vals.zscore), '\n'];
                str = [str, '\t     percentil: \t', sprintf('%30.4f',obj.qb_vals.percent), '\n'];
                str = [str, '\t Qr value:      \t', sprintf('%30.4f',obj.qr_vals.value), '\n'];
                str = [str, '\t     mean:      \t', sprintf('%30.4f',obj.qr_vals.mean), '\n'];
                str = [str, '\t     std:       \t', sprintf('%30.4f',obj.qr_vals.std), '\n'];
                str = [str, '\t     z-score:   \t', sprintf('%30.4f',obj.qr_vals.zscore), '\n'];
                str = [str, '\t     percentil: \t', sprintf('%30.4f',obj.qr_vals.percent), '\n'];
            end
            
            if(obj.nest_done == 1)
                str = [str, 'NODF Nestedness\n'];
                str = [str, '\t Null model:    \t', sprintf('%30s',func2str(obj.nodf_vals.model)), '\n'];
                str = [str, '\t Replicates:    \t', sprintf('%30i',obj.nodf_vals.replicates), '\n'];
                str = [str, '\t NODF value:    \t', sprintf('%30.4f',obj.nodf_vals.value), '\n'];
                str = [str, '\t     mean:      \t', sprintf('%30.4f',obj.nodf_vals.mean), '\n'];
                str = [str, '\t     std:       \t', sprintf('%30.4f',obj.nodf_vals.std), '\n'];
                str = [str, '\t     z-score:   \t', sprintf('%30.4f',obj.nodf_vals.zscore), '\n'];
                str = [str, '\t     percentil: \t', sprintf('%30.4f',obj.nodf_vals.percent), '\n'];
            end
            
            if(obj.temp_done == 1)
                str = [str, 'NTC Nestedness\n'];
                str = [str, '\t Null model:    \t', sprintf('%30s',func2str(obj.tempvals.model)), '\n'];
                str = [str, '\t Replicates:    \t', sprintf('%30i',obj.tempvals.replicates), '\n'];
                str = [str, '\t NTC value:     \t', sprintf('%30.4f',obj.tempvals.value), '\n'];
                str = [str, '\t     mean:      \t', sprintf('%30.4f',obj.tempvals.mean), '\n'];
                str = [str, '\t     std:       \t', sprintf('%30.4f',obj.tempvals.std), '\n'];
                str = [str, '\t     z-score:   \t', sprintf('%30.4f',obj.tempvals.zscore), '\n'];
                str = [str, '\t     percentil: \t', sprintf('%30.4f',obj.tempvals.percent), '\n'];
            end
            
            fprintf(str);  
            
            if(nargin==2)
                Printer.PRINT_TO_FILE(str,filename);
            end
            
        end
        
    end
    
    methods
        
               
        
       
        
        function [out] = GET_DEV_EIG(webbip,rmatrices)
        % NOT TESTED!!!!!!
        % GET_DEV_EIG - Perform a spectral radius Statistical Analysis
        %
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
        
        function [c_rows, c_cols] = GET_NEST_CONTRIBUTIONS(matrix, rmatrices)
        % NOT TESTED!!!!!!
        % GET_NEST_CONTRIBUTIONS - Get NODF nestedness contributions of a
        % matrix
        %
        %   [c_rows c_cols] = GET_NEST_CONTRIBUTIONS(matrix, rmatrices) Get
        %   the NODF contributions for both rows and columns. The
        %   contribution is calculated as the z-score of each individually
        %   randomly permuted row and columns
            nodf = NODF(matrix);
            N = nodf.nodf;
            [n_rows, n_cols] = size(matrix);
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
        
        
    end
    
    methods(Static)
        
        function out = GET_STATISTICAL_VALUES(real_value, random_values)
        % GET_STATISTICAL_VALUES - Get some statistics of the real value
        % with respect to random values
        %
        %   out = GET_STATISTICAL_VALUES(real_value, random_values) Get
        %   some statistics for real_value using the vector random_values. It returns
        %   an structure out with the folloging variables:
        %
        %     value: Just real_value
        %     mean: mean(random_values)
        %     std: std(random_values)
        %     ci: confidence intervals of the t-test.
        %     zscore: the z-score of real_value in the distribution of
        %     random_values: the sorted version of random_values
        %     
            n = length(random_values);
            me = mean(random_values);
            zscore = (real_value - mean(random_values))/std(random_values);
            stad = std(random_values);
            percent = 100.0*sum(real_value>random_values)/n;
            
            out.value = real_value; out.mean = me; out.std = stad;
            %out.p = p; out.ci = ci; 
            out.zscore = zscore; out.percentil = percent;
            out.random_values = sort(random_values);
            
        end
        
        function stest = MODULARITY(matrix,ntrials,null_model,modul_algorithm)
            
            if(nargin==1)
                ntrials = Options.REPLICATES;
                null_model = Options.DEFAULT_NULL_MODEL;
                modul_algorithm = Options.MODULARITY_ALGORITHM;
            elseif(nargin==2)
                null_model = Options.DEFAULT_NULL_MODEL;
                modul_algorithm = Options.MODULARITY_ALGORITHM;
            end
            
            bp = Bipartite(matrix);
            bp.modules = modul_algorithm(bp.matrix);
            stest = StatisticalTest(bp);
            stest.DoNulls(ntrials,null_model);
            stest.TestModularity();
            
            stest.Print();
        end
        
        function stest = NTC(matrix,ntrials,null_model)
            
            if(nargin==1)
                ntrials = Options.REPLICATES;
                null_model = Options.DEFAULT_NULL_MODEL;
            elseif(nargin==2)
                null_model = Options.DEFAULT_NULL_MODEL;
            end
            
            bp = Bipartite(matrix);
            stest = StatisticalTest(bp);
            stest.DoNulls(ntrials,null_model);
            stest.TestNTC();
            
            stest.Print();
            
        end
        
        function stest = NODF(matrix,ntrials,null_model)
            
            if(nargin==1)
                ntrials = Options.REPLICATES;
                null_model = Options.DEFAULT_NULL_MODEL;
            elseif(nargin==2)
                null_model = Options.DEFAULT_NULL_MODEL;
            end
            
            bp = Bipartite(matrix);
            stest = StatisticalTest(bp);
            stest.DoNulls(ntrials,null_model);
            stest.TestNODF();
            
            stest.Print();
            
        end
        
    end

end
