% Printer - Class used for printing output files on screen and files.
%
% Printer Methods:
%    Printer - Main Constructor
%    CreateCytoscapeData - Create input files for Cytoscape software
%
% See also:
%    Reading
classdef Printer < handle
    
    properties
        
        delimiter = ',';
        bipweb = {};
        
    end
    
    methods
        
        function obj = Printer(webbip)
        % PRINTER - Main Constructor
        %   pt = Printer(webbip) Create a Printer object pt wich will be
        %   related to Bipartite object webbip.
        
            obj.bipweb = webbip;
            
        end
        
%         function obj = SpeciesLevel(obj, toscreen, tofile)
%            
%             webbip = obj.bipweb;
%             
%             aindex = webbip.name==' ';
%             namebip = webbip.name;
%             namebip(aindex) = '_';
%             
%             if(nargin == 1)
%                 toscreen = 1;
%                 tofile = 1;
%                 filename = [namebip, '_sp.txt'];
%             elseif(nargin == 2)
%                 tofile = 1;
%                 filename = [namebip, '_sp.txt'];
%             elseif (nargin == 3)
%                 filename = [namebip, '_sp.txt'];
%             end
%             
%             header = '';
%             colwitdth = 9;
%             formatstring = ['%-',num2str(colwitdth),'s'];
%             formatint = ['%-',num2str(colwitdth),'d'];
%             formatfloat = ['%-',num2str(colwitdth),'.3f'];
%             estring = '----';
%             
%             header = [header, sprintf(repmat(formatstring,1,6), 'Specie','Level','Degree','SPE','SSI','RR')];
%             if(webbip.modules.done); header = [header, sprintf(repmat(formatstring,1,1), 'Module')]; end;
%             header = [header, '\n'];
%             
%             info = '';
%             for i = 1:webbip.n_rows
%                 aindex = webbip.row_labels{i}==' ';
%                 label_name = webbip.row_labels{i};
%                 label_name(aindex) = '_';
%                 info = [info, sprintf(formatstring,label_name(1:min(colwitdth,length(label_name))))];
%                 info = [info, sprintf(formatstring, 'top')];
%                 info = [info, sprintf(formatint, webbip.row_degrees(i))];
%                 info = [info, sprintf(formatfloat, webbip.specificity(i))];
%                 info = [info, sprintf(formatfloat, webbip.ssi(i))];
%                 info = [info, sprintf(formatfloat, webbip.rr(i))];
%                 if(webbip.modules.done); info = [info, sprintf(formatint, webbip.modules.row_modules(i))]; end;
%                 info = [info, '\n'];
%             end
%             
%             for i = 1:webbip.n_cols
%                 aindex = webbip.col_labels{i}==' ';
%                 label_name = webbip.col_labels{i};
%                 label_name(aindex) = '_';
%                 info = [info, sprintf(formatstring,label_name(1:min(colwitdth,length(label_name))))];
%                 info = [info, sprintf(formatstring, 'bottom')];
%                 info = [info, sprintf(formatint, webbip.col_degrees(i))];
%                 info = [info, sprintf(repmat(formatstring,1,3),estring,estring,estring)];
%                 if(webbip.modules.done); info = [info, sprintf(formatint, webbip.modules.col_modules(i))]; end;
%                 info = [info, '\n'];
%             end
%             
%             if(toscreen)
%                 fprintf(header);
%                 fprintf(info);
%             end
%             
%             if(tofile)
%                 fid = fopen([namebip,'-sp.txt'],'w');
%                 fprintf(fid,header);
%                 fprintf(fid,info);
%                 fclose(fid);
%             end
%             
%         end
%         
%         function obj = NetworkLevelSingle(obj, toscreen, tofile, filename)
%             To complete, maybe have to delete at the end of the day
%             webbip = obj.bipweb;
%             
%             if(nargin == 1)
%                 toscreen = 1;
%                 tofile = 1;
%                 filename = [webbip.name, '_networklevel.txt'];
%             elseif(nargin == 2)
%                 tofile = 1;
%                 filename = [webbip.name, '_networklevel.txt'];
%             elseif (nargin == 3)
%                 filename = [webbip.name, '_networklevel.txt'];
%             end
%             format = '%5.2f';
%             formatint = '%5d';
%             fileexist = exist(filename,'file');
%             header = '';
%             if(~fileexist)
%                 header = [header, '\nGeneral properties\n'];
%                 header = [header,   '------------------\n'];
%                 header = [header, '\tNetwork Name:      ', webbip.name,'\n'];
%                 header = [header, '\tRows:              ', num2str(webbip.n_rows,formatint),'\n'];
%                 header = [header, '\tCols:              ', num2str(webbip.n_cols,formatint),'\n'];
%                 header = [header, '\tSize:              ', num2str(webbip.size_webmatrix,formatint),'\n'];
%                 header = [header, '\tInteractions:      ', num2str(webbip.n_edges,formatint),'\n'];
%                 header = [header, '\tConnectance:       ', num2str(webbip.connectance,format),'\n'];
%                 header = [header, '\t<Specificity>:     ', num2str(mean(webbip.specificity),format),'\n'];
%                 header = [header, '\t<Res. Range.>:     ', num2str(mean(webbip.rr),format),'\n'];
%                 header = [header, '\tNODF:              ', num2str(webbip.nestedness.nodf,format),'\n'];
%                 header = [header, '\tNODF Rows:         ', num2str(webbip.nestedness.nodf_rows,format),'\n'];
%                 header = [header, '\tNODF Cols:         ', num2str(webbip.nestedness.nodf_cols,format),'\n'];                
%             end
%             
%             if(toscreen)
%                 fprintf(header);
%             end
%             
%         end
%         
%         function obj = NetworkLevel(obj, toscreen, tofile)
%             webbip = obj.bipweb;
%             filename = './networklevel.txt';
%             
%             if(nargin == 1)
%                 toscreen = 1;
%                 tofile = 1;
%             else(nargin == 2)
%                 tofile = 1;
%             end
%             
%             header = '';
%             colwitdth = 9;
%             formatstring = ['%-',num2str(colwitdth),'s'];
%             formatint = ['%-',num2str(colwitdth),'d'];
%             formatfloat = ['%-',num2str(colwitdth),'.3f'];
%             estring = '----';
%             
%             header = [header, sprintf(repmat(formatstring,1,5), 'Name','Conn','Size','n_rows','n_cols')];
%             header = [header, sprintf(repmat(formatstring,1,5), 'nodf','nodf_r','nodf_c','aspe','arr')];
%             header = [header, sprintf(repmat(formatstring,1,5), 'resp','inc','mN','mQb','mQr')];
%             header = [header, sprintf(repmat(formatstring,1,6), 'n_sim','n_pval','nicLow','nicUp','n_z','nco_sim')];
%             header = [header, sprintf(repmat(formatstring,1,5), 'nco_pval','ncoicLow','ncoicUp','nco_z','nro_sim')];
%             header = [header, sprintf(repmat(formatstring,1,4), 'nro_pval','nroicLow','nroicUp','nro_z')];
%             header = [header, sprintf(repmat(formatstring,1,4), 'mQr_sim','mQr_val','mQrIclow','mQrIcup','mQr_z')];
%             header = [header, sprintf(repmat(formatstring,1,4), 'mQb_sim','mQb_val','mQbIclow','mQbIcup','mQb_z')];
%             header = [header, sprintf(repmat(formatstring,1,3), 'tReps','null_name')];
%             header = [header, '\n'];
%             header = [header, '-----------------------------------------------------------------\n'];
%             
%             [R I] = SpeFunc.IR(webbip.webmatrix);
%             aindex = webbip.name==' ';
%             namebip = webbip.name;
%             namebip(aindex) = '_';
%             info = sprintf(formatstring,namebip(1:min(colwitdth-1,length(namebip))));
%             info = [info, sprintf(formatfloat, webbip.connectance)];
%             info = [info, sprintf(formatint, webbip.size_webmatrix)];
%             info = [info, sprintf(formatint, webbip.n_rows)];
%             info = [info, sprintf(formatint, webbip.n_cols)];
%             info = [info, sprintf(formatfloat, webbip.nestedness.nodf)];
%             info = [info, sprintf(formatfloat, webbip.nestedness.nodf_rows)];
%             info = [info, sprintf(formatfloat, webbip.nestedness.nodf_cols)];
%             info = [info, sprintf(formatfloat, mean(webbip.specificity))];
%             info = [info, sprintf(formatfloat, mean(webbip.rr))];
%             info = [info, sprintf(formatfloat, R)];
%             info = [info, sprintf(formatfloat, I)];
%             
%             if(webbip.modules.done)
%                 info = [info, sprintf(formatint, webbip.modules.N)];
%                 info = [info, sprintf(formatfloat, webbip.modules.Qb)];
%                 info = [info, sprintf(formatfloat, webbip.modules.Qr)];
%             else
%                 info = [info, sprintf(repmat(formatstring,1,3),estring,estring,estring)];
%             end
%             
%             if(webbip.statistics.nest_done)
%                 info = [info, sprintf(formatfloat, webbip.statistics.nestvals.ci(1))];
%                 info = [info, sprintf(formatfloat, webbip.statistics.nestvals.p)];
%                 info = [info, sprintf(formatfloat, webbip.statistics.nestvals.ci(2))];
%                 info = [info, sprintf(formatfloat, webbip.statistics.nestvals.ci(3))];
%                 info = [info, sprintf(formatfloat, webbip.statistics.nestvals.zscore)];
%                 info = [info, sprintf(formatfloat, webbip.statistics.nestvals_cols.ci(1))];
%                 info = [info, sprintf(formatfloat, webbip.statistics.nestvals_cols.p)];
%                 info = [info, sprintf(formatfloat, webbip.statistics.nestvals_cols.ci(2))];
%                 info = [info, sprintf(formatfloat, webbip.statistics.nestvals_cols.ci(3))];
%                 info = [info, sprintf(formatfloat, webbip.statistics.nestvals_cols.zscore)];
%                 info = [info, sprintf(formatfloat, webbip.statistics.nestvals_rows.ci(1))];
%                 info = [info, sprintf(formatfloat, webbip.statistics.nestvals_rows.p)];
%                 info = [info, sprintf(formatfloat, webbip.statistics.nestvals_rows.ci(2))];
%                 info = [info, sprintf(formatfloat, webbip.statistics.nestvals_rows.ci(3))];
%                 info = [info, sprintf(formatfloat, webbip.statistics.nestvals_rows.zscore)];
%             else
%                 info = [info, sprintf(repmat(formatstring,1,12),estring,estring,estring,estring,estring,estring, ...
%                     estring,estring,estring,estring,estring,estring)];
%             end
%             
%             if(webbip.statistics.modul_done)
%                 info = [info, sprintf(formatfloat, webbip.statistics.qr_vals.ci(1))];
%                 info = [info, sprintf(formatfloat, webbip.statistics.qr_vals.p)];
%                 info = [info, sprintf(formatfloat, webbip.statistics.qr_vals.ci(2))];
%                 info = [info, sprintf(formatfloat, webbip.statistics.qr_vals.ci(3))];
%                 info = [info, sprintf(formatfloat, webbip.statistics.qr_vals.zscore)];
% 
%                 info = [info, sprintf(formatfloat, webbip.statistics.qb_vals.ci(1))];
%                 info = [info, sprintf(formatfloat, webbip.statistics.qb_vals.p)];
%                 info = [info, sprintf(formatfloat, webbip.statistics.qb_vals.ci(2))];
%                 info = [info, sprintf(formatfloat, webbip.statistics.qb_vals.ci(3))];
%                 info = [info, sprintf(formatfloat, webbip.statistics.qb_vals.zscore)];
%             else
%                 info = [info, sprintf(repmat(formatstring,1,8),estring,estring,estring, ... 
%                     estring,estring,estring,estring,estring)];
%             end
%             
%             if(webbip.statistics.modul_done || webbip.statistics.nest_done)
%                 info = [info, sprintf(formatint, webbip.statistics.replicates)];
%                 info = [info, sprintf(formatstring, func2str(webbip.statistics.model))];
%             else
%                 info = [info, sprintf(repmat(formatstring,1,2),estring,estring)];
%             end
%             
%             info = [info, '\n'];
%             
%             if(toscreen)
%                 fprintf(header);
%                 fprintf(info);
%             end
%             
%             if(tofile)
%                 if(~exist(filename,'file')) 
%                     fid = fopen(filename,'w');
%                     fprintf(fid, header);
%                 else
%                     fid = fopen(filename,'a');
%                 end
%                 fprintf(fid,info);
%                 fclose(fid);
%             end
%             
%         end
        
        function obj = CreateCytoscapeData(obj,filename)
        % CreateCytoscapeData - Create input files for Cytoscape software
        %   obj = CreateCytoscapeData(obj,filename) Create two input files
        %   filename_edges.csv and filename_nodes.csv for Cytoscape
        %   software. The first one is indispensable for Cytoscape and
        %   contains the list of edges of the bipartite network. The second
        %   contain additional information of the nodes such as module id,
        %   name label, etc. This additional information can be used in
        %   Cytoscape for creating layouts with more information.
            
            n_rows = obj.bipweb.n_rows;
            n_cols = obj.bipweb.n_cols;
            
            %inter = sum(sum(obj.bipweb.matrix));
            %adlist = zeros(inter,3);
            
            matrix = obj.bipweb.webmatrix;
            
            %nn = 1;

            fidedges = fopen([filename,'_edges.csv'],'w');
            fidnodes = fopen([filename,'_nodes.csv'],'w');
            
            for i = 1:n_rows; fprintf(fidedges,'%i\n',i); end;
            for j = 1:n_cols; fprintf(fidedges,'%i\n',j+n_cols); end;
            
            modules_done = obj.bipweb.modules.done;
            
            row_modul = mod(find(obj.bipweb.modules.rr'==1),obj.bipweb.modules.N);
            col_modul = mod(find(obj.bipweb.modules.tt'==1),obj.bipweb.modules.N);
            row_modul = row_modul.*(row_modul>0)+(row_modul==0).*obj.bipweb.modules.N;
            col_modul = col_modul.*(col_modul>0)+(col_modul==0).*obj.bipweb.modules.N;
            
            for i = 1:n_rows
                for j = 1:n_cols
                    if(matrix(i,j) > 0)
                        if(~modules_done)
                            fprintf(fidedges,'%i %i %i\n',i,j+n_rows,matrix(i,j));
                        else
                            if(row_modul(i)==col_modul(j))
                                fprintf(fidedges,'%i %i %i %i\n',i,j+n_rows,matrix(i,j),row_modul(i));
                            else
                                fprintf(fidedges,'%i %i %i %i\n',i,j+n_rows,matrix(i,j),0);
                            end
                        end
                    end
                end
            end
            
            fclose(fidedges);
               
            fprintf(fidnodes,'ID,NAME,TYPE,MODULE\n');
            row_labels = obj.bipweb.row_labels;
            col_labels = obj.bipweb.col_labels;
            
            for i = 1:n_rows
                fprintf(fidnodes,'%i,%s,%i,%i\n',i,row_labels{i},1,row_modul(i));
            end
            for j = 1:n_cols
                fprintf(fidnodes,'%i,%s,%i,%i\n',j+n_rows,col_labels{j},2,col_modul(j));
            end
            
            fclose(fidnodes);
%                
%             for j = 1:n_cols
%                 fprintf(fidnodes,'%i %i 2\n',j+n_rows,find(tt_final{cid}(j,:)));
%             end
%             
%             for i = 1:nrows 
%                 for j = 1:n_cols
%                     if(matrix(i,j) > 0)
%                         adlist(nn,:) = [i matrix(i,j) (j+nrows)];
%                         nn = nn+1;
%                     end
%                 end
%             end
            
            %dlmwrite([filename,'.sif'], adlist, 'delimiter', ' ', 'precision', 6);
            %nodeatribs = fopen([filename,'.pvals'],'w');
            %fprintf(nodeatribs,'ID\tDEGREE\tTYPE\n');
            %for i = 1:obj.nRows
            %    fprintf(nodeatribs,'%i\t%i\t%i\n',i,sum(obj.Matrix(i,:)), 1);
            %end
            %for i = 1:obj.nCols
            %    fprintf(nodeatribs,'%i\t%i\t%i\n',i+obj.nRows,sum(obj.Matrix(:,i)), 2);
            %end
        end
        
    end
    
    methods
       
        function PrintGeneralProperties(obj,filename)
            
            str = 'General Properties\n';
            str = [str, '\t Number of species:       \t', sprintf('%6i',obj.bipweb.n_rows+obj.bipweb.n_cols), '\n'];
            str = [str, '\t Number of row species:   \t', sprintf('%6i',obj.bipweb.n_rows), '\n'];
            str = [str, '\t Number of column species:\t', sprintf('%6i',obj.bipweb.n_cols), '\n'];
            str = [str, '\t Number of Interactions:  \t', sprintf('%6i',obj.bipweb.n_edges), '\n'];
            str = [str, '\t Size:                    \t', sprintf('%6i',obj.bipweb.size_webmatrix), '\n'];
            str = [str, '\t Connectance or fill:     \t', sprintf('%6.3f',obj.bipweb.connectance), '\n'];
            
            fprintf(str);
            
            if(nargin==2)
                Printer.PRINT_TO_FILE(str,filename);
            end
            
        end
        
        function PrintStructureValues(obj,filename)
           
            if(obj.bipweb.modules.done == 0)
                obj.bipweb.modules.Detect();
            end
            
            if(obj.bipweb.ntc.done == 0)
                obj.bipweb.ntc.CalculateNestedness;
            end
            
            str = 'Modularity\n';
            str = [str, '\tUsed algorithm:             \t', sprintf('%16s',class(obj.bipweb.modules)), '\n'];
            str = [str, '\tN (Number of modules):      \t', sprintf('%16i',obj.bipweb.modules.N), '\n'];
            str = [str, '\tQb (Standard metric):       \t', sprintf('%16.4f',obj.bipweb.modules.Qb), '\n'];
            str = [str, '\tQr (Ratio of int/ext inter):\t', sprintf('%16.4f',obj.bipweb.modules.Qr), '\n'];
            
            str = [str, 'Nestedness\n'];
            str = [str, '\tNODF metric value:          \t', sprintf('%16.4f',obj.bipweb.nodf.nodf), '\n'];
            str = [str, '\tNTC metric value:           \t', sprintf('%16.4f',obj.bipweb.ntc.N), '\n'];
           
            fprintf(str);  
            
            if(nargin==2)
                Printer.PRINT_TO_FILE(str,filename);
            end
            
        end
        
        function PrintStructureStatistics(obj,filename)
           
            if(obj.bipweb.statistics.modul_done == 0)
                obj.bipweb.statistics.Modularity();
            end
            
            if(obj.bipweb.statistics.nest_done == 0)
                obj.bipweb.statistics.Nestedness();
            end
            
            if(obj.bipweb.statistics.temp_done == 0)
                obj.bipweb.statistics.Temperature();
            end
            
            str = 'Modularity\n';
            str = [str, '\t Used algorithm:\t', sprintf('%30s',obj.bipweb.statistics.qb_vals.algorithm), '\n'];
            str = [str, '\t Null model:    \t', sprintf('%30s',func2str(obj.bipweb.statistics.qb_vals.model)), '\n'];
            str = [str, '\t Replicates:    \t', sprintf('%30i',obj.bipweb.statistics.qb_vals.replicates), '\n'];
            str = [str, '\t Qb value:      \t', sprintf('%30.4f',obj.bipweb.statistics.qb_vals.value), '\n'];
            str = [str, '\t z-score:       \t', sprintf('%30.4f',obj.bipweb.statistics.qb_vals.zscore), '\n'];
            str = [str, '\t percent:       \t', sprintf('%30.4f',obj.bipweb.statistics.qb_vals.percent), '\n'];
            
            str = [str, 'NODF Nestedness\n'];
            str = [str, '\t Null model:    \t', sprintf('%30s',func2str(obj.bipweb.statistics.nestvals.model)), '\n'];
            str = [str, '\t Replicates:    \t', sprintf('%30i',obj.bipweb.statistics.nestvals.replicates), '\n'];
            str = [str, '\t NODF value:    \t', sprintf('%30.4f',obj.bipweb.statistics.nestvals.value), '\n'];
            str = [str, '\t z-score:       \t', sprintf('%30.4f',obj.bipweb.statistics.nestvals.zscore), '\n'];
            str = [str, '\t percent:       \t', sprintf('%30.4f',obj.bipweb.statistics.nestvals.percent), '\n'];
            
            str = [str, 'NTC Nestedness\n'];
            str = [str, '\t Null model:    \t', sprintf('%30s',func2str(obj.bipweb.statistics.tempvals.model)), '\n'];
            str = [str, '\t Replicates:    \t', sprintf('%30i',obj.bipweb.statistics.tempvals.replicates), '\n'];
            str = [str, '\t NTC value:     \t', sprintf('%30.4f',obj.bipweb.statistics.tempvals.value), '\n'];
            str = [str, '\t z-score:       \t', sprintf('%30.4f',obj.bipweb.statistics.tempvals.zscore), '\n'];
            str = [str, '\t percent:       \t', sprintf('%30.4f',obj.bipweb.statistics.tempvals.percent), '\n'];

            fprintf(str);  
            
            if(nargin==2)
                Printer.PRINT_TO_FILE(str,filename);
            end
            
        end
        
        function PrintStructureStatisticsOfModules(obj,filename)
            
            if(nargin == 2)
                obj.bipweb.internal_statistics.gtesting.PrintResults(filename);
            else
                obj.bipweb.internal_statistics.gtesting.PrintResults();
            end
            
        end
        
        function PrintRowModuleDiversity(obj, filename)
           
            headers{1} = 'Module';
            headers{2} = 'index value';
            headers{3} = 'zscore';
            headers{4} = 'percent';
            
            columns = (1:length(obj.bipweb.internal_statistics.row_diversity.value))';
            columns = [columns obj.bipweb.internal_statistics.row_diversity.value'];
            columns = [columns obj.bipweb.internal_statistics.row_diversity.zscore'];
            columns = [columns obj.bipweb.internal_statistics.row_diversity.percent'];
            
            str = Printer.CREATE_FORMATED_STRING(headers,columns,',');
            
            str = ['Random permutations:\t', sprintf('%25i',obj.bipweb.internal_statistics.row_diversity.n_permutations), '\n', str];
            str = ['Diversity index:    \t', sprintf('%25s',func2str(obj.bipweb.internal_statistics.row_diversity.diversity_index)), '\n', str];
            
            fprintf(str);
            
            if(nargin==2)
                Printer.PRINT_TO_FILE(str,filename);
            end
        end
        
        function PrintColumnModuleDiversity(obj, filename)
           
            headers{1} = 'Module';
            headers{2} = 'index value';
            headers{3} = 'zscore';
            headers{4} = 'percent';
            
            columns = (1:length(obj.bipweb.internal_statistics.col_diversity.value))';
            columns = [columns obj.bipweb.internal_statistics.col_diversity.value'];
            columns = [columns obj.bipweb.internal_statistics.col_diversity.zscore'];
            columns = [columns obj.bipweb.internal_statistics.col_diversity.percent'];
            
            str = Printer.CREATE_FORMATED_STRING(headers,columns,',');
            
            str = ['Random permutations:\t', sprintf('%25i',obj.bipweb.internal_statistics.col_diversity.n_permutations), '\n', str];
            str = ['Diversity index:    \t', sprintf('%25s',func2str(obj.bipweb.internal_statistics.col_diversity.diversity_index)), '\n', str];
            
            fprintf(str);
            
            if(nargin==2)
                Printer.PRINT_TO_FILE(str,filename);
            end
            
        end
        
        function PrintGeneralPropertiesOfModules(obj, filename)
            
            fprintf('No. \t H \t P \t S \t I \t M \t C \t Lh \t Lp\n');
            
            headers{1} = 'Module';
            headers{2} = 'm';
            headers{3} = 'n';
            headers{4} = 'S';
            headers{5} = 'I';
            headers{6} = 'M';
            headers{7} = 'C';
            headers{8} = 'Lr';
            headers{9} = 'Lc';
            
            columns = (1:length(obj.bipweb.internal_statistics.module_networks))';
            
            columns = [columns cellfun(@(x) x.n_rows, obj.bipweb.internal_statistics.module_networks)];
            columns = [columns cellfun(@(x) x.n_cols, obj.bipweb.internal_statistics.module_networks)];
            columns = [columns cellfun(@(x) x.n_cols+x.n_rows, obj.bipweb.internal_statistics.module_networks)];
            columns = [columns cellfun(@(x) x.n_edges, obj.bipweb.internal_statistics.module_networks)];
            columns = [columns cellfun(@(x) x.size_webmatrix, obj.bipweb.internal_statistics.module_networks)];
            columns = [columns cellfun(@(x) x.connectance, obj.bipweb.internal_statistics.module_networks)];
            columns = [columns cellfun(@(x) x.n_edges/x.n_rows, obj.bipweb.internal_statistics.module_networks)];
            columns = [columns cellfun(@(x) x.n_edges/x.n_cols, obj.bipweb.internal_statistics.module_networks)];
            
            str = Printer.CREATE_FORMATED_STRING(headers,columns,obj.delimiter);
            
            fprintf(str);
            
            if(nargin==2)
                Printer.PRINT_TO_FILE(str,filename);
            end
            
        end
        
    end
    
    methods(Static)
        
        function PRINT_TO_FILE(str,filename)
            fid = fopen(filename,'w');
            fprintf(fid,str);            
            fclose(fid);
        end
        
        function formated_string = CREATE_FORMATED_STRING(headers,columns,delimiter)
           
            assert((isa(columns,'cell') && length(headers) == length(columns)) ...
                || size(columns,2)==length(headers));
            
            n_cols = length(headers);
            
            if(nargin==2)
                delimiter = ',';
            end
            
            
            
            if(isa(columns,'double'))
                columns_str = arrayfun(@num2str, columns, 'unif', 0);
                [n_r n_c] = size(columns_str);
                spacing = zeros(n_cols,1);
                for i = 1:n_cols
                    spacing(i) = max(max(cellfun('length',columns_str(:,i))),length(headers{i}));
                end
                
                
                row_size = (sum(spacing)+n_cols+1)*(n_r+1);
                %allocate memory for char
                formated_string(row_size) = char(0);
                i_char = 1;
                for i = 1:n_cols
                    format_data = ['%',num2str(spacing(i)),'s'];
                    formated_string(i_char:i_char+spacing(i)-1) = sprintf(format_data, headers{i});
                    if(i ~= n_cols)
                        formated_string(i_char+spacing(i)) = delimiter;
                    end
                    i_char = i_char + spacing(i)+1;
                end
                formated_string(i_char-1:i_char) = '\n';
                i_char = i_char+1;
                for j = 1:n_r
                    for i = 1:n_c
                        format_data = ['%',num2str(spacing(i)),'s'];
                        formated_string(i_char:i_char+spacing(i)-1) = ...
                            sprintf(format_data, columns_str{j,i});
                        if(i ~= n_cols)
                            formated_string(i_char+spacing(i)) = delimiter;
                        end
                        i_char = i_char + spacing(i)+1;
                    end
                    formated_string(i_char-1:i_char) = '\n';
                    i_char = i_char+1;
                end
            end
            
            
        end
        
    end
    
end