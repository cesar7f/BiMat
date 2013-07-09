%==========================================================================
% Name        : NestednessBINMATNEST.m
% Author      : C�sar Flores
% Created     : 22/Jun/2010
% Updated     : 22/Jun/2010
% Description : Represents the algorithm developped by Rodriguez and Santamaria
%==========================================================================

classdef NestednessNODF < Nestedness

    properties(GetAccess = 'public', SetAccess = 'private')
        Np     = 0;
    end
    
    %CONSTRUCTOR AND MAIN PROCEDURE ALGORITHM
    methods
        function obj = NestednessNODF(bipNetwork)
            obj.Matrix = bipNetwork.Matrix > 0; %Normalize the matrix
            obj.nRows = bipNetwork.nRows;  %Number of Rows
            obj.nCols = bipNetwork.nCols; %Number of Columns
            obj.IndexRow = 1:obj.nRows;
            obj.IndexCol = 1:obj.nCols;
            obj.Fill = sum(sum(obj.Matrix>0)) / (bipNetwork.nRows*bipNetwork.nCols);
            obj.Name = bipNetwork.Name;
            obj.NetworkBipartite = bipNetwork;
        end
        
        
        
        function obj = CalculateNestedness(obj)
           
            m = obj.nRows;
            n = obj.nCols;
            denom = n*(n-1)/2 + m*(m-1)/2;
            
            if(m == 1 || n == 1)
                obj.N = 0;
                return;
            end
            
            obj.SortMatrix();
            obj.CalculateNpaired();
            obj.N = obj.Np / (100*denom);
               
        end
           
        function obj = CalculateNestednessWithoutSorting(obj)
           
            m = obj.nRows;
            n = obj.nCols;
            denom = n*(n-1)/2 + m*(m-1)/2;
            
            obj.CalculateNpairedInternalCommunities();
            obj.N = obj.Np / (100*denom);
            
        end
        
    end

    methods
       
        function obj = CalculateNpaired(obj)

            sumrows = sum(obj.Matrix,2);
            sumcols = sum(obj.Matrix,1);
            
            obj.Np = 0;
            %Fill for columns
            for i = 1:obj.nRows
                for j = i+1:obj.nRows
                    if( sumrows(j) < sumrows(i) && sumrows(j) > 0 )
                        obj.Np = obj.Np + 100*sum(obj.Matrix(i, obj.Matrix(j,:)==1))/sumrows(j);
                    end
                end
            end
            
            for k = 1:obj.nCols
                for l = k+1:obj.nCols
                    if( sumcols(l) < sumcols(k) && sumcols(l) > 0 )
                        obj.Np = obj.Np + 100*sum(obj.Matrix(obj.Matrix(:,l)==1,k))/sumcols(l);
                    end
                end
            end
            
        end
        
        
        function obj = CalculateNpairedExternalCommunities(obj)
            
            bn = obj.NetworkBipartite;
            obj.Matrix = bn.Modularity.SortedMatrix ~= 0;
            
            ssRows = bn.Modularity.SSRow;
            ssCols = bn.Modularity.SSCol;
            
            sumrows = sum(obj.Matrix,2);
            sumcols = sum(obj.Matrix,1);
            
            obj.Np = 0;
            %nr = 0;
            for j = 1:obj.nRows
                for k = j+1:obj.nRows
                    if(ssRows(j) ~= ssRows(k))
                        if(sumrows(j) > sumrows(k) && sumrows(k) > 0)
                            obj.Np = obj.Np + 100*sum(obj.Matrix(j, obj.Matrix(k,:)==1))/sumrows(k);
                        elseif(sumrows(j) < sumrows(k) && sumrows(j) > 0)
                            obj.Np = obj.Np + 100*sum(obj.Matrix(k, obj.Matrix(j,:)==1))/sumrows(j);
                        end
                    
                    end
                   
                end
            end
            
            for j = 1:obj.nCols
                for k = j+1:obj.nCols
                    if(ssCols(j) ~= ssCols(k))
                        if(sumcols(j) > sumcols(k) && sumcols(k) > 0)
                            obj.Np = obj.Np + 100*sum(obj.Matrix(obj.Matrix(:,k)==1,j))/sumcols(k);
                        elseif(sumcols(j) < sumcols(k) && sumcols(j) > 0)
                            obj.Np = obj.Np + 100*sum(obj.Matrix(obj.Matrix(:,j)==1,k))/sumcols(j);
                        end
                    end
                end
            end
            
        end
        
        function obj = CalculateNpairedInternalCommunities(obj)

            bn = obj.NetworkBipartite;
            obj.Matrix = bn.Modularity.SortedMatrix ~= 0;
            
            colComSize = bn.Modularity.ColComSize{bn.Modularity.FinalLevel};
            rowComSize = bn.Modularity.RowComSize{bn.Modularity.FinalLevel};
            ssRow = bn.Modularity.SSRow;
            ssCol = bn.Modularity.SSCol;
            ncom = bn.Modularity.CommunityQuantity(bn.Modularity.FinalLevel);
            
            sumrows = sum(obj.Matrix,2);
            sumcols = sum(obj.Matrix,1);
            
            obj.Np = 0;
            for i = 1:ncom
               
                rowcomsize = rowComSize(i);
                colcomsize = colComSize(i);

                irow = find(ssRow==i,1);
                icol = find(ssCol==i,1);
                
                for j = irow:(irow+rowcomsize-1)
                    for k = j+1:(irow+rowcomsize-1)
                        if(sumrows(j) > sumrows(k) && sumrows(k) > 0)
                            obj.Np = obj.Np + 100*sum(obj.Matrix(j, obj.Matrix(k,:)==1))/sumrows(k);
                        elseif(sumrows(j) < sumrows(k) && sumrows(j) > 0)
                            obj.Np = obj.Np + 100*sum(obj.Matrix(k, obj.Matrix(j,:)==1))/sumrows(j);
                        end
                    end
                end
                
                for j = icol:(icol+colcomsize-1)
                    for k = j+1:(icol+colcomsize-1)
                        if(sumcols(j) > sumcols(k) && sumcols(k) > 0)
                            obj.Np = obj.Np + 100*sum(obj.Matrix(obj.Matrix(:,k)==1,j))/sumcols(k);
                        elseif(sumcols(j) < sumcols(k) && sumcols(j) > 0)
                            obj.Np = obj.Np + 100*sum(obj.Matrix(obj.Matrix(:,j)==1,k))/sumcols(j);
                        end
                    end
                end
                
            end
            
        end
        
    end
end
