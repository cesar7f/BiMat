classdef MatrixGenerator < handle
    
    methods (Access = private)
    %private so that you can't instatiate.
        function out = MatrixGenerator

        end
    end 
   
    methods(Static)
       
        function matrix = CreateCheckBoardMatrix(sizeRows, sizeCols)
            %Creates cheack board distribution.
            
            if(nargin == 1)
                [sizeRows sizeCols] = size(sizeRows);
            end
            
            A = zeros(sizeRows, sizeCols);
            row1 = zeros(1,sizeCols);
            row2 = zeros(1,sizeCols);
            
            m = 1;
            for j = 1:sizeCols
                row1(j) = m;
                m = (m == 0) * 1;
                row2(j) = m;
            end 
            
            m = 1;
            for i = 1:sizeRows
                if(m==1)
                    A(i,:) = row1;
                else
                    A(i,:) = row2;
                end
                m = (m==0) * 1;
            end
 
            matrix = A;
            
        end
        
        function matrix = BernoulliRandomMatrix(MatrixOrSizeRows, sizeCols, p)
            %Create a random network where there will be 1's with
            %probability p in every cell            
             if(nargin == 1)
               
                [sizeRows sizeCols] = size(MatrixOrSizeRows);
                p = sum(sum(MatrixOrSizeRows>0))/numel(MatrixOrSizeRows);
                
            else 
                 sizeRows = MatrixOrSizeRows;
            end
            
            matrix = 1.0*(rand(sizeRows,sizeCols) < p);
            
        end
        
        function matrix = BernoulliConstrainedRandomMatrix(MatrixOrSizeRows, sizeCols, ones)
            % function x = bern_rand_matrix(R,C)
            % Creates a Bernoulli random matrix of size R, C with I 1-s 
            if(nargin == 1)
               
                [sizeRows sizeCols] = size(MatrixOrSizeRows);
                ones = sum(sum(MatrixOrSizeRows>0));
                
            else 
                 sizeRows = MatrixOrSizeRows;
            end
            
            tmpr=rand(sizeRows*sizeCols,1);
            [tmpy, tmps]=sort(tmpr);
            x=zeros(sizeRows,sizeCols);
            x(tmps(1:ones))=1;
            
            matrix = x;
        end
        
        function matrix = ColumnRandomMatrix(sizeRows, sizeCols, probCols)
            %Return a random network in wich colums filling according to
            %probCols
            
            if(sizeCols ~= length(probCols))
                error('The column size is different from the assigned probabilities');
            end
            
            randmatrix = rand(sizeRows, sizeCols);
            matrix = zeros(sizeRows, sizeCols);
            
            for i = 1:sizeCols
                
                matrix(:,i) = randmatrix(:,i) < probCols(i);
                
            end
        end
        
        function matrix = RowRandomMatrix(sizeRows, sizeCols, probRows)
            %Return a random network in wich colums filling according to
            %probCols
            
            if(sizeRows ~= length(probRows))
                error('The row size is different from the assigned probabilities');
            end
            
            randmatrix = rand(sizeRows, sizeCols);
            matrix = zeros(sizeRows, sizeCols);
            
            for i = 1:sizeRows
                
                matrix(i,:) = randmatrix(i,:) < probRows(i);
                
            end
        end
        
        function matrix = AverageRandomMatrix(sizeRows, sizeCols, probRows,probCols)
            
            if(nargin ==1)
                matrix = sizeRows;
                [sizeRows sizeCols] = size(matrix);
                probCols = sum(matrix)/sizeRows;
                probRows = sum(matrix,2)/sizeCols;
            end
            
            if(sizeRows ~= length(probRows))
                error('The row size is different from the assigned probabilities');
            end
            if(sizeCols ~= length(probCols))
                error('The column size is different from the assigned probabilities');
            end
            
            ss = size(probRows);
            if(ss(1) == 1); probRows = probRows'; end;
            
            pr = repmat(probRows,1,length(probCols));
            pc = repmat(probCols,length(probRows),1);
            
            %display(pr);
            %display(pc);
            
            mat = (pr+pc)/2;
            
            matrix = rand(sizeRows,sizeCols) <= mat;
        end
        
        function matrix = UnitMatrix(sizerow,sizecol)
           
            if(nargin == 1) sizecol = sizerow; end;
                
            matrix = eye(sizerow,sizecol);
        end
            
        function matrix = UnionMatrix(matrix1, matrix2)
                    
            s1 = size(matrix1);
            s2 = size(matrix2);

            matrix1((s1(1)+1):(s1(1)+s2(1)),(s1(2)+1):(s1(2)+s2(2))) = matrix2;

            matrix = matrix1;
        end
        
        function matrix = ModularMatrix(sizerow,sizecol,sizecomrow,sizecomcol)
        %Return a network with a modular matrix of size sizerow x sizecol and
        %communities sizes of sizecomrow x sizecomcol.
            matrix = zeros(sizerow,sizecol);
            
            i = 1;j = 1;
            while(i <= sizerow && j <= sizecol)
               
                matrix(i:min(sizecomrow+i-1,sizerow),j:min(sizecomcol+j-1,sizecol)) = 1;
                
                i = sizecomrow + i;
                j = sizecomcol + j;
                
            end
            
        end
        
        function matrix = ModularChessMatrix(sizerow,sizecol,sizecomrow,sizecomcol)
        %Return a network with a modular matrix of size sizerow x sizecol and
        %communities sizes of sizecomrow x sizecomcol.
            matrix = zeros(sizerow,sizecol);
            
            i = 1;j = 1;
            while(i <= sizerow && j <= sizecol)
               
                matrix(i:min(sizecomrow+i-1,sizerow),j:min(sizecomcol+j-1,sizecol)) = 1;
                
                i = sizecomrow + i;
                j = sizecomcol + j;
                
            end
            
        end
        
        function matrix = PerfectNestedMatrix(nrows,ncols,fill)
        %Return a perfect square nested network of fill = 0.5
               
            matrix = NestednessBINMATNEST.PerfectNested(nrows,ncols,fill);
            
        end
        
        function matrix = HierarchicalMatrix(nlevels,minsize)
           
            nsize = minsize * 2^(nlevels-1);
            
            matrix = zeros(nsize, nsize);
            
            for i = 1:nlevels
               
                for j = 1:2^(i-1)
                   
                    localsize = nsize/(2^(i-1));
                    matrix( (j-1)*localsize+1 : j*localsize , (j-1)*localsize+1 : j*localsize ) = i^2;
                    matrix
                end
                
            end
            
        end
        
        function matrix = NestedCommunities(nrows, ncols, n, fill)
        %Return a square of size nsize matrices with n communities perfectly nested with
        %a fill of 0.5 in each community
            
            if(nargin==3) fill = 0.5; end; 
        
            nrealrow = floor(nrows/n);
            resrow = nrows - nrealrow*n;
            nrealcol = floor(ncols/n);
            rescol = ncols - nrealcol*n;
            
            matrix = zeros(nrows,ncols);
            
            rowbig = 0;
            colbig = 0;
            indexrow = 1;
            indexcol = 1;
            for i = 1:n
               
                sizerow = (rowbig<resrow)*(nrealrow+1) + (rowbig>=resrow)*nrealrow;
                sizecol = (colbig<rescol)*(nrealcol+1) + (colbig>=rescol)*nrealcol;
                
                matrix(indexrow:indexrow+sizerow-1,indexcol:indexcol+sizecol-1) = MatrixGenerator.PerfectNestedMatrix(sizerow,sizecol,fill);
                indexrow = indexrow+sizerow;
                indexcol = indexcol+sizecol;
                
                rowbig = rowbig + 1;
                colbig = colbig + 1;
            end
            
        end
        
        function matrix = NestedCommunitiesOld(nsize, n)
        %Return a square of size nsize matrices with n communities perfectly nested with
        %a fill of 0.5 in each community
            
            nreal = floor(nsize/n);
            res = nsize - nreal*n;
            
            matrix = zeros(nsize,nsize);
            
            for i = 1:res
               
                start = (nreal+1)*(i-1) + 1;
                endc = (nreal+1)*i;
                matrix(start:endc,start:endc) = flipud(tril(ones(nreal+1,nreal+1)));
                
            end
            
            cont = (nreal+1)*res;
            for i = 1:(n-res)
               
                start = cont + nreal*(i-1) + 1;
                endc = cont + nreal*i;
                matrix(start:endc,start:endc) = flipud(tril(ones(nreal,nreal)));
                
            end
        end
    end
end