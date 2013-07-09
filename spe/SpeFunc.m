classdef SpeFunc
% UNTESTED_CLASS - BE AWARE   
    methods(Static)
        
        function sp_pdi = SP_PDI(vec)
           
            sp_pdi = 0;
            vec = sort(vec,'descend');
            ma = vec(1);
            R = length(vec);
            for i = 2:R
                sp_pdi = sp_pdi + (vec(1)/ma - vec(i)/ma) / (R-1);
            end
            
        end
        
        function spec = SPECIFICITY(webmatrix)

            spec = zeros(size(webmatrix,1),1);
            n_rows = size(webmatrix,1);
            for i = 1:n_rows
                spec(i) = SpeFunc.SP_PDI(webmatrix(i,:));
            end
            
        end
        
        function rr = RESOURCE_RANGE(webmatrix)
            
            row_degrees = sum(webmatrix>0,2);
            [n R] = size(webmatrix);
            rr = zeros(n,1);
            for i = 1:n
                rr(i) = (R-row_degrees(i))/(R-1);
            end
            
        end
        
        function ssi = SPECIES_SPECIFICITY_INDEX(webmatrix)
           
            [n_rows R] = size(webmatrix);
            norm_factor = R * sqrt( (R-1)/R );
            ssi = zeros(n_rows,1);
            for i = 1:n_rows
                
                tspe = 0;
                Pbar = mean(webmatrix(i,:));
                
                for j = 1:R
                    tspe = tspe + (webmatrix(i,j) - Pbar)^2;
                end
                ssi(i) = sqrt(tspe)/(norm_factor * Pbar);
            end
            
        end
        
        function [R,I] = IR(webmatrix)
           
            n_rows = size(webmatrix,1);
            %Correlation matrix of responses
            rho = corrcoef(webmatrix');
            %Standard deviations of responses
            sig = std(webmatrix,0,2);
            %Correction factor
            cfac = n_rows*(n_rows-1);
            %Measures
            tR = 0;
            tI = 0;
            for i = 1:n_rows-1
                for j = i+1:n_rows
                    tR = tR + (sig(i)-sig(j))^2;
                    tI = tI + sig(i)*sig(j)*(1-rho(i,j));
                end
            end

            R = tR / (2 * cfac);
            I = tI / cfac;
            
            VarW = std(webmatrix(:),1)*std(webmatrix(:),1);

            if (VarW == 0)
                R = 0;
                I = 0;
            else
                R = R/VarW;
                I = I/VarW;
            end
            
        end
        
    end

end