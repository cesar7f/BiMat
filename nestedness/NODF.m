function [nodf,nodf_rows,nodf_cols] = NODF(matrix, nodf_strict)
% NODF - Calculate normalized nodf value of a matrix. To know how this
% nestedness metric works you can consult the following paper:
%
%   Almeida-Neto, Mario and Guimaraes, Paulo and Guimaraes, Paulo R and 
%   Loyola, Rafael D and Ulrich, Werner. A consistent metric for nestedness
%   analysis in ecological systems: reconciling concept and measurement.
%   Oikos 2008
%
%   [nodf nodf_rows nodf_cols] = NODF(matrix) - Calculate the normalized nodf
%   value for rows+columns (nodf), rows (nodf_rows), and columns (nodf_cols)
%   of matrix forcing an strict decreasing filling.
%
%   [nodf nodf_rows nodf_cols] = NODF(matrix,nodf_strict) - Calculate the normalized nodf
%   value for rows+columns (nodf), rows (nodf_rows), and columns (nodf_cols)
%   of matrix forcing an strict decreasing filling nodf_stric=1 or counting
%   equal degree rows/columns as nested nodf_strict=0.
%
% See also:
%   NestednessBINMATNEST
    if(isempty(matrix))
        if(nargin == 1)
            nodf.nodf = NaN;
            nodf.nodf_rows = NaN;
            nodf.nodf_cols = NaN;
        end
        if(nargout==2)
            nodf_rows = NaN;
            nodf = NaN;
        elseif(nargout == 3)
            nodf_rows = NaN;
            nodf_cols = NaN;
            nodf = NaN;
        end
        return;
    end

    if(nargin>=2 && nodf_strict==0)

        [n_rows, n_cols] = size(matrix);
        nest_rows = 0;
        nest_cols = 0;
        deg_rows = sum(matrix,2);
        deg_cols = sum(matrix,1);
        for ii = 1:n_rows
            for jj = ii+1:n_rows
                if(deg_rows(ii)==0 || deg_rows(jj)==0); continue; end;
                nest_rows = nest_rows + sum(matrix(ii,:).*matrix(jj,:))/min(deg_rows(ii),deg_rows(jj));
            end
        end
        for ii = 1:n_cols
            for jj = ii+1:n_cols
                if(deg_cols(ii)==0 || deg_cols(jj)==0); continue; end;
                nest_cols = nest_cols + sum(matrix(:,ii).*matrix(:,jj))/min(deg_cols(ii),deg_cols(jj));
            end
        end

        nest_rows(isnan(nest_rows)) = 0;
        nest_cols(isnan(nest_cols)) = 0;

        denom = max(n_rows*(n_rows-1)/2 + n_cols*(n_cols-1)/2,1);
        nodf.nodf_rows = 1.0*nest_rows/max(n_rows*(n_rows-1)/2,1);
        nodf.nodf_cols = 1.0*nest_cols/max(n_cols*(n_cols-1)/2,1);   
        nodf.nodf = (nest_rows+nest_cols)/denom;

    else
        if(sum(sum(matrix)) == numel(matrix))
            nr = 0;
            nc = 0;
            nvalue = 0;
        else
            [n_rows, n_cols] = size(matrix);

            C1 = (1.0*matrix)*matrix';
            D1 = repmat(diag(C1),1,n_rows);
            nr = (D1~=D1');%
            nr = C1.*nr;
            nr = nr./(min(D1,D1'));
            %nr = (C1-eye(size(C1)))./(min(D1,D1'));

            C2 = (1.0*matrix')*matrix;
            D2 = repmat(diag(C2),1,n_cols);
            nc = C2.*(D2~=D2');
            nc = nc./(min(D2,D2'));

            nr = nansum(nr(:))/2;
            nc = nansum(nc(:))/2;
            nvalue = nr+nc;

            denom = max(n_rows*(n_rows-1)/2 + n_cols*(n_cols-1)/2,1);
            nr = nr/max(n_rows*(n_rows-1)/2,1);
            nc = nc/max(n_cols*(n_cols-1)/2,1);   
            nvalue = nvalue/denom;
        end

        nodf.nodf_rows = nr;
        nodf.nodf_cols = nc;
        nodf.nodf = nvalue;

    end


    if(nargout==2)
        nodf_rows = nodf.nodf_rows;
        nodf = nodf.nodf;
    elseif(nargout == 3)
        nodf_rows = nodf.nodf_rows;
        nodf_cols = nodf.nodf_cols;
        nodf = nodf.nodf;
    end
end