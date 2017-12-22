function percthr = PercolationThr(adjmat, varargin)
    %PERCOLATIONTHR Calculate the percolation threshold for an adjacency matrix.
    %   P = PERCOLATIONTHR(M)
    %   P = PERCOLATIONTHR(M, method)
    %       M is an adjacency matrix
    %       method can be either 'cutoff' or 'iterativeEdge', if left blank, function will default to 'cutoff'
    %               cuttoff: a heuristic approximation of the percolation threshold calculated to the 0.001 accuracy.
    %               iterativeEdge: the true percolation threshold calculated by iteratively removing the weakest edges
    %                              until the threshold is found. This can be very time consuming.
    
    if nargin > 2
        error('PercolationThr:TooManyInputs','Too many inputs')
    end
    switch(nargin)
        case 1
            method = 'cutoff';
        case 2
            method = char(varargin(1));
    end
    
    adjmat(isnan(adjmat)) = 0;
    switch(method)
        case 'cutoff'
            percthr = PercolationThreshold_cutoff(adjmat);
        case 'iterativeEdge'
            percthr = PercolationThreshold_iterativeEdge(adjmat);
        otherwise
            error('PercolationThr:InvalidMethod','Invalid Percolation Method. Try ''cutoff'' or ''iterativeEdge''')
    end
    disp(percthr)

end

function percthr = PercolationThreshold_cutoff(adjmat)
   for thr = 1:-0.01:0
       thrmat = threshold_proportional(adjmat, thr);
       [comps, compsize] = get_components(thrmat);
       if(length(unique(comps))>1)
          %Coarse threshold found, now compute to 0.001 accuracy
          if(thr >= 0.01)
            for newthr = (thr+0.01):-0.001:(thr-0.01)
                thrmat = threshold_proportional(adjmat, newthr);
                [comps, compsize] = get_components(thrmat);
                if(length(unique(comps))>1)
                    %Percolation threshold found
                    return
                end
                percthr = newthr;
            end
          else 
             for newthr = 0.01:-0.001:0
                thrmat = threshold_proportional(adjmat, newthr);
                [comps, compsize] = get_components(thrmat);
                if(length(unique(comps))>1)
                    %Percolation threshold found
                    return
                end
                percthr = newthr;
              end    
          end
       end
       percthr = thr;
   end
end

function percthr = PercolationThreshold_iterativeEdge(adjmat)
    edges = adj2edge(adjmat);
    edges = sortrows(edges, 3);
    for ed = 1:length(edges)
        [comps, compsize] = get_components(adjmat);
        if(length(unique(comps))>1)
            return
        end
        row = edges(ed,1);
        col = edges(ed,2);
        %percthr = edges(ed,3); %absolute value percthr
        percthr = 1 - ed./length(edges); %relative percthr
        adjmat(row,col) = 0;
    end
end
