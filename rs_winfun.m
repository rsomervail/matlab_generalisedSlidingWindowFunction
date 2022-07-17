%   
%   Computes any sliding window with a chosen input function
% 
%   Usage:    out = rs_winfun( data, cfg)
%   
%   data - first dimension is the one you slide through (usually time)
% 
%   cfg - structure with the following fields: 
% 
%        ** essential **
%          fs - sampling rate of data (Hz)
%          win - sliding window dimensions in seconds: [timeBefore timeAfter]
%          func - handle of function to be applied, e.g. @mean
% 
%        ** optional **
%          func_params - parameters for the function as a cell array
%          toi - sample numbers to compute value for in first dimension of data (if empty computes for all timepoints)
%          edge_mode - 'exclude' - ignore edge values which do not have a complete window (default)
%                      'include' - compute these values anyway with limited values in the window
%          verbosity - 'quiet'   - use evalc to supress command window output of the function
%                      'verbose' - just use feval (maybe more efficient)
% 
%% 
function out = rs_winfun( data, cfg)


if ~isfield(cfg, 'toi')
    cfg.toi = 1:size(data,1);
    warning 'rs_winfun: no timepoints of interest specified; using all timepoints of first dimension of data'
end
if ~isfield(cfg, 'edge_mode')
    cfg.edge_mode = 'exclude';
    warning 'rs_winfun: no edge_mode specified, using default mode of excluding timepoints in toi with incomplete windows'
end
if ~isfield(cfg, 'verbosity')
    cfg.verbosity = 'verbose';
end
if ~isfield(cfg, 'win')
    error 'rs_winfun: you need specify a window size in seconds:  [ pre_seconds , post_seconds ]'
end
if ~isfield(cfg, 'func')
    error 'rs_winfun: not function handle specified (cfg.func); need to choose a function to apply';
end
if ~isfield(cfg, 'func_params')
    cfg.func_params = [];
end
if ~isa( cfg.func, 'function_handle' )
    error 'rs_winfun: not function handle specified (cfg.func); need to choose a function to apply';
end

% compute window samples
win = round( cfg.win * cfg.fs );
if win(1) < 0 
    win(1) = -1*(win(1));
end

% get dimensionality of input data (for later indexing)
input_dim = ndims(data);
input_size = size(data);
input_size = input_size(2:end);
% make comma-separated list index of other input dimensions
for k = 1:input_dim-1
    input_inds{k} = 1:input_size(k);
end

% loop through timepoints of interest and compute sliding window function
for t = 1 :length( cfg.toi )

    % get timepoints in the window centered on this timepoint
    tempwin = [cfg.toi(t)-win(1) : cfg.toi(t)+win(2)];
    switch cfg.edge_mode
        case 'exclude'
            if any( tempwin <= 0 ) || any( tempwin > size(data,1) )
                skip_flag = true; % skip this iteration
            else 
                skip_flag = false; % don't skip
            end
        case 'include'
            skip_flag = false; % don't skip
            tempwin( tempwin <= 0 | tempwin > size(data,1) ) = []; % compute anyway with limited timepoints in this window
    end

    if ~skip_flag

        % make comma-separated lists index from the temporary window indices
        windy_in = [ {tempwin}, input_inds ];
    
        % run function on window
        switch cfg.verbosity
            case 'verbose'
                if ~isempty(cfg.func_params)
                    tempout = feval( cfg.func, data(windy_in{:}), cfg.func_params{:}  );
                else
                    tempout = feval( cfg.func, data(windy_in{:}) );
                end
            case 'quiet'
                if ~isempty(cfg.func_params)
                    evalc( 'tempout = feval( cfg.func, data(windy_in{:}), cfg.func_params{:}  ); ' )
                else
                    evalc( 'tempout = feval( cfg.func, data(windy_in{:}) );' )
                end
        end
    
        % if first timepoint, get dimensionality of output
        if ~exist('out','var')
            output_dims = [length(cfg.toi), size(tempout) ];
            output_dims(output_dims==1) = [];
            if length(output_dims)==1, output_dims(2) = 1; end
            out = nan(  output_dims );
    
            % make initial output index
            output_size = size(tempout); 
            output_size = output_size(2:end); 
            if any(output_size>1) % if multiple dimensions to account for then make comma-seperated list of these other dims
                for k = 1:length(output_size)
                    output_inds{k} = 1:output_size(k);
                end
            end
        end
    
        % store temporary out in main output matrix
        if exist('output_inds','var') % if multiple dimensions to account for
            windy_out = [ {t} , output_inds{:}   ]; % make comma-seperated list of indices for output
            out(windy_out{:}) = tempout;
        else
            out(t) = tempout;
        end
    
    end % end if ~ skip_flag

    % display output
    if t>1
        fprintf(repmat(char(8), [1 numel(msg)])); % clear previous line
    end
    msg = sprintf( 'rs_winfun: completed t = %d/%d\n' ,t, length( cfg.toi ) );
    fprintf(msg) 
    
end % end toi loop
fprintf('rs_winfun: finished!')


end