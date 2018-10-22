function [data, varargout ] = extractDataInputs( data, substruct, names, vals )
% Tim Whalen June 2017
% Checks if expected fields specifiying analysis parameteers exist in data
% substruct and sets them to default otherwise
% Should be used at the start of every analysis script using the data
% struct format

% Inputs:
% data: data struct (i.e. from open_data)
% substruct: string name of substruct in which analysis parameters are specified
% names: cell array string names of expected fields in substruct
% vals: cell array of default values for each element of names

% Outputs:
% data: same as input data, with expected fields filled in as eeeded
% varargout: variable numebr of arguments, the value to set each of names
%   to (either the specified value or default)

% Example usage, where the sub-struct is called "sync"
% names = {'nstd', 'bin', 'len', 'cutoff'};
% defaults = {3, .01, 120, 10};
% [data, nstd, bin, len, cutoff] = extractDataInputs(data,'sync',names,defaults);

varargout = cell(size(names));

if isfield(data,substruct)
    for i = 1:length(names) % check each expected name in substruct
        if ~isfield(data.(substruct),names{i})
            data.(substruct).(names{i}) = vals{i}; % if none, set to default
            disp([names{i},' is undefined, setting to default ', num2str(vals{i})])
        end
        varargout{i} = data.(substruct).(names{i}); % set to val - either user-specified or now default
    end

else % no substruct
    disp([substruct, ' is undefined, setting defaults:'])
    for i = 1:length(names)
        disp([names{i}, ' = ', num2str(vals{i})])
        data.(substruct).(names{i}) = vals{i};
        varargout{i} = data.(substruct).(names{i});
    end
end

end

