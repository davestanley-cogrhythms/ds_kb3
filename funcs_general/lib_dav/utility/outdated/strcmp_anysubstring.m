function varargout = strcmp_anysubstring(varargin)
    % Outdated; Alias for strcmp_substr

    [varargout{1:nargout}] = strcmp_substr(varargin{:});
end