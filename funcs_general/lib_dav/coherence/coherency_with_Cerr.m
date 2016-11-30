

function [varargout] = coherency_with_Cerr(varargin)
%     Usage: [varargout] = coherency_with_Cerr(fname,varargin)
%     Implements [varargout] = fname(varargin). If this function
%     fails, assumes it is because of the Cerr output being requested.
%     In this case, calls fname again only without requesting the final
%     argument (Cerr) and instead generates some fake data for Cerr.
%     Forms:
%         [varargout] = coherency_with_Cerr(fname,varargin)
%     Input:
%         fname - function handle
%         varargin - inputs to fname
%         varargout - outputs requested from fname. If fails, will contain a
%             final entry for Cerr

    fname = varargin{1};
    args = varargin(2:end);
    try
        [varargout{1:nargout}] = fname(args{:});
    catch
        [varargout{1:nargout-1}] = fname(args{:});
        Cerr1_fake =  zeros(1,length(varargout{1}(:)));
        Cerr2_fake =  ones(1,length(varargout{1}(:)));          % Needs to match the dimensionality of Cave for cases where we do skip_Jackknife_for_ctg_alltrials17
        varargout{nargout} = [Cerr1_fake; Cerr2_fake];
    end


end