


function [h] = imagescy(varargin)
%     Calls imagesc and then flips y-axis
%     Same inputs and outputs as imagesc.
%     CONTACT: David Stanley, Boston University (stanleyd@bu.edu, https://github.com/davestanley)
% 


    h = imagesc(varargin{:});
    set(gca,'YDir','normal');
end