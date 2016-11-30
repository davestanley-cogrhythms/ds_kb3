

function im_out = calc_similarity(im,imm,func)
    % im - image vectors dim_x_Nimages
    % imm - morphlines/prototype vectors that images will be compared; dim X Nmorphs
    % against
    
    if nargin < 3
        func = @(x,y) sum((x-y).^2);
%         func = @mycorrcoef;
    end
    
    N_images = size(im,2);
    N_bases = size(imm,2);  % Number of basis vectors

    im_out = size(N_bases,N_images);
    for j = 1:N_images
        for i = 1:N_bases
            im_out(i,j) = func(im(:,j),imm(:,i));
        end
    end
    
end

function out = mycorrcoef(x,y)

    out = corrcoef(x,y);
    out=out(1,2);
end