

function imc = reshape_image(imc)


rs = @(x) x(:);
imc = cellfunu(rs,imc); imc = horzcat(imc{:}); imc=double(imc);



end