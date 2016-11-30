function z = pls_complex2angle_cent0(z)

    if ~isreal(z)
        z = angle(z);            % Gets the angle of the data
        z(z>pi) = z(z>pi) - 2*pi;     % Shift from being 0 to 2pi -> -pi to pi
        %z = abs(z);
    end

end