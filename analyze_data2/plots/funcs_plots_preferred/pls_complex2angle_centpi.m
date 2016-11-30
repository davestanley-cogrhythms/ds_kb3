function z = pls_complex2angle_centpi(z)

    if ~isreal(z)
        z = angle(z);            % Gets the angle of the data
        z(z<0) = z(z<0) + 2*pi;     % Shift from being -pi to pi -> 0 to 2pi, since things tend to cluster around 180
    end

end