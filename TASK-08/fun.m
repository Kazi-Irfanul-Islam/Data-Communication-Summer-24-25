function dn = fun(txt)
    dec = double(txt);              % Text to ASCII (decimal)
    p2  = 2.^(0:-1:-7);             % 2^0, 2^-1, ..., 2^-7
    B   = mod(floor(p2' * dec), 2); % Decimal -> binary bits (columns per char)
    dn  = reshape(B, 1, numel(B));  % Bytes to serial conversion
end
