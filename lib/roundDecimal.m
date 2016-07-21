function d = roundDecimal(n, digits)
    d = round(n.*10.^digits)./10^digits;
end