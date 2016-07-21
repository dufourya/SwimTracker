function [wmean, wsd] = weighted(v, w)
    all_nan = isnan(v) | isnan(w);
    v(all_nan) = [];
    w(all_nan) = [];
    %http://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance
    V1 = sum(w);
    V2 = sum(w.^2);

    wmean = sum(v .* w)/V1;
    wsd = sqrt(V1/(V1.^2-V2) .* sum(w .* (v - wmean).^2));
end