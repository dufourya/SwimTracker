function C = replicate(A, B)
    %replicate.m -- replicate each A(i) B(i) times.
    % http://stackoverflow.com/a/2382382
    is_nan = isnan(A) | isnan(B);
    A(is_nan) = [];
    B(is_nan) = [];
    index = zeros(1,sum(B));
    index(cumsum([1 B(1:end-1)])) = 1;
    index = cumsum(index);

    C = A(index);
end

