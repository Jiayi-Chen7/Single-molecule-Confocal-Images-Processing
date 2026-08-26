function [r0, c0, r1, c1] = clipPoints(x, y, M, N)

m = diff(y) / diff(x);
b = y(1) - m*x(1);
if diff(x) < 0, x = fliplr(x); y = fliplr(y); end
if x(1) < 1, x(1) = 1; y(1) = m + b; end
if x(2) > N, x(2) = N; y(2) = m*N + b; end
if diff(y) < 0, x = fliplr(x); y = fliplr(y); end
if y(1) < 1, y(1) = 1; x(1) = (1-b) / m; end
if y(2) > M, y(2) = M; x(2) = (M-b) / m; end
x = round(x); y = round(y);
r0 = y(1); c0 = x(1); r1 = y(2); c1 = x(2);
end