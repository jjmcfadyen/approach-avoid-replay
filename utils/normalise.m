function Y = normalise(X,thismin,thismax)

thismin = min(X(:));
thismax = max(X(:));

Y = (X - thismin) / ( thismax - thismin );

end