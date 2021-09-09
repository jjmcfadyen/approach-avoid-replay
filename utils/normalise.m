function Y = normalise(X,thismin,thismax)

Y = (X - min(X(:))) / ( max(X(:)) - min(X(:)) );
if nargin>1
   thislim = thismax-thismin;
   Y = Y*thislim + thismin;
end

end