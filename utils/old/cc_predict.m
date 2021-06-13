function [pred] = cc_predict(X,classifier,useLogistic,useIntercept)

B = classifier.B';  % betas
I = classifier.I;  % intercept

% Scale
X = X ./ prctile(abs(X(:)),95);

% Predict
if useLogistic
    if useIntercept
        pred = 1 ./ (1 + exp(-(X*B + repmat(I, [size(X,1) 1]))));
    else
        pred = 1 ./ (1 + exp(-(X*B)));
    end
else
    if useIntercept
        pred = X*B + I;
    else
        pred = X*B;
    end
end

end