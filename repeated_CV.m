function accuracy = repeated_CV_matlab(A,kFolds,Nrun)

% This function calculates cross-validation accuracy
% using Support Vector Machines with linear kernel.
% The data is randomly split into "kFolds" partitions several times ("Nrun" times).
%
%
% A : input matrix [m n+1]; m: subjects (cases/observations), n: features
% last column must be targets (for e.g., Controls=1/MCI=2)
% kFolds : Number of folds
% Nrun : Number of runs to reapeate CV
% accuracy: holds cross-validation accuracies for each run
%
% Henson R.N 2020, Vaghari D 2020
% -------------------------------------------------------------------------

allData = A(:,1:end-1);
targets = A(:,end);

parfor rr = 1:Nrun
    
    kIdx = crossvalind('Kfold', length(targets), kFolds);
    acc = zeros(1, kFolds);
    
    for k = 1:kFolds
        
        xtrain = allData(kIdx~=k, :);
        ytrain = targets(kIdx~=k);
        xtest = allData(kIdx==k, :);
        ytest = targets(kIdx==k);
        
          SVMmodel = fitcsvm(xtrain, ytrain, ...
              'KernelFunction', 'linear', 'KernelScale', 'auto', ...
              'BoxConstraint', 1,'Standardize', true,'solver','ISDA');
%         SVMmodel = fitclinear(xtrain, ytrain,'Solver',{'dual','bfgs'});
%        Score(k) = loss(SVMmodel,xtest, ytest); % default: 'LossFun','classiferror'
        ts_pred = predict(SVMmodel,xtest);
        acc(k) = ((sum(ts_pred==1 & ytest==1)/sum(ytest==1))+...
                    (sum(ts_pred==2 & ytest==2)/sum(ytest==2)))/2*100;
        
    end
    
    accuracy(rr) = (mean(acc));
end
