function accuracy = nestedCVSP (A,kFolds)

% A multi-site, multi-participant magnetoencephalography resting-state dataset
% to study dementia: The BioFIND dataset
%
% This function calculates nested cross-validation
% using a well-designed Support Vector Machine with linear kernel
% which in inner loop finds best C hyper-parameter
%
% A : input matrix [m n+1]; m: subjects (cases), n: features
% last column must be targets (Controls=1/MCI=2)
% kFolds : Number of folds
% accuracy: nested cross-validation accuracy
%
% Henson R.N 2020, Vaghari D 2020
% -------------------------------------------------------------------------

if nargin < 2
    kFolds = 10;
end

allData = A(:,1:end-1);
targets = A(:,end);
    
%CVratio=[(kFolds-1)/kFolds 1-((kFolds-1)/kFolds)];

%parfor r = 1:100
kIdx = crossvalind('Kfold', length(targets), kFolds);
bestScore = zeros(1, kFolds);

parfor k = 1:kFolds
    
    xtrain = allData(kIdx~=k, :);
    ytrain = targets(kIdx~=k);
    xtest = allData(kIdx==k, :);
    ytest = targets(kIdx==k);
    
    % Best Hyper-parameter
    bestCScore = inf;
    bestC = NaN;
    %gridC = 2.^(-10:0.1:15);
    gridC = [0.01 0.1 1 2 5 10 100];
    
    for C = gridC
        % cross validation for parameter C
        rng(1) % For reproducibility
        kIdxC = crossvalind('Kfold', length(ytrain), kFolds);
        L = zeros(1, kFolds);
        
        for kC = 1:kFolds
            
            xtrainC = xtrain(kIdxC~=kC, :);
            ytrainC = ytrain(kIdxC~=kC);
            xtestC = xtrain(kIdxC==kC, :);
            ytestC = ytrain(kIdxC==kC);
            
            anSVMModel = fitcsvm(xtrainC, ytrainC, ...
                'KernelFunction', 'linear', 'KernelScale', 'auto', ...
                'BoxConstraint', C,'Standardize', true,'solver','SMO',...
                'PolynomialOrder',[] );
            L(kC) = loss(anSVMModel,xtestC, ytestC);
%             [y_test_predicted] = predict(anSVMModel,xtestC); % y_test_predicted: predicted class labels
%             L(kC) = sum(y_test_predicted == ytestC)/length(ytestC);  % Note unbalanced, but assume equal 1,2
            
        end
        L = mean(L);
        if L < bestCScore
            bestCScore = L;
            bestC = C;
        end
    end
    % we need to retrain here and save the SVM for the best C
    bestCSVM = fitcsvm(xtrain, ytrain,'KernelScale', 'auto', ...
        'KernelFunction', 'linear','Standardize', true,  ...
        'BoxConstraint',bestC,'solver','SMO',...
        'PolynomialOrder',[]);
     bestScore(k) = loss(bestCSVM,xtest, ytest); % default: 'LossFun','classiferror'
    K(k)=bestC;
%     [y_test_predicted] = predict(bestCSVM,xtest); % y_test_predicted: predicted class labels
%     acc(k) = sum(y_test_predicted == ytest)/length(ytest);  % Note unbalanced, but assume equal 1,2
%     [y_test_predicted] = predict(bestCSVM,xtrain); % y_test_predicted: predicted class labels
%     acc1(k) = sum(y_test_predicted == ytrain)/length(ytrain);  % Note unbalanced, but assume equal 1,2
    
end

accuracy = mean(1-bestScore);

