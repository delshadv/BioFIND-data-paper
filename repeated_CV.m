function accuracy = repeated_CV(A,CVratio,kFolds,Nrun)
%
%
% This function calculates cross-validation
% using Support Vector Machine with linear kernel
%
%
% A : input matrix [m n+1]; m: subjects (cases), n: features
% last column must be targets (Controls=1/MCI=2)
% kFolds : Number of folds
% accuracy: cross-validation accuracy
%
% Henson R.N 2020, Vaghari D 2020
% -------------------------------------------------------------------------
allData = A(:,1:end-1);
targets = A(:,end);

for c = 1:2
    N(c) = length(find(targets==c));
    Ntrain(c)  = floor(CVratio(1)*N(c)); % was "ceil" - needs to be "floor"
    Ntest(c) = N(c) - Ntrain(c);
end

% Ensure trained on equal number of each category
balN = min(Ntrain); fldN = balN*CVratio(2);
for c = 1:2
    if Ntrain(c)>balN
        Ntest(c) = Ntest(c)+(Ntrain(c)-balN);
        Ntrain(c) = balN;
    end
end
accuracy = [];

parfor rr = 1:Nrun
    
    rp = {};
    for c = 1:2
        rp{c} = randperm(N(c));
    end
    %tr = cell(1,5); te=tr;
   Score = zeros(1, kFolds);
    
    for k = 1:kFolds
        
        xtrain = []; xtest = []; ytrain = []; ytest = [];
        for c = 1:2
            ii = setdiff([1:(balN+fldN)],[1:fldN]+(k-1)*fldN);
            ri = rp{c}(ii);
            rj = setdiff([1:N(c)],ri);
            
            ii = find(targets==c);
            %tr{f} = [tr{f}; ii(ri)]; te{f} = [te{f}; ii(rj)];
            
            xtrain = [xtrain; allData(ii(ri),:)];
            xtest  = [xtest;  allData(ii(rj),:)];
            ytrain = [ytrain; targets(ii(ri))];
            ytest  = [ytest;  targets(ii(rj))];
        end
        
        
        C = 1; %literature
        % we need to retrain here and save the SVM for the best C
        SVMModel = fitcsvm(xtrain, ytrain, ...
            'KernelFunction', 'linear', 'KernelScale', 'auto', ...
            'BoxConstraint', C,'Standardize', true,'solver','ISDA');
%       SVMModel = fitclinear(xtrain, ytrain)
        
        
        Score(k) = loss(SVMModel,xtest, ytest); % default: 'LossFun','classiferror'
        
    end
    accuracy(rr) = (mean(1-Score))*100;
    
    
end


