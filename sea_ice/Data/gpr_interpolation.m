load('ProcessedData.mat')
tMin = tMinutesArray(1);
tMax = tMinutesArray(end);
deltaT = 5;
tInterpArray = tMin:deltaT:tMax;
nTime = length(tInterpArray);
obsTimesIdx = round(tMinutesArray/5)+1;

xSigmaArray = nan(nFloes);
ySigmaArray = nan(nFloes);

xInterpArray = nan(nFloes, nTime);
yInterpArray = nan(nFloes, nTime);
for iFloe = 1:nFloes
    tData = tMinutesArray;
    xData = xMetersArray(iFloe, :);
    yData = yMetersArray(iFloe, :);
    hasObsIdx = find(xData~=0);
    tData = tData(hasObsIdx);
    xData = xData(hasObsIdx);
    yData = yData(hasObsIdx);
    interpIdx = obsTimesIdx(hasObsIdx(1)):obsTimesIdx(hasObsIdx(end));
    tInterpValues = tInterpArray(interpIdx);

    % gprMdlX = fitrgp(tData.', xData.', "KernelFunction", "squaredexponential", "KernelParameters", [100e3, 5*60*24] , "Sigma",250, 'ConstantSigma', true, 'SigmaLowerBound', 1);
    % gprMdlY = fitrgp(tData.', yData.', "KernelFunction", "squaredexponential", "KernelParameters", [100e3, 5*60*24] , "Sigma",250, 'ConstantSigma', true, 'SigmaLowerBound', 1);   

    gprMdlX = fitrgp(tData.', xData.');
    gprMdlY = fitrgp(tData.', yData.');   


    xSigmaArray(iFloe) = gprMdlX.Sigma;
    ySigmaArray(iFloe) = gprMdlY.Sigma;

    xInterpValues = predict(gprMdlX,tInterpValues.').';
    yInterpValues = predict(gprMdlY,tInterpValues.').';
    xInterpArray(iFloe, interpIdx) = xInterpValues;
    yInterpArray(iFloe, interpIdx) = yInterpValues;   
end

save('GprInterpData', "tInterpArray", "xInterpArray", "yInterpArray", "nTime", "obsTimesIdx", 'xSigmaArray', 'ySigmaArray')