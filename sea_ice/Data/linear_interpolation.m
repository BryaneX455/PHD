load('ProcessedData.mat')
tMin = tMinutesArray(1);
tMax = tMinutesArray(end);
deltaT = 5;
tInterpArray = tMin:deltaT:tMax;
nTime = length(tInterpArray);
obsTimesIdx = round(tMinutesArray/5)+1;


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
    xInterpValues = interp1(tData, xData, tInterpValues, "linear");
    yInterpValues = interp1(tData, yData, tInterpValues, "linear");
    xInterpArray(iFloe, interpIdx) = xInterpValues;
    yInterpArray(iFloe, interpIdx) = yInterpValues;   
end

save('LinearInterpData', "tInterpArray", "xInterpArray", "yInterpArray", "nTime", "obsTimesIdx")