function weightedCorr = ReplayQualityNormalized(probDecode, decodedPos)
%%

decodedPosTemp = decodedPos(decodedPos ~= 0);

maxPos = max(decodedPosTemp);
minPos = min(decodedPosTemp);

probDecode = probDecode(minPos:maxPos,:);
[nBins, nTimeBins] = size(probDecode);
binVector = 1:nBins;

% if maxPos ~= minPos
%     binVectorTemp = interp1([minPos maxPos], [1 100], minPos:maxPos);
%     diffBinVecTemp = diff(binVectorTemp);
%     diffBinVecTemp = diffBinVecTemp(1);
%     nBinsUsed = maxPos - minPos + 1;
%     nBinsStart = minPos - 1;
%     nBinsEnd = nBins - maxPos;
%     if nBinsStart ~= 0
%         binVectorStart = 1 - (nBinsStart*diffBinVecTemp):diffBinVecTemp:1 - diffBinVecTemp;
%     else
%         binVectorStart = [];
%     end
%     if nBinsEnd ~=0
%         binVectorEnd = 100 + diffBinVecTemp:diffBinVecTemp:100 + nBinsEnd*diffBinVecTemp;
%     else
%         binVectorEnd = [];
%     end
%     
%     binVector = [binVectorStart, binVectorTemp, binVectorEnd];
% else
%     binVector = 1:nBins;
% end

numTerm = bsxfun(@times, probDecode, binVector');
numTerm = sum(numTerm(:));
denTerm = sum(probDecode(:));
weightedMeanPos = numTerm/denTerm;

numTerm = bsxfun(@times, probDecode, 1:nTimeBins);
numTerm = sum(numTerm(:));
weightedMeanTime = numTerm/denTerm;

weightedCovNum = bsxfun(@times, probDecode, (binVector - weightedMeanPos)');
weightedCovNum = bsxfun(@times, weightedCovNum, (1:nTimeBins) - weightedMeanTime);
weightedCov = sum(weightedCovNum(:))/denTerm;

numTerm = bsxfun(@times, probDecode, ((binVector - weightedMeanPos).^2)');
posCov = sum(numTerm(:))/denTerm;

numTerm = bsxfun(@times, probDecode, (((1:nTimeBins) - weightedMeanTime).^2));
timeCov = sum(numTerm(:))/denTerm;

weightedCorr = weightedCov/(sqrt(posCov * timeCov));