%% Position Shuffle (firing rate template shuffle) of each place cell
% Left and right direction was performed separately

% left
numberofshuffles = 1000;
allpossleft = repmat(lefttuningcurve, numberofshuffles, 1);
shiftby=randi([0 size(allpossleft,2)],size(allpossleft,1),1);
[mb nb] = size(allpossleft);
[R C]=ndgrid(1:mb,1:nb);
C=mod(bsxfun(@plus,C,shiftby(:)-1),nb)+1;
shuffledout = allpossleft;
shuffledout(:) = allpossleft(sub2ind([mb nb], R, C));
allpossleft=mat2cell(shuffledout,[ones(numberofshuffles,1)*size(lefttuningcurve,1)],size(lefttuningcurve,2));

% right...


%% Monte carlo p value calculation, 
% store weighted correlation value
wts1sig=[];
for n=1:length(realleftcorrs)
    actualvalleft=abs(realleftcorrs(n));
    actualvalright=abs(realrightcorrs(n));
    
    shuffvalleft=shuffledrvalsleft(:,n);
    shuffvalright=shuffledrvalsright(:,n);
    
    shuffvalleft=abs(shuffvalleft);
    shuffvalright=abs(shuffvalright);
    
    shuffvalleft(shuffvalleft<actualvalleft)=[];
    shuffvalright(shuffvalright<actualvalright)=[];
    
    shuffvalleft=numel(shuffvalleft);
    shuffvalright=numel(shuffvalright);
    
    pvalleft=(shuffvalleft+1)/(1000+1);
    pvalright=(shuffvalright+1)/(1000+1);
    
    wts1sig=[[wts1sig];[pvalleft,pvalright]];  % store two results in one array, left, right
end
