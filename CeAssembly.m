% Cell assemblies were detected from left and right direction separately

%% left direction
leftcorrelationMatrix = (leftspikeTrains * leftspikeTrains')./nWindowsTotal;
 
leftpfCorrMatrix = zeros(size(leftcorrelationMatrix));
leftassemblyMembership = zeros(size(leftcorrelationMatrix));
 
[prncplComps,~,vars] = pca(leftspikeTrains');
minVar = (1 + sqrt(nCells/nWindowsTotal))^2;
nPatterns = length(vars(vars > minVar));

prncplComps = prncplComps(:, 1:nPatterns);

pcProj = prncplComps' * leftspikeTrains;
[icasig, mixMatrix, unmixMatrix] = fastica(pcProj);
assmblPtrnWgts = prncplComps * unmixMatrix';% each cell's WEIGHT within each assembly
    
% define if it is a cell assembly
cellID=[1:length(lineartrackcells)]';
leftassmblPtrnCellIDs=cell(1,nPatterns);
for p = 1:nPatterns
    assmblPtrnWgts(:,p) = assmblPtrnWgts(:,p)./norm(assmblPtrnWgts(:,p));
    if max(assmblPtrnWgts(:,p)) < abs(min(assmblPtrnWgts(:,p)))
        assmblPtrnWgts(:,p) = -assmblPtrnWgts(:,p);
    end
    memberInds = find(assmblPtrnWgts(:,p) > ...
        (mean(assmblPtrnWgts(:,p)) + 2*std(assmblPtrnWgts(:,p))));
    leftassmblPtrnCellIDs{p} = cellID(memberInds, :); 
    if length(leftassmblPtrnCellIDs{p})>1
    pairCombs = nchoosek(memberInds, 2);
    linearInd = sub2ind(size(leftassemblyMembership), pairCombs(:,1), pairCombs(:,2));
    leftspikeTrains(linearInd) = leftassemblyMembership(linearInd(:)) + 1;
    end

    if length(leftassmblPtrnCellIDs{p})<2
    leftassmblPtrnCellIDs{p}=[];
    end
end
 
assmblPtrnWgts_left=assmblPtrnWgts;  
  
eachcell_contri_left=[];
for n=1:size(assmblPtrnWgts_left,1)
 temp=mean(assmblPtrnWgts_left(n,:).^2);
 eachcell_contri_left=[eachcell_contri_left;temp];
end

leftassmblPtrnCellIDs=leftassmblPtrnCellIDs(~cellfun('isempty',leftassmblPtrnCellIDs));      

leftprojmtrx={};
for m=1:size(assmblPtrnWgts,2)
    tempvar=assmblPtrnWgts(:,m);
    tempvar=kron(tempvar,tempvar');
    tempvar(logical(eye(size(tempvar)))) = 0;
    leftprojmtrx{m}=tempvar;
end



%% right direction
rightcorrelationMatrix = (rightspikeTrains * rightspikeTrains')./nWindowsTotal;
rightpfCorrMatrix = zeros(size(rightcorrelationMatrix));
rightassemblyMembership = zeros(size(rightcorrelationMatrix));


[prncplComps,~,vars] = pca(rightspikeTrains');
minVar = (1 + sqrt(nCells/nWindowsTotal))^2;
nPatterns = length(vars(vars > minVar));
prncplComps = prncplComps(:, 1:nPatterns);
pcProj = prncplComps' * rightspikeTrains;
[icasig, mixMatrix, unmixMatrix] = fastica(pcProj);
assmblPtrnWgts = prncplComps * unmixMatrix';

cellID=[1:length(lineartrackcells)]';
rightassmblPtrnCellIDs=cell(1,nPatterns);
    for p = 1:nPatterns
       assmblPtrnWgts(:,p) = assmblPtrnWgts(:,p)./norm(assmblPtrnWgts(:,p));
        if max(assmblPtrnWgts(:,p)) < abs(min(assmblPtrnWgts(:,p)))
            assmblPtrnWgts(:,p) = -assmblPtrnWgts(:,p);
        end
        memberInds = find(assmblPtrnWgts(:,p) > ...
            (mean(assmblPtrnWgts(:,p)) + 2*std(assmblPtrnWgts(:,p))));
        rightassmblPtrnCellIDs{p} = cellID(memberInds, :);
        if length(rightassmblPtrnCellIDs{p})>1
        pairCombs = nchoosek(memberInds, 2);
        linearInd = sub2ind(size(rightassemblyMembership), pairCombs(:,1), pairCombs(:,2));
        rightspikeTrains(linearInd) = rightassemblyMembership(linearInd(:)) + 1;
        end
        
        if length(rightassmblPtrnCellIDs{p})<2
        rightassmblPtrnCellIDs{p}=[];
        end
    end
 
assmblPtrnWgts_right=assmblPtrnWgts;  

eachcell_contri_right=[];
for n=1:size(assmblPtrnWgts_right,1)
 temp=mean(assmblPtrnWgts_right(n,:).^2);
 eachcell_contri_right=[eachcell_contri_right;temp];
end

rightassmblPtrnCellIDs=rightassmblPtrnCellIDs(~cellfun('isempty',rightassmblPtrnCellIDs));       
    
rightprojmtrx={};
for m=1:size(assmblPtrnWgts,2)
    tempvar=assmblPtrnWgts(:,m);
    tempvar=kron(tempvar,tempvar');
    tempvar(logical(eye(size(tempvar)))) = 0;
    rightprojmtrx{m}=tempvar;
end    