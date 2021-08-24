function [ post ] = decode_calcPosterior( rm, spkT, tBinSz, event_length )
% ARGS
% rm            ratemaps [nPosBin x nCells]. The smoothed expected firing
%               rate (Hz) for every cell in every spatial bin [spatial bins
%               in y-axis and cells on x-axis)
%
% spkT          spike times [cell array vector nCells - each cell has all 
%               spike times for that cell]. Spike times are just for the 
%               period to be decoded (e.g. event). In ms
%
% tBinSz        temporal bin size [scalar] in ms [typically around 20ms]

%--- MAIN -----------------------------------------------------------------

% ---- HOUSE KEEPING
nPBin           =size(rm,1); %N positional bin
nCell           =size(rm,2); %N cells in total


% ---- 1. Bin spikes and determine how many were fired in each temporal bin
spkCells        =find(~cellfun(@isempty, spkT)); %just cells that spk
spkT2           =spkT(spkCells); %just cells that spk
fstSpk          =min(cellfun(@min, spkT2));
lstSpk          =max(cellfun(@max, spkT2));

nTBin           =ceil((lstSpk-fstSpk)/tBinSz) +1; %N time bins needed
binEdge         =fstSpk + (0:tBinSz:nTBin*tBinSz);
binEdge         =binEdge(1:end-1); % previous step adds a time bin
clear spkT2 lstSpk fstSpk 

%Create and populate a nCell (just spiking ones) x nTBin mat with number of
%spikes emited by each cell in that time bin
spkPerBin       =zeros(1,nTBin, nCell); %Preallocate
for nn          =1:length(spkCells) %just loop over cells that fired
    cellInd     =spkCells(nn); %current cell
    spkPerBin(1,:,cellInd)  =histc(spkT{cellInd}, binEdge); %[1 x nTBin x nCell]
end
nSpkPerTBin     =squeeze(sum(spkPerBin,3)); %[nTBin x 1] number of spikes in tBin
clear  nn binEdge

% ---- 2. Now start to calculate the posterior code 
rm              =rm+ (eps.^8); %Add a small number so there are no zeros
expecSpk        =rm.*tBinSz./1000; %[nPos x nCells] Expected number of spikes per bin
expecSpk        =reshape(expecSpk, [nPBin,1, nCell]); %nPos x 1 x nCell]


expon           =exp(-expecSpk); %Exponent of equation.
factSpkPerBin   =factorial(spkPerBin); %Factorial to divide by

wrking          =bsxfun(@power, expecSpk, spkPerBin); %[nPos x nTbin x nCell]
wrking          =bsxfun(@rdivide, wrking, factSpkPerBin); %[nPos x nTbin x nCell]
wrking          =bsxfun(@times,wrking, expon); %[nPos x nTbin x nCell]
post            =prod(wrking,3); %Non normalised prob [nPos x Tbin]
clear wrking factSpkPerBin expon expecSpk


% ---- 3. Do some normalisation and mark timebins that had no spikes and so
% shouldn't be decoded
%Normalise each column i.e. time bin to sum to 1
post            =bsxfun(@rdivide, post, sum(post)); %Each column sums to 1



%Set time bins with no spikes to be all nan - don't try and use these to
%decode.
% post(:,nSpkPerTBin==0)  =NaN;
post(:,nSpkPerTBin==0)  =0; 




end

