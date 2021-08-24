% function [CPH, TimeLags] = CondProbHist(EVr, EVt, bin_sz_sec, max_lag_sec, norm_t)
% Calculate what neuroscientists call 'cross(auto)-correlogram', which
% strictly speaking is a "conditional probability histogram", where each
% bin shows a conditional probability of an event (spike, ripple, etc.) 
% at time t0+t on the condition that there is a reference event at time t0.
% Input arguments:
% EVr - COLUMN vector of reference event timestamps. Must be sorted. 
% EVt - COLUMN vector of tested event timestamps. Must be sorted. 
% If EVr and EVt are the same then auto-correllogram will be calculated.
% bin_sz_sec - scalar, gives the size of a bin in second.
% max_lag_sec - scalar, output time lag value
% will be in units of EVr (or EVt) and (+/-)max_lag_sec in size.
% norm_t - string, one of 'cntbin', 'evtsec' or 'prob' describe
% normalization method used for output histogram values.
% So for example if your spike timestamps in microseconds and you want the
% histogram to be +/- 0.1 sec width and 1 ms bin size, then use:
% >> CondProbHist(EVr/1e6, EVt/1e6, 0.001, 0.1, 'evtsec');
% Output arguments:
% CPH - COLUMN vector, histogram in units according to norm_t.
% TimeLags - COLUMN vector, same length as CPH, - time lag values, this
% is useful for plotting, i.e.: bar(TimeLags, CPH); axis tight;
% If no output arguments specified, then output histogram will be plotted.
% WARNING: Note, that because improper choise of the bin size will lead to
% wrong calculation of probability values this function will warn you if
% bin size provided as input is too big (for norm_t == 'prob').

% *  Copyright (C) 2012,2013 Denis Polygalov,
% *  Laboratory for Neural Circuit and Behavioral Physiology,
% *  RIKEN Brain Science Institute, Saitama, Japan.
% *
% *  This program is free software; you can redistribute it and/or modify
% *  it under the terms of the GNU General Public License as published by
% *  the Free Software Foundation; either version 3 of the License, or
% *  (at your option) any later version.
% *
% *  This program is distributed in the hope that it will be useful,
% *  but WITHOUT ANY WARRANTY; without even the implied warranty of
% *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% *  GNU General Public License for more details.
% *
% *  You should have received a copy of the GNU General Public License
% *  along with this program; if not, a copy is available at
% *  http://www.fsf.org/

function [CPH, TimeLags] = CondProbHist(EVr, EVt, bin_sz_sec, max_lag_sec, norm_t)

if ~isvector(EVr)
   error('CondProbHist: EVr must be a vector');
end

if ~isvector(EVt)
   error('CondProbHist: EVt must be a vector');
end

if isvector(EVr) && size(EVr,1) == 1 % row vector is given
   EVr = EVr';
   fprintf('CondProbHist: row vector supplied as EVr. Check your code.\n');
end

if isvector(EVt) && size(EVt,1) == 1 % row vector is given
   EVt = EVt';
   fprintf('CondProbHist: row vector supplied as EVt. Check your code.\n');
end

if ~isempty(find(diff(EVr) < 0, 1))
   error('CondProbHist: EVr must be a vector of monotonically increasing values');
end

if ~isempty(find(diff(EVt) < 0, 1))
   error('CondProbHist: EVt must be a vector of monotonically increasing values');
end

if     strcmp(norm_t, 'cntbin')
   norm_t_code = 1;
elseif strcmp(norm_t, 'evtsec')
   norm_t_code = 2;
elseif strcmp(norm_t, 'prob')
   norm_t_code = 3;
else
   error('CondProbHist: unknown normalization method requested.');
end

nbins_half = floor(max_lag_sec/bin_sz_sec);
CPH = zeros(2 * nbins_half + 1, 1);

% TimeLags = (-nbins_half:nbins_half) * bin_sz_sec;
TimeLags = [(-nbins_half:nbins_half) * bin_sz_sec, nbins_half * bin_sz_sec + bin_sz_sec]; % ADJUST
if isvector(TimeLags) && size(TimeLags,1) == 1 % row vector is obtained
   TimeLags = TimeLags'; % make sure we return this as a colummn vector
end

for ev_idx = 1:numel(EVr) % for each 'reference' event...
   
   edges = EVr(ev_idx) + TimeLags;
   X = EVt(EVt >= edges(1) & EVt < edges(end));
   phist = histc(X, edges);
   
   if isvector(phist) && size(phist,1) == 1 % row vector is obtained
      phist = phist'; % convert to column
   end
   
   if norm_t_code == 3 && (~isempty(find(phist > 1, 1)))
      warning('MATLAB:paramAmbiguous', ...
      'CondProbHist: bin size is too big, event probability is incorrect');
   end
   
   if size(CPH, 1) ~= size(phist,1) - 1  || size(CPH, 2) ~= size(phist, 2) % ADJUST
      error('CondProbHist: wrong histogram size, ambiguous data?');
   end
   
   % CPH = CPH + phist; % add partial histogram to the main one
   CPH = CPH + phist(1:end-1); % add partial histogram to the main one % ADJUST
end

% apply requested normalization
if norm_t_code == 2
   % convert bin counts to number of events per second
   CPH = CPH / (numel(EVr) * bin_sz_sec);
elseif norm_t_code == 3
   % covert bin counts to conditional probability
   CPH = CPH / numel(EVr);
end

TimeLags = TimeLags(1:end-1); % ADJUST

% plot output histogram if no output raguments requested
if (nargout < 1)
   bar(TimeLags, CPH);
   xlabel('Time, sec');
   if norm_t_code == 1
      ylabel('Number of events');
   elseif norm_t_code == 2
      ylabel('Number of events per time bin (event/second)');
   elseif norm_t_code == 3
      ylabel('Conditional probability of event');
   end
   axis tight;
end


end


