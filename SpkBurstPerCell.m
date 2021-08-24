function [ BURST_OUT ] = SpkBurstPerCell( BURST_P, CELL_TS_USEC, TS_TRIAL_USEC )
% function [ BURST_OUT ] = SpkBurstPerCell( BURST_P, CELL_TS_USEC, TS_TRIAL_USEC )

% *  Copyright (C) 2016 Denis Polygalov,
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
% *  http://fsf.org/

BURST_IN = struct();
BURST_OUT = struct();

for cid = 1:length(CELL_TS_USEC)
   
   BURST_IN.T = CELL_TS_USEC{cid}/1e6;
   B_OUT  = BurstDetectISIn(BURST_IN, BURST_P.N, BURST_P.ISI_N);
   B_TS = [B_OUT.T_start; B_OUT.T_end]';
   
   spk_per_burst = zeros(numel(B_OUT.T_start),1);
   for bid = 1:numel(spk_per_burst)
      spk_per_burst(bid) = numel( find( ...
         (BURST_IN.T <= B_TS(bid,2)) & ...
         (BURST_IN.T >= B_TS(bid,1)) ) ...
      );
   end
   
   BURST_OUT(cid).num_total = numel(B_OUT.T_start);
   BURST_OUT(cid).num_per_min = BURST_OUT(cid).num_total / (diff(TS_TRIAL_USEC)/60e6);
   BURST_OUT(cid).mean_ibi_sec = mean(diff(B_OUT.T_start));
   BURST_OUT(cid).mean_dur_ms = mean(B_OUT.T_end - B_OUT.T_start) * 1e3;
   BURST_OUT(cid).burst_spikes = SpkFilterByTrials( BURST_IN.T, B_TS, 'rm_out' );
   BURST_OUT(cid).burst_spk_total_spk = numel(BURST_OUT(cid).burst_spikes) / numel(BURST_IN.T);
   BURST_OUT(cid).mean_spk_per_burst = mean(spk_per_burst);
end

end

