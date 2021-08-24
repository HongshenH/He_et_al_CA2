clc;
clear;

nSD_frmap_smth = 3;
bin_sz_cm = 2;
dsp_delay = 984;
s_fname_bhv = 'bhv.mat';
s_fname_suff_plf = '_plf_fine.mat';

s_fname_sites = 'sites.txt';
s_fname_csv = 'CA2Cre_DREADD_plf.csv';
dset_list = LoadTargets('trials_lt_fine.txt');

% 2 spikes within 10 ms considered as a burst on a single cell basis.
BURST_P = struct();
BURST_P.N = 2;
BURST_P.ISI_N = 0.01;

lcnt = 1;
list_out = cell(1,1);

for did = 1:numel(dset_list)
   
   C = textscan(dset_list{did},'%s %s','Delimiter','\t');
   dset_list{did} = C{1}{1};
   s_trial_fine   = C{2}{1};
   fprintf('Process dataset: %s\n', dset_list{did});
   
   [s, s_trial] = fileparts(dset_list{did});
   [s, s_dataset] = fileparts(s);
   [s, s_mouse] = fileparts(s);
   [s, s_group] = fileparts(s);
   clearvars s;
   
   ntt_file_list = LoadTargetFlist(dset_list{did},'*_TT?.NTT');
   
   % load behavior data
   BHV = load(fullfile(dset_list{did},s_fname_bhv));
   TS_TRIAL = [BHV.BEHAV.pos_ts_usec(1) BHV.BEHAV.pos_ts_usec(end)];
   pos_bins = BHV.BEHAV.outline(1):bin_sz_cm:BHV.BEHAV.outline(2);
   
   % 1D occupancy map 
   OCCMAP = BehavLtOccMap(BHV.BEHAV, BHV.LAPS, pos_bins);
   
   for fid = 1:numel(ntt_file_list)
      fprintf('\tProcess file: %s\n', ntt_file_list{fid});
      
      s_site = GetSite(fullfile(s_group,s_mouse,s_fname_sites), strcat('CSC',ntt_file_list{fid}(end-4),'.ncs'));
      if isempty(s_site)
         error('unable to extract site name');
      end
      
      % Load spike data from every ntt file. Same code will work for *.nse files too.
      % It is assumed that the input file contain data recorded from only single
      % trial (run1 or run2 etc.) and spikes are already sorted.
      [all_ts, all_cids, all_wfs, all_fets] = NlxGetSpikesAll( ntt_file_list{fid} );
      all_ts = all_ts - dsp_delay;
      
      % Split spike features by cell
      cell_fet = SpkSelectFet(all_fets, all_cids, false);
      
      % Split all spike timestamps by cell
      [cell_ts, cell_ids] = SpkSelectTs(all_ts, all_cids, false);
      
      % Split all spike waveforms by cell
      cell_wfs = SpkSelectWf(all_wfs, all_cids, false);
      %
      bad_cells = false(length(cell_ids),1);
      for ii = 1:length(cell_ids)
         if numel(cell_ts{ii}) < 50
            bad_cells(ii) = true;
         end
      end
      cell_fet(bad_cells) = [];
      cell_ids(bad_cells) = [];
      cell_ts(bad_cells)  = [];
      cell_wfs(bad_cells) = [];
      %
      % Calculate waveform properties
      WFP = SpkWFormProp2(cell_ids, cell_wfs, 1, true);
      WFP_BEST_CH = [WFP(1).best_ch; vertcat(WFP.best_ch)];
      
      % Calculate spike train properties
      STP = SpkTrainProp2(cell_ids, cell_ts, cell_wfs, cell_fet, WFP_BEST_CH, TS_TRIAL, 1);
      
      % remove cluster zero before place field calculation
      cell_fet = cell_fet(2:end);
      cell_ids = cell_ids(2:end);
      cell_ts  = cell_ts(2:end);
      cell_wfs = cell_wfs(2:end);
      
      % Calculate place fields etc.
      PLF = CalcPlaceFields(BHV.PARAM, BHV.BEHAV, cell_ts);
      
      % Cell types
      CT  = CalcCellTypeStrict2( WFP, STP, PLF );
      
      % Bursting properties on a per-cell basis
      BURST = SpkBurstPerCell(BURST_P, cell_ts, TS_TRIAL);
      
      %
      save( strcat(ntt_file_list{fid}(1:end-4),s_fname_suff_plf), ...
         'cell_ts', 'cell_ids', 'WFP', 'STP', 'PLF', 's_*', 'CT', ...
         'pos_bins', 'nSD_frmap_smth', 'bin_sz_cm', 'BURST*');
      %
      for cid = 1:length(PLF)
         fprintf('\tProcess cell: %i\n', cell_ids{cid});
         %
         % additional cell selection
         if isnan(PLF(cid).max_frf_x) || isnan(PLF(cid).max_frf_y)
            p2eot_proxy = NaN;
         else
            % calculate firing map peak-to-end-of-track proximity. This value
            % is bounded between 0 (peak at either end of the track)
            % and 0.5 (peak in the middle of the track).
            p2eot_proxy = min( [ PLF(cid).max_frf_x, abs(size(PLF(cid).rmap_smooth,2) - PLF(cid).max_frf_x) ] );
            p2eot_proxy = p2eot_proxy / size(PLF(cid).rmap_smooth,2);
         end
         %
         MUA_SPKD = BehavLtSpkDistr(BHV.BEHAV, BHV.LAPS, pos_bins, {cell_ts{cid}});
         
         Nspk_L = numel(vertcat(MUA_SPKD.spk_ts_L{:}));
         Nspk_R = numel(vertcat(MUA_SPKD.spk_ts_R{:}));
         
         if (Nspk_L + Nspk_R) < 10
            CT{cid} = strcat(CT{cid}, '_LOW_LpSp');
         end
         
         FRMAP_L = smoothn(MUA_SPKD.nspk_L ./ OCCMAP.Locmap, nSD_frmap_smth, bin_sz_cm, 'nanexcl', 1, 'correct', 1);
         FRMAP_R = smoothn(MUA_SPKD.nspk_R ./ OCCMAP.Rocmap, nSD_frmap_smth, bin_sz_cm, 'nanexcl', 1, 'correct', 1);
         
         FRMAP_L(isnan(FRMAP_L)) = 0;
         FRMAP_R(isnan(FRMAP_R)) = 0;
         
         if isrow(FRMAP_L) || isrow(FRMAP_R)
            error('wrong shape of FRMAP_L and/or FRMAP_R');
         end
         
         % 1D spatial cross-correlation (same distance)
         spacorr_same_dist = FRmapSpatXcorr(FRMAP_L, FRMAP_R);
         % 1D spatial cross-correlation (same position)
         spacorr_same_pos  = FRmapSpatXcorr(FRMAP_L, flipud(FRMAP_R));
         
         %
         s_out = sprintf('%s,%s,%s,%s,%s,%s,%s,%s', s_group, s_dataset, s_trial, s_trial_fine, ...
            ntt_file_list{fid}, num2str(cell_ids{cid}), ...
            s_mouse, s_site );
         s_out = sprintf('%s,%s,', s_out, CT{cid} ); % CELL TYPE
         
         % GCID
         [~, s_tt] = fileparts(ntt_file_list{fid});
         s_gcid = strcat(s_mouse(5:end),s_dataset(4:5),num2str(100*str2num(s_tt(end))+cell_ids{cid}));
         s_out = sprintf('%s%s,', s_out, s_gcid );
         
         % WFP...
         s_out = sprintf('%s%.3f,', s_out, WFP(cid).wf_bad_prop);
         s_out = sprintf('%s%.3f,', s_out, WFP(cid).wf_peak);
         s_out = sprintf('%s%.3f,', s_out, WFP(cid).wf_swing);
         s_out = sprintf('%s%.3f,', s_out, WFP(cid).wf_width);
         s_out = sprintf('%s%.3f,', s_out, WFP(cid).wf_amp_ass);
         s_out = sprintf('%s%.3f,', s_out, WFP(cid).rms);
         
         % STP...
         s_out = sprintf('%s%i,', s_out, STP(cid).num_spk); % whole train!
         s_out = sprintf('%s%.4f,', s_out, STP(cid).frate_peak);
         s_out = sprintf('%s%.4f,', s_out, STP(cid).frate_mean);
         s_out = sprintf('%s%.4f,', s_out, STP(cid).perc_isi_u2ms);
         s_out = sprintf('%s%.4f,', s_out, STP(cid).csi_swing);
         s_out = sprintf('%s%.4f,', s_out, STP(cid).csi_peaks);
         s_out = sprintf('%s%.6f,', s_out, STP(cid).lratio);
         s_out = sprintf('%s%.4f,', s_out, STP(cid).isold);
         
         % PLF...
         s_out = sprintf('%s%i,', s_out, numel(PLF(cid).spk_pos_x)); % PLF!
         s_out = sprintf('%s%.4f,', s_out, PLF(cid).frate_peak);
         s_out = sprintf('%s%.4f,', s_out, PLF(cid).frate_mean);
         s_out = sprintf('%s%.4f,', s_out, PLF(cid).frate_std);
         s_out = sprintf('%s%.4f,', s_out, PLF(cid).frate_thr);
         s_out = sprintf('%s%.4f,', s_out, PLF(cid).spatinf_bspk);
         s_out = sprintf('%s%.4f,', s_out, PLF(cid).spatinf_bsec);
         s_out = sprintf('%s%.4f,', s_out, PLF(cid).max_frf_x);
         s_out = sprintf('%s%.4f,', s_out, PLF(cid).max_frf_y);
         s_out = sprintf('%s%.4f,', s_out, PLF(cid).max_frf_size);
         s_out = sprintf('%s%.4f,', s_out, PLF(cid).coherence);
         s_out = sprintf('%s%.4f,', s_out, PLF(cid).snr);
         s_out = sprintf('%s%.4f,', s_out, PLF(cid).in_field_frate);
         s_out = sprintf('%s%.4f,', s_out, PLF(cid).out_field_frate);
         s_out = sprintf('%s%.4f,', s_out, PLF(cid).perc_whole_sz);
         s_out = sprintf('%s%.4f,', s_out, PLF(cid).perc_pix_above_thr);
         s_out = sprintf('%s%.4f,', s_out, PLF(cid).sparsity_2d);
         
         s_out = sprintf('%s%i,', s_out, Nspk_L );
         s_out = sprintf('%s%i,', s_out, Nspk_R );
         s_out = sprintf('%s%.4f,', s_out, mean(MUA_SPKD.frate_L));
         s_out = sprintf('%s%.4f,', s_out, mean(MUA_SPKD.frate_R));
         s_out = sprintf('%s%.4f,', s_out, MUA_SPKD.frate_D);
         s_out = sprintf('%s%.4f,', s_out, MUA_SPKD.frate_DI1);
         s_out = sprintf('%s%.4f,', s_out, MUA_SPKD.frate_DI2);
         s_out = sprintf('%s%.4f,', s_out, spacorr_same_dist);
         s_out = sprintf('%s%.4f,', s_out, spacorr_same_pos);
         s_out = sprintf('%s%.4f,', s_out, p2eot_proxy);
         
         % BURST
         s_out = sprintf('%s%.4f,', s_out, BURST(cid).num_total);
         s_out = sprintf('%s%.4f,', s_out, BURST(cid).num_per_min);
         s_out = sprintf('%s%.4f,', s_out, BURST(cid).mean_ibi_sec);
         s_out = sprintf('%s%.4f,', s_out, BURST(cid).mean_dur_ms);
         s_out = sprintf('%s%.4f,', s_out, BURST(cid).burst_spk_total_spk);
         s_out = sprintf('%s%.4f,', s_out, BURST(cid).mean_spk_per_burst);
         
         s_out = sprintf('%s\n',s_out); % NEW LINE
         list_out{lcnt,1} = s_out;
         lcnt = lcnt + 1;
         %

      end
      fclose('all');
   end
end

%
f_out = fopen(s_fname_csv,'w');
fprintf(f_out, '%s', 'Group,Dataset,Trial,Trial_Fine,File,Cell_ID,Mouse,Site,Cell_TYPE,GCID,');
fprintf(f_out, '%s', 'Perc_Bad_WF,');
fprintf(f_out, '%s', 'WF_Peak_mV,');
fprintf(f_out, '%s', 'WF_Swing_mV,');
fprintf(f_out, '%s', 'WF_Width_usec,');
fprintf(f_out, '%s', 'WF_Amp_Assym,');
fprintf(f_out, '%s', 'WF_RMS,');
fprintf(f_out, '%s', 'train_Nspk,');
fprintf(f_out, '%s', 'train_FRpeak,');
fprintf(f_out, '%s', 'train_FRmean,');
fprintf(f_out, '%s', 'train_perc_ISI_u2ms,');
fprintf(f_out, '%s', 'train_CSI_swing,');
fprintf(f_out, '%s', 'train_CSI_peaks,');
fprintf(f_out, '%s', 'train_Lratio,');
fprintf(f_out, '%s', 'train_IsolD,');
fprintf(f_out, '%s', 'PLF_Nspk,');
fprintf(f_out, '%s', 'PLF_FRpeak,');
fprintf(f_out, '%s', 'PLF_FRmean,');
fprintf(f_out, '%s', 'PLF_SD,');
fprintf(f_out, '%s', 'PLF_Thr,');
fprintf(f_out, '%s', 'Spat_inf_bspk,');
fprintf(f_out, '%s', 'Spat_inf_bsec,');
fprintf(f_out, '%s', 'Max_FRF_X,');
fprintf(f_out, '%s', 'Max_FRF_Y,');
fprintf(f_out, '%s', 'Max_FRF_Size,');
fprintf(f_out, '%s', 'Coherence,');
fprintf(f_out, '%s', 'SNR,');
fprintf(f_out, '%s', 'In_field_FR,');
fprintf(f_out, '%s', 'Out_field_FR,');
fprintf(f_out, '%s', 'Perc_sampled_size,');
fprintf(f_out, '%s', 'Perc_bins_above_thr,');
fprintf(f_out, '%s', 'Sparsity_2D,');

fprintf(f_out, '%s', 'Nspk_L,');
fprintf(f_out, '%s', 'Nspk_R,');
fprintf(f_out, '%s', 'FRate_L,');
fprintf(f_out, '%s', 'FRate_R,');
fprintf(f_out, '%s', 'FRate_D,');
fprintf(f_out, '%s', 'FRate_DI1,');
fprintf(f_out, '%s', 'FRate_DI2,');
fprintf(f_out, '%s', 'Spat_Corr_same_dist,');
fprintf(f_out, '%s', 'Spat_Corr_same_pos,');
fprintf(f_out, '%s', 'Peak2EoT_proximity,');

fprintf(f_out, '%s', 'B_num_total,');
fprintf(f_out, '%s', 'B_num_per_min,');
fprintf(f_out, '%s', 'B_ibi_sec,');
fprintf(f_out, '%s', 'B_dur_msec,');
fprintf(f_out, '%s', 'Nspk_B_Nspk_train_ratio,');
fprintf(f_out, '%s', 'Spk_per_burst,');

fprintf(f_out, '\n'); % NEW LINE
for ii = 1:length(list_out)
   fprintf(f_out, '%s', list_out{ii});
end
%
fclose('all');


