
clc;
clear;

s_fname_sites = 'sites.txt';
dset_list = LoadTargets('trials_sleep.txt');
s_fname_suff_rpl = '*_rpl.mat';
s_fname_csv = 'CA2Cre_DREADD_rpl.csv';

lcnt = 1;
list_out = cell(1,1);

for did = 1 : numel(dset_list) % % ---------------------------------------------------
   fprintf('Process dataset: %s\n', dset_list{did});
   
   [s, s_trial] = fileparts(dset_list{did});
   [s, s_dataset] = fileparts(s);
   [s, s_mouse] = fileparts(s);
   [s, s_group] = fileparts(s);
   clearvars s;
   
   mat_file_list = LoadTargetFlist(dset_list{did},s_fname_suff_rpl);
      
   for fid = 1:numel(mat_file_list)
      fprintf('\tProcess file: %s\n', mat_file_list{fid});
      
      R = load(mat_file_list{fid});
      
      if R.RPL_DATA.num_detected == 0
         fprintf('NO RIPPLES DETECTED. SKIP TARGET\n');
         continue;
      end
      
      thr_amp_high = R.RPL_DATA.smoenv_mean + R.RPL_PAR.thr_nSD(2) * R.RPL_DATA.smoenv_std;
      r_amp_raw_mV = mean(R.RPL_DATA.amp_norm * thr_amp_high);
      
      RPL_TS = R.RPL_DATA.ts(:,2);
      r_num_left = R.RPL_DATA.num_detected;
      r_num_per_min = r_num_left / ((R.RPL_DATA.ts_end_usec-R.RPL_DATA.ts_beg_usec)/(1e6*60));
      r_mean_iri_sec = mean(diff( RPL_TS/1e6 ));
      r_ts = R.RPL_DATA.ts ./ 1e6;
      r_mean_dur_msec = 1000 * mean(r_ts(:,3) - r_ts(:,1));
      r_peak_psd_freq_raw = mean(R.RPL_DATA.peak_psd_freq_raw);
      r_peak_psd_freq_norm = mean(R.RPL_DATA.peak_psd_freq_norm);
      r_amp_norm = mean(R.RPL_DATA.amp_norm);
      r_type = R.RPL_DATA.type;
      r_num_slow = numel(find(r_type == 1));
      r_num_fast = numel(find(r_type == 2));
      r_num_na   = numel(find(r_type == 3));
      r_perc_slow = (100 * r_num_slow) / r_num_left;
      r_perc_fast = (100 * r_num_fast) / r_num_left;
      r_rate_per_min_slow = r_num_slow / ((R.RPL_DATA.ts_end_usec-R.RPL_DATA.ts_beg_usec)/(1e6*60));
      r_rate_per_min_fast = r_num_fast / ((R.RPL_DATA.ts_end_usec-R.RPL_DATA.ts_beg_usec)/(1e6*60));
      r_nspk_within_rpl = mean(R.RPL_DATA.nspk_within_rpl);
      
      s_out = sprintf('%s,%s,%s,%s,%s,%s,', R.s_group, R.s_mouse, R.s_dataset, R.s_trial, R.s_site, R.s_ncs);
      
      s_out = sprintf('%s%i,',   s_out, r_num_left);
      s_out = sprintf('%s%.0f,', s_out, r_num_per_min);
      s_out = sprintf('%s%.3f,', s_out, r_mean_iri_sec);
      s_out = sprintf('%s%.3f,', s_out, r_mean_dur_msec);
      s_out = sprintf('%s%.3f,', s_out, r_peak_psd_freq_raw);
      s_out = sprintf('%s%.3f,', s_out, r_peak_psd_freq_norm);
      s_out = sprintf('%s%.3f,', s_out, r_amp_raw_mV);
      s_out = sprintf('%s%.3f,', s_out, r_amp_norm);
      
      s_out = sprintf('%s%i,', s_out, r_num_fast);
      s_out = sprintf('%s%i,', s_out, r_num_slow);
      s_out = sprintf('%s%i,', s_out, r_num_na);
      
      s_out = sprintf('%s%0.2f,', s_out, r_perc_fast );
      s_out = sprintf('%s%0.2f,', s_out, r_perc_slow );
      
      s_out = sprintf('%s%0.2f,', s_out, r_rate_per_min_fast );
      s_out = sprintf('%s%0.2f,', s_out, r_rate_per_min_slow );
      
      s_out = sprintf('%s%.3f,',  s_out, r_nspk_within_rpl);
      
      s_out = sprintf('%s\n', s_out);     % NEW LINE
      
      list_out{lcnt,1} = s_out;
      lcnt = lcnt + 1;
   end
   fclose('all');
end

f_out = fopen(s_fname_csv,'w');
fprintf(f_out, '%s', 'Group,Mouse,Dataset,Trial,Site,File,');
fprintf(f_out, '%s', 'Num total,');
fprintf(f_out, '%s', 'Num. per min.,');
fprintf(f_out, '%s', 'Mean IRI sec,');
fprintf(f_out, '%s', 'Mean dur. ms,');
fprintf(f_out, '%s', 'Peak PSD(raw) freq. Hz,');
fprintf(f_out, '%s', 'Peak PSD(norm) freq. Hz,');
fprintf(f_out, '%s', 'Mean ampl. mV,');
fprintf(f_out, '%s', 'Mean ampl. nSD above thr.,');
fprintf(f_out, '%s', 'Num. fast,');
fprintf(f_out, '%s', 'Num. slow,');
fprintf(f_out, '%s', 'Num. NA,');
fprintf(f_out, '%s', '% fast,');
fprintf(f_out, '%s', '% slow,');
fprintf(f_out, '%s', 'Rate of fast r/min,');
fprintf(f_out, '%s', 'Rate of slow r/min,');
fprintf(f_out, '%s', 'Nspk within rpl.,');
fprintf(f_out, '\n');

for ii = 1:length(list_out)
   fprintf(f_out, '%s', list_out{ii});
end

fclose('all');






