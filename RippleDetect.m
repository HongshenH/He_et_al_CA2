clc;
clear;

s_fname_sites = 'sites.txt';
dset_list = LoadTargets('trials_sleep.txt');
s_fname_suff_rpl = '_rpl.mat';            

RPL_PAR = struct();
RPL_PAR.smooth_sec = 0.005; % 5 ms
RPL_PAR.thr_nSD = [0 3.0];  % number of SDs used as thr. for detection
RPL_PAR.w_lim = 0.03; % shortest ripple width limit, in seconds
RPL_PAR.nw = 1.5;     % see 'doc pmtm'
RPL_PAR.nfft = 2^12;  % length of fft 4096
RPL_PAR.num_trg = 2000;
RPL_PAR.slowrange = [100 130]; % in Herz
RPL_PAR.fastrange = [140 180]; % in Herz
RPL_PAR.t_prev_us = 2e5;
RPL_PAR.t_post_us = 2e5;

for did = 1 : numel(dset_list) % % ---------------------------------------------------
   fprintf('Process dataset: %s\n', dset_list{did});
   
   [s, s_trial] = fileparts(dset_list{did});
   [s, s_dataset] = fileparts(s);
   [s, s_mouse] = fileparts(s);
   [~, s_group] = fileparts(s);
   clearvars s;
   
   NTT_file_list = LoadTargetFlist(dset_list{did},'*_TT?.NTT');
   
   ncs_file_list = cell(size(NTT_file_list));
   
   if isempty(dir(fullfile(dset_list{did},'*HR_FSI20.ncs')))
      s_tmp = '_FSI20.ncs';
   else
      s_tmp = 'HR_FSI20.ncs';
   end
   
   for ii = 1:length(ncs_file_list)
      [s1, s2] = fileparts(NTT_file_list{ii});
      ncs_file_list{ii,1} = fullfile(s1, strcat('CSC',s2(end),s_tmp));
   end
   
   clearvars s s1 s2 s_tmp NTT_file_list;
   
   for fid = 1:numel(ncs_file_list)
      fprintf('\tProcess file: %s\n', ncs_file_list{fid});
      %
      [~, s_ncs] = fileparts(ncs_file_list{fid});
      s_site = GetSite(fullfile(s_group,s_mouse,s_fname_sites), strcat('CSC',s_ncs(4),'.ncs'));
      if isempty(s_site)
         error('unable to extract site name');
      end
      %
      % search for matched NTT file.
      ntt_file_list = LoadTargetFlist(dset_list{did},strcat('*_TT',s_ncs(4),'.NTT'));
      if length(ntt_file_list) ~= 1
         error('unable to find matched NTT file');
      end
      %
      s_mat = strcat(ncs_file_list{fid}(1:end-4), s_fname_suff_rpl );
      %
      if exist(s_mat, 'file')
         fprintf('*** SKIP FILE ***\n');
         continue;
      end
      %
      % count number of units in the matched NTT file
      [all_ts, all_cids] = NlxGetSpikesAll( ntt_file_list{1} );
      cell_ids = unique(all_cids);
      cell_ids(cell_ids == 0) = [];
      if isempty(cell_ids)
         error('unsorted file?');
      end
      % remove unsorted spikes
      all_ts(all_cids == 0) = [];
      % detect ripples
      RPL_DATA = RipplesDetectWideBand(ncs_file_list{fid}, RPL_PAR, all_ts);
      if RPL_DATA.num_detected == 0
         % error('NO RIPPLES DETECTED');
         fprintf('NO RIPPLES DETECTED\n');
      end
      save(s_mat, 'RPL_DATA', 'RPL_PAR', 's_*');
      %
   end
   fclose('all');
end
fclose('all');






