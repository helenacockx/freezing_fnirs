function lsl_events=test_delay(REC)
 
%% load events
events{1}=[REC.trig_events.timestamp];
idx_events=find(endsWith(REC.lsldert0_events.Data, {'TTL_on', 'TTL_off', 'start_run', 'stop_run', 'start_block', 'stop_block', 'sync'}));
lsldert0_events=lsl_correct_lsl_timestamps(REC.lsldert0_events);
events{2}=lsldert0_events(idx_events);
idx_events=find(endsWith(REC.lsldert1_events.Data, {'TTL_on', 'TTL_off', 'start_run', 'stop_run','start_block', 'stop_block', 'sync'}));
lsldert1_events=lsl_correct_lsl_timestamps(REC.lsldert1_events);
events{3}=lsldert1_events(idx_events);
idx_events=find(endsWith(REC.lsldert2_events.Data, {'TTL_on', 'TTL_off', 'start_run', 'stop_run','start_block', 'stop_block', 'sync'}));
lsldert2_events=lsl_correct_lsl_timestamps(REC.lsldert2_events);
events{4}=lsldert2_events(idx_events);
idx_events=find(endsWith(REC.lsldert3_events.Data, {'TTL_on', 'TTL_off', 'start_run', 'stop_run','start_block', 'stop_block', 'sync'}));
lsldert3_events=lsl_correct_lsl_timestamps(REC.lsldert3_events);
events{5}=lsldert3_events(idx_events);
idx_events=find(endsWith(REC.lsldert4_events.Data, {'TTL_on', 'TTL_off', 'start_run', 'stop_run','start_block', 'stop_block', 'sync'}));
lsldert4_events=lsl_correct_lsl_timestamps(REC.lsldert4_events);
events{6}=lsldert4_events(idx_events);
if isfield(REC, 'lsldert5_events')
  idx_events=find(endsWith(REC.lsldert5_events.Data, {'TTL_on', 'TTL_off', 'start_run', 'stop_run','start_block', 'stop_block', 'sync'}));
  lsldert5_events=lsl_correct_lsl_timestamps(REC.lsldert5_events);
  events{7}=lsldert5_events(idx_events);
end
 
%% plot delays 
figure;
subplot(3,1,1); title('timestamps for all events'); hold on
for i=1:7
  try
  plot(events{i}, ones(size(events{i}))*i, '+')
  end
end
ylim([-3 10])
if isfield(REC, 'lsldert5_events')
  legend({'trigger', 'lsldert0', 'lsldert1', 'lsldert2', 'lsldert3', 'lsldert4', 'lsldert5'}, 'Location', 'eastoutside')
else
  legend({'trigger', 'lsldert0', 'lsldert1', 'lsldert2', 'lsldert3', 'lsldert4'}, 'Location', 'eastoutside')
end
% figure; title('delay between creating the trigger and sending it via lsldert0'); hold on;
% plot(events{1}, events{2}-events{1})
% ylabel('time delay (sec.)'); xlabel('timestamps');
% ylim([0 0.1]);
 
subplot(3,1,2); title('delay of timestamp relative to when the trigger was sent via lsldert0'); hold on;
for i=3:7
  try
  plot(events{2}, events{i}-events{2})
  end
end
if isfield(REC, 'lsldert5_events')
  legend({'lsldert1', 'lsldert2', 'lsldert3', 'lsldert4', 'lsldert5'},'Location', 'eastoutside')
else
  legend({'lsldert1', 'lsldert2', 'lsldert3', 'lsldert4'},'Location', 'eastoutside')
end
ylabel('time delay (sec.)'); xlabel('timestamps');
ylim([-0.005 0.05]);

subplot(3,1,3); title('delay between raspberry pies (relative to lsldert4)'); hold on;
for i=[3:5 7]
  try
  plot(events{6}, events{i}-events{6})
  end
end
if isfield(REC, 'lsldert5_events')
  legend({'lsldert1', 'lsldert2', 'lsldert3', 'lsldert5'},'Location', 'eastoutside')
else
  legend({'lsldert1', 'lsldert2', 'lsldert3'},'Location', 'eastoutside')
end
ylabel('time delay (sec.)'); xlabel('timestamps');
ylim([-0.005 0.05]);

%% calculate corrections relative to lsldert4
% lsldert04
idx_events4=find(endsWith(REC.lsldert4_events.Data, {'TTL_on', 'TTL_off', 'start_run', 'stop_run', 'start_block', 'stop_block', 'sync'}));
lsl_timestamps4=lsl_correct_lsl_timestamps(REC.lsldert4_events);
lsl_timestamps4=lsl_timestamps4(idx_events4)';
correction=zeros(length(lsl_timestamps4), 1);
duration=zeros(length(lsl_timestamps4), 1);
type=repmat({'lsldert4_events'},length(lsl_timestamps4),1);
% start making a table
lsl_events=table(lsl_timestamps4, correction, duration, type, {REC.lsldert4_events.Data{idx_events4}}', 'VariableNames', {'onset', 'correction', 'duration', 'type', 'value'});
% loop over other lsl streams and calculate corrections relative to
% lsldert04
lsl_streams={'lsldert1_events', 'lsldert2_events', 'lsldert3_events'};
for i=1:3
  idx_eventsi=find(endsWith(REC.(lsl_streams{i}).Data, {'TTL_on', 'TTL_off', 'start_run', 'stop_run', 'start_block', 'stop_block', 'sync'}));
  lsl_timestampsi=lsl_correct_lsl_timestamps(REC.(lsl_streams{i}));
  lsl_timestampsi=lsl_timestampsi(idx_eventsi)';
  if length(lsl_timestampsi)~= length(lsl_timestamps4)
    warning('%s has not the same number of events as lsldert4_events. Not calculating corrections', lsl_streams{i})
    continue
  end
  lsl_corri=lsl_timestampsi-lsl_timestamps4;
  duration=zeros(length(lsl_timestampsi),1);
  type=repmat({lsl_streams{i}},length(lsl_timestampsi),1);
  % make table
  tablei=table(lsl_timestampsi, lsl_corri, duration, type, {REC.(lsl_streams{i}).Data{idx_eventsi}}', 'VariableNames', {'onset', 'correction', 'duration', 'type', 'value'});
  if any(lsl_corri>0.02)
    warning('some events of %s had a delay >20 msec:', lsl_streams{i})
    idx=find(lsl_corri>0.02);
    display(tablei(idx,:))
  end
  % concatenate with lsl_events
  lsl_events=[lsl_events; tablei];
end

% lsldert05: compare to TTLon en TTLoff for lsldert04(FIX ME: maybe use lsldert0 for all lsl streams?, but lsldert0 and 4 are not so different...)
if isfield(REC, 'lsldert5_events')
  idx_events4=find(endsWith(REC.lsldert4_events.Data, {'TTL_on', 'TTL_off'}));
  lsl_timestamps4=lsl_correct_lsl_timestamps(REC.lsldert4_events);
  lsl_timestamps4=lsl_timestamps4(idx_events4)';
  idx_events5=find(endsWith(REC.lsldert5_events.Data, {'PULSEIR_3'}));
  lsl_timestamps5=lsl_correct_lsl_timestamps(REC.lsldert5_events);
  lsl_timestamps5=lsl_timestamps5(idx_events5)';
  if length(lsl_timestamps5)~= length(lsl_timestamps4)
    warning('lsldert5 has not the same number of events as lsldert4_events. Not calculating corrections')
  else
    lsl_corr5=lsl_timestamps5-lsl_timestamps4; % remark: this is quite long because lsldert5 first produces a beep of 0.5 s (start_run) or 0.2 s (stop_run), before creating the IR pulse
    duration=zeros(length(lsl_timestamps5),1);
    type=repmat({'lsldert5_events'},length(lsl_timestamps5),1);
    table5=table(lsl_timestamps5, lsl_corr5, duration, type, {REC.lsldert5_events.Data{idx_events5}}', 'VariableNames', {'onset', 'correction', 'duration', 'type', 'value'});
    lsl_events=[lsl_events; table5];
  end
end
