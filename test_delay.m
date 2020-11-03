function test_delay(REC)
%% 
events{1}=[REC.trig_events.timestamp];
idx_events=find(endsWith(REC.lsldert0_events.Data, {'start_run', 'stop_run', 'start_block', 'stop_block', 'sync'}));
lsldert0_events=lsl_correct_lsl_timestamps(REC.lsldert0_events);
events{2}=lsldert0_events(idx_events);
idx_events=find(endsWith(REC.lsldert1_events.Data, {'start_run', 'stop_run','start_block', 'stop_block', 'sync'}));
lsldert1_events=lsl_correct_lsl_timestamps(REC.lsldert1_events);
events{3}=lsldert1_events(idx_events);
idx_events=find(endsWith(REC.lsldert2_events.Data, {'start_run', 'stop_run','start_block', 'stop_block', 'sync'}));
lsldert2_events=lsl_correct_lsl_timestamps(REC.lsldert2_events);
events{4}=lsldert2_events(idx_events);
idx_events=find(endsWith(REC.lsldert3_events.Data, {'start_run', 'stop_run','start_block', 'stop_block', 'sync'}));
lsldert3_events=lsl_correct_lsl_timestamps(REC.lsldert3_events);
events{5}=lsldert3_events(idx_events);
idx_events=find(endsWith(REC.lsldert4_events.Data, {'start_run', 'stop_run','start_block', 'stop_block', 'sync'}));
lsldert4_events=lsl_correct_lsl_timestamps(REC.lsldert4_events);
events{6}=lsldert4_events(idx_events);
 
 
figure; title('timestamps for all events'); hold on
for i=1:6
  plot(events{i}, ones(size(events{i}))*i, '+')
end
ylim([-23 30])
legend({'trigger', 'lsldert0', 'lsldert1', 'lsldert2', 'lsldert3', 'lsldert4'})

figure; title('delay between creating the trigger and sending it via lsldert0'); hold on;
plot(events{1}, events{2}-events{1})
ylabel('time delay (sec.)'); xlabel('timestamps');
ylim([0 0.1]);
 
figure; title('delay of timestamp relative to when the trigger was sent via lsldert0'); hold on;
for i=3:6
  plot(events{2}, events{i}-events{2})
end
legend({'lsldert1', 'lsldert2', 'lsldert3', 'lsldert4'})
ylabel('time delay (sec.)'); xlabel('timestamps');
ylim([0 0.1]);

figure; title('delay between raspberry pies (relative to lsldert1)'); hold on;
for i=4:6
  plot(events{3}, events{i}-events{3})
end
legend({'lsldert2', 'lsldert3', 'lsldert4'})
ylabel('time delay (sec.)'); xlabel('timestamps');
ylim([-0.05 0.05]);



