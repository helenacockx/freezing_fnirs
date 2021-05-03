function lsl_events=test_xdfdelay(xdfevent)

xdfevent=xdfevent(contains({xdfevent.value}, {'TTL_on', 'TTL_off', 'start_run', 'stop_run', 'start_block', 'stop_block', 'sync'})); % only select relevant events
lsldert00_events=ft_filter_event(xdfevent, 'type', 'Digital Triggers @ lsldert00');
events{1}=[lsldert00_events.timestamp];
values{1}={lsldert00_events.value};
lsldert01_events=ft_filter_event(xdfevent, 'type', 'Digital Triggers @ lsldert01');
events{2}=[lsldert01_events.timestamp];
values{2}={lsldert01_events.value};
lsldert02_events=ft_filter_event(xdfevent, 'type', 'Digital Triggers @ lsldert02');
events{3}=[lsldert02_events.timestamp];
values{3}={lsldert02_events.value};
lsldert03_events=ft_filter_event(xdfevent, 'type', 'Digital Triggers @ lsldert03');
events{4}=[lsldert03_events.timestamp];
values{4}={lsldert03_events.value};
lsldert04_events=ft_filter_event(xdfevent, 'type', 'Digital Triggers @ lsldert04');
events{5}=[lsldert04_events.timestamp];
values{5}={lsldert04_events.value};
lsldert05_events=ft_filter_event(xdfevent, 'type', 'Digital Triggers @ lsldert05');
events{6}=[lsldert05_events.timestamp];
values{6}={lsldert05_events.value};

%% plot delays 
figure;
subplot(3,1,1); title('timestamps for all events'); hold on
for i=1:6
  try
  plot(events{i}, ones(size(events{i}))*i, '+')
  end
end
ylim([-3 10])
if ~isempty(lsldert05_events)
  legend({'lsldert0', 'lsldert1', 'lsldert2', 'lsldert3', 'lsldert4', 'lsldert5'}, 'Location', 'eastoutside')
else
  legend({'lsldert0', 'lsldert1', 'lsldert2', 'lsldert3', 'lsldert4'}, 'Location', 'eastoutside')
end
% figure; title('delay between creating the trigger and sending it via lsldert0'); hold on;
% plot(events{1}, events{2}-events{1})
% ylabel('time delay (sec.)'); xlabel('timestamps');
% ylim([0 0.1]);
 
subplot(3,1,2); title('delay of timestamp relative to when the trigger was sent via lsldert0'); hold on;
for i=2:6
  try
  plot(events{1}, events{i}-events{1})
  end
end
if ~isempty(lsldert05_events)
  legend({'lsldert1', 'lsldert2', 'lsldert3', 'lsldert4', 'lsldert5'},'Location', 'eastoutside')
else
  legend({'lsldert1', 'lsldert2', 'lsldert3', 'lsldert4'},'Location', 'eastoutside')
end
ylabel('time delay (sec.)'); xlabel('timestamps');
ylim([-0.005 0.05]);

subplot(3,1,3); title('delay between raspberry pies (relative to lsldert4)'); hold on;
for i=[2:4 6]
  try
  plot(events{5}, events{i}-events{5})
  end
end
if ~isempty(lsldert05_events)
  legend({'lsldert1', 'lsldert2', 'lsldert3', 'lsldert5'},'Location', 'eastoutside')
else
  legend({'lsldert1', 'lsldert2', 'lsldert3'},'Location', 'eastoutside')
end
ylabel('time delay (sec.)'); xlabel('timestamps');
ylim([-0.005 0.05]);

%% calculate corrections relative to lsldert4
% lsldert04
lsl_timestamps4=[lsldert04_events.timestamp]';
correction=zeros(length(lsl_timestamps4), 1);
duration=zeros(length(lsl_timestamps4), 1);
type=repmat({'lsldert4_events'},length(lsl_timestamps4),1);
% start making a table
lsl_events=table(lsl_timestamps4, correction, duration, type, {lsldert04_events.value}', 'VariableNames', {'onset', 'correction', 'duration', 'type', 'value'});
% loop over other lsl streams and calculate corrections relative to
% lsldert04
lsl_streams={'lsldert0_events' 'lsldert1_events', 'lsldert2_events', 'lsldert3_events'};
for i=2:4
  lsl_timestampsi=events{i}';
  if length(lsl_timestampsi)~= length(lsl_timestamps4)
    warning('%s has not the same number of events as lsldert4_events. Not calculating corrections', lsl_streams{i})
    continue
  end
  lsl_corri=lsl_timestampsi-lsl_timestamps4;
  duration=zeros(length(lsl_timestampsi),1);
  type=repmat({lsl_streams{i}},length(lsl_timestampsi),1);
  % make table
  tablei=table(lsl_timestampsi, lsl_corri, duration, type, values{i}', 'VariableNames', {'onset', 'correction', 'duration', 'type', 'value'});
  if any(lsl_corri>0.02)
    warning('some events of %s had a delay >20 msec:', lsl_streams{i})
    idx=find(lsl_corri>0.02);
    display(tablei(idx,:))
  end
  % concatenate with lsl_events
  lsl_events=[lsl_events; tablei];
end

% lsldert05: compare to TTLon en TTLoff for lsldert04(FIX ME: maybe use lsldert0 for all lsl streams?, but lsldert0 and 4 are not so different...)
if ~isempty(lsldert05_events)
  TTL_events=xdfevent(endsWith({xdfevent.value},  {'TTL_on', 'TTL_off'}));
  events_lsl4=xdfevent(strcmp({xdfevent.type}, 'Digital Triggers @ lsldert04') & endsWith({xdfevent.value},  {'TTL_on', 'TTL_off'}));
  events_lsl5=xdfevent(strcmp({xdfevent.type}, 'Digital Triggers @ lsldert05') & endsWith({xdfevent.value},  {'PULSEIR_3'}));
  lsl_timestamps5=[events_lsl5.timestamp]';
  if length(events_lsl5)~= length(events_lsl4)
    warning('lsldert5 has not the same number of events as lsldert4_events. Not calculating corrections')
  else
    lsl_corr5=[events_lsl5.timestamp]'-[events_lsl4.timestamp]'; % remark: this is quite long because lsldert5 first produces a beep of 0.5 s (start_run) or 0.2 s (stop_run), before creating the IR pulse
    duration=zeros(length(lsl_timestamps5),1);
    type=repmat({'lsldert5_events'},length(lsl_timestamps5),1);
    table5=table(lsl_timestamps5, lsl_corr5, duration, type, {events_lsl5.value}', 'VariableNames', {'onset', 'correction', 'duration', 'type', 'value'});
    lsl_events=[lsl_events; table5];
  end
end
