function test_xdfdelay(xdfevent)

xdfevent=xdfevent(contains({xdfevent.value}, {'TTL_on', 'TTL_off', 'start_run', 'stop_run', 'start_block', 'stop_block', 'sync'})); % only select relevant events
lsldert00_events=ft_filter_event(xdfevent, 'type', 'Digital Triggers @ lsldert00');
events{1}=[lsldert00_events.timestamp];
lsldert01_events=ft_filter_event(xdfevent, 'type', 'Digital Triggers @ lsldert01');
events{2}=[lsldert01_events.timestamp];
lsldert02_events=ft_filter_event(xdfevent, 'type', 'Digital Triggers @ lsldert02');
events{3}=[lsldert02_events.timestamp];
lsldert03_events=ft_filter_event(xdfevent, 'type', 'Digital Triggers @ lsldert03');
events{4}=[lsldert03_events.timestamp];
lsldert04_events=ft_filter_event(xdfevent, 'type', 'Digital Triggers @ lsldert04');
events{5}=[lsldert04_events.timestamp];

figure
subplot(3,1,1); title('timestamps for all events'); hold on
for i=1:5
  plot(events{i}, ones(size(events{i}))*i, '+')
end
ylim([-3 10])
legend({'lsldert00', 'lsldert01', 'lsldert02', 'lsldert03', 'lsldert04'}, 'Location', 'eastoutside')

subplot(3,1,2); title('delay of timestamp relative to when the trigger was sent via lsldert00'); hold on;
for i=2:5
  plot(events{1}, (events{i}-events{1})*1000)
end
legend({'lsldert01', 'lsldert02', 'lsldert03', 'lsldert04'}, 'Location', 'eastoutside')
ylabel('time delay (msec.)'); xlabel('timestamps');yticks([0:10:50]); grid on
ylim([0 50]);

subplot(3,1,3); title('delay between raspberry pies (relative to lsldert04)'); hold on;
for i=2:4
  plot(events{5}, (events{i}-events{5})*1000)
  if any([(events{i}-events{5})*1000])>100
    warning('delay between some of the events was more than 100 msec')
  end
end
legend({'lsldert01', 'lsldert02', 'lsldert03'}, 'Location', 'eastoutside')
ylabel('time delay (msec.)'); xlabel('timestamps'); yticks([0:10:50]); grid on
ylim([0 50]);

