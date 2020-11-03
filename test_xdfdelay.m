function test_xdfdelay(xdfevent)

events{1}=[xdfevent(strcmp({xdfevent.type}, 'Digital Triggers @ lsldert00')).timestamp];
events{1}=events{1}(1:end-1); % last event is empty
events{2}=[xdfevent(strcmp({xdfevent.type}, 'Digital Triggers @ lsldert01')).timestamp];
events{3}=[xdfevent(strcmp({xdfevent.type}, 'Digital Triggers @ lsldert02')).timestamp];
events{4}=[xdfevent(strcmp({xdfevent.type}, 'Digital Triggers @ lsldert03')).timestamp];
events{5}=[xdfevent(strcmp({xdfevent.type}, 'Digital Triggers @ lsldert04')).timestamp];

figure;
subplot(3,1,1); title('timestamps for all events'); hold on
for i=1:5
  plot(events{i}, ones(size(events{i}))*i, '+')
end
ylim([-23 30])
legend({'lsldert00', 'lsldert01', 'lsldert02', 'lsldert03', 'lsldert04'}, 'Location', 'eastoutside')

subplot(3,1,2); title('delay of timestamp relative to when the trigger was sent via lsldert00'); hold on;
for i=2:5
  plot(events{1}, (events{i}-events{1})*1000)
end
legend({'lsldert01', 'lsldert02', 'lsldert03', 'lsldert04'}, 'Location', 'eastoutside')
ylabel('time delay (msec.)'); xlabel('timestamps');
ylim([0 100]);

subplot(3,1,3); title('delay between raspberry pies (relative to lsldert04)'); hold on;
for i=2:4
  plot(events{5}, (events{i}-events{5})*1000)
  if any([(events{i}-events{5})*1000])>100
    warning('delay between some of the events was more than 100 msec')
  end
end
legend({'lsldert01', 'lsldert02', 'lsldert03'}, 'Location', 'eastoutside')
ylabel('time delay (msec.)'); xlabel('timestamps');
ylim([0 100]);

