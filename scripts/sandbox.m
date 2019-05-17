

stop = 4000:4000:max(lhdata.time{:})
start = 1:4000:max(lhdata.time{:})
start = start(1:end-1)
offset = repmat(-3000, size(start))

trl = [start; stop; offset]'
cfg = [];
cfg.trl = trl;
tmlck = ft_redefinetrial(cfg, lhdata)


