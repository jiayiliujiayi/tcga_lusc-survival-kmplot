combine <- inner_join(ROBO1,EPHA2)
combine$diff <- combine$ROBO1 - combine$EPHA2

pcombine <- ggplot(combine, aes(ROBO1, EPHA2))
pcombine + geom_point()


#qqplot for robo1
ggplot(data = combine, aes(sample = ROBO1)) + 
  geom_qq() +
  stat_qq_line()
#qqplot for epha2
ggplot(data = combine, aes(sample = EPHA2)) + 
  geom_qq() +
  stat_qq_line()
#qqplot for diff
ggplot(data = combine, aes(sample = diff)) + 
  geom_qq() +
  stat_qq_line()
