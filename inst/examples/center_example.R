# Less data for shorter example
repdata2 <- repdata[1:10,c(1:3,6:10,51)]

### One variable ----

# grand-mean
grand.mean(repdata2$item01)
grand.mean(repdata2$item01,wt = repdata2$wt)

# group-mean
group.mean(repdata2$item01,group = repdata2$GROUP)
group.mean(repdata2$item01,group = repdata2$GROUP,wt = repdata2$wt)

### More than one variable with the same rule ----

# grand-mean
grand.mean(repdata2[,4:8])
grand.mean(repdata2[,4:8],wt = repdata2$wt)

# group-mean
group.mean(repdata2[,4:8],group = repdata2$GROUP)
group.mean(repdata2[,4:8],group = repdata2$GROUP,wt = repdata2$wt)

### More than one variable with different rules ----
center(repdata2, group = repdata2$GROUP, grandmean = 4:5, groupmean = 6:8, wt = repdata2$wt)
center(repdata2, group = repdata2$GROUP, grandmean = 6:8, groupmean = 4:5, wt = repdata2$wt)

center(repdata2, group = repdata2$GROUP, wt = repdata2$wt,
       grandmean = paste0("item0",1:3), groupmean = paste0("item0",4:5))
center(repdata2, group = repdata2$GROUP, wt = repdata2$wt,
       grandmean = paste0("item0",4:5), groupmean = paste0("item0",1:3))
