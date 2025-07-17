library(R.matlab)

baro_table <- readMat('data/msftbarot_NorESM2-LM_ssp126_detection_data_Atlantic.mat')

baro_data <- data.frame(lons = baro_table$LONS[1,], lats = baro_table$LATS[1,], max_detect = baro_table$detection.index[1,], when_detect = baro_table$tip.time[1,])

baro_data_na_rm <- baro_data[-which(is.na(baro_data$max_detect)),]

baro_data_na_rm$max_detect_scaled <- (baro_data_na_rm$max_detect+1)/2
baro_data_na_rm$when_detect_scaled <- (baro_data_na_rm$when_detect-1850)/250

baro_wss <- rep(NA, 15)

for (i in 1:15) {
	print(i)
	baro_wss[i] <- kmeans(cbind(baro_data_na_rm$max_detect_scaled, baro_data_na_rm$when_detect_scaled), i)$tot.withinss
}

ratio_diffs <- array(NA, dim=c(1,14))
ratio_diffs[1,] <- diff(baro_wss)/baro_wss[1:14]

cluster_nums <- rep(3, 1)

baro_clusters <- kmeans(cbind(baro_data_na_rm$max_detect_scaled, baro_data_na_rm$when_detect_scaled), cluster_nums[1])

baro_data_na_rm$cluster <- baro_clusters$cluster

write.csv(baro_data_na_rm, file='msftbarot_NorESM2-LM_ssp126_clusters_scaled.csv')
