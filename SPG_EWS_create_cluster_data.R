library(R.matlab)

tas_table <- readMat('Data/tas_NorESM2-LM_ssp126_detection_data_Atlantic.mat')
baro_table <- readMat('Data/msftbarot_NorESM2-LM_ssp126_detection_data_Atlantic.mat')
lotst_table <- readMat('Data/mlotst_NorESM2-LM_ssp126_detection_data_Atlantic.mat')

tas_data <- data.frame(lons = tas_table$LONS[1,], lats = tas_table$LATS[1,], max_detect = tas_table$detection.index[1,], when_detect = tas_table$tip.time[1,])
baro_data <- data.frame(lons = baro_table$LONS[1,], lats = baro_table$LATS[1,], max_detect = baro_table$detection.index[1,], when_detect = baro_table$tip.time[1,])
lotst_data <- data.frame(lons = lotst_table$LONS[1,], lats = lotst_table$LATS[1,], max_detect = lotst_table$detection.index[1,], when_detect = lotst_table$tip.time[1,])

tas_data_na_rm <- tas_data
baro_data_na_rm <- baro_data[-which(is.na(baro_data$max_detect)),]
lotst_data_na_rm <- lotst_data[-which(is.na(lotst_data$max_detect)),]

tas_data_na_rm$max_detect_scaled <- (tas_data_na_rm$max_detect+1)/2
tas_data_na_rm$when_detect_scaled <- (tas_data_na_rm$when_detect-1850)/250
baro_data_na_rm$max_detect_scaled <- (baro_data_na_rm$max_detect+1)/2
baro_data_na_rm$when_detect_scaled <- (baro_data_na_rm$when_detect-1850)/250
lotst_data_na_rm$max_detect_scaled <- (lotst_data_na_rm$max_detect+1)/2
lotst_data_na_rm$when_detect_scaled <- (lotst_data_na_rm$when_detect-1850)/250

tas_wss <- rep(NA, 15)
baro_wss <- rep(NA, 15)
lotst_wss <- rep(NA, 15)

for (i in 1:15) {
	print(i)
	tas_wss[i] <- kmeans(cbind(tas_data_na_rm$max_detect_scaled, tas_data_na_rm$when_detect_scaled), i)$tot.withinss
	baro_wss[i] <- kmeans(cbind(baro_data_na_rm$max_detect_scaled, baro_data_na_rm$when_detect_scaled), i)$tot.withinss
	lotst_wss[i] <- kmeans(cbind(lotst_data_na_rm$max_detect_scaled, lotst_data_na_rm$when_detect_scaled), i)$tot.withinss
}

ratio_diffs <- array(NA, dim=c(3,14))
ratio_diffs[1,] <- diff(tas_wss)/tas_wss[1:14]
ratio_diffs[2,] <- diff(baro_wss)/baro_wss[1:14]
ratio_diffs[3,] <- diff(lotst_wss)/lotst_wss[1:14]


cluster_nums <- rep(3, 3)

tas_clusters <- kmeans(cbind(tas_data_na_rm$max_detect_scaled, tas_data_na_rm$when_detect_scaled), cluster_nums[1])
baro_clusters <- kmeans(cbind(baro_data_na_rm$max_detect_scaled, baro_data_na_rm$when_detect_scaled), cluster_nums[2])
lotst_clusters <- kmeans(cbind(lotst_data_na_rm$max_detect_scaled, lotst_data_na_rm$when_detect_scaled), cluster_nums[3])

tas_data_na_rm$cluster <- tas_clusters$cluster
baro_data_na_rm$cluster <- baro_clusters$cluster
lotst_data_na_rm$cluster <- lotst_clusters$cluster

write.csv(tas_data_na_rm, file='tas_NorESM2-LM_ssp126_clusters_scaled.csv')
write.csv(baro_data_na_rm, file='msftbarot_NorESM2-LM_ssp126_clusters_scaled.csv')
write.csv(lotst_data_na_rm, file='mlotst_NorESM2-LM_ssp126_clusters_scaled.csv')









load('spg_clusters_140224.RData')

pdf('cluster_choice_spg.pdf', width=10, height=5)
par(mfrow=c(1,2))
plot(tas_wss, type='b', main='tas-GC31-LL', xlab='Clusters', ylab='Total WSS')
points(cluster_nums[1], tas_wss[cluster_nums[1]], col='red', pch=16)
plot(baro_wss, type='b', main='CAMS-CSM1-0', xlab='Clusters', ylab='Total WSS')
points(cluster_nums[2], baro_wss[cluster_nums[2]], col='red', pch=16)
dev.off()

#pdf('tas_cluster_plot_scaled.pdf')
#plot(tas_data_na_rm$max_detect, tas_data_na_rm$when_detect, pch=16, col=rainbow(6)[tas_data_na_rm$cluster], main='tas-GC31-LL', xlab='Max Detection', ylab='Time of Tipping')
#points(tas_clusters$centers[,1]*2-1, tas_clusters$centers[,2]*250+1850, cex=3)
#dev.off()

load('spg_clusters_080224.RData')





