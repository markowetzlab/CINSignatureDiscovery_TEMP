x=readRDS("2_combined_mixmodels.rds")

# Replace components in segment size
s=x$segsize
s[14,] = c(8333403, 932916, 0.0334)
s=s[-15,]
s[16,]=c(13435586, 1439442, 0.0414)
s=s[-c(17:18),]
x$segsize=s

# Replace components in changepoint
s=x$changepoint
s[4,]=c(1.008185883, 0.2528096, 0.37043868)
s=s[-c(5:9),]
s[6,]=c(1.447916491, 0.4018533, 0.04827426775)
s=s[-c(7:9),]
s[7,]=c(2.20266882, 0.4656683, 0.2117234981)
s=s[-c(8:11),]
x$changepoint=s

saveRDS(x, "2_combined_mixmodels_merged_components.rds")
