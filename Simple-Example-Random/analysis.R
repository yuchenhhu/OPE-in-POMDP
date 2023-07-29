library(ggplot2)
theme_set(theme_bw())

point_estimate0 <- read.csv("result1.csv",header = F)
sd_estimate0 <- read.csv("result2.csv",header = F)
lepski_mse0 <- read.csv("result3.csv",header = F)
best_k0 <- read.csv("result4.csv",header = F)

k_all = -1:8
t = 300
n_all = seq(10,30,5)^2

point_estimate <- as.matrix(point_estimate0)
point_estimate <- array(as.numeric(point_estimate),dim=c(length(k_all),length(n_all),10000))
sd_estimate <- as.matrix(sd_estimate0)
sd_estimate <- array(sd_estimate,dim=c(length(k_all),length(n_all),10000))

# find the empirical mean
mean_point_estimate <- apply(point_estimate,c(1,2),mean)
mean_point_estimate <- replicate(10000, mean_point_estimate, simplify="array")

# find the standardized empirical distribution
standardized_point_estimate <- (point_estimate-mean_point_estimate)/sd_estimate

# find the p-value from a ks test
ks_statistic <- array(dim=c(length(n_all),length(k_all)))
for (i in 1:length(n_all)) {
    for (k in 1:length(k_all)) {
      ks_statistic[i,k] <- ks.test(scale(standardized_point_estimate[k,i,]),"pnorm")$statistic
  }
}
plot_ks <- data.frame(ks_statistic=as.numeric(ks_statistic),
                      n=(rep((n_all),10)),
                      k=rep(k_all,each=5))
plot_ks$k <- as.character(plot_ks$k)
ggplot(plot_ks) +
  geom_line(aes(k,ks_statistic,group=n,colour=n),size=1.1) +
  geom_point(aes(k,ks_statistic,group=n,colour=n))+
  scale_colour_gradientn(colours = c("orange","darkred"))+
  ylab('Kolmogorov–Smirnov statistics') +
  theme(
    #legend.position=c(.4, .7),
    legend.direction="horizontal",
    legend.position="none",
    panel.grid.minor = element_blank(),
    axis.text=element_text(size=20),
    axis.title=element_text(size=20,face="bold")
  )

# do the mse vs k plot
lepski_mse <- as.matrix(lepski_mse0)
lepski_mse <- array(as.numeric(lepski_mse),dim=c(length(k_all)+2,length(n_all),10000))
lepski_mse <- apply(lepski_mse,c(1,2),mean)

plot_mse <- data.frame(MSE=as.numeric(lepski_mse[1:length(k_all),]),
                       n=(rep((n_all),each=length(k_all))),
                       k=rep(paste0(k_all),5))

ggplot(plot_mse) +
  geom_line(aes(k,MSE,group=n,colour=n),size=1.1) +
  geom_point(aes(k,MSE,group=n,colour=n))+
  scale_colour_gradientn(colours = c("orange","darkred"))+
  ylab('Mean Squared Error') +
  theme(
    #legend.position=c(.4, .7),
    legend.direction="horizontal",
    legend.position="none",
    panel.grid.minor = element_blank(),
    axis.text=element_text(size=20),
    axis.title=element_text(size=20,face="bold")
  )

# do the mse lepski's plot
plot_lepski <- data.frame(MSE=as.numeric(lepski_mse),
                          n=(rep((n_all),each=(length(k_all)+2))),
                          k=rep((-1:length(k_all)),5))

ggplot(plot_lepski[!(plot_lepski$k%in%c(9,10)),]) +
  geom_line(aes(n,MSE,group=k,colour=k),size=1) +
  geom_point(aes(n,MSE,group=k,colour=k)) +
  #scale_color_brewer(palette="Greens")+
  geom_line(aes(n,MSE,group=k),data=plot_lepski[plot_lepski$k==9,],size=1.4,color="darkorchid",linetype="dashed") +
  geom_point(aes(n,MSE),data=plot_lepski[plot_lepski$k==9,],color="darkorchid") +
  geom_line(aes(n,MSE,group=k),data=plot_lepski[plot_lepski$k==10,],size=1.4,color="cornflowerblue",linetype="dotdash") +
  geom_point(aes(n,MSE),data=plot_lepski[plot_lepski$k==10,],color="cornflowerblue") +
  #scale_y_continuous(trans='sqrt')+
  #scale_x_discrete(expand=c(0.05,0.02)) +
  scale_colour_gradientn(colours = c("orange","darkred"))+
  ylab('Mean Squared Error') +
  theme(
    #legend.position=c(.4, .7),
    legend.direction="horizontal",
    legend.position="none",
    panel.grid.minor = element_blank(),
    axis.text=element_text(size=20),
    axis.title=element_text(size=20,face="bold")
  )

# lepski's share plot
best_k <- as.matrix(best_k0)
best_k <- array(as.numeric(best_k),dim=c(length(n_all),10000))
k_all = -1:8

percent <- NULL
for (i in 1:length(n_all)) {
  for (k in k_all) {
    percent <- c(percent,sum(best_k[i,]==k)/10000)
  }
}

plot_df <- data.frame(percentage=percent,
                      n=(rep(n_all,each=10)),
                      k=rep(-1:8,5))
library(RColorBrewer)
blues_fun <- colorRampPalette(brewer.pal(9,"Reds"))
plot_df$k <- as.factor(plot_df$k)
ggplot(plot_df[plot_df$k!=-1,], aes(fill=k, y=percentage, x=n,group=k)) + 
  scale_fill_manual("Fancy title",values=blues_fun(10))+
  geom_area(position="fill")+
  annotate("text", x = 215, y = 0.68, label = "k = 2", size = 6) +
  annotate("text", x = 540, y = 0.45, label = "k = 3", size = 6) +
  annotate("text", x = 750, y = 0.17, label = "k = 4", size = 6) +
  theme(
    #legend.position=c(.4, .7),
    legend.position="none",
    panel.grid.minor = element_blank(),
    axis.text=element_text(size=20),
    axis.title=element_text(size=20,face="bold")
  )
