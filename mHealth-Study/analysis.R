library(ggplot2)

point_estimate0 <- read.csv("/home/users/yuchenhu/POMDP/mhealth_n/all_info/result1.csv",header = F)
sd_estimate0 <- read.csv("/home/users/yuchenhu/POMDP/mhealth_n/all_info/result2.csv",header = F)
lepski_mse0 <- read.csv("/home/users/yuchenhu/POMDP/mhealth_n/all_info/result3.csv",header = F)
best_k0 <- read.csv("/home/users/yuchenhu/POMDP/mhealth_n/all_info/result4.csv",header = F)

k_all = -1:12
t = 240
n_all = seq(10,30,5)^2
b = 100000

point_estimate <- as.matrix(point_estimate0)
point_estimate <- array(as.numeric(point_estimate),dim=c(length(k_all),length(n_all),b))
sd_estimate <- as.matrix(sd_estimate0)
sd_estimate <- array(sd_estimate,dim=c(length(k_all),length(n_all),b))

# do the mse vs k plot
lepski_mse <- as.matrix(lepski_mse0)
lepski_mse <- array(as.numeric(lepski_mse),dim=c(length(k_all)+1,length(n_all),b))
lepski_mse <- apply(lepski_mse,c(1,2),mean)

plot_mse <- data.frame(MSE=as.numeric(lepski_mse[1:length(k_all),]),
                       n=(rep((n_all),each=length(k_all))),
                       k=rep(k_all,length(n_all)))
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
                          n=(rep((n_all),each=(length(k_all)+1))),
                          k=rep(-1:(k_all[length(k_all)]+1),length(n_all)))

ggplot(plot_lepski[plot_lepski$k!=k_all[length(k_all)]+1,]) +
  geom_line(aes(n,MSE,group=k,colour=k),size=1) +
  geom_point(aes(n,MSE,group=k,colour=k)) +
  #scale_color_brewer(palette="Greens")+
  geom_line(aes(n,MSE,group=k),data=plot_lepski[plot_lepski$k==k_all[length(k_all)]+1,],size=1.4,color="darkorchid",linetype="dashed") +
  geom_point(aes(n,MSE),data=plot_lepski[plot_lepski$k==k_all[length(k_all)]+1,],color="darkorchid") +
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

