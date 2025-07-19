#################
# Plot the density of psi under different sampling design
rm(list = ls())
library(dplyr)
library(ggplot2)
library(randomForest)
library(MASS)
library(foreach)
require(doSNOW)
library(caret)
library(ggpubr)
source("somefunction17.R")
set.seed(1)
# method_set <- c("aipw", "ipw", "reg")
# N_set = c(10000)
# n_m_set = c(500, 1000)
# n_set = c(1000, 1500, 2000, 2500)
# K_set = c(1:5)*2
# alpha_set = c(0.05, 0.1, 0.2, 0.3, 0.5)
# 
# expand.grid(method_set, N_set, n_m_set, n_set,K_set, alpha_set)->settings
# colnames(settings) <- c("method", "N", "n_m", "n", "K", "alpha")
# 
# settings[order(settings$method, settings$N, 
#                settings$n_m, settings$n, 
#                settings$K, settings$alpha), ] -> settings
# rownames(settings) <- c(1:nrow(settings))
# 
# j=186
# method = settings$method[j]
# N = settings$N[j]
# n_m = settings$n_m[j]
# n = settings$n[j]
# K = settings$K[j]
# alpha = settings$alpha[j]

method = "aipw"
N = 20000
n_m = 3000
n = 3000
K = 5
alpha = 0.1
cat(" Method=", method, ", N=", N,", n_m=", n_m,", n=", n, ", K=", K,  ", alpha=", alpha,"\n")


Generatedt(N) -> dt1 -> dt
#######
# Outcome dependent group
te <- dt1 %>% dplyr::select(id, Y)
te[order(te$Y), ] -> te
R_O = rep(1:K, each = ceiling(N/K)) %>% as.factor()
te$R_O =  R_O[1:N]
dt1$R_O = te[order(te$id), "R_O"]
######
# Influence function based group
x = cbind(dt1$X1, dt1$U1); y = dt1$Y; A = dt1$Z
esATE(x, y, A, method) -> fit1
esATE(dt1$X1, y, A, method) -> fit1_c
dt1$phi = fit1$infl

c(fit1$est - 2* sqrt(fit1$ve), fit1$est + 2* sqrt(fit1$ve))
c(fit1_c$est - 2* sqrt(fit1_c$ve), fit1_c$est + 2* sqrt(fit1_c$ve))

dt1$phi_m = fit1_c$infl
phi_bar = mean(dt1$phi)
phi_bar_m = mean(dt1$phi_m)
dt1$phi_a = abs(dt1$phi-phi_bar)
dt1$phi_m_a = abs(dt1$phi_m-phi_bar_m)
dt1 <- dt1[order(dt1$phi_m_a), ]; 
R_M = rep(1:K, each = ceiling(N/K)) %>% as.factor()
dt1$R_M = R_M[1:N]
dt1 <- dt1[order(dt1$phi_a), ];
R_T = rep(1:K, each = ceiling(N/K)) %>% as.factor()
dt1$R_T = R_T[1:N]

#####
# uniform
dt1$pi_u <- rep((n)/N, N)

#####
# best optimal
getpi(dt1$phi_a, (n), alpha)*(n) -> dt1$pi_b

############
# the middle phase data
midphase_ind <- sample(1:N, size = n_m, replace = F)
#### The interndata
dt2 <- dt1[midphase_ind,]
#### The otherdata
dt1_2 <- dt1[-midphase_ind,]
N1_2 = nrow(dt1_2)

x2 = cbind(dt2$X1, dt2$U1); y2 = dt2$Y; A2 = dt2$Z
esATE(x2, y2, A2, method) -> fit2
dt2$phi_hat <- fit2$infl
phi_hat_bar <- mean(dt2$phi_hat)
dt2$phi_hat_a <- abs(dt2$phi_hat - phi_hat_bar)
dt2 <- dt2[order(dt2$phi_hat_a), ]
R_L = rep(1:K, each = ceiling(n_m/K)) %>% as.factor()
dt2$R_L = R_L[1:n_m]
te <- c()
for(k in 1:K){
  te[k] <- mean((dt2[dt2$R_L==k,"phi_hat"]-phi_hat_bar)^2)
}
te2 <- te
getpi_stra(sqrt(te2), N1_2, n, K, p = rep(1/K, K), alpha) -> pi_m
dt2$pi_m = NA
for(k in 1:K){
  dt2[dt2$R_L==k, "pi_m"] = pi_m[k]
}

##############
# The group optimal design
rep(pi_m, each = ceiling(N1_2/K)) -> pi_1_0
pi_1_0[1:N1_2] -> dt1_2$pi_1_0


##############
#Outcome dependent group 
te <- c()
for(k in 1:K){
  te[k] <- mean((dt2[dt2$R_O==k,"phi_hat"]-phi_hat_bar)^2)
}
te2 <- te
getpi_stra(sqrt(te2), N1_2, n, K, p = rep(1/K, K), alpha) -> pi_O

pr_O <- data.frame(R_O = 1:K, p_O = pi_O)
merge(dt1_2, pr_O, by="R_O") -> dt1_2_O
merge(dt2, pr_O, by="R_O") -> dt2_O

#######
# M1 model
te <- c()
for(k in 1:K){
  te[k] <- mean((dt2[dt2$R_M==k,"phi_hat"]-phi_hat_bar)^2)
}
te2 <- te
getpi_stra(sqrt(te2), N1_2, n, K, p = rep(1/K, K), alpha) -> p_M
pr_M <- data.frame(R_M = 1:K, p_M = p_M)
merge(dt1_2, pr_M, by="R_M") -> dt1_2_M1

#######
# M2 model
pr_M <- data.frame(R_M = 1:K, p_M = pi_m)
merge(dt1_2, pr_M, by="R_M") -> dt1_2_M2
###########
# Z model
digit.x = dt2[, c("X1", "Y", "Z")]
digit.y = dt2[, "R_L"]
digit.ml <- train(x=digit.x, y=digit.y,
                  method="nnet",
                  tuneGrid=expand.grid(
                    # 5个隐藏神经元
                    .size=10,
                    # 衰变率
                    .decay=0.1
                  ),
                  trControl=trainControl(method="none"),
                  # 最大权重数量
                  MaxNWts=10000,
                  # 最大迭代次数
                  maxit=500)
dt1_2$R_Z <- predict(digit.ml, newdata = dt1_2)
pr_Z <- data.frame(R_Z = 1:K, p_Z = pi_m)
merge(dt1_2, pr_Z, by="R_Z") -> dt1_2_Z

##########
# ZM model
digit.x = dt2[, c("X1", "Y", "Z", "phi_m")]
digit.y = dt2[, "R_L"]
digit.zm <- train(x=digit.x, y=digit.y,
                  method="nnet",
                  tuneGrid=expand.grid(
                    # 5个隐藏神经元
                    .size=10,
                    # 衰变率
                    .decay=0.1
                  ),
                  trControl=trainControl(method="none"),
                  # 最大权重数量
                  MaxNWts=10000,
                  # 最大迭代次数
                  maxit=500)
dt1_2$R_ZM <- predict(digit.zm, newdata = dt1_2)
pr_ZM <- data.frame(R_ZM = 1:K, p_ZM = pi_m)
merge(dt1_2, pr_ZM, by="R_ZM") -> dt1_2_ZM


#####################
#The random foest method for classification: F method
fit_f <- randomForest(R_L~X1+Y+Z , data = dt2, proximity = T)
dt1_2$R_F <- predict(fit_f, newdata = dt1_2)
pr_F <- data.frame(R_F = 1:K, p_F = pi_m)
merge(dt1_2, pr_F, by="R_F") -> dt1_2_F

######################
##The random foest method for classification: FM method
fit_fm <- randomForest(R_L~X1+Y+Z+phi_m, data = dt2, proximity = T)
dt1_2$R_FM <- predict(fit_fm, newdata = dt1_2)
pr_FM <- data.frame(R_FM = 1:K, p_FM = pi_m)
merge(dt1_2, pr_FM, by="R_FM") -> dt1_2_FM

# ######
sum((dt1_2$phi_a)^2/(n/(N-n_m)))/N/N 
sum((dt1_2_FM$phi_a)^2/dt1_2_FM$pi_b)/N/N
sum((dt1_2_O$phi_a)^2/dt1_2_O$p_O)/N/N
sum((dt1_2_FM$phi_a)^2/dt1_2_FM$p_FM)/N/N
# sum((dt1_2_F$phi_a)^2/dt1_2_F$p_F)/N/N

low = min(dt1_2_O$phi)-1
high = max(dt1_2_O$phi)+1
breaks = seq(low, high, by = (high-low)/50)
dt1_2_O$V_b = rbinom(length(dt1_2_O$pi_b), 1, prob = dt1_2_O$pi_b)
dt1_2_O$V_u = rbinom(length(dt1_2_O$pi_u), 1, prob = dt1_2_O$pi_u)
dt1_2_O$V_O = rbinom(length(dt1_2_O$p_O), 1, prob = dt1_2_O$p_O)
dt1_2_O$V_1_0 = rbinom(length(dt1_2_O$pi_1_0), 1, prob = dt1_2_O$pi_1_0)
dt1_2_FM$V_FM = rbinom(length(dt1_2_FM$p_FM), 1, prob = dt1_2_FM$p_FM)

library(ggplot2)
theme_clean <- theme_minimal() +
  theme(
    plot.background = element_blank(),
    axis.line = element_line(colour = "gray")
  )
# my_theme <- theme_minimal(base_size = 12) +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold"),
#     legend.position = "bottom",
#     panel.grid.minor = element_blank()
#   )
# 定义一组新的颜色
fill_color <- "#AED6F1"
hist_border_color <- "#3498DB"
density_line_color <- "IndianRed"

gg_1 <- ggplot(dt1_2_O[dt1_2_O$V_u==1,], aes(x=phi)) +
  geom_histogram(aes(y=..density..), color=hist_border_color, fill=fill_color, bins=length(breaks)) +
  geom_density(color=density_line_color, lwd=0.7) + 
  xlim(c(-100, 100)) +ylim(c(0, 0.027))+
  labs(title="SimRan", x = expression(psi(D, theta[0])), y = expression(paste("Density")))+theme_clean+
  annotate("text", x = -55, 
           y = 0.023, label = expression(Lambda[n]*(bold(pi)) == 0.114), 
           parse = TRUE)

gg_2 <- ggplot(dt1_2_O[dt1_2_O$V_b==1,], aes(x=phi)) +
  geom_histogram(aes(y=..density..), color=hist_border_color, fill=fill_color, bins=length(breaks)) +
  geom_density(color=density_line_color, lwd=0.7) + 
  xlim(c(-100, 100)) +ylim(c(0, 0.027))+
  labs(title="Oracle", x = expression(psi(D, theta[0])), y = expression(paste("Density")))+theme_clean+
  annotate("text", x = -55, 
           y = 0.023, label = expression(Lambda[n]*(bold(pi)) == 0.077), 
           parse = TRUE) 

gg_3 <- ggplot(dt1_2_O[dt1_2_O$V_O==1,], aes(x=phi)) +
  geom_histogram(aes(y=..density..), color=hist_border_color, fill=fill_color, bins=length(breaks)) +
  geom_density(color=density_line_color, lwd=0.7) + 
  xlim(c(-100, 100)) + ylim(c(0, 0.027))+
  labs(title="FixStrat", x = expression(psi(D, theta[0])), y = expression(paste("Density")))+theme_clean+
  annotate("text", x = -55, 
           y = 0.023, label = expression(Lambda[n]*(bold(pi)) == 0.101), 
           parse = TRUE)

gg_4 <- ggplot(dt1_2_O[dt1_2_O$V_1_0==1,], aes(x=phi)) +
  geom_histogram(aes(y=..density..), color=hist_border_color, fill=fill_color, bins=length(breaks)) +
  geom_density(color=density_line_color, lwd=0.7) +
  xlim(c(-100, 100)) +ylim(c(0, 0.027))+
  labs(title="AdaStrat", x = expression(psi(D, theta[0])), y = expression(paste("Density")))+theme_clean+
  annotate("text", x = -55, 
           y = 0.023, label = expression(Lambda[n]*(bold(pi)) == 0.089), 
           parse = TRUE)

# 使用gridExtra包整合各个图形输出为一个2x2的图表
library(gridExtra)
pdf(file = paste0("result_figure/illustrate_p1/density.pdf"),
    onefile = F, width = 7, height = 6)
grid.arrange(gg_1, gg_2, gg_3, gg_4, ncol=2)
dev.off()



data1 <- dt1_2_O$phi[dt1_2_O$V_u==1]
data2 <- dt1_2_O$phi[dt1_2_O$V_b==1]
data3 <- dt1_2_O$phi[dt1_2_O$V_O==1]
data4 <- dt1_2_FM$phi[dt1_2_FM$V_FM==1]

# 创建数据框
df <- data.frame(Value = c(data1, data2, data3, data4),
                 Type = factor(c(rep('SRS', length(data1)), rep('ORACALE', length(data2)),
                                 rep('FS-ODS', length(data3)), rep('AS-CC', length(data4)))))

# 加载ggplot2库
library(ggplot2)

# 创建颜色向量
colors <- c('SRS' = 'red', 'ORACALE' = 'blue', 'FS-ODS' = 'green', 'AS-CC' = 'purple')

# 生成图像
ggplot(df, aes(x = Value, color = Type)) +
  geom_density() +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = 'Value', y = 'Density', color = 'Type')+xlim(-100, 100)





# 
# # install.packages("rpart")
# library(rpart)
# fit_tr <- rpart(R_L~X1+Y+Z, data = dt2, method = "class")
# dt1_2$R_tr <- predict(fit_tr, newdata = dt1_2, type = "class")
# pr_tr <- data.frame(R_tr = 1:K, p_tr = pi_m)
# merge(dt1_2, pr_tr, by="R_tr") -> dt1_2_tr
# table(dt1_2_tr$R_tr, dt1_2_tr$R_T)
# sum((dt1_2_tr$phi_a)^2/dt1_2_tr$p_tr)/N/N
# par(cex = 0.8)  # 调整文本大小
# plot(fit_tr, uniform = TRUE, compress = TRUE, box.palette = c("green", "blue"), branch = 0.8)
# text(fit_tr, use.n = TRUE, cex = 0.8, font = 2) 
# 
# 
# fit_trm <- rpart(R_L~X1+Y+Z+phi_m, data = dt2, method = "class")
# dt1_2$R_trm <- predict(fit_trm, newdata = dt1_2, type = "class")
# pr_trm <- data.frame(R_trm = 1:K, p_trm = pi_m)
# merge(dt1_2, pr_trm, by="R_trm") -> dt1_2_trm
# table(dt1_2_trm$R_trm, dt1_2_trm$R_T)
# sum((dt1_2_trm$phi_a)^2/dt1_2_trm$p_trm)/N/N
# par(cex = 0.8)  # 调整文本大小
# plot(fit_trm, uniform = TRUE, compress = TRUE, box.palette = c("green", "blue"), branch = 0.8)
# text(fit_trm, use.n = TRUE, cex = 0.8, font = 2) 
# 
# # pdf(file = "result_figure/illustrate_p1/phi.pdf", width = 12, height = 4)
# # par(mfrow=c(1,4))
# # hist(dt1_2_O$phi, breaks = breaks, freq = F, xlim = c(-50, 50), xlab = "phi", main = "USD")
# # hist(dt1_2_O$phi[dt1_2_O$V_u==1], freq = F, breaks = breaks, add = T, col = "red")
# # hist(dt1_2_O$phi, breaks = breaks, freq = F, xlim = c(-50, 50), xlab = "phi", main = "Orcale")
# # hist(dt1_2_O$phi[dt1_2_O$V_b==1], freq = F, breaks = breaks, add = T, col = "red")
# # hist(dt1_2_O$phi, breaks = breaks, freq = F, xlim = c(-50, 50), xlab = "phi", main = "Ks-fs-ODS")
# # hist(dt1_2_O$phi[dt1_2_O$V_O==1], freq = F,  breaks = breaks, add = T, col = "red")
# # hist(dt1_2_F$phi, breaks = breaks, freq = F, xlim = c(-50, 50), xlab = "phi", main = "Ks-as-CC")
# # hist(dt1_2_FM$phi[dt1_2_FM$V_FM==1], freq = F,  breaks = breaks, add = T, col = "red")
# 
# # 随机生成一些数据
# par(mfrow=c(2,2))
# # 生成直方图
# hist(dt1_2_O$phi[dt1_2_O$V_u==1], freq=FALSE, main="SRS", xlab="Data", breaks = breaks,
#      ylab="Density", border="black", col="lightblue", xlim = c(-100, 100))
# # 生成密度曲线
# lines(density(dt1_2_O$phi[dt1_2_O$V_u==1]), col="red", lwd=2)
# 
# hist(dt1_2_O$phi[dt1_2_O$V_b==1], freq=FALSE, main="ORCALE", xlab="Data", breaks = breaks,
#      ylab="Density", border="black", col="lightblue", xlim = c(-100, 100))
# # 生成密度曲线
# lines(density(dt1_2_O$phi[dt1_2_O$V_b==1]), col="red", lwd=2)
# 
# hist(dt1_2_O$phi[dt1_2_O$V_O==1], freq=FALSE, main="FS-ODS", xlab="Data", breaks = breaks,
#      ylab="Density", border="black", col="lightblue", xlim = c(-100, 100))
# # 生成密度曲线
# lines(density(dt1_2_O$phi[dt1_2_O$V_O==1]), col="red", lwd=2)
# 
# hist(dt1_2_FM$phi[dt1_2_FM$V_FM==1], freq=FALSE, main="AS-CC", xlab="Data", breaks = breaks,
#      ylab="Density", border="black", col="lightblue", xlim = c(-100, 100))
# # 生成密度曲线
# lines(density(dt1_2_FM$phi[dt1_2_FM$V_FM==1]), col="red", lwd=2)





library(ggplot2)
library(ggpubr)

# 设置统一主题

# 颜色设置
strata_colors <- scale_color_brewer(palette = "Set1", name = "Stratum")

tt <- sample(nrow(dt1_2), 2000, replace = F)
# 第一组图形：散点图
p1 <- ggplot(dt1_2_FM[tt,], aes(x = phi, y = Y, color = factor(R_FM))) +
  geom_point(size = 3, alpha = 0.7) +
  labs(
    title = "AdaStrat",
    x = expression(psi(D, theta[0])), 
    y = "Y"
  ) +
  strata_colors +
  theme_clean +
  guides(color = guide_legend(nrow = 1))+
  annotate("text", x = -100, 
           y = 20, label = expression(Lambda[n]*(bold(pi)) == 0.089), 
           parse = TRUE)

p2 <- ggplot(dt1_2_FM[tt,], aes(x = phi, y = Y, color = factor(R_O))) +
  geom_point(size = 3, alpha = 0.7) +
  labs(
    title = "FixStrat",
    x = expression(psi(D, theta[0])), 
    y = "Y"
  ) +
  strata_colors +
  theme_clean +
  guides(color = guide_legend(nrow = 1)) + 
  annotate("text", x = -100, 
  y = 20, label = expression(Lambda[n]*(bold(pi)) == 0.101), 
  parse = TRUE)

# 第二组图形：分层分布图
p3 <- ggplot(dt1_2_FM[tt,], aes(x = factor(R_FM), y = phi)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.4, size = 2, aes(color = factor(R_FM))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  coord_cartesian(ylim = c(-200, 200)) +  # 比ylim更友好，不丢弃数据
  labs(
    title = "AdaStrat",
    x = expression(psi(D, theta[0])), 
    y = "Y"
  ) +
  strata_colors +
  theme_clean +
  theme(legend.position = "none")+
  annotate("text", x = 1, 
           y = 150, label = expression(Lambda[n]*(bold(pi)) == 0.089), 
           parse = TRUE)

p4 <- ggplot(dt1_2_O[tt,], aes(x = factor(R_O), y = phi)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.4, size = 2, aes(color = factor(R_O))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  coord_cartesian(ylim = c(-200, 200)) +
  labs(
    title = "FixStrat",
    x = expression(psi(D, theta[0])), 
    y = "Y"
  ) +
  strata_colors +
  theme_clean +
  theme(legend.position = "none") + 
  annotate("text", x = 1, 
           y = 150, label = expression(Lambda[n]*(bold(pi)) == 0.101), 
           parse = TRUE)

# 输出图形
pdf("result_figure/illustrate_p1/points.pdf", width = 10, height = 5)
ggarrange(p1, p2, nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom")
dev.off()

pdf("result_figure/illustrate_p1/strata.pdf", width = 10, height = 5)
ggarrange(p3, p4, nrow = 1, ncol = 2)
dev.off()




p5 <- ggplot(dt1_2_FM[tt,], aes(x = X1, y = Y, color = factor(R_T))) +
  geom_point(size = 3, alpha = 0.7) +
  labs(
    title = "AdaStrat",
    x = expression(X[e]), 
    y = "Y"
  ) +
  strata_colors +
  theme_clean +
  guides(color = guide_legend(nrow = 1))

p6 <- ggplot(dt1_2_FM[tt,], aes(x = X1, y = Y, color = factor(R_O))) +
  geom_point(size = 3, alpha = 0.7) +
  labs(
    title = "FixStrat",
    x = expression(X[e]), 
    y = "Y"
  ) +
  strata_colors +
  theme_clean +
  guides(color = guide_legend(nrow = 1)) 

pdf("result_figure/illustrate_p1/XY.pdf", width = 10, height = 5)
ggarrange(p5, p6, nrow = 1, ncol = 2)
dev.off()



# 合并数据集
merge(dt1_2_FM, dt1_2_O[, c("id", "p_O")], by = "id") -> aa

# 设置统一的颜色范围（基于所有颜色变量的范围）
color_min <- min(c(aa$p_FM, aa$pi_b, aa$p_O, aa$pi_u), na.rm = TRUE)
color_max <- max(c(aa$p_FM, aa$pi_b, aa$p_O, aa$pi_u), na.rm = TRUE)

# 使用改进的渐变色条图例
common_color_scale <- scale_color_gradient2(
  low = "darkblue",
  mid = "gray90",
  high = "darkred",
  midpoint = (color_min + 0.36)/2,
  limits = c(color_min, 0.36),
  name = expression(pi[i]),  # 数学符号表示倾向得分
  oob = scales::squish,
  guide = guide_colorbar(  # 优化后的色条图例
    direction = "vertical",  # 水平方向
    barwidth = unit(0.5, "cm"),  # 增加色条长度
    barheight = unit(5, "cm"),  # 增加高度提升可见性
    title.position = "top",  # 标题在色条上方
    title.hjust = 0.5,  # 标题居中
    ticks = TRUE,  # 显示刻度
    ticks.colour = "black",  # 黑色刻度线增强对比度
    frame.colour = "black"  # 添加边框增强区分度
  )
)

# 第一组图形：散点图（保持原样）
p1 <- ggplot(aa[tt, ], aes(x = phi, y = Y, color = pi_1_0)) +
  geom_point(size = 0.1, alpha = 0.8, shape = 6) +
  labs(
    title = "AdaStrat",
    x = expression(psi[i]), 
    y = expression(Y[i])
  ) +
  common_color_scale +
  theme_clean +
  annotate("text", 
           x = -55, 
           y = 20, 
           label = expression(Lambda[n]*(bold(pi)) == 0.089), 
           parse = TRUE)

p2 <- ggplot(aa[tt, ], aes(x = phi, y = Y, color = pi_b)) +
  geom_point(size = 0.1, alpha = 0.8, shape = 6) +
  labs(
    title = "Oracle",
    x = expression(psi[i]), 
    y = expression(Y[i])
  ) +
  common_color_scale +
  theme_clean +
  annotate("text", 
           x = -55, 
           y = 20, 
           label = expression(Lambda[n]*(bold(pi)) == 0.077), 
           parse = TRUE)

p3 <- ggplot(aa[tt, ], aes(x = phi, y = Y, color = p_O)) +
  geom_point(size = 0.1, alpha = 0.8, shape = 6) +
  labs(
    title = "FixStra",
    x = expression(psi[i]), 
    y = expression(Y[i])
  ) +
  common_color_scale +
  theme_clean +
  annotate("text", 
           x = -55, 
           y = 20, 
           label = expression(Lambda[n]*(bold(pi)) == 0.101), 
           parse = TRUE)

p4 <- ggplot(aa[tt, ], aes(x = phi, y = Y, color = pi_u)) +
  geom_point(size = 0.1, alpha = 0.8, shape = 6) +
  labs(
    title = "SimRan",
    x = expression(psi[i]), 
    y = expression(Y[i])
  ) +
  common_color_scale +
  theme_clean +
  annotate("text", 
           x = -55, 
           y = 20, 
           label = expression(Lambda[n]*(bold(pi)) == 0.114), 
           parse = TRUE)

# 输出图形（确保图例位置正确）
pdf("result_figure/illustrate_p1/points_p_color.pdf", width = 7, height = 6)
ggarrange(p4, p2, p3, p1, 
          nrow = 2, ncol = 2, 
          common.legend = TRUE, 
          legend = "right")  # 明确指定图例在底部
dev.off()
