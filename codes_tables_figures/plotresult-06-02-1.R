## Plot results 
# Simulations: EDR methods,
# Appendix A.3 (Figure S1, Table S1)

rm(list = ls())
library(dplyr)
library(ggplot2)
name = "06-02-1test"
if(! dir.exists(paste0("result/",name))){
  dir.create(paste0("result/",name))
}


N_set = c(10000, 20000)
n_m_set = c(1000)
n_set = c(1000, 1200, 1500, 2000, 4000)
K_set = c(6, 8, 10)
alpha_set = c(0.3)

expand.grid(N_set, n_m_set, n_set,K_set, alpha_set)->settings
colnames(settings) <- c("N", "n_m", "n", "K", "alpha")

settings[order(settings$N, settings$n_m, settings$n, 
               settings$K, settings$alpha), ] -> settings
rownames(settings) <- c(1:nrow(settings))


library(stringr)
isalpha <- function(str){
  str <- ifelse(grepl("[a-zA-Z]", str), NA, str) %>% unlist %>% as.numeric()
  return(str)
}
result_td <- function(re){
  re = cbind((apply(re, 2, mean, na.rm = T)+0.61)/(-0.61) *100,
             apply(re, 2, var, na.rm = T),
             (apply(re, 2, mean, na.rm = T)+0.61)^2 +apply(re, 2, var, na.rm = T) )
  colnames(re) <- c("bias", "var", "mse")
  return(re)
}




###################
if(! dir.exists(paste0("result_figure/",name))){
  dir.create(paste0("result_figure/",name))
}
r_td <- c()
for(j in 1:nrow(settings)){
  N = settings$N[j]
  n_m = settings$n_m[j]
  n = settings$n[j]
  K = settings$K[j]
  alpha = settings$alpha[j]
  rrr <- c()
  path = paste0("result/", name, "/", paste(N, n_m, n, K, alpha, sep = "_"), ".rdata")
  # print(path)
  load(file = path)
  apply(result, 1, isalpha)%>%t() -> result
  result_td(result[, 1:4]) %>% as.data.frame() -> rr
  cbind(rr,
        var_pot=colMeans(result[ ,5:8], na.rm = T),
        N = rep(N, 4),
        n_m = rep(n_m,4),
        n = rep(n, 4),
        K = rep(K, 4),
        alpha = rep(alpha, 4),
        des = c("u", "z", "f", "o"))->rrr
  colnames(rrr) <- c("Bias", "Var", "MSE",
                     "Var_est", "N",
                     "n_m", "n", "K", "alpha","des")
  rbind(r_td, rrr) -> r_td
}







########################
# Plot 
#######################
n_m_i = 1000; K_i = 6; N_i = 10000
r_td %>% filter(n_m == n_m_i, K==K_i, N == N_i,
                des %in%c("u", "f", "o")) %>%
  ggplot(aes(x = n, y = MSE, color = des, shape = des, linetype = des)) +
  geom_line()+
  geom_point(size = 3) +
  scale_color_manual(values = c("u" = "#66c2a5", "f" = "#8da0cb", "o" = "#e78ac3"),
                     labels = c("u" ="SimRan"  , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_shape_manual(values = c("u" = 15, "f" = 17, "o" = 18),
                     labels = c("u" ="SimRan"  , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_linetype_manual(values = c("u" = 15, "f" = 17, "o" = 18),
                        labels = c("u" ="SimRan"  , "f" ="AdaStrat" , "o" ="FixStrat"))+
  labs(color = "Sampling design", shape = "Sampling design", linetype = "Sampling design", 
       x = expression(tilde(n)[S]), 
       y = expression("MSE of "~hat(mu[1])))+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"))+
  ggtitle("(a)")-> p1

n_m_i = 1000; K_i = 8; N_i = 10000
r_td %>% filter(n_m == n_m_i, K==K_i, N == N_i,
                des %in%c("u", "f", "o")) %>%
  ggplot(aes(x = n, y = MSE, color = des, shape = des, linetype = des)) +
  geom_line()+
  geom_point(size = 3) +
  scale_color_manual(values = c("u" = "#66c2a5", "f" = "#8da0cb", "o" = "#e78ac3"),
                     labels = c("u" ="SimRan"  , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_shape_manual(values = c("u" = 15, "f" = 17, "o" = 18),
                     labels = c("u" ="SimRan"  , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_linetype_manual(values = c("u" = 15, "f" = 17, "o" = 18),
                        labels = c("u" ="SimRan"  , "f" ="AdaStrat" , "o" ="FixStrat"))+
  labs(color = "Sampling design", shape = "Sampling design", linetype = "Sampling design", 
       x = expression(tilde(n)[S]), 
       y = expression("MSE of "~hat(mu[1])))+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"))+
  ggtitle("(b)")-> p2


n_m_i = 1000; K_i = 10; N_i = 10000
r_td %>% filter(n_m == n_m_i, K==K_i, N == N_i,
                des %in%c("u", "f", "o")) %>%
  ggplot(aes(x = n, y = MSE, color = des, shape = des, linetype = des)) +
  geom_line()+
  geom_point(size = 3) +
  scale_color_manual(values = c("u" = "#66c2a5", "f" = "#8da0cb", "o" = "#e78ac3"),
                     labels = c("u" ="SimRan"  , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_shape_manual(values = c("u" = 15, "f" = 17, "o" = 18),
                     labels = c("u" ="SimRan"  , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_linetype_manual(values = c("u" = 15, "f" = 17, "o" = 18),
                        labels = c("u" ="SimRan"  , "f" ="AdaStrat" , "o" ="FixStrat"))+
  labs(color = "Sampling design", shape = "Sampling design", linetype = "Sampling design", 
       x = expression(tilde(n)[S]), 
       y = expression("MSE of "~hat(mu[1])))+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"))+
  ggtitle("(c)")-> p3

n_m_i = 1000; K_i = 6; N_i = 20000
r_td %>% filter(n_m == n_m_i, K==K_i, N == N_i,
                des %in%c("u", "f", "o")) %>%
  ggplot(aes(x = n, y = MSE, color = des, shape = des, linetype = des)) +
  geom_line()+
  geom_point(size = 3) +
  scale_color_manual(values = c("u" = "#66c2a5", "f" = "#8da0cb", "o" = "#e78ac3"),
                     labels = c("u" ="SimRan"  , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_shape_manual(values = c("u" = 15, "f" = 17, "o" = 18),
                     labels = c("u" ="SimRan"  , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_linetype_manual(values = c("u" = 15, "f" = 17, "o" = 18),
                        labels = c("u" ="SimRan"  , "f" ="AdaStrat" , "o" ="FixStrat"))+
  labs(color = "Sampling design", shape = "Sampling design", linetype = "Sampling design", 
       x = expression(tilde(n)[S]), 
       y = expression("MSE of "~hat(mu[1])))+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"))+
  ggtitle("(d)")-> p4


n_m_i = 1000; K_i = 8; N_i = 20000
r_td %>% filter(n_m == n_m_i, K==K_i, N == N_i,
                des %in%c("u", "f", "o")) %>%
  ggplot(aes(x = n, y = MSE, color = des, shape = des, linetype = des)) +
  geom_line()+
  geom_point(size = 3) +
  scale_color_manual(values = c("u" = "#66c2a5", "f" = "#8da0cb", "o" = "#e78ac3"),
                     labels = c("u" ="SimRan"  , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_shape_manual(values = c("u" = 15, "f" = 17, "o" = 18),
                     labels = c("u" ="SimRan"  , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_linetype_manual(values = c("u" = 15, "f" = 17, "o" = 18),
                        labels = c("u" ="SimRan"  , "f" ="AdaStrat" , "o" ="FixStrat"))+
  labs(color = "Sampling design", shape = "Sampling design", linetype = "Sampling design", 
       x = expression(tilde(n)[S]), 
       y = expression("MSE of "~hat(mu[1])))+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"))+
  ggtitle("(e)")-> p5

n_m_i = 1000; K_i = 10; N_i = 20000
r_td %>% filter(n_m == n_m_i, K==K_i, N == N_i,
                des %in%c("u", "f", "o")) %>%
  ggplot(aes(x = n, y = MSE, color = des, shape = des, linetype = des)) +
  geom_line()+
  geom_point(size = 3) +
  scale_color_manual(values = c("u" = "#66c2a5", "f" = "#8da0cb", "o" = "#e78ac3"),
                     labels = c("u" ="SimRan"  , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_shape_manual(values = c("u" = 15, "f" = 17, "o" = 18),
                     labels = c("u" ="SimRan"  , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_linetype_manual(values = c("u" = 15, "f" = 17, "o" = 18),
                        labels = c("u" ="SimRan"  , "f" ="AdaStrat" , "o" ="FixStrat"))+
  labs(color = "Sampling design", shape = "Sampling design", linetype = "Sampling design", 
       x = expression(tilde(n)[S]), 
       y = expression("MSE of "~hat(mu[1])))+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"))+
  ggtitle("(f)") -> p6



plots <- list(p1, p2, p3, p4, p5, p6)
library(ggpubr)
pdf(file = paste0("result_figure/", name, "/", name,".pdf"),
    onefile = F, width = 7.5, height = 5)
ggarrange(p1, p2, p3, p4, p5, p6, 
          nrow = 2, ncol = 3, common.legend = TRUE)
dev.off()
# 使用 annotate_figure 添加字母标签


#############################