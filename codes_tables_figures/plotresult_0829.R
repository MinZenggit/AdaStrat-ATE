###############
# simulation: main results in Figure 3
# ## Plot results 
# ## Simulation results, but with legends revised
# # 修改method_set, 跑多次; "ipw", "reg"的ylim都要调整的
# 
# rm(list = ls())
# library(dplyr)
# library(ggplot2)
# name = "05-28-1test"
# method_set <- c("ipw", "reg", "aipw")
# N_set = c(10000)
# n_m_set = c(500, 1000)
# n_set = c(1000, 1200, 1500, 2000, 4000)
# K_set = c(2, 4, 6, 8, 10)
# alpha_set = c(0.05, 0.1, 0.2, 0.3)
# 
# expand.grid(method_set, N_set, n_m_set, n_set,K_set, alpha_set)->settings
# colnames(settings) <- c("method", "N", "n_m", "n", "K", "alpha")
# settings[order(settings$method, settings$N, 
#                settings$n_m, settings$n, 
#                settings$K, settings$alpha), ] -> settings
# rownames(settings) <- c(1:nrow(settings))
# 
# 
# library(stringr)
# isalpha <- function(str){
#   str <- ifelse(grepl("[a-zA-Z]", str), NA, str) %>% unlist %>% as.numeric()
#   return(str)
# }
# result_td <- function(re){
#   re = cbind((apply(re, 2, mean, na.rm = T)-0.50)/0.50 *100,
#              apply(re, 2, var, na.rm = T),
#              (apply(re, 2, mean, na.rm = T)-0.50)^2 +apply(re, 2, var, na.rm = T) )
#   colnames(re) <- c("Bias", "var", "mse")
#   return(re)
# }
# 
# 
# 
# 
# if(! dir.exists(paste0("result_figure/",name))){
#   dir.create(paste0("result_figure/",name))
# }
# r_td <- c()
# for(j in 1:nrow(settings)){
#   method = settings$method[j]
#   N = settings$N[j]
#   n_m = settings$n_m[j]
#   n = settings$n[j]
#   K = settings$K[j]
#   alpha = settings$alpha[j]
#   rrr <- c()
#   path = paste0("result/", name, "/", paste(N, n_m, n, K, alpha, method, sep = "_"), ".rdata")
#   # print(path)
#   load(file = path)
#   apply(result, 1, isalpha)%>%t() -> result
#   result_td(result[,  c(1,3,5,7,9,11, 13, 15, 17)]) %>% as.data.frame() -> rr
#   cbind(rr,
#         var_pot=colMeans(result[,c(1,3,5,7,9,11, 13, 15, 17)+1], na.rm = T),
#         method = rep(method, 9),
#         N = rep(N, 9),
#         n_m = rep(n_m,9),
#         n = rep(n, 9),
#         K = rep(K, 9),
#         alpha = rep(alpha, 9),
#         des = c("u", "b", "m1", "m2","z", "zm","f", "fm", "o"))->rrr
#   colnames(rrr) <- c("Bias", "Var", "MSE",
#                      "Var_est",
#                      "method", "N", "n_m", "n", "K", "alpha","des")
#   rbind(r_td, rrr) -> r_td
# }
# 
# r_td$des %>% factor(level = c("u", "b", "o", "fm", "f", "m1", "m2", "z", "zm")) -> r_td$des
# save(r_td, file = paste0("result_figure/", name, "/r_td.rdata"))


########################
# Plot 
#######################
library(dplyr)
library(ggplot2)
name = "05-28-1test"
alpha_i = 0.3; n_m_i = 1000; 
load(file = paste0("result_figure/", name, "/r_td.rdata"))

###########
# AIPW
K_i = 2; method_i = "aipw"
r_td %>% filter(method == method_i,n_m == n_m_i, K==K_i, alpha == alpha_i,
                des %in%c("u", "b", "o", "fm")) %>%
  ggplot(aes(x = n, y = MSE, color = des, shape = des, linetype = des)) +
  geom_line()+
  geom_point(size = 2.5) +
  scale_color_manual(values = c("u" = "#66c2a5", "b" = "#fc8d62", "fm" = "#8da0cb", "o" = "#e78ac3"),
                     labels = c("u" ="SimRan"  , "b" = "Oracle" , "fm" ="AdaStrat" , "o" ="FixStrat")) +
  scale_shape_manual(values = c("u" = 15, "b" = 16, "fm" = 17, "o" = 18),
                     labels = c("u" ="SimRan"  , "b" = "Oracle" , "fm" ="AdaStrat" , "o" ="FixStrat")) +
  scale_linetype_manual(values = c("u" = 15, "b" = 16, "fm" = 17, "o" = 18),
                        labels = c("u" ="SimRan"  , "b" = "Oracle" , "fm" ="AdaStrat" , "o" ="FixStrat"))+
  ylim(c(0.02, 0.18))+
  labs(color = "", shape = "", linetype = "",
       x = expression(tilde(n)[S]), 
       y = "MSE")+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"),
        legend.text = element_text(size = 12))+
  ggtitle("(a1)")-> p1

K_i = 6; method_i = "aipw"
r_td %>% filter(method == method_i,n_m == n_m_i, K==K_i, alpha == alpha_i,
                des %in%c("u", "b", "o", "fm")) %>%
  ggplot(aes(x = n, y = MSE, color = des, shape = des, linetype = des)) +
  geom_line()+
  geom_point(size = 2.5) +
  scale_color_manual(values = c("u" = "#66c2a5", "b" = "#fc8d62", "fm" = "#8da0cb", "o" = "#e78ac3"),
                     labels = c("u" ="SimRan"  , "b" = "Oracle" , "fm" ="AdaStrat" , "o" ="FixStrat")) +
  scale_shape_manual(values = c("u" = 15, "b" = 16, "fm" = 17, "o" = 18),
                     labels = c("u" ="SimRan"  , "b" = "Oracle" , "fm" ="AdaStrat" , "o" ="FixStrat")) +
  scale_linetype_manual(values = c("u" = 15, "b" = 16, "fm" = 17, "o" = 18),
                        labels = c("u" ="SimRan"  , "b" = "Oracle" , "fm" ="AdaStrat" , "o" ="FixStrat"))+
  ylim(c(0.02, 0.18))+
  labs(color = "", shape = "", linetype = "", 
       x = expression(tilde(n)[S]), 
       y = "MSE")+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"),
        legend.text = element_text(size = 12))+
  ggtitle("(a2)")-> p2

K_i = 10; method_i = "aipw"
r_td %>% filter(method == method_i,n_m == n_m_i, K==K_i, alpha == alpha_i,
                des %in%c("u", "b", "o", "fm")) %>%
  ggplot(aes(x = n, y = MSE, color = des, shape = des, linetype = des)) +
  geom_line()+
  geom_point(size = 2.5) +
  scale_color_manual(values = c("u" = "#66c2a5", "b" = "#fc8d62", "fm" = "#8da0cb", "o" = "#e78ac3"),
                     labels = c("u" ="SimRan"  , "b" = "Oracle" , "fm" ="AdaStrat" , "o" ="FixStrat")) +
  scale_shape_manual(values = c("u" = 15, "b" = 16, "fm" = 17, "o" = 18),
                     labels = c("u" ="SimRan"  , "b" = "Oracle" , "fm" ="AdaStrat" , "o" ="FixStrat")) +
  scale_linetype_manual(values = c("u" = 15, "b" = 16, "fm" = 17, "o" = 18),
                        labels = c("u" ="SimRan"  , "b" = "Oracle" , "fm" ="AdaStrat" , "o" ="FixStrat"))+
  ylim(c(0.02, 0.18))+
  labs(color = "", shape = "", linetype = "", 
       x = expression(tilde(n)[S]), 
       y = "MSE")+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"),
        legend.text = element_text(size = 12))+
  ggtitle("(a3)")-> p3





###########
# IPW
K_i = 2; method_i = "ipw"
r_td %>% filter(method == method_i,n_m == n_m_i, K==K_i, alpha == alpha_i,
                des %in%c("u", "b", "o", "fm")) %>%
  ggplot(aes(x = n, y = MSE, color = des, shape = des, linetype = des)) +
  geom_line()+
  geom_point(size = 2.5) +
  scale_color_manual(values = c("u" = "#66c2a5", "b" = "#fc8d62", "fm" = "#8da0cb", "o" = "#e78ac3"),
                     labels = c("u" ="SimRan"  , "b" = "Oracle" , "fm" ="AdaStrat" , "o" ="FixStrat")) +
  scale_shape_manual(values = c("u" = 15, "b" = 16, "fm" = 17, "o" = 18),
                     labels = c("u" ="SimRan"  , "b" = "Oracle" , "fm" ="AdaStrat" , "o" ="FixStrat")) +
  scale_linetype_manual(values = c("u" = 15, "b" = 16, "fm" = 17, "o" = 18),
                        labels = c("u" ="SimRan"  , "b" = "Oracle" , "fm" ="AdaStrat" , "o" ="FixStrat"))+
  ylim(c(0.05, 0.4))+
  labs(color = "", shape = "", linetype = "",
       x = expression(tilde(n)[S]), 
       y = "MSE")+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"),
        legend.text = element_text(size = 12))+
  ggtitle("(b1)")-> p4

K_i = 6; method_i = "ipw"
r_td %>% filter(method == method_i,n_m == n_m_i, K==K_i, alpha == alpha_i,
                des %in%c("u", "b", "o", "fm")) %>%
  ggplot(aes(x = n, y = MSE, color = des, shape = des, linetype = des)) +
  geom_line()+
  geom_point(size = 2.5) +
  scale_color_manual(values = c("u" = "#66c2a5", "b" = "#fc8d62", "fm" = "#8da0cb", "o" = "#e78ac3"),
                     labels = c("u" ="SimRan"  , "b" = "Oracle" , "fm" ="AdaStrat" , "o" ="FixStrat")) +
  scale_shape_manual(values = c("u" = 15, "b" = 16, "fm" = 17, "o" = 18),
                     labels = c("u" ="SimRan"  , "b" = "Oracle" , "fm" ="AdaStrat" , "o" ="FixStrat")) +
  scale_linetype_manual(values = c("u" = 15, "b" = 16, "fm" = 17, "o" = 18),
                        labels = c("u" ="SimRan"  , "b" = "Oracle" , "fm" ="AdaStrat" , "o" ="FixStrat"))+
  ylim(c(0.05, 0.4))+
  labs(color = "", shape = "", linetype = "", 
       x = expression(tilde(n)[S]), 
       y = "MSE")+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"),
        legend.text = element_text(size = 12))+
  ggtitle("(b1)")-> p5

K_i = 10; method_i = "ipw"
r_td %>% filter(method == method_i,n_m == n_m_i, K==K_i, alpha == alpha_i,
                des %in%c("u", "b", "o", "fm")) %>%
  ggplot(aes(x = n, y = MSE, color = des, shape = des, linetype = des)) +
  geom_line()+
  geom_point(size = 2.5) +
  scale_color_manual(values = c("u" = "#66c2a5", "b" = "#fc8d62", "fm" = "#8da0cb", "o" = "#e78ac3"),
                     labels = c("u" ="SimRan"  , "b" = "Oracle" , "fm" ="AdaStrat" , "o" ="FixStrat")) +
  scale_shape_manual(values = c("u" = 15, "b" = 16, "fm" = 17, "o" = 18),
                     labels = c("u" ="SimRan"  , "b" = "Oracle" , "fm" ="AdaStrat" , "o" ="FixStrat")) +
  scale_linetype_manual(values = c("u" = 15, "b" = 16, "fm" = 17, "o" = 18),
                        labels = c("u" ="SimRan"  , "b" = "Oracle" , "fm" ="AdaStrat" , "o" ="FixStrat"))+
  ylim(c(0.05, 0.4))+
  labs(color = "", shape = "", linetype = "", 
       x = expression(tilde(n)[S]), 
       y = "MSE")+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"),
        legend.text = element_text(size = 12))+
  ggtitle("(b3)")-> p6






###########
# reg
K_i = 2; method_i = "reg"
r_td %>% filter(method == method_i,n_m == n_m_i, K==K_i, alpha == alpha_i,
                des %in%c("u", "b", "o", "fm")) %>%
  ggplot(aes(x = n, y = MSE, color = des, shape = des, linetype = des)) +
  geom_line()+
  geom_point(size = 2.5) +
  scale_color_manual(values = c("u" = "#66c2a5", "b" = "#fc8d62", "fm" = "#8da0cb", "o" = "#e78ac3"),
                     labels = c("u" ="SimRan"  , "b" = "Oracle" , "fm" ="AdaStrat" , "o" ="FixStrat")) +
  scale_shape_manual(values = c("u" = 15, "b" = 16, "fm" = 17, "o" = 18),
                     labels = c("u" ="SimRan"  , "b" = "Oracle" , "fm" ="AdaStrat" , "o" ="FixStrat")) +
  scale_linetype_manual(values = c("u" = 15, "b" = 16, "fm" = 17, "o" = 18),
                        labels = c("u" ="SimRan"  , "b" = "Oracle" , "fm" ="AdaStrat" , "o" ="FixStrat"))+
  ylim(c(0.02, 0.18))+
  labs(color = "", shape = "", linetype = "",
       x = expression(tilde(n)[S]), 
       y = "MSE")+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"),
        legend.text = element_text(size = 12))+
  ggtitle("(c1)")-> p7

K_i = 6; method_i = "reg"
r_td %>% filter(method == method_i,n_m == n_m_i, K==K_i, alpha == alpha_i,
                des %in%c("u", "b", "o", "fm")) %>%
  ggplot(aes(x = n, y = MSE, color = des, shape = des, linetype = des)) +
  geom_line()+
  geom_point(size = 2.5) +
  scale_color_manual(values = c("u" = "#66c2a5", "b" = "#fc8d62", "fm" = "#8da0cb", "o" = "#e78ac3"),
                     labels = c("u" ="SimRan"  , "b" = "Oracle" , "fm" ="AdaStrat" , "o" ="FixStrat")) +
  scale_shape_manual(values = c("u" = 15, "b" = 16, "fm" = 17, "o" = 18),
                     labels = c("u" ="SimRan"  , "b" = "Oracle" , "fm" ="AdaStrat" , "o" ="FixStrat")) +
  scale_linetype_manual(values = c("u" = 15, "b" = 16, "fm" = 17, "o" = 18),
                        labels = c("u" ="SimRan"  , "b" = "Oracle" , "fm" ="AdaStrat" , "o" ="FixStrat"))+
  ylim(c(0.02, 0.18))+
  labs(color = "", shape = "", linetype = "", 
       x = expression(tilde(n)[S]), 
       y = "MSE")+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"),
        legend.text = element_text(size = 12))+
  ggtitle("(c2)")-> p8

K_i = 10; method_i = "reg"
r_td %>% filter(method == method_i,n_m == n_m_i, K==K_i, alpha == alpha_i,
                des %in%c("u", "b", "o", "fm")) %>%
  ggplot(aes(x = n, y = MSE, color = des, shape = des, linetype = des)) +
  geom_line()+
  geom_point(size = 2.5) +
  scale_color_manual(values = c("u" = "#66c2a5", "b" = "#fc8d62", "fm" = "#8da0cb", "o" = "#e78ac3"),
                     labels = c("u" ="SimRan"  , "b" = "Oracle" , "fm" ="AdaStrat" , "o" ="FixStrat")) +
  scale_shape_manual(values = c("u" = 15, "b" = 16, "fm" = 17, "o" = 18),
                     labels = c("u" ="SimRan"  , "b" = "Oracle" , "fm" ="AdaStrat" , "o" ="FixStrat")) +
  scale_linetype_manual(values = c("u" = 15, "b" = 16, "fm" = 17, "o" = 18),
                        labels = c("u" ="SimRan"  , "b" = "Oracle" , "fm" ="AdaStrat" , "o" ="FixStrat"))+
  ylim(c(0.02, 0.18))+
  labs(color = "", shape = "", linetype = "", 
       x = expression(tilde(n)[S]), 
       y = "MSE")+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"),
        legend.text = element_text(size = 12))+
  ggtitle("(c3)")-> p9





###############
# plots <- list(p1, p3, p5, p6, p8, p10)
library(ggpubr)
pdf(file = paste0("result_figure/", name, "/", "np=", n_m_i, "_a=",alpha_i, "_0829.pdf"),
    onefile = F, width = 9, height = 9)
plot(ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, 
          nrow = 3, ncol = 3, common.legend = TRUE))
dev.off()
# 使用 annotate_figure 添加字母标签


#############################