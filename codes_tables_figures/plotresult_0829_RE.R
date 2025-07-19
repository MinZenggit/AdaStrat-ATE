# Evaluate the efficiency gain estimation
# Appendix A.4 (Table S2, Figure S2)

## Plot the estimated relative efficiency.
rm(list = ls())
library(dplyr)
library(ggplot2)
library(randomForest)
library(MASS)
library(foreach)
require(doSNOW)
library(caret)

name = "05-28-2test"
if(! dir.exists(paste0("result_figure/",name))){
  dir.create(paste0("result_figure/",name))
}

result_td_r <- function(x){
  mean(result[, 6]) -> r_hat
  mean(result[, 5]) -> r
  var(result[, 1])/var(result[, 3]) -> r_em
  return(c(r, r_hat, r_em))
}
method_set <- c("aipw", "ipw", "reg")
N_set = c(10000)
n_m_set = c(500, 1000, 2000)

n_set = c(2000)
K_set = c(2, 4, 6, 8, 10)

expand.grid(method_set, N_set, n_m_set, n_set,K_set)->settings
colnames(settings) <- c("method", "N", "n_m", "n", "K")

settings[order(settings$method, settings$N, 
               settings$n_m, settings$n, 
               settings$K), ] -> settings
rownames(settings) <- c(1:nrow(settings))

settings$r = settings$r_hat = settings$r_em = 0
# result_td_r(result)
for(j in 1:nrow(settings)){
  method = settings$method[j]
  N = settings$N[j]
  n_m = settings$n_m[j]
  n = settings$n[j]
  K = settings$K[j]
  path = paste0("result/", name, "/", paste(N, n_m, n, K, method, sep = "_"), ".rdata")
  load(file = path)
  result_td_r(result) -> rr
  settings$r[j] = rr[1]
  settings$r_hat[j] = rr[2]
  settings$r_em[j] = rr[3]
}

# ####################################
# # rm(p)
# alpha = 0.1
# for(method_i in "aipw"){
#   for(n_m_i in n_m_set){
#     for(K_i in K_set){
#       pdf(file = paste0("result_figure/", name, "/","n_m=",n_m_i, "_", "K=", K_i,"_", method_i, ".pdf"),
#           onefile = F)
#       r_td %>% filter(method == method_i,n_m == n_m_i, K==K_i, alpha == 0.3,
#                       des %in%c("u", "b", "o", "fm")) %>%
#         ggplot(aes(x = n, y = MSE, color = des, shape = des)) +
#         geom_line()+
#         geom_point(size = 2.5) +
#         scale_color_manual(values = c("u" = "#66c2a5", "b" = "#fc8d62", "fm" = "#8da0cb", "o" = "#e78ac3"),
#                            labels = c("u" ="RS"  , "b" = "Oracle" , "fm" ="K-as-opt" , "o" ="K-fs-opt")) +
#         scale_shape_manual(values = c("u" = 16, "b" = 17, "fm" = 18, "o" = 3),
#                            labels = c("u" ="RS"  , "b" = "Oracle" , "fm" ="K-as-opt" , "o" ="K-fs-opt")) +
#         theme_minimal() -> p
#       plot(p)
#       dev.off()
#       rm(p)
#     }
#   }
# }

method_i = "aipw"; n_m_i = 500; n_i = 2000; 
settings %>% filter(method == method_i, n_m == n_m_i, n== n_i) -> re
re <- data.frame(x = rep(re$K, 3), 
                 r = c(re$r_em, re$r_hat, re$r), 
                 type = c(rep("Empirical RE", 5), rep("Estimated RE", 5), rep("True RE", 5)))
ggplot(re, aes(x = x, y = r, color = type, shape = type, linetype = type)) +
  geom_line()+
  geom_point(size = 2.5)+
  ylim(c(0.5, 2))+
  labs(color = "", shape = "", linetype = "", 
       x = "K", 
       y = expression("RE"))+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"),
        legend.text = element_text(size = 12))+
  ggtitle("(a1)") -> p1

method_i = "aipw"; n_m_i = 1000; n_i = 2000; 
settings %>% filter(method == method_i, n_m == n_m_i, n== n_i) -> re
re <- data.frame(x = rep(re$K, 3), 
                 r = c(re$r_em, re$r_hat, re$r), 
                 type = c(rep("Empirical RE", 5), rep("Estimated RE", 5), rep("True RE", 5)))
ggplot(re, aes(x = x, y = r, color = type, shape = type, linetype = type)) +
  geom_line()+
  geom_point(size = 2.5)+
  ylim(c(0.5, 2))+
  labs(color = "", shape = "", linetype = "", 
       x = "K", 
       y = expression("RE"))+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"),
        legend.text = element_text(size = 12))+
  ggtitle("(a2)") -> p2

method_i = "aipw"; n_m_i = 2000; n_i = 2000; 
settings %>% filter(method == method_i, n_m == n_m_i, n== n_i) -> re
re <- data.frame(x = rep(re$K, 3), 
                 r = c(re$r_em, re$r_hat, re$r), 
                 type = c(rep("Empirical RE", 5), rep("Estimated RE", 5), rep("True RE", 5)))
ggplot(re, aes(x = x, y = r, color = type, shape = type, linetype = type)) +
  geom_line()+
  geom_point(size = 2.5)+
  ylim(c(0.5, 2))+
  labs(color = "", shape = "", linetype = "", 
       x = "K", 
       y = expression("RE"))+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"),
        legend.text = element_text(size = 12))+
  ggtitle("(a3)") -> p3



method_i = "ipw"; n_m_i = 500; n_i = 2000; 
settings %>% filter(method == method_i, n_m == n_m_i, n== n_i) -> re
re <- data.frame(x = rep(re$K, 3), 
                 r = c(re$r_em, re$r_hat, re$r), 
                 type = c(rep("Empirical RE", 5), rep("Estimated RE", 5), rep("True RE", 5)))
ggplot(re, aes(x = x, y = r, color = type, shape = type, linetype = type)) +
  geom_line()+
  geom_point(size = 2.5)+
  ylim(c(0.5, 2))+
  labs(color = "", shape = "", linetype = "", 
       x = "K", 
       y = expression("RE"))+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"),
        legend.text = element_text(size = 12))+
  ggtitle("(b1)") -> p4

method_i = "ipw"; n_m_i = 1000; n_i = 2000; 
settings %>% filter(method == method_i, n_m == n_m_i, n== n_i) -> re
re <- data.frame(x = rep(re$K, 3), 
                 r = c(re$r_em, re$r_hat, re$r), 
                 type = c(rep("Empirical RE", 5), rep("Estimated RE", 5), rep("True RE", 5)))
ggplot(re, aes(x = x, y = r, color = type, shape = type, linetype = type)) +
  geom_line()+
  geom_point(size = 2.5)+
  ylim(c(0.5, 2))+
  labs(color = "", shape = "", linetype = "", 
       x = "K", 
       y = expression("RE"))+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"),
        legend.text = element_text(size = 12))+
  ggtitle("(b2)") -> p5

method_i = "ipw"; n_m_i = 2000; n_i = 2000; 
settings %>% filter(method == method_i, n_m == n_m_i, n== n_i) -> re
re <- data.frame(x = rep(re$K, 3), 
                 r = c(re$r_em, re$r_hat, re$r), 
                 type = c(rep("Empirical RE", 5), rep("Estimated RE", 5), rep("True RE", 5)))
ggplot(re, aes(x = x, y = r, color = type, shape = type, linetype = type)) +
  geom_line()+
  geom_point(size = 2.5)+
  ylim(c(0.5, 2))+
  labs(color = "", shape = "", linetype = "", 
       x = "K", 
       y = expression("RE"))+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"),
        legend.text = element_text(size = 12))+
  ggtitle("(b3)") -> p6


method_i = "reg"; n_m_i = 500; n_i = 2000; 
settings %>% filter(method == method_i, n_m == n_m_i, n== n_i) -> re
re <- data.frame(x = rep(re$K, 3), 
                 r = c(re$r_em, re$r_hat, re$r), 
                 type = c(rep("Empirical RE", 5), rep("Estimated RE", 5), rep("True RE", 5)))
ggplot(re, aes(x = x, y = r, color = type, shape = type, linetype = type)) +
  geom_line()+
  geom_point(size = 2.5)+
  ylim(c(0.5, 2))+
  labs(color = "", shape = "", linetype = "", 
       x = "K", 
       y = expression("RE"))+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"),
        legend.text = element_text(size = 12))+
  ggtitle("(c1)") -> p7

method_i = "reg"; n_m_i = 1000; n_i = 2000; 
settings %>% filter(method == method_i, n_m == n_m_i, n== n_i) -> re
re <- data.frame(x = rep(re$K, 3), 
                 r = c(re$r_em, re$r_hat, re$r), 
                 type = c(rep("Empirical RE", 5), rep("Estimated RE", 5), rep("True RE", 5)))
ggplot(re, aes(x = x, y = r, color = type, shape = type, linetype = type)) +
  geom_line()+
  geom_point(size = 2.5)+
  ylim(c(0.5, 2))+
  labs(color = "", shape = "", linetype = "", 
       x = "K", 
       y = expression("RE"))+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"),
        legend.text = element_text(size = 12))+
  ggtitle("(c2)") -> p8

method_i = "reg"; n_m_i = 2000; n_i = 2000; 
settings %>% filter(method == method_i, n_m == n_m_i, n== n_i) -> re
re <- data.frame(x = rep(re$K, 3), 
                 r = c(re$r_em, re$r_hat, re$r), 
                 type = c(rep("Empirical RE", 5), rep("Estimated RE", 5), rep("True RE", 5)))
ggplot(re, aes(x = x, y = r, color = type, shape = type, linetype = type)) +
  geom_line()+
  geom_point(size = 2.5)+
  ylim(c(0.5, 2))+
  labs(color = "", shape = "", linetype = "", 
       x = "K", 
       y = expression("RE"))+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"),
        legend.text = element_text(size = 12))+
  ggtitle("(c3)") -> p9





library(ggpubr)
pdf(file = paste0("result_figure/", name, "/",  "rhat_0829.pdf"),
    onefile = F, width = 7.5, height = 7.5)
ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, 
          nrow = 3, ncol = 3, common.legend = TRUE)
dev.off()

