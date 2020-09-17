data_trt <- unique(apply(MTCdata[ , c("treat1", "treat2")], 1, FUN = function(x){paste0(x[1],":",x[2])}))
MTCdata$comp <- paste0(MTCdata$treat1, ":", MTCdata$treat2)
comp_trt <- sig_comp_all_df[rownames(sig_comp_all_df) %in% data_trt,]
comp_trt <- comp_trt %>% tibble::rownames_to_column() %>% 
  left_join(., 
            MTCdata %>% group_by(comp) %>% summarise(n=n()) %>% mutate(fullname=paste0(comp, " (", n, ")")),
            c("rowname"="comp")) %>% 
  tibble::column_to_rownames("fullname") %>% 
  dplyr::select(-rowname)
comp <- data.frame(comparisons = rep(rownames(comp_trt),2), lower = c(comp_trt$MLE_single_lower, comp_trt$MLE_multiple_lower),upper = c(comp_trt$MLE_single_upper, comp_trt$MLE_multiple_upper), type = rep(c(1, 2), each = nrow(comp_trt)))




library(ggplot2)
comp$comparisons <- factor(comp$comparisons, levels = c(rownames(comp_trt)[grep("^N",comp_trt %>% rownames())], rownames(comp_trt)[setdiff(1:24, grep("^N",comp_trt %>% rownames()))]))
ggplot(comp, aes(ymin = lower, ymax = upper, x = comparisons)) +
  geom_linerange(aes(color = factor(comp$type, levels = c(1,2))),position = position_dodge(width = -0.3), size = 1) + 
  
  guides(color=guide_legend(title= expression(paste("number of ", tau^{2}, "s", sep = "")) )) +
  
  ylab("log odds ratio") +
  theme(axis.text.y = element_text(colour = c(rep("darkblue", 11), rep("black", 13)))) +
  coord_flip()




