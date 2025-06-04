library(tidyverse)
library(ggplot2)


df <- read.csv('./dados/base_completa.csv') %>%
    filter(Texto %in% c(0:5, 9))

df.pivoted <- df %>%
    select(Equipamento, Data, `Denominação`, Texto) %>%
    pivot_wider(
        names_from = `Denominação`,
        values_from = `Texto`,
        values_fn = function(x) mean(as.double(x), na.rm = T),
        values_fill = NA
    )


reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1 - cormat)/2)
    hc <- hclust(dd)
    cormat <- cormat[hc$order, hc$order]
}

cormat <- round(cor(df.pivoted %>%
                        select(!Equipamento:Data), use = 'pairwise.complete.obs'), 2) %>%
    # reorder_cormat() %>%
    as.data.frame()
View(cormat)
cormat <- cormat[, (colSums(is.na(cormat)) - 1) < nrow(cormat)]
cormat <- cormat[(rowSums(is.na(cormat)) - 1) < ncol(cormat), ]
cormat$Var1 <- rownames(cormat)
vars <- rownames(cormat)


melted_cormat <- cormat %>%
    pivot_longer(
        cols = !Var1,
        names_to = 'Var2'
    ) %>% mutate(
        Var1 = factor(Var1, levels = vars),
        Var2 = factor(Var2, levels = rev(vars))
    )


length(vars)
v1 <- c(vars[1], vars[75])
v2 <- c(vars[2], vars[149])
vars[vars %>% startsWith('Vazamento')]


melted_cormat %>%
    ggplot(aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "red", high = "blue", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name = "Correlação de Pearson") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank())

cormat %>%
    arrange(desc(`Vazamento Ar Sistema Pneumático`)) %>%
    select(`Vazamento Ar Sistema Pneumático`) %>%
    View()


melted_cormat %>%
    ggplot(aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "red", high = "blue", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name = "Correlação de Pearson") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 12, hjust = 1)) +
    coord_fixed() +
    xlim(vars[150 - 12:24]) +
    ylim(rev(vars[150 - 12:24]))

df %>%
    filter(Equipamento == 188988,
           Denominação == 'Comp/tes Elétricos-Instrum Risco de Avar',
           Data == '2024-05-14') %>%
    select(Equipamento, Data, Cód.valorização, Denominação, Texto) %>% knitr::kable()


df %>%
    filter(ValMed.PosTCont > LimSupIntMed.) %>% View()


df %>%
   filter(ValMed.PosTCont > LimSupIntMed.,
          Denominação == 'Tempo de Abertura Camara 1 Fase B',
          Texto %in% c(3, 5)) %>%
    select(Equipamento, Data, ValMed.PosTCont, LimSupIntMed., Texto, Denominação) %>%
    View()
