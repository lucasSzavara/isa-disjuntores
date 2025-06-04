library(dplyr)
library(cmdstanr)
library(ggplot2)
library(bayesplot)
library(tidyr)
library(qqplotr)


df <- read.csv('./dados/base_completa.csv')
df$avaliacao <- df$Texto %>% as.numeric()
df <- df %>%
    filter(!is.na(LimInfIntMed.),
           !is.na(LimSupIntMed.),
           !is.na(avaliacao),
           avaliacao <= 5)

df$y <- df$avaliacao >= 4
df$anos_operacao <- (as.Date(df$Data) - as.Date(df$Dt.entr.serviço)) / 365
mod <- cmdstan_model('./reg_logistica.stan')

df <- df[complete.cases(df[, c(
    'anos_operacao',
    'Código.ABC',
    'Tensão',
    'Tipo.Instalação',
    'Classe',
    'Valor.teórico',
    'y',
    'Equipamento',
    'Fabricante',
    'Modelo',
    'Unid.caracter.',
    'Denominação'
)]), ]

df <- df %>%
    filter(startsWith(Denominação, 'Tempo de Abertura'))

fixed_effects <- model.matrix(~ anos_operacao + Classe + Código.ABC + Tensão + Tipo.Instalação + Valor.teórico, df)

X_p <- model.matrix(~ Classe - 1, df)

nrow(fixed_effects)
nrow(df)
nrow(X_p)
data <- list(
    N = nrow(fixed_effects),
    K = ncol(fixed_effects),
    # P = ncol(X_p),
    medicao = log(df$ValMed.PosTCont),
    fixed_effects = fixed_effects,
    # efeitos_aleatorios = X_p,
    y = df$y,
    pred_cutoff = 0.5,
    beta_f = 1,
    dose_letal = 0.5
)
start <- Sys.time()
opt <- mod$optimize(
    data = data,
    seed = 42,
    jacobian = T,
    # init = list(list(
    #     beta_prob = array(0, dim = ncol(fixed_effects)),
    #     theta_prob = array(0, dim = ncol(X_p)),
    #     beta_medicao_prob = 0,
    #     beta_sqrd_medicao_prob = 0,
    #     log_desvio_theta_prob = 0
    # )),
    algorithm = 'lbfgs',
    iter = 1e5
    # tol_param = 1e-20
)
fit <- mod$laplace(
    data = data,
    seed = 42,
    # parallel_chains = 4
    init = opt,
    mode = opt
)
# fit <- mod$sample(
#     data = data,
#     seed = 42,
#     # init = opt,
#     chains = 4,
#     parallel_chains = 4
# )
end <- Sys.time()
print(end - start)

resultados <- fit$summary(variables = c(
    'beta_prob',
    'beta_medicao_prob',
    'beta_sqrd_medicao_prob'
    # 'acuracia',
    # 'f_beta_score',
    # 'kappa'
)) %>% as.data.frame()
resultados$mean[startsWith(resultados$variable, 'beta')] <- exp(resultados$mean[startsWith(resultados$variable, 'beta')])
resultados$q5[startsWith(resultados$variable, 'beta')] <- exp(resultados$q5[startsWith(resultados$variable, 'beta')])
resultados$q95[startsWith(resultados$variable, 'beta')] <- exp(resultados$q95[startsWith(resultados$variable, 'beta')])
#
resultados$Info[startsWith(resultados$variable, 'beta_prob')] <- colnames(fixed_effects)
resultados$variable[startsWith(resultados$variable, 'beta_prob')] <- paste('efeito fixo prob', colnames(fixed_effects))
resultados$variable[startsWith(resultados$variable, 'beta_medicao_prob')] <- paste('efeito fixo prob medição')

# Validação p
y_rep <- fit$draws('b', format = 'matrix')
y <- as.double(df$y)
ppc_bars_grouped(y, y_rep, group = df$Código.ABC)

ppc_stat_grouped(y, y_rep, stat = 'sum', group = df$Código.ABC) +
    labs(title = 'Quantidade de medições positivas por criticidade')


n_sample <- 16
resid <- fit$draws('residuo_quantil_normalizado', format = 'matrix')
linha <- sample(1:nrow(resid), size = n_sample)

r <- t(resid[linha, ]) %>% as.data.frame()
colnames(r) <- 1:n_sample
r$delta <- df$y
r$valor <- df$ValMed.PosTCont
r <- r %>%
    pivot_longer(
        cols = 1:n_sample,
        names_to = 'observacao',
        values_to = 'residuo'
    )
r %>%
    ggplot(aes(sample = residuo)) +
    stat_qq_band(qprobs = c(0.025, 0.975),
                 distribution = 'norm',
                 dparams = list(mean = 0, sd = 1),
                 bandType = 'ell') +
    stat_qq_line(identity=T) +
    stat_qq_point() +
    facet_wrap(observacao~.)

r %>%
    ggplot(aes(x = log(valor), y = residuo)) +
    geom_point() +
    facet_wrap(observacao~.)


metricas <- fit$draws(variables = c(
    'falso_positivo',
    'falso_negativo',
    'verdadeiro_positivo',
    'verdadeiro_negativo',
    'recall',
    'precisao',
    'acuracia',
    'f_beta_score',
    'kappa'
), format = 'matrix')
metricas.pontual <- metricas %>% apply(2, mean)
metricas.pontual

classe_real <- c('Acordo', 'Acordo', 'Sem acordo', 'Sem acordo')
classe_prevista <- c('Acordo', 'Sem acordo', 'Acordo', 'Sem acordo')
y.metricas <- metricas.pontual[c('verdadeiro_positivo', 'falso_negativo', 'falso_positivo', 'verdadeiro_negativo')]
matriz.confusao <- data.frame(classe_real, classe_prevista, y.metricas)
ggplot(data =  matriz.confusao, mapping = aes(x = classe_real, y = classe_prevista)) +
    geom_tile(aes(fill = y.metricas), colour = "white") +
    geom_text(aes(label = sprintf("%1.0f", y.metricas)), vjust = 1) +
    xlab('Resultado real') +
    ylab('Resultado predito') +
    # scale_fill_gradient(low = "blue", high = "red") +
    theme_bw() +
    theme(legend.position = "none")

p_rep <- fit$draws('p', format = 'matrix')
ic_p <- p_rep %>%
    apply(2, function(x) quantile(x, probs=c(0.025, 0.25, 0.5, 0.75, 0.975)))
p_aproximado <- round(ic_p["50%",] * 10)
df['p_aproximado'] <-  p_aproximado
df['p_estimado'] <-  ic_p["50%",]



df %>%
    group_by(p_aproximado) %>%
    summarise(media = mean(y),
              Qtd = n()) %>%
    mutate(erro = qnorm(0.025)*sqrt(media * (1 - media) / Qtd)) %>%
    ggplot(aes(x = p_aproximado / 10, y = media_acordos)) +
    geom_point(aes(size = Qtd)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    geom_segment(aes(y = media_acordos - erro, yend = media_acordos + erro)) +
    geom_abline(slope = 1, intercept = 0) +
    xlab('Probabilidade estimada') +
    ylab('Probabilidade real')

df %>%
    ggplot(aes(x = p_estimado, y = as.double(y))) +
    geom_smooth(method = "lm", se = T) +
    geom_abline(slope = 1, intercept = 0) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    xlab('Probabilidade estimada') +
    ylab('Probabilidade real')

dose_letal_50.1 <- fit$draws('dose_letal_1', format = 'matrix')
dose_letal_50.2 <- fit$draws('dose_letal_2', format = 'matrix')
ic_y <- dose_letal_50.1 %>%
    exp() %>%
    apply(2, function(x) quantile(x[!is.infinite(x)], probs=c(0.025, 0.25, 0.5, 0.75, 0.975)))

df[, 'Dose Letal 1'] <-  round(ic_y["50%",])
df[, 'IC (50%) para DL1'] <-  paste(round(ic_y["25%",], 2), '-', round(ic_y["75%",], 2))
df[, 'IC (95%) para DL1'] <-  paste(round(ic_y["2.5%",], 2), '-', round(ic_y["97.5%",], 2))

ic_y <- dose_letal_50.2 %>%
    exp() %>%
    apply(2, function(x) quantile(x[!is.infinite(x)], probs=c(0.025, 0.25, 0.5, 0.75, 0.975)))

df[, 'Dose Letal 2'] <-  round(ic_y["50%",])
df[, 'IC (50%) para DL2'] <-  paste(round(ic_y["25%",], 2), '-', round(ic_y["75%",], 2))
df[, 'IC (95%) para DL2'] <-  paste(round(ic_y["2.5%",], 2), '-', round(ic_y["97.5%",], 2))
