data {
    int<lower=0> N;   // number of data items
    int<lower=0> K;   // number of predictors
    // int<lower=0> P;   // numero de efeitos aleatorios
    // int<lower=0> M;   // numero de itens na base de teste
    matrix[N, K] fixed_effects;
    // matrix[N, P] efeitos_aleatorios;
    vector[N] medicao;
    array[N] int<lower=0, upper=1> y;
    real<lower=0, upper=1> pred_cutoff;
    real beta_f;
    real<lower=0, upper=1> dose_letal;
}

parameters {
    vector[K] beta_prob;
    // vector[P] theta_prob;
    real beta_medicao_prob;
    real beta_sqrd_medicao_prob;
    // real<lower=0> desvio_theta_prob;
}

transformed parameters {
    // real<lower=0> desvio_theta_prob = exp(log_desvio_theta_prob);
    vector[N] p = inv_logit(
        fixed_effects * beta_prob +
        medicao * beta_medicao_prob +
        medicao^2 * beta_sqrd_medicao_prob
    );
}

model {
    // print(max(p));
    // print(min(p));
    y ~ bernoulli_logit(
        fixed_effects * beta_prob +
        medicao * beta_medicao_prob +
        medicao^2 * beta_sqrd_medicao_prob
    );
    // beta_sqrd_medicao_prob ~ normal(0, 1);
    // beta_medicao_prob ~ normal(0, 1);
    beta_prob ~ normal(0, 3);
    // theta_prob ~ normal(0, desvio_theta_prob);
    // desvio_theta_prob ~ gamma(2, 2);
}

generated quantities {
    array[N] int b = bernoulli_logit_rng(fixed_effects * beta_prob +
                                    medicao * beta_medicao_prob +
                                    medicao^2 * beta_sqrd_medicao_prob);
    array[N] int<lower=0> predicoes;
    real<lower=0> falso_positivo = 0;
    real<lower=0> falso_negativo = 0;
    real<lower=0> verdadeiro_positivo = 0;
    real<lower=0> verdadeiro_negativo = 0;

    vector[N] a;
    vector[N] bi;

    for (i in 1:N) {
        if (y[i] == 1) {
            a[i] = 1 - inv_logit(fixed_effects[i] * beta_prob +
                                    medicao[i] * beta_medicao_prob +
                                    medicao[i]^2 * beta_sqrd_medicao_prob);
            bi[i] = 1;
        } else {
            a[i] = 0;
            bi[i] = 1 - inv_logit(fixed_effects[i] * beta_prob +
                                    medicao[i] * beta_medicao_prob +
                                    medicao[i]^2 * beta_sqrd_medicao_prob);
        }
        predicoes[i] = inv_logit(fixed_effects[i] * beta_prob +
                                    medicao[i] * beta_medicao_prob +
                                    medicao[i]^2 * beta_sqrd_medicao_prob) > pred_cutoff;
        if (y[i] == 0 && predicoes[i] == 1) {
            falso_positivo += 1;
        } else if (y[i] == 1 && predicoes[i] == 0) {
            falso_negativo += 1;
        } else if (y[i] == 1 && predicoes[i] == 1) {
            verdadeiro_positivo += 1;
        } else {
            verdadeiro_negativo += 1;
        }
    }

    array[N] real u = uniform_rng(a, bi);
    array[N] real residuo_quantil_normalizado = std_normal_qf(u);

    real recall = verdadeiro_positivo / (verdadeiro_positivo + falso_negativo);
    real precisao = verdadeiro_positivo / (verdadeiro_positivo + falso_positivo);
    real acuracia = (verdadeiro_positivo + verdadeiro_negativo) /  N;
    real f_beta_score = (1 + beta_f^2) * (precisao * recall) / ((beta_f^2) * precisao + recall);
    real kappa = (acuracia - mean(y)) / (1 - mean(y));
    vector[N] c = fixed_effects * beta_prob - logit(dose_letal);
    vector[N] sqrt_delta = sqrt(beta_medicao_prob^2 - 4 * beta_sqrd_medicao_prob * c);
    vector[N] dose_letal_1 = -(beta_medicao_prob + sqrt_delta) / (2 * beta_sqrd_medicao_prob);
    vector[N] dose_letal_2 = -(beta_medicao_prob - sqrt_delta) / (2 * beta_sqrd_medicao_prob);
}
