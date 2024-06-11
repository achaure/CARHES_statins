pacman::p_load(
  dplyr,
  tidyr,
  rio,
  lubridate,
  smd,
  ggplot2,
  boot,
  splitstackshape
)

base_analisis <- base_larga_trabajo
limite <- .01
n_muestras <- 5

load("datos/ensayos_itt")
load("datos/base_prueba")
load("original/Indice_HL.RDAta")

ft <- 
  "iniciador ~ sexo_fem + 
  edad + I(edad*edad) + 
  ldl_t + I(ldl_t*ldl_t) +
  hdl_t + I(hdl_t*hdl_t) +
  colest_t + I(colest_t*colest_t) +
  glucosa_t + I(glucosa_t*glucosa_t) +
  n_analiticas_t + I(n_analiticas_t * n_analiticas_t) +
  diabetes_t + hta_t + as.factor(cohorte)"

fc <- 
  "cens == 0 ~ sexo_fem + 
  edad + I(edad*edad) + 
  ldl_t + I(ldl_t*ldl_t) +
  hdl_t + I(hdl_t*hdl_t) +
  colest_t + I(colest_t*colest_t) +
  glucosa_t + I(glucosa_t*glucosa_t) +
  n_analiticas_t + I(n_analiticas_t * n_analiticas_t) +
  diabetes_t + hta_t + as.factor(cohorte)"

analisis_prot <- function(
    base_analisis, 
    formula_trat,
    formula_cens,
    limite = 0, 
    seguimiento, 
    n_muestras = 500,
    formula_grafico = length(formulas),
    tipo_grafico = "risk"
    ){
  analisis_descriptivo_prot <- function(base_analisis){
    
    # RESTO DESCRIPTIVO
    # Eventos totales
    et_1 <- base_analisis %>%
      filter(evento == 1) %>%
      filter(iniciador == 1) %>%
      nrow()
    
    et_0 <- base_analisis %>%
      filter(evento == 1) %>%
      filter(iniciador == 0) %>%
      nrow()
    
    # eventos reales
    er_1 <- base_analisis %>%
      filter(evento == 1) %>%
      filter(iniciador == 1) %>%
      distinct(patient_id) %>%
      nrow()
    
    er_0 <- base_analisis %>%
      filter(evento == 1) %>%
      filter(iniciador == 0) %>%
      distinct(patient_id) %>%
      nrow()
    
    # Tiempo de seguimiento
    t_1 <- base_analisis %>%
      distinct(patient_id, fecha_inicio, .keep_all = TRUE) %>%
      filter(iniciador == 1) %>%
      {paste0(median(.$seguim), " (", quantile(.$seguim, .25), "-", quantile(.$seguim, .75), ")")}
    
    t_0 <- base_analisis %>%
      distinct(patient_id, fecha_inicio, .keep_all = TRUE) %>%
      filter(iniciador == 0) %>%
      {paste0(median(.$seguim), " (", quantile(.$seguim, .25), "-", quantile(.$seguim, .75), ")")}
    

    resultados <- c("Eventos totales",
                    "Eventos reales",
                    "Seguimiento")
    
    iniciadores <- c(et_1, er_1, t_1)
    no_iniciadores <- c(et_0, er_0, t_0)
    
    resto <- data.frame(iniciadores, no_iniciadores)
    row.names(resto) <- resultados
    
    return(resto)
    
  }
  
  descriptivo <- analisis_descriptivo_prot(base_analisis)
  
  analisis_prot_tte <- function(base_analisis, n_muestras){
    
    funcion_boot <- function(base, indices) {
      
      patient_ids_unicos <- unique(base$patient_id)
      muestra_patient_ids <- patient_ids_unicos[indices]
      b <- base %>% filter(patient_id %in% muestra_patient_ids)
      
      b <- b %>%
        mutate(cohorte = year(fecha_inicio))


      # Calculo y añado los pesos de tratamiento
      base_corta <- filter(b, tiempo == 0)
      
      p_denom <- glm(formula_trat, data = base_corta, family = binomial())
      p_num <- glm(iniciador ~ 1, data = base_corta, family = binomial())
      
      base_corta <- base_corta %>%
        {mutate(., denom = predict(p_denom, ., type = "response"),
                num = predict(p_num, ., type = "response"))} %>%
        mutate(peso_trat = if_else(iniciador == 1, 
                                   num/denom, 
                                   (1 - num)/(1 - denom))) %>%
        select(patient_id, fecha_inicio, peso_trat)
      
      base_larga <- b %>%
        left_join(base_corta, by = c("patient_id", "fecha_inicio"))
      
      # Calculo y añado los pesos de la censura
      
      # 1. Identifico las censuras
      base_larga <- base_larga %>%
        mutate(cens = if_else(tiempo == seguim - 1 & 
                                evento == 0 & 
                                fecha_fin != max(fecha_fin),
                              1, 0))
      
      # 2. Calculo los pesos 
      p_denom_c <- glm(formula_cens, data = base_larga, family = binomial())
      p_num_c <- glm(cens == 0 ~ 1, data = base_larga, family = binomial())
      
      base_larga <- base_larga %>%
        {mutate(., denom = predict(p_denom_c, ., type = "response"),
                num = predict(p_num_c, ., type = "response"))} %>%
        mutate(peso_cens = if_else(cens == 0, 
                                num/denom, 
                                (1-num)/(1-denom))) %>%
        mutate(peso_cens = cumprod(peso_cens), .by = c("patient_id", "fecha_inicio"),
               peso = peso_trat * peso_cens) 
      
      # recorto los pesos extremos
      base_larga <- base_larga %>%
        filter(peso <= quantile(peso, 1-limite),
               quantile(peso, limite) <= peso)
      
      modelo_hazards <- glm(
        evento == 0 ~ iniciador + I(iniciador*tiempo) + I(iniciador*tiemposq) + tiempo + tiemposq, 
        family = binomial(),
        weights = peso,
        data = base_larga
      )
      
      trat0 <- data.frame(cbind(seq(0, seguimiento),0,(seq(0, seguimiento))^2))
      trat1 <- data.frame(cbind(seq(0, seguimiento),1,(seq(0, seguimiento))^2))
      
      colnames(trat0) <- c("tiempo", "iniciador", "tiemposq")
      colnames(trat1) <- c("tiempo", "iniciador", "tiemposq")
      
      trat0 <- trat0 %>%
        {mutate(., p_noevent0 = predict(modelo_hazards, ., type = "response"))} %>%
        mutate(surv0 = cumprod(p_noevent0),
               risk0 = 1 - surv0)
      
      trat1 <- trat1 %>%
        {mutate(., p_noevent1 = predict(modelo_hazards, ., type = "response"))} %>%
        mutate(surv1 = cumprod(p_noevent1),
               risk1 = 1 - surv1)
      
      graf_hazards <- merge(trat0, trat1, by = c("tiempo", "tiemposq"))
      
      
      
      # Riesgos
      # Riesgo no iniciadores
      r0 <- round(max(trat0$risk0)*100, 2)
      
      # Riesgo iniciadores
      r1 <- round(max(trat1$risk1)*100, 2)
      
      # Razón de riesgos
      rr <- round(r1 / r0, 2)
      
      # Diferencia de riesgos
      rd <- r1 - r0
      
      v1 <- c("Riesgo en no iniciadores:",
              "Riesgo en iniciadores:",
              "Razón de riesgos:",
              "Diferencia de riesgos:")
      
      return(c(r0, r1, rr, rd))
      }
    set.seed(1234)
    boot_results <- boot(data = base_analisis, statistic = funcion_boot, R = n_muestras)
    
    v1 <- c("Analisis por protocolo")
    v2 <- boot_results$t0[1]
    v3 <- boot_results$t0[2]
    v4 <- boot_results$t0[3]
    v5 <- boot_results$t0[4]
    
    v2.2 <- c(paste0( "(",
                      boot.ci(boot_results, type = "perc", index = 1)$percent[4],
                      "%-",
                      boot.ci(boot_results, type = "perc", index = 1)$percent[5],
                      "%)"
    ))
    
    v3.2 <- c(paste0( "(",
                      boot.ci(boot_results, type = "perc", index = 2)$percent[4],
                      "%-",
                      boot.ci(boot_results, type = "perc", index = 2)$percent[5],
                      "%)"
    ))
    
    v4.2 <- c(paste0( "(",
                      boot.ci(boot_results, type = "perc", index = 3)$percent[4],
                      "-",
                      boot.ci(boot_results, type = "perc", index = 3)$percent[5],
                      ")"
    ))
    
    v5.2 <- c(paste0( "(",
                      boot.ci(boot_results, type = "perc", index = 4)$percent[4],
                      "%-",
                      boot.ci(boot_results, type = "perc", index = 4)$percent[5],
                      ")"))
    
    v2 <- paste0(round(v2,2), "% ", v2.2)
    v3 <- paste0(round(v3,2), "% ", v3.2)
    v4 <- paste0(round(v4,2), " ", v4.2)
    v5 <- paste0(round(v5,2), "% ", v5.2)
    
    tabla <- data.frame(v2, v3, v4, v5)
    colnames(tabla) <- c("Risk in non-initiators", "Risk in initiators", "RR", "RD")
    row.names(tabla) <- v1
    tabla
    
  }
  prot <- analisis_prot_tte(base_analisis, n_muestras)

  return(list(descriptivo = descriptivo, 
              Protocolo = prot))
  
}

resultados <- analisis_prot(
  base_analisis = base_larga_trabajo,
  formula_trat = ft,
  formula_cens = fc,
  limite = 0.01,
  seguimiento = 60,
  n_muestras = 5
  )
