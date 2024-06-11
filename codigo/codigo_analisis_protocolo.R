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

base <- base_larga_trabajo

limite <- .01
muestras <- 5

load("datos/ensayos_itt")
load("datos/base_prueba")
load("original/Indice_HL.RDAta")
analisis_prot <- function(base, muestras = 500){
  analisis_descriptivo_prot <- function(base){
    
    # RESTO DESCRIPTIVO
    # Eventos totales
    et_1 <- base %>%
      filter(evento == 1) %>%
      filter(iniciador == 1) %>%
      nrow()
    
    et_0 <- base %>%
      filter(evento == 1) %>%
      filter(iniciador == 0) %>%
      nrow()
    
    # eventos reales
    er_1 <- base %>%
      filter(evento == 1) %>%
      filter(iniciador == 1) %>%
      distinct(patient_id) %>%
      nrow()
    
    er_0 <- base %>%
      filter(evento == 1) %>%
      filter(iniciador == 0) %>%
      distinct(patient_id) %>%
      nrow()
    
    # Tiempo de seguimiento
    t_1 <- base %>%
      distinct(patient_id, fecha_inicio, .keep_all = TRUE) %>%
      filter(iniciador == 1) %>%
      {paste0(median(.$seguim), " (", quantile(.$seguim, .25), "-", quantile(.$seguim, .75), ")")}
    
    t_0 <- base %>%
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
  
  descriptivo <- analisis_descriptivo_prot(base)
  
  analisis_prot_tte <- function(base, muestras){
    
    funcion_boot <- function(base, indices) {
      
      # Obtener los identificadores únicos de los pacientes
      patient_ids_unicos <- unique(base$patient_id)
      
      # Convertir índices en identificadores de pacientes seleccionados
      muestra_patient_ids <- patient_ids_unicos[indices]
      
      # Filtrar los datos para incluir solo los sujetos seleccionados
      b <- base %>% filter(patient_id %in% muestra_patient_ids)
      
      b <- b %>%
        mutate(cohorte = year(fecha_inicio))
      
      # CALCULO LOS PESOS DEL TRATAMIENTO
      
      base_corta <- filter(b, tiempo == 0)
      
      # estimo el denominador 
      p_denom_2 <- glm(iniciador ~ sexo_fem + 
                         edad + I(edad*edad) + 
                         ldl_t + I(ldl_t*ldl_t) +
                         hdl_t + I(hdl_t*hdl_t) +
                         colest_t + I(colest_t*colest_t) +
                         nhdl_t + I(nhdl_t*nhdl_t) +
                         glucosa_t + I(glucosa_t*glucosa_t) +
                         n_analiticas_t + I(n_analiticas_t * n_analiticas_t) +
                         diabetes_t + 
                         as.factor(cohorte) +
                         hta_t,
                       data = base_corta, family = binomial())
      
      # estimo el numerador
      p_num <- glm(iniciador ~ 1, data = base_corta, family = binomial())
      
      # Calculo los pesos
      base_corta <- base_corta %>%
        {mutate(., pd_trat_2 = predict(p_denom_2, ., type = "response"),
                pn_trat = predict(p_num, ., type = "response"))} %>%
        mutate( pes_2 = if_else(iniciador == 1, 
                               pn_trat/pd_trat_2, 
                               (1-pn_trat)/(1-pd_trat_2))) %>%
        select(patient_id, fecha_inicio, pes_2)
      
      base_larga <- b %>%
        left_join(base_corta, by = c("patient_id", "fecha_inicio"))
      
      
      # Calculo los pesos de la censura
      # Identifico las censuras
      
      base_larga <- base_larga %>%
        mutate(cens = if_else(tiempo == seguim - 1 & evento == 0 & fecha_fin != max(fecha_fin),
                              1, 0))
      
      # estimo el denominador 
      p_denom_c <- glm(cens == 0 ~ sexo_fem + 
                         edad + I(edad*edad) + 
                         ldl_t + I(ldl_t*ldl_t) +
                         hdl_t + I(hdl_t*hdl_t) +
                         colest_t + I(colest_t*colest_t) +
                         nhdl_t + I(nhdl_t*nhdl_t) +
                         glucosa_t + I(glucosa_t*glucosa_t) +
                         n_analiticas_t + I(n_analiticas_t * n_analiticas_t) +
                         diabetes_t + 
                         as.factor(cohorte) +
                         hta_t,
                       data = base_larga, family = binomial())
      
      # estimo el numerador
      p_num_c <- glm(cens == 0 ~ 1, data = base_larga, family = binomial())
      
      # Calculo los pesos
      base_larga <- base_larga %>%
        {mutate(., pd_c = predict(p_denom_c, ., type = "response"),
                pn_c = predict(p_num_c, ., type = "response"))} %>%
        mutate( pes_c = if_else(cens == 0, 
                                pn_c/pd_c, 
                                (1-pn_c)/(1-pd_c))) %>%
        mutate(pes_c = cumprod(pes_c), .by = c("patient_id", "fecha_inicio"),
               peso = pes_c * pes_2) 
      
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
      
      trat0 <- data.frame(cbind(seq(0, 60),0,(seq(0, 60))^2))
      trat1 <- data.frame(cbind(seq(0, 60),1,(seq(0, 60))^2))
      
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
    boot_results <- boot(data = base, statistic = funcion_boot, R = muestras)
    
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
  prot <- analisis_prot_tte(base, muestras)
  
  
  # grafico_ajustado <- function(base){
  #   
  #   # La única diferencia entre ajustado o no ajustado son los pesos. Así que meto todos los
  #   # pesos en una única base y luego analizo con o sin pesos.
  #   base <- base %>%
  #     mutate(cohorte = year(fecha_inicio))  
  #   
  #   
  #   # estimo el denominador del ajuste
  #   p_denom_2 <- glm(iniciador ~ sexo_fem + 
  #                      edad + I(edad*edad) + 
  #                      ldl_t + I(ldl_t*ldl_t) +
  #                      hdl_t + I(hdl_t*hdl_t) +
  #                      colest_t + I(colest_t*colest_t) +
  #                      nhdl_t + I(nhdl_t*nhdl_t) +
  #                      glucosa_t + I(glucosa_t*glucosa_t) +
  #                      n_analiticas_t + I(n_analiticas_t * n_analiticas_t) +
  #                      diabetes_t + 
  #                      hta_t +
  #                      as.factor(cohorte),
  #                    data = base, family = binomial())
  #   
  #   # estimo el numerador
  #   p_num <- glm(iniciador ~ 1, data = base, family = binomial())
  #   
  #   # Calculo los pesos
  #   base <- base %>%
  #     {mutate(., pd_trat_2 = predict(p_denom_2, ., type = "response"),
  #             pn_trat = predict(p_num, ., type = "response"))} %>%
  #     mutate(pes_2 = if_else(iniciador == 1, 
  #                            pn_trat/pd_trat_2, 
  #                            (1-pn_trat)/(1-pd_trat_2)))
  #   
  #   base_larga <- base %>%
  #     mutate(seguim = as.numeric(round((fecha_fin - fecha_inicio)/30.4375)),
  #            id = paste0(patient_id, fecha_inicio)) %>%
  #     expandRows("seguim", drop = FALSE) %>%
  #     {mutate(., tiempo = sequence(rle(.$id)$lengths)-1)} %>%
  #     mutate(evento = if_else(tiempo == seguim - 1 & (iam == 1 | ictus == 1),
  #                             1, 0),
  #            tiemposq = tiempo^2)
  #   
  #   modelo_hazards <- glm(
  #     evento == 0 ~ iniciador + I(iniciador*tiempo) + I(iniciador*tiemposq) + tiempo + tiemposq, 
  #     family = binomial(),
  #     weights = pes_2,
  #     data = base_larga
  #   )
  #   
  #   trat0 <- data.frame(cbind(seq(0, 60),0,(seq(0, 60))^2))
  #   trat1 <- data.frame(cbind(seq(0, 60),1,(seq(0, 60))^2))
  #   
  #   colnames(trat0) <- c("tiempo", "iniciador", "tiemposq")
  #   colnames(trat1) <- c("tiempo", "iniciador", "tiemposq")
  #   
  #   trat0 <- trat0 %>%
  #     {mutate(., p_noevent0 = predict(modelo_hazards, ., type = "response"))} %>%
  #     mutate(surv0 = cumprod(p_noevent0),
  #            risk0 = 1 - surv0)
  #   
  #   trat1 <- trat1 %>%
  #     {mutate(., p_noevent1 = predict(modelo_hazards, ., type = "response"))} %>%
  #     mutate(surv1 = cumprod(p_noevent1),
  #            risk1 = 1 - surv1)
  #   
  #   graf_hazards <- merge(trat0, trat1, by = c("tiempo", "tiemposq"))
  #   
  #   ggplot(graf_hazards, aes(x = tiempo, y = risk)) +
  #     geom_line(aes(y = risk0, colour = "Non-initiators")) +
  #     geom_line(aes(y = risk1, colour = "Initiators")) +
  #     xlab("Months") +
  #     scale_x_continuous(limits = c(0, 60), breaks=seq(0,60,12)) +
  #     scale_y_continuous(limits=c(0, max(c(graf_hazards$risk0, graf_hazards$risk1)) + 0.1),
  #                        breaks=seq(0, max(c(graf_hazards$risk0, graf_hazards$risk1)) + 0.1, 0.1)) +
  #     ylab("Risk") +
  #     ggtitle("Cumulative risk from hazards model") +
  #     labs(colour="A:") +
  #     theme_bw() +
  #     theme(legend.position="bottom")
  # }
  # graf <- grafico_ajustado(base)
  
  return(list(descriptivo = descriptivo, 
              Protocolo = prot))
  
}
resultados <- analisis_prot(base_larga_trabajo, 10)
