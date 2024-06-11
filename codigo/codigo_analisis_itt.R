# Análisis por intención de tratar

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

set.seed(1234)
# 1. Cargo y modifico las bases y las entradas a la función ---------------
load("datos/ensayos_itt")
load("datos/base_prueba")
load("datos/basescarhes.rdata")

base_analisis <- base_trabajo
n_muestras <- 10

confusores <- c("hdl", "ldl", "colest", "nhdl", "glucosa", "sexo_fem", "edad", "n_analiticas", "diabetes", "hta")
titulos_tabla <- list(
  sexo_fem = c("Sex", "Male", "Female"),
  edad = "Median age, in years (IQR)",
  colest = "Total cholesterol, in mg/dl (IQR)",
  hdl = "HDLc, in mg/dl (IQR)",
  ldl = "LDLc, in mg/dl (IQR)",
  nhdl = "Non-HDLc, in mg/dl (IQR)",
  glucosa = "Glucose, in mg/dl (IQR)",
  n_analiticas = "Number of blood test in the previous 6 months (IQR)",  
  diabetes = c("Diabetes", "No", "Yes"),
  hta = c("Hypertension", "No", "Yes")
)

base_perdidas <- Indice %>%
  filter(!is.na(fecha_salida)) %>%
  mutate(muerte = if_else(motivo_baja == "FA", 1, 0)) %>%
  select(patient_id, fecha_salida, muerte)

formulas <- c(
  "iniciador ~ sexo_fem + edad + I(edad*edad) + cohorte",
  "iniciador ~ sexo_fem + 
  edad + I(edad*edad) + 
  ldl + I(ldl*ldl) +
  hdl + I(hdl*hdl) +
  colest + I(colest*colest) +
  glucosa + I(glucosa*glucosa) +
  n_analiticas + I(n_analiticas * n_analiticas) +
  diabetes + hta + as.factor(cohorte)")


# 2. Defino la función ----------------------------------------------------
analisis_tte <- function(
    base_analisis, 
    base_perdidas, 
    confusores,
    titulos_tabla,
    formulas,
    n_muestras = 500,
    formula_grafico = length(formulas),
    tipo_grafico = "risk"
    ){
  
  base <- base_analisis
  
  #############################################################################
  #################### F1. ANALISIS DESCRIPTIVO ###############################
  #############################################################################
    
  # F1.1 Defino la función que lo realiza
  analisis_descriptivo_tte <- function(base){
    
    
    
    
    
    # F1.1.1 --- Figura---

    # Número de sujetos elegibles:
    r_1 <- base %>%
      distinct(patient_id) %>%
      nrow()
    
    # Número de sujetos asignados al tratamiento
    r_2 <- base %>%
      filter(iniciador == 1) %>%
      distinct(patient_id) %>%
      nrow()
    
    # Número de sujetos no asignados al tratamiento
    r_3 <- base %>%
      filter(iniciador == 0) %>%
      distinct(patient_id) %>%
      nrow()
    
    # Número de eventos reales
    # Identifico los eventos
    fechas_no_eventos <- c("fecha_inicio", "fecha_fin")
    fechas_totales <- names(base)[which(grepl("fecha_", names(base)))]
    fechas_eventos <- setdiff(fechas_totales, fechas_no_eventos)
    eventos <- substring(fechas_eventos, 7)

    # Para cada evento sumo los que hay
    
    # En iniciadores
    r_4 <- 0
    for(e in eventos){
      ne <- base[which(base[["iniciador"]] == 1 & base[[e]] == 1),]
      ne <- distinct(ne, patient_id) %>%
        nrow()
      r_4 <- r_4 + ne
    }    
  
    # En no iniciadores
    r_5 <- 0
    for(e in eventos){
      ne <- base[which(base[["iniciador"]] == 0 & base[[e]] == 1),]
      ne <- distinct(ne, patient_id) %>%
        nrow()
      r_5 <- r_5 + ne
    }    
    
    # Número de sujetos-ensayo que suponen los iniciadores
    r_6 <- base %>%
      filter(iniciador == 1) %>%
      nrow()
    
    # Número de sujetos-ensayo que suponen los no iniciadores
    r_7 <- base %>%
      filter(iniciador == 0) %>%
      nrow()
    
    # Creo una tabla con los resultados
    v1 <- c("Nº de sujetos elegibles:",
            "Nº de sujetos asignados al tratamiento:",
            "Nº de sujetos no asignados al tratamiento:",
            "Nº de eventos reales en iniciadores:",
            "Nº de eventos reales en no iniciadores:",
            "Nº de sujetos-ensayo iniciadores:",
            "Nº de sujetos-ensayo no iniciadores"
    )
    v2 <- c(r_1, r_2, r_3, r_4, r_5, r_6, r_7)
    figura_1 <- data.frame(v2, row.names = v1)
    names(figura_1) <- NULL
    # --- Fin de creación de la figura ---
    
    
    
    
    
    # F1.1.2 --- Tabla ---
    tabla_base <- base %>%
      mutate(cohorte = year(fecha_inicio))
    
    # identifico las variables categóricas binarias
    bin <- c()
    for(c in confusores){
      
      min <- min(base[[c]], na.rm = TRUE)
      max <- max(base[[c]], na.rm = TRUE)
      
      if(min == 0 & max == 1) {
        bin <- c(bin, c)
      }
    }
    
    # Defino funciones para:
      # Variables categóricas
    tabla_categorica <- function(var){
      if(var %in% bin){
        ncats <- 2
        cats <- c(0,1)
      } else{
        ncats <- max(tabla_base[[var]], na.rm = TRUE) - min(tabla_base[[var]], na.rm = TRUE) + 1
        cats <- min(tabla_base[[var]], na.rm = TRUE):max(tabla_base[[var]], na.rm = TRUE)
      }
      
      v1 <- titulos_tabla[[var]]
      
      linea_cat <- function(n, i){
        a1 <- tabla_base[which(tabla_base[["iniciador"]] == i & tabla_base[[var]] == n),]
        nrow(a1)
      }
      
      v2 <- c(NA, sapply(cats, linea_cat, 1))
      v3 <- paste(round(v2/sum(v2, na.rm = TRUE)*100, 1), "%")
      v4 <- c(NA, sapply(cats, linea_cat, 0))
      v5 <- paste(round(v4/sum(v4, na.rm = TRUE)*100, 1), "%")
      smd <- round(abs(smd(as.factor(tabla_base[[var]]), as.factor(tabla_base$iniciador ))$estimate),2)
      v6 <- c(smd, rep(NA, ncats))
      data.frame(v1, v2, v3, v4, v5, v6) 
    }
    
      # Variables numéricas (Mediana e IQR)
    tabla_mediana <- function(var){
      
      var <- as.symbol(var)
      
      t <- tabla_base %>%
        group_by(iniciador) %>%
        summarize(median = median({{ var }}),
                  q25  = quantile({{ var }}, 0.25),
                  q75 = quantile({{ var }}, 0.75)) %>%
        ungroup() %>%
        dplyr::select(iniciador , median, q25, q75) 
      
      v1 <- c(titulos_tabla[[var]])
      v2 <- c(round(t$median[which(t$iniciador == 1)]))
      v3 <- c(paste0(round(t$q25[which(t$iniciador == 1)]), 
                     "-",
                     round(t$q75[which(t$iniciador == 1)])))
      v4 <- c(round(t$median[which(t$iniciador == 0)]))
      v5 <- c(paste0(round(t$q25[which(t$iniciador == 0)]), 
                     "-",
                     round(t$q75[which(t$iniciador == 0)])))
      
      v6 <- c(abs(round(smd(tabla_base[[paste(substitute(var))]], as.factor(tabla_base$iniciador ))$estimate, 2)))
      data.frame(v1, v2, v3, v4, v5, v6)
    }
    
    # para saber que función aplicar, busco el número de valores distintos:
    # Si más de 10, lo considero como numérica.
    cats <- c()
    num <- c()
    
    for(c in confusores){
      l <- length(unique(base[[c]]))
      if(l <= 10){
        cats <- c(cats, c)
      } else {
        num <- c(num, c)
      }
    }
    
    # A cada tipo de variable le aplico la función que le corresponde
    tablas_cat <- lapply(cats, tabla_categorica) %>%
      bind_rows()
    
    tablas_num <- lapply(num, tabla_mediana) %>%
      bind_rows()

    # Uno las bases y depuro la tabla
    tabla_1 <- bind_rows(tablas_cat,
                         tablas_num) %>%
      mutate(across(.cols = c("v3", "v5"), 
                    .fns = ~if_else(.x == "NA %", NA_character_, .x)))
    
    names(tabla_1) <- c("Characteristic", "Initiators", 
                        paste0("N = ", r_6, " (", r_2, " unique individuals)"), 
                        "Non-initiators", 
                        paste0("N = ", r_7, " (", r_3, " unique individuals)"), 
                        "SMD")
    
    # --- Fin de creación de la tabla ---
    
    
    
    
    
    # F1.1.3 --- Otros resultados ---
    # Número de eventos totales
    et_1 <- 0
    for(e in eventos){
      ne <- base[which(base[["iniciador"]] == 1 & base[[e]] == 1),]
      ne <- nrow(ne)
      et_1 <- et_1 + ne
    }    
    
    et_0 <- 0
    for(e in eventos){
      ne <- base[which(base[["iniciador"]] == 0 & base[[e]] == 1),]
      ne <- nrow(ne)
      et_0 <- et_0 + ne
    }   
    
    # Número de eventos reales
    er_1 <- r_4
    er_0 <- r_5
    
    # Tiempo de seguimiento (mediana e IQR)
    t_1 <- base %>%
      mutate(seguim = round((fecha_fin - fecha_inicio)/30.4375)) %>%
      filter(iniciador == 1) %>%
      {paste0(median(.$seguim), " (", quantile(.$seguim, .25), "-", quantile(.$seguim, .75), ")")}
    
    t_0 <- base %>%
      mutate(seguim = round((fecha_fin - fecha_inicio)/30.4375)) %>%
      filter(iniciador == 0) %>%
      {paste0(median(.$seguim), " (", quantile(.$seguim, .25), "-", quantile(.$seguim, .75), ")")}

    # Tipo de evento:
    # Con un bucle, voy añadiendo cada tipo de evento y su porcentaje con respecto al 
    # total de eventos, para iniciadores y no iniciadores.
    vector_ini <- c()
    for(e in eventos){
      indi <- base[which(base[["iniciador"]] == 1 & base[[e]] == 1), ]
      indi <- distinct(indi, patient_id) %>% nrow()
      indi <- paste0(indi, " (", round(indi/er_1*100), "%)")
      vector_ini <- c(vector_ini, indi)
    }
    
    vector_no_ini <- c()
    for(e in eventos){
      indi <- base[which(base[["iniciador"]] == 0 & base[[e]] == 1), ]
      indi <- distinct(indi, patient_id) %>% nrow()
      indi <- paste0(indi, " (", round(indi/er_0*100), "%)")
      vector_no_ini <- c(vector_no_ini, indi)
    }
    
    # Creo una tabla que agrupe los resultados
    resultados <- c("Eventos totales",
                    "Eventos reales",
                    "Seguimiento",
                    eventos)
    
    iniciadores <- c(et_1, er_1, t_1, vector_ini)
    no_iniciadores <- c(et_0, er_0, t_0, vector_no_ini)
    
    resto <- data.frame(iniciadores, no_iniciadores)
    row.names(resto) <- resultados
    
    # --- Fin de "otros resultados" ---  
    
    
    
    
    # F1.1.4 Pérdidas en el seguimiento 

    # identifico perdidos en el seguimiento
    ids <- base[which(base$fecha_fin != max(base$fecha_fin)),]
    for(e in eventos){
      ids <- ids[which(ids[[e]] == 0),]
    }
    ids <- distinct(ids, patient_id)
    
    # Identifico muertes durante el seguimiento
    muertes <- base_perdidas %>%
      filter(patient_id %in% ids$patient_id,
             muerte == 1) 
    
    # Creo un texto final que resuma el total de pérdidas y las que son por fallecimiento
    perdidas_fin <- paste0("Total de pérdidas: ", nrow(ids), 
                           " (", round(nrow(ids) / r_1 * 100, 1), "%); ",
                           "por fallecimiento: ", nrow(muertes), 
                           " (", round(nrow(muertes) / r_1 * 100, 1), "%)")
    
    # --- Fin de "pérdidas" ---  
    
    
    
    
    # F1.1.5 Devuelvo en una lista los 4 resultados calculados
    return(list(figura = figura_1, 
                tabla = tabla_1, 
                otros = resto, 
                perdidas = perdidas_fin))
    
  }
  
  # F1.2 Aplico la función que acabo de definir a la base
  descriptivo <- analisis_descriptivo_tte(base)
  
  
  #############################################################################
  ############## F2. ANALISIS POR INTENCION DE TRATAR #########################
  #############################################################################
  
  # F2.1 Defino la función
  analisis_itt_tte <- function(base, n_muestras){
    
    
    # F2.1.2 Defino dentro de la funcion_boot el análisis a realizar, que se repetirá
    #        para calcular los intervalos de confianza.
    funcion_boot <- function(base, indices) {
      
      
      
      # F2.1.2.1 Muestra aleatoria de pacientes (para el bootstrap)
      patient_ids_unicos <- unique(base$patient_id)
      muestra_patient_ids <- patient_ids_unicos[indices]
      b <- base %>% filter(patient_id %in% muestra_patient_ids)
      
      
      
      
      # F2.1.2.2 Calculo los pesos IP para los análisis ajustados
      # La función principal admite varias fórmulas, así que por cada fórmula calcularemos un denominador
      
      # Primero estimo el numerador (común a todos)
      p_num <- glm(iniciador ~ 1, data = b, family = binomial())

      b <- b %>%
        {mutate(., num = as.numeric(predict(p_num, ., type = "response")))}      
      
      # Después los denominadores
      b <- b %>%
        mutate(cohorte = year(fecha_inicio))
      
      for (n in 1:length(formulas)){
        var <- paste0("peso_", n)
        f <- formulas[n]
        p_denom <- glm(f, data = b, family = binomial())
        b$denom <- predict(p_denom, b, type = "response")
        b[[var]] <- if_else(b$iniciador == 1, 
                            b$num / b$denom, 
                            (1-b$num)/(1-b$denom))
      }
      
      
      
      # F2.1.2.3 Creo la "base larga" para el análisis de regresión
      base_larga <- b %>%
        mutate(seguim = as.numeric(round((fecha_fin - fecha_inicio)/30.4375)),
               id = paste0(patient_id, fecha_inicio)) %>%
        expandRows("seguim", drop = FALSE) %>%
        {mutate(., tiempo = sequence(rle(.$id)$lengths)-1)} %>%
        mutate(tiemposq = tiempo^2,
               fecha_cal = fecha_inicio + months(tiempo),
               evento = 0) 

      # Para cada evento lo cuento sólo al final del seguimiento
      fechas_no_eventos <- c("fecha_inicio", "fecha_fin", "fecha_cal")
      fechas_totales <- names(b)[which(grepl("fecha_", names(b)))]
      fechas_eventos <- setdiff(fechas_totales, fechas_no_eventos)
      eventos <- substring(fechas_eventos, 7)
      
      for(e in eventos){
        base_larga[["evento"]] <- if_else(base_larga[["tiempo"]] == base_larga[["seguim"]] - 1 & 
                                            base_larga[[e]] == 1, 
                                          1, 
                                          base_larga[["evento"]])
      }
      
      
      
      # F2.1.2.4 Hago los análisis de regresión y los cálculos de riesgos a 5 años
      # Para cada uno de los pesos, y para ningún peso, se hace un modelo glm
      
      # Empiezo por ningún peso
      # Regresión
      modelo_hazards <- glm(
        evento == 0 ~ iniciador + I(iniciador*tiempo) + I(iniciador*tiemposq) + tiempo + tiemposq, 
        family = binomial(), 
        data = base_larga
      )
      
      # Predigo resultados en base en blanco
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
      
      # Calculo Riesgos
      r0 <- round(max(trat0$risk0)*100, 2)
      r1 <- round(max(trat1$risk1)*100, 2)
      rr <- round(r1 / r0, 2)
      rd <- r1 - r0
      
      # Guardo los resultados en tabla
      v1 <- c("Riesgo en no iniciadores:",
              "Riesgo en iniciadores:",
              "Razón de riesgos:",
              "Diferencia de riesgos:")
      tabla <- data.frame(v1 = v1, v2 = c(r0, r1, rr, rd))
      
      
      # Para cada peso (formula) uso una función para crear la tabla de riesgos
      resultados_plr <- function(n){
        var_peso <- paste0("peso_", n)
        
        # Regresión
        modelo_hazards <- glm(
          evento == 0 ~ iniciador + I(iniciador*tiempo) + I(iniciador*tiemposq) + tiempo + tiemposq, 
          family = binomial(),
          weights = base_larga[[var_peso]],
          data = base_larga
        )
        
        # Predigo resultados en base en blanco
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
        
        # Calculo riesgos
        r0 <- round(max(trat0$risk0)*100, 2)
        r1 <- round(max(trat1$risk1)*100, 2)
        rr <- round(r1 / r0, 2)
        rd <- r1 - r0

        # Almaceno los resultados en una tabla
        v1 <- c("Riesgo en no iniciadores:",
                "Riesgo en iniciadores:",
                "Razón de riesgos:",
                "Diferencia de riesgos:")
        
        return(data.frame(v1 = v1, v2 = c(r0, r1, rr, rd)))
        
      }
      
      # Aplico la función a los pesos o fórmulas
      resultados <- lapply(1:length(formulas), resultados_plr)
      
      
      
      
      # F2.1.2.5 Mediante un bucle por las fómrulas, almaceno los resultados en 
      # un único vector (r), que devuelvo con la funcion_boot
      
      r <- c(tabla$v2)
      for(n in 1:length(formulas)){
        r <- c(r, resultados[[n]]$v2)
      }
      
      return(r)
      
    }
    
    
    
    # F2.1.3 Aplico boot
    boot_results <- boot(data = base, statistic = funcion_boot, R = n_muestras)
    
    
    
    # F2.1.4 Creo una tabla con los resultados, sacándolos de boot_results
    texto_formulas <- paste0("formula ", 1:length(formulas))
    v1 <- c("No ajustado", texto_formulas)

    # Resultados originales
    v2 <- boot_results$t0[seq(from = 1, by = 4, length.out = length(formulas) + 1)]
    v3 <- boot_results$t0[seq(from = 2, by = 4, length.out = length(formulas) + 1)]
    v4 <- boot_results$t0[seq(from = 3, by = 4, length.out = length(formulas) + 1)]
    v5 <- boot_results$t0[seq(from = 4, by = 4, length.out = length(formulas) + 1)]
    
    # Para extraer los intervalos de confianza, creo una función
    extr_ci_inf <- function(i){
      boot.ci(boot_results, type = "perc", index = i)$percent[4]
    }
    
    extr_ci_sup <- function(i){
      boot.ci(boot_results, type = "perc", index = i)$percent[5]
    }
    
    v2.1 <- sapply(seq(from = 1, by = 4, length.out = length(formulas) + 1), extr_ci_inf)
    v2.2 <- sapply(seq(from = 1, by = 4, length.out = length(formulas) + 1), extr_ci_sup)
    
    v3.1 <- sapply(seq(from = 2, by = 4, length.out = length(formulas) + 1), extr_ci_inf)
    v3.2 <- sapply(seq(from = 2, by = 4, length.out = length(formulas) + 1), extr_ci_sup)
    
    v4.1 <- sapply(seq(from = 3, by = 4, length.out = length(formulas) + 1), extr_ci_inf)
    v4.2 <- sapply(seq(from = 3, by = 4, length.out = length(formulas) + 1), extr_ci_sup)
    
    v5.1 <- sapply(seq(from = 4, by = 4, length.out = length(formulas) + 1), extr_ci_inf)
    v5.2 <- sapply(seq(from = 4, by = 4, length.out = length(formulas) + 1), extr_ci_sup)

    tabla <- data.frame(v2, v2.1, v2.2, v3, v3.1, v3.2, v4, v4.1, v4.2, v5, v5.1, v5.2) %>%
      mutate(v2 = paste0(round(v2, 2), "% (", round(v2.1, 2), "-", round(v2.2, 2), "%)"),
             v3 = paste0(round(v3, 2), "% (", round(v3.1, 2), "-", round(v3.2, 2), "%)"),
             v4 = paste0(round(v4, 2), " (", round(v4.1, 2), "-", round(v4.2, 2), ")"),
             v5 = paste0(round(v5, 2), "% (", round(v5.1, 2), "-", round(v5.2, 2), "%)")) %>%
      select(v2, v3, v4, v5)
    
    colnames(tabla) <- c("Risk in non-initiators", "Risk in initiators", "RR", "RD")
    row.names(tabla) <- v1
    
    
    
    
    # F2.1.5 Devuelvo la tabla con los resultados
    return(tabla)
    
  }
  
  
  # F2.2 Aplico la función
  itt <- analisis_itt_tte(base, n_muestras)
  
  
  
  #############################################################################
  ############## F3. GRÁFICO POR INTENCIÓN DE TRATAR AJUSTADO #################
  #############################################################################
  
  # F3. 1 Defino la función
  grafico_ajustado <- function(base){
    
    formula <- formulas[formula_grafico]
    base <- base %>%
      mutate(cohorte = year(fecha_inicio))  
    
    
    # Calculo los pesos IP
    p_num <- glm(iniciador ~ 1, data = base, family = binomial())
    p_denom <- glm(formula, data = base, family = binomial())
    base <- base %>%
      {mutate(., denom = predict(p_denom, ., type = "response"),
              num = predict(p_num, ., type = "response"))} %>%
      mutate(peso = if_else(iniciador == 1, 
                             num/denom, 
                             (1-num)/(1-denom)))

    
    # Creo la "base larga" para el análisis de regresión
    base_larga <- base %>%
      mutate(seguim = as.numeric(round((fecha_fin - fecha_inicio)/30.4375)),
             id = paste0(patient_id, fecha_inicio)) %>%
      expandRows("seguim", drop = FALSE) %>%
      {mutate(., tiempo = sequence(rle(.$id)$lengths)-1)} %>%
      mutate(tiemposq = tiempo^2,
             fecha_cal = fecha_inicio + months(tiempo),
             evento = 0) 
    
    fechas_no_eventos <- c("fecha_inicio", "fecha_fin", "fecha_cal")
    fechas_totales <- names(base_larga)[which(grepl("fecha_", names(base_larga)))]
    fechas_eventos <- setdiff(fechas_totales, fechas_no_eventos)
    eventos <- substring(fechas_eventos, 7)
    
    for(e in eventos){
      base_larga[["evento"]] <- if_else(base_larga[["tiempo"]] == base_larga[["seguim"]] - 1 & 
                                          base_larga[[e]] == 1, 
                                        1, 
                                        base_larga[["evento"]])
    }
    
    
    # Ajusto el modelo de regresión
    modelo_hazards <- glm(
      evento == 0 ~ iniciador + I(iniciador*tiempo) + I(iniciador*tiemposq) + tiempo + tiemposq, 
      family = binomial(),
      weights = peso,
      data = base_larga
    )
    
    
    # Aplico los resultados a la base en blanco
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
    
    
    # Dibujo el gráfico de riesgo o de supervivencia
    graf_hazards <- merge(trat0, trat1, by = c("tiempo", "tiemposq"))
    
    if(tipo_grafico == "risk"){
      ggplot(graf_hazards, aes(x = tiempo, y = risk)) +
        geom_line(aes(y = risk0, colour = "Non-initiators")) +
        geom_line(aes(y = risk1, colour = "Initiators")) +
        xlab("Months") +
        scale_x_continuous(limits = c(0, 60), breaks=seq(0,60,12)) +
        scale_y_continuous(limits=c(0, max(c(graf_hazards$risk0, graf_hazards$risk1)) + 0.1),
                           breaks=seq(0, max(c(graf_hazards$risk0, graf_hazards$risk1)) + 0.1, 0.1)) +
        ylab("Risk") +
        ggtitle("Cumulative risk from hazards model") +
        labs(colour="A:") +
        theme_bw() +
        theme(legend.position="bottom")
    } else if(tipo_grafico == "surv"){
      
      ggplot(graf_hazards, aes(x = tiempo, y = surv)) +
        geom_line(aes(y = surv0, colour = "Non-initiators")) +
        geom_line(aes(y = surv1, colour = "Initiators")) +
        xlab("Months") +
        scale_x_continuous(limits = c(0, 60), breaks=seq(0,60,12)) +
        scale_y_continuous(limits=c(min(c(graf_hazards$surv0, graf_hazards$surv1)) - 0.1, 1),
                           breaks=seq(min(c(graf_hazards$surv0, graf_hazards$surv1)) - 0.1, 1, 0.1)) +
        ylab("Survival") +
        ggtitle("Survival from hazards model") +
        labs(colour="A:") +
        theme_bw() +
        theme(legend.position="bottom")
    }

  }
  
  # F3.2 Aplico la función
  graf <- grafico_ajustado(base)
  
  
  
  # F4. Devuelvo los resultados
  return(list(descriptivo = descriptivo, 
              ITT = itt,
              grafico = graf))
  
}

lista <- ls()

ba <- base_trabajo
bp <- base_perdidas
c <- confusores
tt <- titulos_tabla
f <- formulas
n <- n_muestras
fg <- 2

rm(list = lista)

# 3. Aplico la función ----------------------------------------------------
resultados <- analisis_tte(
  base_analisis = ba, 
  base_perdidas = bp, 
  confusores = c,
  titulos_tabla = tt,
  formulas = f,
  n_muestras = 5,
  formula_grafico = fg,
  tipo_grafico = "risk"
)
