pacman::p_load(
  dplyr,
  tidyr,
  rio,
  skimr,
  lubridate,
  miceadds,
  splitstackshape
)


load("datos/basesCARHES.Rdata")
load("datos/ensayos_itt")

crear_base_protocolo <- function(
    base_itt,
    base_farma,
    base_confusores_1,
    base_confusores_2,
    gap_max = 61,
    lavado = 180){
  
  n_pasos <- length(base_confusores_1) + length(base_confusores_2) + 1 + 1
  l_paso <- round(100/n_pasos)
  prog <- 0
  
  base <- base_itt
  
  # Identifico los tratados que dejan el tratamiento
  interrupciones <- base_farma %>%
    filter(patient_id %in% base$patient_id,
           !is.na(fecha_dispensacion)) %>%
    arrange(patient_id, fecha_dispensacion) %>%
    mutate(gap = as.numeric(lead(fecha_dispensacion) - fecha_dispensacion),
           gap = if_else(lead(patient_id) != patient_id, NA, gap)) %>%
    filter(gap_max < gap | is.na(gap) ) %>%
    mutate(fecha_interrupcion = floor_date(fecha_dispensacion + gap_max, "months"),
           iniciador = 1) %>%
    select(patient_id, fecha_interrupcion, iniciador)
  
  # Identifico las fechas de inicio de tratamiento de los no tratados
  inicios <- base_farma %>%
    arrange(patient_id) %>%
    filter(patient_id %in% base$patient_id,
           !is.na(fecha_dispensacion)) %>%
    mutate(gap = as.numeric( fecha_dispensacion - lag(fecha_dispensacion)),
           gap = if_else(lag(patient_id) != patient_id, NA, gap)) %>%
    filter(lavado < gap | is.na(gap) ) %>%
    mutate(fecha_interrupcion = ceiling_date(fecha_dispensacion, "months"),
           iniciador = 0) %>%
    select(patient_id, fecha_interrupcion, iniciador)
  
  interrupciones <- bind_rows(interrupciones, inicios)
  
  base <- base %>%
    left_join(interrupciones, by = c("patient_id", "iniciador"), relationship = "many-to-many") %>%
    mutate(fecha_interrupcion = if_else(fecha_interrupcion < fecha_inicio, NA, fecha_interrupcion)) %>%
    arrange(fecha_interrupcion) %>%
    distinct(patient_id, fecha_inicio, .keep_all = TRUE) %>%
    # Ajustar la fecha de fin según la fecha de interrupción
    mutate(fecha_fin = pmin(fecha_fin, fecha_interrupcion, na.rm = TRUE))
  
  fechas_no_eventos <- c("fecha_inicio", "fecha_fin", "fecha_interrupcion")
  fechas_totales <- names(base_itt)[which(grepl("fecha_", names(base_itt)))]
  fechas_eventos <- setdiff(fechas_totales, fechas_no_eventos)
  eventos <- substring(fechas_eventos, 7)
  
  prog <- prog + l_paso
  print(paste0(prog, "%...........fecha de fin modificada"))
  
  # Según fecha de fin de seguimiento, modifico evento
  for(e in eventos){
    fecha_e <- paste0("fecha_", e)
    base[[e]] <- if_else(base[["fecha_fin"]] != base[[fecha_e]] & !is.na(base[[fecha_e]]), 0, base[[e]])
  }
  
  # Una vez modificada la fecha de fin de seguimiento, se crea la "base larga"
  base_larga <- base %>%
    mutate(seguim = as.numeric(round((fecha_fin - fecha_inicio)/30.4375)),
           id = paste0(patient_id, fecha_inicio)) %>%
    expandRows("seguim", drop = FALSE) %>%
    {mutate(., tiempo = sequence(rle(.$id)$lengths)-1)} %>%
    mutate(tiemposq = tiempo^2,
           fecha_cal = fecha_inicio + months(tiempo),
           evento = 0) 
  
  for(e in eventos){
    base_larga[["evento"]] <- if_else(base_larga[["tiempo"]] == base_larga[["seguim"]] - 1 & 
                                        base_larga[[e]] == 1, 
                                      1, 
                                      base_larga[["evento"]])
  }
  
  prog <- prog + l_paso
  print(paste0(prog, "%...........base larga creada"))
  # 2. Incluir los confusores (los tratamientos no hacen falta porque el 
  # seguimiento ya se ha puesto en base a la continuidad)
  
  # Hay que seleccionar las variables previamente en la lista "base_confusores_1"
  confusores <- names(base_confusores_1)
  
  
  for(confusor_fun in confusores){
    confusor_fun_t <- paste0(confusor_fun, "_t")
    confusor <- base_confusores_1[[confusor_fun]] %>%
      filter(patient_id %in% base$patient_id) %>%
      mutate(fecha_cal = floor_date(fecha, "months")) %>%
      distinct(patient_id, fecha_cal, .keep_all = TRUE) %>%
      dplyr::select(patient_id, valor, fecha_cal)
    
    names(confusor) <- c("patient_id", paste0(confusor_fun, "_t"), "fecha_cal")
    
    base_larga <- base_larga %>% 
      left_join(confusor, by = c("patient_id", "fecha_cal"))
    
    base_larga[[confusor_fun_t]] <- if_else(base_larga[["tiempo"]] == 0,
                                            base_larga[[confusor_fun]],
                                            base_larga[[confusor_fun_t]])
    base_larga <- base_larga %>%
      fill(all_of(confusor_fun_t))
    
    prog <- prog + l_paso
    print(paste0(prog, "%...........variable ",confusor_fun_t ," incorporada"))
    
    
  }
  confusores <- names(base_confusores_2)
  for(confusor_fun in confusores){
    confusor_fun_t <- paste0(confusor_fun, "_t")
    base_larga <- base_larga %>% 
      left_join(base_confusores_2[[confusor_fun]], by = c("patient_id", "fecha_cal" = "fecha")) %>%
      mutate(valor = if_else(is.na(valor), 0, valor))
    
    names(base_larga)[which(names(base_larga) == "valor")] <- confusor_fun_t
    
    base_larga[[confusor_fun_t]] <- if_else(base_larga[["tiempo"]] == 0,
                                            base_larga[[confusor_fun]],
                                            base_larga[[confusor_fun_t]])
    
    prog <- prog + l_paso
    print(paste0(prog, "%...........variable ",confusor_fun_t ," incorporada"))
  }
  
  
  return(base_larga)
  
}

confusores_1_prot <- confusores_1[c("ldl", "hdl", "colest", "nhdl", "glucosa", "diabetes", "hta")]

base_larga_trabajo <- crear_base_protocolo(
  base_itt = base_trabajo,
  base_farma = DispensacionC10_HL,
  base_confusores_1 = confusores_1_prot_2,
  base_confusores_2 = confusores_2,
  gap_max = 61)

save(base_larga_trabajo, file = "datos/base_prueba_prot")  
