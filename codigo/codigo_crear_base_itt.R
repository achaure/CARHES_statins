pacman::p_load(
  dplyr,
  tidyr,
  rio,
  skimr,
  lubridate,
  miceadds,
  splitstackshape
)

# 1. Lectura de bases de datos -------------------------------------------
load("datos/basesCARHES.Rdata")


# 2. Función y creación de base para ITT -------------------------------------
crear_bases <- function(
    base_indice,
    base_farma,
    base_evento,
    base_confusores_1,
    base_confusores_2,
    fecha_inicio, 
    fecha_fin,
    fecha_fin_segui,
    lavado = 180
    ){
  fechas <- seq.Date(from = dmy(fecha_inicio), to = dmy(fecha_fin), by = "month")
  crear_base_ensayo <- function(fecha){
    
    print(format(fecha, "%B de %Y"))
    fecha_ini_seg <- floor_date(fecha + 32, "month")
    
    # 1. SELECCION DE LOS SUJETOS 
    # A. No prescripción de estatinas en los 6 meses previos.
    excluidos_farma <- base_farma %>%
      filter(fecha_dispensacion >= fecha - lavado & fecha_dispensacion < fecha) %>%
      distinct(patient_id) %>%
      {.$patient_id}
    
    # B. Sujetos sin evento cardiovascular previo
    excluidos_mace <- base_evento %>%
      filter(fecha_evento < fecha_ini_seg) %>%
      distinct(patient_id) %>%
      {.$patient_id}
    
    # Aplico las exclusiones
    base <- base_indice %>%
      filter(fecha_entrada < fecha,
             fecha_salida > fecha | is.na(fecha_salida),
             !(patient_id %in% excluidos_farma),
             !(patient_id %in% excluidos_mace)) %>%
      mutate(fecha_inicio = fecha)
    
    # Reduzco las otras bases
    farma <- base_farma %>%
      filter(patient_id %in% base$patient_id,
             fecha_dispensacion >= fecha)
    
    # Me quedo con un solo evento (el primero pero posterior al inicio) por paciente
    mace <- base_evento %>%
      filter(patient_id %in% base$patient_id,
             fecha < fecha_evento) %>%
      arrange(fecha_evento) %>%
      distinct(patient_id, .keep_all = TRUE)
    
    
    # 2. ASIGNACION DE TRATAMIENTO
    # Identifico iniciadores
    iniciadores <- farma %>%
      filter(month(fecha_dispensacion) == month(fecha),
             year(fecha_dispensacion) == year(fecha)) %>%
      distinct(patient_id, .keep_all = TRUE)
    
    base <- base %>%
      mutate(iniciador = if_else(patient_id %in% iniciadores$patient_id, 1, 0))
    
    
    # 3. RESULTADO
    # Identifico resultado
    eventos <- unique(base_evento$evento)
    
    for (evento_fun in eventos){
      solo_evento <- mace %>%
        mutate(fecha_evento =  floor_date(fecha_evento, unit = "months"),
               evento = if_else(evento == evento_fun, 1, 0)) %>%
        filter(evento == 1) %>%
        dplyr::select(patient_id, evento, fecha_evento)
      
      names(solo_evento) <- c("patient_id", evento_fun, paste0("fecha_", evento_fun))
      base <- base %>%
        left_join(solo_evento, by = "patient_id")
      
      base[[evento_fun]] <- if_else(is.na(base[[evento_fun]]), 0, base[[evento_fun]])
      
    }
    
    
    # 4.1 Confusores que rellenan valores futuros (p.ej. colesterol, si un mes no hay analítica
    # se asume que el valor es el del mes anterior)
    confusores <- names(base_confusores_1)
    
    for(confusor_fun in confusores){
      confusor <- base_confusores_1[[confusor_fun]] %>%
        filter(patient_id %in% base$patient_id,
               .$fecha < fecha_ini_seg) %>%
        arrange(desc(fecha)) %>%
        distinct(patient_id, .keep_all = TRUE) %>%
        dplyr::select(patient_id, valor)
      
      names(confusor) <- c("patient_id", confusor_fun)
      
      base <- base %>% 
        left_join(confusor, by = "patient_id")
    }

    # 4.2 Confusores que tienen valores fijos (p.ej. analíticas en los 6 últimos meses, si no están en la base
    # es porque el resultado es 0)
    
    confusores <- names(base_confusores_2)
    
    for(confusor_fun in confusores){
      base <- base %>% 
        left_join(base_confusores_2[[confusor_fun]], by = c("patient_id", "fecha_inicio" = "fecha")) %>%
        mutate(valor = if_else(is.na(valor), 0, valor))
      
      names(base)[which(names(base) == "valor")] <- confusor_fun
    }
    
    
    # 5. Fecha de fin de seguimiento
    fechas_fin <- c(paste0("fecha_", eventos))
    
    base <- base
    base$fecha_fin <- base$fecha_salida
    
    for (f in fechas_fin){
      base <- base %>%
        mutate(fecha_fin = pmin(.[["fecha_fin"]], .[[f]], na.rm = TRUE))
    }
    
    base <- base %>%
      mutate(fecha_fin = pmin(fecha_fin, dmy(fecha_fin_segui), na.rm = TRUE))

    # 6. según fecha de fin de seguimiento, modifico evento
    
    for(e in eventos){
      fecha_e <- paste0("fecha_", e)
      base[[e]] <- if_else(base[["fecha_fin"]] != base[[fecha_e]] & !is.na(base[[fecha_e]]), 0, base[[e]])
    }
    return(base)
  }

  base_total <- lapply(fechas, crear_base_ensayo) %>%
    bind_rows()
  return(base_total)
}


base_total <- crear_bases(
  base_indice = Indice,
  base_farma = DispensacionC10_HL,
  base_evento = base_MACE,
  base_confusores_1 = confusores_1,
  base_confusores_2 = confusores_2,
  fecha_inicio = "1/2/2018", 
  fecha_fin = "1/12/2022",
  fecha_fin_segui = "1/1/2023"
  ) 

save(base_total, file = "datos/base_ensayos_cruda")
# Exploro la base para ver con qué variables me quedo (me quedo con las que tengan más de un 70%)
# skim(base_total)

load("datos/base_ensayos_cruda")

# Me quedo con la variable LDL y elimino los que no la tienen
base_total <- base_total %>%
  filter(!is.na(ldl))

# Me quedo con la variable HDL y elimino los que no la tienen
base_total <- base_total %>%
  filter(!is.na(hdl))

# Me quedo con la variable glucosa y elimino los que no la tienen
base_total <- base_total %>%
  filter(!is.na(glucosa))

# Me quedo con la variable nhdl y elimino los que no la tienen
base_total <- base_total %>%
  filter(!is.na(nhdl))

# Selecciono las variables finales
base_total <- base_total %>%
  mutate(edad = as.numeric((fecha_inicio - fecha_nac )/365.25),
         sexo_fem = if_else(sexo == "MUJER", 1, 0),
         diabetes = if_else(is.na(diabetes), 0, 1),
         hta = if_else(is.na(hta), 0, 1)) %>%
  select(patient_id, fecha_inicio, fecha_fin, iniciador, iam, ictus, edad, sexo_fem, 
         ldl, hdl, colest, nhdl, glucosa, diabetes, hta, n_analiticas,
         fecha_iam, fecha_ictus)

save(base_total, file = "datos/ensayos_itt")




