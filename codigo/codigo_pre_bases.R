pacman::p_load(
  dplyr,
  tidyr,
  rio,
  skimr,
  lubridate,
  miceadds,
  splitstackshape
)

# Preparación de las bases de datos para la función de crear base
load("original/DispensacionC10_HL.RDAta")
load("original/Indice_HL.RDAta")
load("original/MACE_HL.RDAta")

load("original/colesterol_HDL.RDAta")
load("original/colesterol_LDL.RDAta")
load("original/colesterol_largo.RData")

load("original/Datos_tabaquismo.RDAta")
load("original/tabaco_situacion.RDAta")

load("original/glucosa.RDAta")
load("original/hemoglobina_glicada.RDAta")

load("original/peso.RDAta")
load("original/talla.RDAta")

load("original/TAD.RDAta")
load("original/TAS.RDAta")

load("original/Indice_HL2.RData")
load("original/Tratamientos.RData")


# 1.3 Formateo de los datos -----------------------------------------------
peso_PHL$valor <- as.numeric(peso_PHL$valor)
talla_PHL$valor <- as.numeric(talla_PHL$valor)

# 1.2 Validación de los datos (aplico límites admitidos) -------------------------------------
TAD_PHL <- filter(TAD_PHL, 30 < valor & valor < 160 )
TAS_PHL <- filter(TAS_PHL, 50 < valor & valor < 280 )
glucosa_PHL <- filter(glucosa_PHL, 50 < valor & valor < 400 )
hemoglobina_glicada_PHL <- filter(hemoglobina_glicada_PHL, 4 < valor & valor < 20 )
colesterol_PHL <- filter(colesterol_PHL, 43 < valor & valor < 970 )
colesterol_HDL_PHL <- filter(colesterol_HDL_PHL, 9 < valor & valor < 200 )
colesterol_LDL_PHL <- filter(colesterol_LDL_PHL, 30 < valor & valor < 500 )
peso_PHL <- filter(peso_PHL, 25 < valor & valor < 180 )
talla_PHL <- filter(talla_PHL, 140 < valor & valor < 216 )

TA_PHL <- merge(TAD_PHL, TAS_PHL, all = TRUE, by = c("patient_id", "fecha")) %>%
  filter(valor.y > valor.x | is.na(valor.y) | is.na(valor.x))

TAD_PHL <- dplyr::select(TA_PHL, patient_id, fecha, valor = valor.x)
TAS_PHL <- dplyr::select(TA_PHL, patient_id, fecha, valor = valor.y)

colesterol_NHDL <- merge(colesterol_PHL, colesterol_HDL_PHL, all = TRUE, by = c("patient_id", "fecha")) %>%
  mutate(valor =  valor.x - valor.y)

# Preparación de las bases para la función de crear ensayos  --------------------------
# Creo base grande de analíticas, como tarda mucho en hacerse, la guardo para futuras lecturas
analiticas <- colesterol_HDL_PHL %>%
  merge(colesterol_LDL_PHL, by = c("patient_id", "fecha"), all = TRUE) %>%
  merge(glucosa_PHL, by = c("patient_id", "fecha"), all = TRUE) %>%
  merge(hemoglobina_glicada_PHL, by = c("patient_id", "fecha"), all = TRUE) %>%
  merge(TAD_PHL, by = c("patient_id", "fecha"), all = TRUE) %>%
  merge(TAS_PHL, by = c("patient_id", "fecha"), all = TRUE)

save(analiticas, file = "datos/analiticas.RData")

# la "base_envento" tiene que tener una variable llamada fecha_evento y otra llamada evento con
# el nombre del evento en texto
base_MACE <- MACE_HL %>% rename(fecha_evento = Fecha_MACE1) %>%
  mutate(evento = if_else(grepl("I21", CIE_MACE1), "iam", "ictus")) %>%
  dplyr::select(patient_id, evento, fecha_evento = fecha_evento)

# la "base_indice" tiene que tener una variable llamada fecha_entrada y otra fecha_salida
Indice <- Indice_HL %>% rename(fecha_entrada = fecha_FRCV, 
                               fecha_salida = fecha_baja_bdu)

# Tenemos que tener un objeto llamado "confusores". Es una lista que incluye varias bases de 
# datos con los confusores. El nombre del elemento en la lista será el del confusor en la base
# Tienen que tener al menos 3 variables: patient_id, valor y fecha

diabeticos <- Indice_HL2 %>%
  filter(DM == "SI") %>%
  mutate(valor = 1) %>%
  rename(fecha = fecha_FRCV_DM) %>%
  select(patient_id, valor, fecha)

hipertensos <- Indice_HL2 %>%
  filter(HTA == "SI") %>%
  mutate(valor = 1) %>%
  rename(fecha = fecha_FRCV_HTA) %>%
  select(patient_id, valor, fecha)

suma_analiticas_6m <- function(fecha_fun, analiticas) {
  analiticas <- analiticas %>%
    select(patient_id, fecha) %>%
    filter(fecha <= fecha_fun,
           fecha_fun - 180 <= fecha) %>%
    summarize(valor = n(), .by = patient_id) %>%
    mutate(fecha = fecha_fun)
  print(fecha_fun)
  return(analiticas)
}

analiticas_6m <- lapply(seq.Date(from = round_date(min(analiticas$fecha) + months(6), "months"), 
                                 to = round_date(max(analiticas$fecha), "months"), 
                                 by = "months"), 
                        suma_analiticas_6m, analiticas) %>%
  bind_rows()

tabaquismo <- Datos_tabaquismo_PHL %>%
  mutate(valor = if_else(Tabaquismo == "N", 0, 1))

confusores_1 <- list(hdl = colesterol_HDL_PHL, 
                     ldl = colesterol_LDL_PHL,
                     nhdl = colesterol_NHDL,
                     colest = colesterol_PHL,
                     glucosa = glucosa_PHL,
                     glicada = hemoglobina_glicada_PHL,
                     altura = talla_PHL,
                     peso = peso_PHL,
                     tas = TAS_PHL,
                     tad = TAD_PHL,
                     diabetes = diabeticos,
                     hta = hipertensos,
                     tabaco = tabaquismo
)

confusores_2 <- list(n_analiticas = analiticas_6m)

save(Indice, DispensacionC10_HL, base_MACE, confusores_1, confusores_2, file = "datos/basesCARHES.Rdata")
