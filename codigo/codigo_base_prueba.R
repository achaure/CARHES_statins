# Creación de la base de prueba
id1 <- filter(base_total, iniciador == 1, iam == 1) %>% {.$patient_id[1:10]}
id2 <- filter(base_total, iniciador == 1, iam == 0) %>% {.$patient_id[1:10]}
id3 <- filter(base_total, iniciador == 0, iam == 1) %>% {.$patient_id[1:10]}
id4 <- filter(base_total, iniciador == 0, iam == 0) %>% {.$patient_id[1:10]}

base_trabajo <- filter(base_total, patient_id %in% c(id1, id2, id3, id4))