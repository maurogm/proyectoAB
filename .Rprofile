
exports_code_in_Rmd_as_script <- function(path, input_file_name, output_file_name = NA, documentation = 1) {
  if (is.na(output_file_name)) output_file_name <- paste0(input_file_name, '_rscript')
  knitr::purl(paste0(path, '/', input_file_name, '.Rmd'),
              paste0(path, '/', output_file_name, '.R'),
              documentation = documentation)
}
append_archivo_secundario_a_archivo_principal <- function(path_archivo_principal, nombre_archivo_principal,
                            path_archivo_secundario, nombre_archivo_secundario,
                            sobreescribir_archivo_principal = TRUE) {
  #Si sobreescribir_archivo_principal es FALSE, escribe la nueva versi?n en el path principal con el nombre nombre_archivo_principal_bis
  archivo_principal <- fread(paste(path_archivo_principal, nombre_archivo_principal, sep = '/'), na.string="NULL", encoding="UTF-8")
  archivo_secundario <- fread(paste(path_archivo_secundario, nombre_archivo_secundario, sep = '/'), na.string="NULL", encoding="UTF-8")
  
  if (sobreescribir_archivo_principal) {
    write.table(rbind(archivo_principal, archivo_secundario),
                paste(path_archivo_principal, nombre_archivo_principal, sep = '/'),
                row.names = FALSE, sep='\t')
    } else {
    write.table(rbind(archivo_principal,archivo_secundario),
                paste(path_archivo_principal, paste('appended', nombre_archivo_principal,sep='_'), sep = '/'),
                row.names = FALSE, sep='\t') }
}
object_sizes_all <- function() {
  size_Mb <- sapply(ls(globalenv()),
                    function(x) { object.size(get(x)) / 1024^2 })
  as.data.frame(size_Mb) %>% 
    tibble::rownames_to_column('object')
}
source_all_files <- function(path_from_wd = "source") {
  for (file in list.files(path_from_wd)) {
    source(paste0(path_from_wd, "/", file), encoding = 'UTF-8')
  }
}

source_all_files()