install_and_check <- function(pkg, bioc = FALSE) {
    if (bioc) BiocManager::install(pkg) else install.packages(pkg)
    if (!requireNamespace(pkg, quietly = TRUE))
        stop("Failed to install: ", pkg)
}
install_and_check("limma", TRUE)
install_and_check("edgeR", TRUE) # prerequisite for sva
install_and_check("sva", TRUE)
install_and_check("readr")
install_and_check("dplyr")
install_and_check("jsonlite")
install_and_check("Rfit")
install_and_check("emmeans")
install_and_check("testthat")