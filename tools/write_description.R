source(here::here("shared", "functions", "libraries.R"))

desc_content <- paste0(
  "Type: project\n",
  "Package: scRNAseqTools\n",
  "Title: Internal analysis pipeline (not a real package)\n",
  "Version: ", sub("-dev$", ".9000", trimws(readLines(here::here("VERSION"), warn = FALSE)[1])), "\n",
  "Imports:\n    ",
  paste(library_registry$package, collapse = ",\n    "),
  "\n"
)
writeLines(desc_content, here::here("DESCRIPTION"))
