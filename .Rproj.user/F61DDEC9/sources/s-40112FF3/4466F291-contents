# building the reference table of element masses
# source: http://www.ciaaw.org/atomic-masses.htm
.build_elements.tsv <- function() {

  ciaaw.df <- utils::read.table(system.file("extdata/reference_tables/vocDB_tables/ciaaw.tsv",
                                            package = "ptairMS"),
                                header = TRUE,
                                quote = "\"",
                                sep = "\t",
                                stringsAsFactors = FALSE)

  element.df <- ciaaw.df[!is.na(ciaaw.df[, "Z"]), ]

  element.df[, "mass"] <- vapply(element.df[, "mass"],
                                 function(mass.c) {

                                   if (mass.c == " 12(exact") {

                                     return(12)

                                   } else {

                                     mass.c <- substr(mass.c, 1, nchar(mass.c) - 3)
                                     mass_split.vc <- unlist(strsplit(mass.c, split = ""))
                                     mass_split.vl <- !is.na(as.numeric(mass_split.vc)) | mass_split.vc == "."
                                     mass_split.vc <- mass_split.vc[mass_split.vl]
                                     mass.n <- floor(as.numeric(paste(mass_split.vc, collapse = "")) * 1e6) / 1e6

                                     return(mass.n)

                                   }

                                 },FUN.VALUE = 1.01)

  # adding supplementary information about electron, proton, and neutron
  supp.df <- element.df[c(1,2,3), ]
  supp.df[, "Z"] <- rep(NA_integer_, 3)
  supp.df[, "symbol"] <- c("electron", "proton", "neutron")
  supp.df[, "element"] <- rep("", 3)
  supp.df[, "A"] <- rep("", 3)
  supp.df[, "mass"] <- c(0.000548, 1.0072765, 1.008665549)

  element.df <- rbind.data.frame(supp.df,
                                 element.df,
                                 stringsAsFactors = FALSE)

  utils::write.table(element.df,
                     file = "elements.tsv",
                     row.names = FALSE,
                     sep = "\t")

}

# source: de_Lacy_Costello_2014
# .build_vocDB.tsv <- function(databases = "de_Lacy_Costello_2014") {
#
#   if ("de_Lacy_Costello_2014" %in% databases) {
#
#     # loading the de Lacy Costello database (2014)
#     # [CAS have been converted to CHEBIs by Pierrick Roger-Mele with the 'biodb' package]
#     cost.df <- utils::read.table(file = system.file("extdata/reference_tables/vocDB_tables/de_Lacy_Costello_2014.tsv",
#                                              package = "ptairMS"),
#                           comment.char = "",
#                           header = TRUE,
#                           quote = "",
#                           sep = "\t",
#                           stringsAsFactors = FALSE)
#
#     message("Initial number annotations: ", nrow(cost.df))
#
#
#     # adding extra information from CTS conversion of CAS ids
#     # https://cts.fiehnlab.ucdavis.edu/batch
#     cts.df <- utils::read.table(file = system.file("extdata/reference_tables/vocDB_tables/de_Lacy_Costello_2014_cts-20191029142351.csv",
#                                             package = "ptairMS"),
#                          comment.char = "",
#                          header = TRUE,
#                          quote = "\"",
#                          sep = ",",
#                          stringsAsFactors = FALSE)
#     cts.df[cts.df[, "ChEBI"] == "No result", "ChEBI"] <- ""
#
#     ## additional ChEBI to those already found by biodb
#     add_chebi.vi <- setdiff(which(cts.df[, "ChEBI"] != ""),
#                             which(cost.df[, "chebi.accession"] != ""))
#     add_chebi.vc <- cts.df[add_chebi.vi, "ChEBI"]
#     add_chebi.vc <- gsub("CHEBI:", "", gsub("\nCHEBI:", "|", add_chebi.vc))
#     add_chebi.vl <- !grepl("|", add_chebi.vc, fixed = TRUE)
#     add_chebi.vi <- add_chebi.vi[add_chebi.vl]
#     add_chebi.vc <- add_chebi.vc[add_chebi.vl]
#
#     ## getting information for those new ChEBIs with biodb
#     mybiodb <- biodb::Biodb()
#
#     chebi <- mybiodb$getFactory()$createConn('chebi')
#
#     select.vl <- sapply(add_chebi.vc, function(chebi.c) {
#       entry.ls <- try(chebi$getEntry(chebi.c))
#       !inherits(entry.ls, "try-error")
#     }) # 4 ChEBI needs to be discarded because of 'UTF-8' errors
#
#     add_chebi.vi <- add_chebi.vi[select.vl]
#     add_chebi.vc <- add_chebi.vc[select.vl]
#
#     entries.ls <- chebi$getEntry(add_chebi.vc)
#     chebi.df <- mybiodb$entriesToDataframe(entries.ls,
#                                            fields = c("chebi.id", "formula", "monoisotopic.mass", "name", "inchi"))
#
#     mybiodb$terminate()
#
#     add_chebi.vi <- add_chebi.vi[-14] # weird formula (XXX)nH2O
#     chebi.df <- chebi.df[-14, ]
#     cost.df[add_chebi.vi, "chebi.accession"] <- chebi.df[, "chebi.id"]
#     cost.df[add_chebi.vi, "chebi.formula"] <- chebi.df[, "formula"]
#     cost.df[add_chebi.vi, "chebi.monoisotopic.mass"] <- chebi.df[, "monoisotopic.mass"]
#     cost.df[add_chebi.vi, "chebi.name"] <- chebi.df[, "name"]
#     cost.df[add_chebi.vi, "chebi.inchi"] <- chebi.df[, "inchi"]
#
#     ## getting mz_Hplus and formula_Hplus for those new ChEBI
#     mz_formula.vn <- formula2mass(chebi.df[, "formula"],
#                                           protonizeL = TRUE)
#     cost.df[add_chebi.vi, "mz_Hplus"] <- mz_formula.vn
#     cost.df[add_chebi.vi, "formula_Hplus"] <- names(mz_formula.vn)
#
#
#     # removing all rows (molecules) without CAS converted ids (ie unkown mz and formula)
#     cost.df <- cost.df[!is.na(cost.df[, "mz_Hplus"]), ]
#
#     message("Number of annotations with converted CAS ids: ", nrow(cost.df))
#
#     # indexing isomers
#
#     isomer.i <- 1
#
#     isomer.vi <- numeric(nrow(cost.df))
#
#     for (k in 1:nrow(cost.df)) {
#
#       if (isomer.vi[k] == 0) {
#
#         isomer.vi[which(cost.df[, "formula_Hplus"] == cost.df[k, "formula_Hplus"])] <- isomer.i
#         isomer.i <- isomer.i + 1
#
#       }
#
#     }
#
#     cost.df[, "isomer.group"] <- isomer.vi
#
#     cas.vi <- suppressWarnings(as.numeric(sapply(cost.df[, "cas.number"],
#                                                  function(cas.c) {
#                                                    gsub("-", "", cas.c)
#                                                  })))
#
#     cost.df <- cost.df[order(cost.df[, "mz_Hplus"],
#                              cas.vi,
#                              cost.df[, "compound.name"]), ]
#
#     # concatenating isomers
#
#     cost.remaining.vi <- 1:nrow(cost.df)
#     cost.discarding.vi <- integer()
#
#     while (length(cost.remaining.vi) > 0) {
#
#       cost.selecting.vi <- which(cost.df[, "isomer.group"] == cost.df[cost.remaining.vi[1], "isomer.group"])
#
#       if (length(cost.selecting.vi) > 1) {
#
#         for (colname.c in setdiff(colnames(cost.df), c("mz_Hplus",
#                                                        "formula_Hplus",
#                                                        "reference",
#                                                        "chebi.formula",
#                                                        "chebi.monoisotopic.mass"))) {
#
#           cost.df[cost.selecting.vi[1], colname.c] <- paste(cost.df[cost.selecting.vi, colname.c],
#                                                             collapse = "|")
#
#         }
#
#         cost.discarding.vi <- c(cost.discarding.vi, cost.selecting.vi[-1])
#
#       }
#
#       cost.remaining.vi <- setdiff(cost.remaining.vi, cost.selecting.vi)
#
#     }
#
#     cost.df <- cost.df[-cost.discarding.vi, ]
#
#     message("Final number of annotations after concatenating isomers: ", nrow(cost.df))
#
#     colnames(cost.df) <- gsub("cas.number", "cas", colnames(cost.df))
#     colnames(cost.df) <- gsub("compound.name", "cas.name", colnames(cost.df))
#     colnames(cost.df) <- gsub("chebi.accession", "chebi", colnames(cost.df))
#     colnames(cost.df) <- gsub("chebi.formula", "formula", colnames(cost.df))
#     colnames(cost.df) <- gsub("chebi.monoisotopic.mass", "monoisotopic.mass", colnames(cost.df))
#     colnames(cost.df) <- gsub("chebi.inchi", "inchi", colnames(cost.df))
#     cost.df[, "isomer.group"] <- NULL
#
#     # building a 'summary' annotation
#
#     summary.vc <- apply(cost.df, 1,
#                         function(voc.vc) {
#
#                           mz.n <- floor(as.numeric(voc.vc["mz_Hplus"]) * 1000)/1000
#                           formula.c <- voc.vc["formula_Hplus"]
#                           name.c <- voc.vc["cas.name"]
#                           if (grepl("|", name.c, fixed = TRUE)) {
#                             name.c <- unlist(strsplit(name.c, "|", fixed = TRUE))[1]
#                             suffix.c <- "..."
#                           } else
#                             suffix.c <- ""
#
#                           if (grepl(", ", name.c, fixed = TRUE)) {
#                             name.c <- unlist(strsplit(name.c, ", ", fixed = TRUE))[1]
#                             suffix.c <- "..."
#                           }
#
#                           paste0(mz.n, ", ",
#                                  formula.c, ", ",
#                                  name.c,
#                                  suffix.c)
#
#                         })
#
#     cost.df[, "summary"] <- summary.vc
#
#     cost.df[, "chebi"] <- paste0("CHEBI:", cost.df[, "chebi"])
#
#     vocDB.df <- cost.df[, c("mz_Hplus",
#                             "formula_Hplus",
#                             "cas.name",
#                             "summary",
#                             "cas",
#                             "blood",
#                             "breath",
#                             "faeces",
#                             "milk",
#                             "saliva",
#                             "skin",
#                             "urine",
#                             "reference",
#                             "chebi",
#                             "monoisotopic.mass",
#                             "formula",
#                             "chebi.name",
#                             "inchi",
#                             "kegg.compound.id")]
#
#
#
#   }
#
#   utils::write.table(vocDB.df,
#               file = "vocDB.tsv",
#               quote = FALSE,
#               row.names = FALSE,
#               sep = "\t")
#
#   # 'vocDB.tsv' file to be moved to the 'ptairMS/extdata/reference_tables/' directory
#
#   invisible(vocDB.df)
#
# }
