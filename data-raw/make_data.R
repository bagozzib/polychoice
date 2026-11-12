## make_data.R -- build the bundled example datasets from the replication
## archive of Bagozzi and Marchetti (2017, PSRM). Run from the package root.
## Source: the Harvard Dataverse deposit doi:10.7910/DVN/3JR2YL (version 1.1,
## released 2015-07-16 under CC0 1.0), whose single file "PSRM Replication
## Files.zip" (md5 43c6af048c93226486ef4a4fb32f96dd) unzips to the folder
## "PSRM Replication Files". The raw .dta files are NOT part of the package;
## RAW points at that folder (environment variable POLYCHOICE_RAW), which holds
##   Replication I/Least_Informed.dta       Arceneaux and Kolodny (2009, AJPS) as used in Replication I
##   Replication II/replicationdata.dta     Campbell and Monson (2008, POQ) as used in Replication II
## and both are checked against their md5 sums in the deposit (verified
## 2026-10-14) before anything is built.
## Only the analysis variables of the published specifications are kept, and
## only the complete cases the published tables use (N = 1,998 and N = 1,341).
RAW <- Sys.getenv("POLYCHOICE_RAW")
if (!nzchar(RAW)) stop("Set the environment variable POLYCHOICE_RAW to the folder \"PSRM Replication Files\" of the deposit")
md5 <- c("Replication I/Least_Informed.dta" = "0a0964fc9183e1506e893695000675a2",
         "Replication II/replicationdata.dta" = "29ee7e0d84871d898fb19e1450b5c9a0")
got <- unname(tools::md5sum(file.path(RAW, names(md5))))
if (!identical(got, unname(md5))) stop("The input files do not match the Dataverse deposit (doi:10.7910/DVN/3JR2YL): ", paste(names(md5)[got != md5], collapse = ", "))
suppressMessages(library(haven))
asdf <- function(x) { x <- as.data.frame(x); x[] <- lapply(x, function(v) { attributes(v) <- NULL; v }); x }

## 1. pa2006: 2006 Pennsylvania state House vote choice (Arceneaux and Kolodny 2009; Replication I)
a <- asdf(read_dta(file.path(RAW, "Replication I", "Least_Informed.dta")))
a <- a[stats::complete.cases(a), ]
stopifnot(nrow(a) == 1998L)
## votechoice_dk: 0 = Republican/other candidate, 1 = abstain (did not vote), 2 = Democratic candidate, 3 = nonresponse
pa2006 <- data.frame(
  vote = factor(a$votechoice_dk, levels = c(1, 2, 0, 3), labels = c("abstain", "democrat", "republican_other", "nonresponse")),
  canvass = as.integer(a$canvass_all), phone = as.integer(a$phone_only),
  democrat = as.integer(a$sid_dem), republican = as.integer(a$sid_rep),
  household = as.integer(a$hhcount), age = as.numeric(a$age), female = as.integer(a$female),
  voted2004 = as.integer(a$vg04), district156 = as.integer(a$d156),
  stringsAsFactors = FALSE)
## the published interactions are exact products of the bundled columns
stopifnot(all(a$cXdem == a$canvass_all * a$sid_dem), all(a$cXrep == a$canvass_all * a$sid_rep),
          all(a$poXdem == a$phone_only * a$sid_dem), all(a$poXrep == a$phone_only * a$sid_rep))
rownames(pa2006) <- NULL
print(table(pa2006$vote))

## 2. pres2004: 2004 presidential vote choice (Campbell and Monson 2008; Replication II)
b <- asdf(read_dta(file.path(RAW, "Replication II", "replicationdata.dta")))
b <- b[stats::complete.cases(b), ]
stopifnot(nrow(b) == 1341L)
pres2004 <- data.frame(
  vote = factor(b$vote_turn, levels = 0:2, labels = c("abstain", "kerry", "bush")),
  educ = as.integer(b$educ), female = as.integer(b$female), black = as.integer(b$black), hispanic = as.integer(b$hispanic),
  pid7 = as.integer(b$party7), bush_values = as.integer(b$w3gbvalu2), battleground = as.integer(b$presbatt),
  south = as.integer(b$south), gmb = as.integer(b$gayban2), agegroup = as.integer(b$agegroup2),
  evangelical = as.integer(b$bornagn_w), catholic = as.integer(b$catholic), other_religion = as.integer(b$other_rel),
  secular = as.integer(b$secular), ideology = as.integer(b$ideo), mobilization = as.integer(b$w3mobidx),
  stringsAsFactors = FALSE)
stopifnot(all(b$brnag_wXgmb2 == b$bornagn_w * b$gayban2), all(b$cathXgmb2 == b$catholic * b$gayban2), all(b$secXgmb2 == b$secular * b$gayban2))
rownames(pres2004) <- NULL
print(table(pres2004$vote))

save(pa2006, file = "data/pa2006.rda", compress = "xz")
save(pres2004, file = "data/pres2004.rda", compress = "xz")
cat("written data/pa2006.rda and data/pres2004.rda\n")
