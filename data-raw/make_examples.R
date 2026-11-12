## make_examples.R -- build the example datasets drawn from published
## replication archives. Run from the package root. The raw files are
## NOT part of the package; RAW points at the folder that holds one subfolder
## per study (environment variable POLYCHOICE_REPLICATION), and every input is
## checked against the md5 sum the Harvard Dataverse lists for it (verified
## 2026-10-20) before anything is built. All seven deposits are released under
## CC0 1.0; bes1987 and bes2010 are extracts of British Election Study surveys
## that the UK Data Service distributes under CC BY 4.0 (SN 2568, SN 7529).
##
##   afghan_leaders  Bhatia and Monroe (2024, BJPS)        doi:10.7910/DVN/OIV11K  Afghan_Data_All_Clean.csv
##   fluvax          Liu et al. (2024, JJPS)               doi:10.7910/DVN/BYAFJR  Data.dta
##   cabinets        Baeck, Debus and Imre (2026, PP)      doi:10.7910/DVN/CJIVXC  baeck_et_al_2025_data.RDATA
##   sanctions       Liou (2024, FPA)                      doi:10.7910/DVN/9QQYPB  DiversionarySanctions_Main.dta
##   bes1987         Glasgow (2001, PA)                    doi:10.7910/DVN/PCPVWB  glasgow.web.zip (glasgow.data.asc)
##   bes2010         Arzheimer and Evans (2012, Pol Geog)  hdl:1902.1/17940        replication.dta
##   defection       Cook, Niehaus and Zuhlke (2018, R&P)  doi:10.7910/DVN/GOYSFK  replication data.dta (CC0 1.0;
##                   added 2026-11-05, md5 as listed by the Dataverse)
RAW <- Sys.getenv("POLYCHOICE_REPLICATION")
if (!nzchar(RAW)) stop("Set the environment variable POLYCHOICE_REPLICATION to the folder that holds one subfolder per study")
md5 <- c("Bhatia Monroe/Afghan_Data_All_Clean.csv"   = "8c9d3d11d81efdad4ded7770a1dd6e77",
         "Liu et al/Data.dta"                        = "40cad902c1e17f2de59b6dac0aacc5da",
         "Back 2025/baeck_et_al_2025_data.RDATA"     = "e046b31d8e83dea2e5f01ed498f0cf32",
         "Liou 2024/DiversionarySanctions_Main.dta"  = "5230d7cda77acc7165ff4e459bbe0508",
         "Glasgow 2001/glasgow.web.zip"              = "ab0c27430d5b281eed1859ebb88b64bc",
         "Arzheimer 2012/replication.dta"            = "d82e986aef3f1389c6ee3babb0a15566",
         "Cook 2018/replication data.dta"            = "c08c38bcb3369ae38940b0a0e1723244")
got <- unname(tools::md5sum(file.path(RAW, names(md5))))
if (!identical(got, unname(md5))) stop("Input files differ from the Dataverse deposits: ", paste(names(md5)[got != md5], collapse = ", "))
suppressMessages(library(haven))
asdf <- function(x) { x <- as.data.frame(x); x[] <- lapply(x, function(v) { attributes(v) <- NULL; v }); x }
is_ <- function(x, vals) !is.na(x) & x %in% vals             # Stata's `if x == c`: a missing x never matches

## 1. afghan_leaders: rankings of six leader attributes in a survey experiment
##    in northern Afghanistan (Bhatia and Monroe 2024). The codebook (Bhatia
##    Codebook Survey BJPolS_Clean.xlsx) gives the codes used below. Kept: the
##    four arms of the vignette experiment (Treatment 1-4) and complete
##    rankings 1-6 without ties; the respondent's village, district and place
##    of birth are not carried over.
b <- utils::read.csv(file.path(RAW, "Bhatia Monroe", "Afghan_Data_All_Clean.csv"), stringsAsFactors = FALSE)
att <- c(religious = "Att_Strong_Rel", educated = "Att_Well_Ed", mujahideen = "Att_Muj",
         government = "Att_Gov_Exp", military = "Att_Afg_Mil", honest = "Att_Hon")
b <- b[b$Treatment %in% 1:4 & apply(b[, att], 1, function(v) !anyNA(v) && all(sort(v) == 1:6)), ]
b <- b[stats::complete.cases(b[, c("Res_Age", "Res_Gender", "Res_Education", "Res_Ethnicity", "Province")]) &
       b$Res_Gender %in% 1:2 & b$Res_Education %in% c(1:6, 8) & b$Res_Ethnicity %in% 1:7, ]
resp <- data.frame(
  id = as.integer(b$ID),
  treatment = factor(b$Treatment, levels = 1:4, labels = c("control", "bribery", "nepotism", "insecurity")),
  female = as.integer(b$Res_Gender == 2), age = as.numeric(b$Res_Age),
  education = factor(b$Res_Education, levels = c(1:6, 8),
                     labels = c("none", "primary", "secondary", "vocational", "some_university", "university", "madrassa")),
  ethnicity = factor(b$Res_Ethnicity, levels = 1:7, labels = c("pashtun", "tajik", "uzbek", "turkmen", "hazara", "baloch", "other")),
  province = factor(b$Province, levels = c("Balkh", "Kundoz", "Sar-e-pull"), labels = c("Balkh", "Kunduz", "Sar-e Pol")),
  stringsAsFactors = FALSE)
stopifnot(!anyNA(resp), !anyDuplicated(resp$id))
afghan_leaders <- do.call(rbind, lapply(names(att), function(a)
  data.frame(id = resp$id, attribute = factor(a, levels = names(att)), rank = as.integer(b[[att[[a]]]]), resp[, -1])))
afghan_leaders <- afghan_leaders[order(afghan_leaders$id, afghan_leaders$attribute), ]
rownames(afghan_leaders) <- NULL
cat("afghan_leaders:", length(unique(afghan_leaders$id)), "respondents\n"); print(table(resp$treatment))

## 2. fluvax: intention to take a flu vaccine developed by the US alone or with
##    China or Japan, a randomized survey experiment (Liu et al. 2024).
##    Variables rebuilt from the raw items as DoFile.do builds them (Table 3,
##    Model 1); the codebook (Codebook.pdf) gives the answer codes.
l <- asdf(read_dta(file.path(RAW, "Liu et al", "Data.dta")))
q14 <- as.character(l$Q14)
race <- rep(NA_character_, nrow(l)); race[is_(q14, "1")] <- "white"; race[is_(q14, "2")] <- "black"
other <- c("1,2", "1,2,3", "1,2,3,4", "1,2,3,4,5,6,7", "1,2,3,6", "1,2,6", "1,3", "1,3,6", "1,3,7", "1,4", "1,5", "1,6", "1,7",
           "2,3", "2,4,5", "2,6", "2,7", "3", "3,6", "4", "4,5", "4,6", "5", "6", "7")
race[is_(q14, other)] <- "other"
p7 <- rep(NA_real_, nrow(l))                                 # the do-file's sequential replaces
p7[is_(l$Q10, 1)] <- 1; p7[is_(l$Q10, 2)] <- 2; p7[is_(l$Q12, 3)] <- 3; p7[is_(l$Q12, 2)] <- 4
p7[is_(l$Q12, 1)] <- 5; p7[is_(l$Q11, 2)] <- 6; p7[is_(l$Q11, 1)] <- 7
party <- ifelse(is_(p7, 1:3), "democrat", ifelse(is_(p7, 4), "independent", ifelse(is_(p7, 5:7), "republican", "other")))
pk <- cbind(ifelse(is_(l$Q112, 2), 1, ifelse(is_(l$Q112, c(1, 3, 4, 5)), 0, NA)), ifelse(is_(l$Q113, 4), 1, ifelse(is_(l$Q113, c(1, 2, 3, 5)), 0, NA)),
            ifelse(is_(l$Q114, 2), 1, ifelse(is_(l$Q114, c(1, 3, 4)), 0, NA)), ifelse(is_(l$Q115, 2), 1, ifelse(is_(l$Q115, c(1, 3)), 0, NA)))
gv <- apply(as.matrix(l[, c("Q75", "Q76", "Q77")]), 1, function(r) if (all(is.na(r))) NA_real_ else max(r, na.rm = TRUE))
fluvax <- data.frame(
  vaccinate = factor(ifelse(is_(gv, 1), "no", ifelse(is_(gv, 3), "dont_know", ifelse(is_(gv, 2), "yes", NA))), levels = c("no", "dont_know", "yes")),
  condition = factor(as.character(l$Group), levels = c("Control", "China", "Japan"), labels = c("us_alone", "us_china", "us_japan")),
  party = factor(party, levels = c("democrat", "independent", "republican", "other")),
  age = as.numeric(l$Q8), female = as.integer(l$Q120 - 1), educ = as.integer(l$Q6), income = as.integer(l$Q7),
  race = factor(race, levels = c("white", "black", "other")), knowledge = rowMeans(pk),
  stringsAsFactors = FALSE)
fluvax <- fluvax[stats::complete.cases(fluvax), ]
stopifnot(nrow(fluvax) == 1334L, all(fluvax$age < 999))
rownames(fluvax) <- NULL
print(table(fluvax$vaccinate, fluvax$condition))

## 3. cabinets: government formation in the German states, 1990-2023 (Baeck,
##    Debus and Imre 2025): one row per formation opportunity and potential
##    cabinet. The residuals pea_res and rej_res of the two-stage residual
##    inclusion are rebuilt as in baeck_et_al_2025_script.R (logit first stages
##    on the rows with s != 1, before any row is dropped). Dropped: Thuringia
##    2019 (two chosen cabinets) and Hesse 2023 (none).
env <- new.env(); load(file.path(RAW, "Back 2025", "baeck_et_al_2025_data.RDATA"), envir = env)
k <- as.data.frame(env$data)
rhs <- "incXno_change + incXpm_change + minder + minwc + bprop + strongest + ccoal + pheteco_s + phetsoc_s"
first <- k[!is.na(k$s) & k$s != 1, , drop = FALSE]
k$pea_res <- k$pea - stats::predict(stats::glm(stats::as.formula(paste("pea ~", rhs)), data = first, family = stats::binomial()), newdata = k, type = "response")
k$rej_res <- k$rej - stats::predict(stats::glm(stats::as.formula(paste("rej ~", rhs)), data = first, family = stats::binomial()), newdata = k, type = "response")
parties <- c("cdu", "spd", "fdp", "gru", "pds", "rep", "dvu", "afb", "sch", "sta", "npd", "ssw", "pir", "afd")
k$cabinet <- apply(as.matrix(k[, parties]), 1, function(r) paste(parties[which(r == 1)], collapse = "+"))
dup <- stats::ave(seq_len(nrow(k)), k$election, k$cabinet, FUN = seq_along)
k$cabinet <- ifelse(dup > 1, paste0(k$cabinet, "#", dup), k$cabinet)
w <- c("incXno_change", "incXpm_change", "minder", "minwc", "bprop", "strongest", "ccoal", "pheteco_s", "phetsoc_s", "pea", "rej", "pea_res", "rej_res")
k <- k[stats::complete.cases(k[, c("election", "cabinet", "choice", w)]), ]
npos <- tapply(k$choice, k$election, sum); keep <- names(npos)[npos == 1]
k <- k[k$election %in% keep, ]
cabinets <- data.frame(election = as.character(k$election), state = as.character(k$Bundesland), year = as.integer(k$year),
                       cabinet = k$cabinet, chosen = as.integer(k$choice), k[, w], stringsAsFactors = FALSE)
names(cabinets)[names(cabinets) == "incXno_change"] <- "incumbent"; names(cabinets)[names(cabinets) == "incXpm_change"] <- "incumbent_newpm"
names(cabinets)[names(cabinets) == "minder"] <- "minority"; names(cabinets)[names(cabinets) == "minwc"] <- "min_winning"
names(cabinets)[names(cabinets) == "bprop"] <- "bargaining"; names(cabinets)[names(cabinets) == "ccoal"] <- "cross_cutting"
names(cabinets)[names(cabinets) == "pheteco_s"] <- "het_economic"; names(cabinets)[names(cabinets) == "phetsoc_s"] <- "het_social"
## ASCII strings (CRAN): the archive's "Baden-Württemberg" becomes
## "Baden-Wuerttemberg" in the state and the election identifiers
de_ascii <- function(s) { for (p in list(c("ä", "ae"), c("ö", "oe"), c("ü", "ue"), c("Ä", "Ae"), c("Ö", "Oe"), c("Ü", "Ue"), c("ß", "ss")))
                            s <- gsub(p[1], p[2], s, fixed = TRUE)
                          s }
cabinets$election <- de_ascii(cabinets$election); cabinets$state <- de_ascii(cabinets$state)
stopifnot(!any(grepl("[^ -~]", c(cabinets$election, cabinets$state, cabinets$cabinet))))
cabinets <- cabinets[order(cabinets$election, cabinets$cabinet), ]
cabinets$cabinet <- factor(cabinets$cabinet)
rownames(cabinets) <- NULL
stopifnot(length(unique(cabinets$election)) == 125L, nrow(cabinets) == 3352L)
cat("cabinets:", nrow(cabinets), "rows,", length(unique(cabinets$election)), "formation opportunities\n")

## 4. sanctions: the choice between no coercion, economic sanctions and
##    military force in country-years 1975-2005 (Liou 2024, Table 1): the
##    complete cases of the published specification (N = 2,985, 137 countries).
s <- asdf(read_dta(file.path(RAW, "Liou 2024", "DiversionarySanctions_Main.dta")))
s <- s[!is.na(s$year) & s$year > 1974 & s$year < 2006, ]
v <- c("diversion", "lbinonvio_campaign", "lbivio_campaign", "loppseatshare", "lgdpgr", "llninflation", "Flegelec", "Fexelec", "llngdppc",
       "ltradeopen", "llnmilex", "lvdem_polyarchy", "lagtheta_mean", "coldwar", "yrssince_divert", "yrssince_divertsq", "yrssince_divertcubic", "ccode", "year")
s <- s[stats::complete.cases(s[, v]), v]
sanctions <- data.frame(
  outcome = factor(s$diversion, levels = 0:2, labels = c("none", "sanctions", "force")),
  ccode = as.integer(s$ccode), year = as.integer(s$year),
  nonviolent = as.integer(s$lbinonvio_campaign), violent = as.integer(s$lbivio_campaign), opposition = s$loppseatshare,
  growth = s$lgdpgr, inflation = s$llninflation, leg_election = as.integer(s$Flegelec), exec_election = as.integer(s$Fexelec),
  gdppc = s$llngdppc, trade = s$ltradeopen, milex = s$llnmilex, democracy = s$lvdem_polyarchy, human_rights = s$lagtheta_mean,
  coldwar = as.integer(s$coldwar), peaceyears = s$yrssince_divert)
stopifnot(nrow(sanctions) == 2985L, all(abs(s$yrssince_divertsq - s$yrssince_divert^2) < 1e-6), all(abs(s$yrssince_divertcubic - s$yrssince_divert^3) < 1e-3))
sanctions <- sanctions[order(sanctions$ccode, sanctions$year), ]
rownames(sanctions) <- NULL
print(table(sanctions$outcome))

## 5. bes1987: vote choice in the 1987 British general election with issue
##    distances to the three parties (Glasgow 2001; BES 1987 cross-section):
##    the 39 columns of glasgow.data.asc in the order of the GAUSS load
##    statement (glasgow.mxlcode.prg lines 59-104), in the long layout.
td <- tempfile(); utils::unzip(file.path(RAW, "Glasgow 2001", "glasgow.web.zip"), exdir = td)
asc <- list.files(td, "glasgow.data.asc$", recursive = TRUE, full.names = TRUE)
stopifnot(length(asc) == 1L, unname(tools::md5sum(asc)) == "eaab3920305e75e2d2e2da23005922af")
g <- as.data.frame(matrix(scan(asc, what = numeric(), quiet = TRUE), ncol = 39L, byrow = TRUE))
issues <- c("defense", "phillips", "taxation", "nationalization", "redistribution", "crime", "welfare")
names(g) <- c("id", "votechc", "respinf", "respunem", "resptax", "south", "mid", "north", "wales", "scot", "union2", "public", "blue", "women",
              "age", "homeown", "faminc", "educage", as.vector(outer(issues, 1:3, paste0)))
stopifnot(nrow(g) == 2131L, all(g$votechc %in% 1:3))
bes1987 <- do.call(rbind, lapply(1:3, function(j) {
  out <- data.frame(id = as.integer(g$id), party = factor(c("conservative", "labour", "alliance")[j], levels = c("conservative", "labour", "alliance")),
                    chosen = as.integer(g$votechc == j))
  for (i in issues) out[[i]] <- g[[paste0(i, j)]]
  cbind(out, data.frame(south = g$south, midlands = g$mid, north = g$north, wales = g$wales, scotland = g$scot, union = g$union2,
                        public = g$public, working_class = g$blue, female = g$women, age = g$age, homeowner = g$homeown,
                        income = g$faminc, education = g$educage, inflation = g$respinf, taxes = g$resptax, unemployment = g$respunem))
}))
bes1987 <- bes1987[order(bes1987$id, bes1987$party), ]
rownames(bes1987) <- NULL
print(table(bes1987$party[bes1987$chosen == 1]))

## 6. bes2010: vote choice in English constituencies in 2010 with
##    candidate-voter driving distances (Arzheimer and Evans 2012; BES 2010
##    face-to-face survey): the sample of the do-file's fourth model
##    ("deprivation"; replication.do line 52), rows dropped as Stata's clogit
##    drops them, and the post-election ratings of the parties (bq17_1-3).
a <- asdf(read_dta(file.path(RAW, "Arzheimer 2012", "replication.dta")))
a <- a[is.na(a$touse) | a$touse != 0, ]
a$party <- ifelse(a$labour == 1, "labour", ifelse(a$conservatives == 1, "conservative", ifelse(a$liberal_democrats == 1, "libdem", NA)))
a$aconsnam[!is.na(a$aconsnam) & trimws(a$aconsnam) == ""] <- NA
a$rating <- ifelse(a$party == "labour", a$bq17_1, ifelse(a$party == "conservative", a$bq17_2, a$bq17_3))
a <- a[stats::complete.cases(a[, c("id", "party", "votechoice", "aq16_", "inccand", "drdistp", "absdepdist", "aconsnam")]), ]
npos <- tapply(a$votechoice, a$id, sum); nrow_ <- tapply(a$votechoice, a$id, length)
a <- a[as.character(a$id) %in% names(npos)[npos == 1 & nrow_ >= 2], ]
bes2010 <- data.frame(id = as.integer(a$id), party = factor(a$party, levels = c("labour", "conservative", "libdem")),
                      chosen = as.integer(a$votechoice), feeling = a$aq16_, incumbent = as.integer(a$inccand), distance = a$drdistp,
                      deprivation = a$absdepdist, rating = a$rating, constituency = a$aconsnam, stringsAsFactors = FALSE)
bes2010 <- bes2010[order(bes2010$id, bes2010$party), ]
rownames(bes2010) <- NULL
stopifnot(length(unique(bes2010$id)) == 674L, nrow(bes2010) == 1788L)
cat("bes2010:", nrow(bes2010), "rows,", length(unique(bes2010$id)), "respondents\n")

## 7. defection: how ex-combatants of Colombian left-wing guerrilla groups
##    left them (captured, individually demobilized, or switched to the
##    paramilitaries), from a survey of the Fundacion Ideas para la Paz
##    analysed by Oppenheim et al. (2015, JCR) and redistributed by Cook,
##    Niehaus and Zuhlke (2018). Value labels of the deposit's .dta:
##    DV_cap_inddemob_switcher 0 captured, 1 individually demobilized,
##    2 switchers; p402g 0 no, 1 si; p101 0 femenino, 1 masculino. Sample as
##    in Cook et al.'s replication script (a department of residence recorded,
##    p504dantes != "") and the complete cases of their models (N = 582).
x <- asdf(read_dta(file.path(RAW, "Cook 2018", "replication data.dta")))
x <- x[x$p504dantes != "", ]
defection <- data.frame(id = as.integer(x$var00002),
                        outcome = factor(x$DV_cap_inddemob_switcher, levels = 0:2, labels = c("captured", "demobilized", "switched")),
                        economic = as.integer(x$join_ec_need), ideological = as.integer(x$join_ideo),
                        indoctrination = as.integer(x$p409dummy), abuse = as.integer(x$p402g),
                        age_joined = as.integer(x$p104a), birth_year = as.integer(x$yrbirth), school_age = as.integer(x$p203),
                        male = as.integer(x$p101), besieged = as.integer(x$p410_recoded))
defection <- defection[stats::complete.cases(defection[, setdiff(names(defection), "besieged")]), ]
rownames(defection) <- NULL
stopifnot(nrow(defection) == 582L, !anyDuplicated(defection$id))
defection$id <- seq_len(nrow(defection))       # a case number: the survey's own identifiers are not distributed
cat("defection:", nrow(defection), "cases\n"); print(table(defection$outcome))

for (nm in c("afghan_leaders", "fluvax", "cabinets", "sanctions", "bes1987", "bes2010", "defection"))
  save(list = nm, file = file.path("data", paste0(nm, ".rda")), compress = "xz")
cat("written", paste0("data/", c("afghan_leaders", "fluvax", "cabinets", "sanctions", "bes1987", "bes2010", "defection"), ".rda"), "\n")
