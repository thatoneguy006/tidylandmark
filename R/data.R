#' Mayo Clinic Primary Biliary Cirrhosis Baseline and Status Data for Randomized
#' Patients
#'
#' Adapted from [survival::pbcseq], this dataset contains baseline and outcome
#' status data among randomized patients in the Mayo Clinic primary biliary
#' cirrhosis trial. Accompanies [pbcseq_visits] which contains follow-up visit
#' data preceding `status`.
#'
#' @format
#' A data frame with 312 observations of six variables:
#' \describe{
#'   \item{id}{case number}
#'   \item{age}{in years}
#'   \item{sex}{m/f}
#'   \item{trt}{1/2 for D-penicillmain, placebo}
#'   \item{futime}{Time at `status`}
#'   \item{status}{status at endpoint, 0/1/2 for censored, transplant, dead}
#' }
#'
#' @usage pbcseq_baseline
#'
#' @source Adapted from [survival::pbcseq].
#'
#' T Therneau and P Grambsch, "Modeling Survival Data: Extending the Cox Model",
#' Springer-Verlag, New York, 2000. ISBN: 0-387-98784-3.
#'
"pbcseq_baseline"

#' Mayo Clinic Primary Biliary Cirrhosis Follow-up Visit Data for Randomized
#' Patients
#'
#' Adapted from [survival::pbcseq], this dataset contains follow-up visit
#' data preceding final study status among randomized patients in the Mayo
#' Clinic primary biliary cirrhosis trial. Accompanies [pbcseq_baseline] which
#' contains baseline and outcome status data.
#'
#' @format
#' A data frame with 1,945 observations of 14 variables:
#' \describe{
#'   \item{id}{case number}
#'   \item{day}{	number of days between enrollment and this visit date
#'   all measurements below refer to this date}
#'   \item{ascites}{presence of ascites}
#'   \item{hepato}{presence of hepatomegaly or enlarged liver}
#'   \item{spiders}{blood vessel malformations in the skin}
#'   \item{edema}{0 no edema, 0.5 untreated or successfully treated
#'   edema despite diuretic therapy}
#'   \item{bili}{serum bilirunbin (mg/dl)}
#'   \item{chol}{serum cholesterol (mg/dl)}
#'   \item{albumin}{serum albumin (mg/dl)}
#'   \item{alk.phos}{alkaline phosphotase (U/liter)}
#'   \item{ast}{aspartate aminotransferase, once called SGOT (U/ml)}
#'   \item{platelet}{platelet count}
#'   \item{protime}{standardised blood clotting time}
#'   \item{stage}{histologic stage of disease (needs biopsy)}
#' }
#'
#' @usage pbcseq_visits
#'
#' @source Adapted from [survival::pbcseq].
#'
#' T Therneau and P Grambsch, "Modeling Survival Data: Extending the Cox Model",
#' Springer-Verlag, New York, 2000. ISBN: 0-387-98784-3.
"pbcseq_visits"
