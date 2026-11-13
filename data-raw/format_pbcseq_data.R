# Select baseline covariates and final status
pbcseq_baseline <- survival::pbcseq |>
  dplyr::select(
    id,
    trt,
    age,
    sex,
    futime,
    status
  ) |>
  dplyr::distinct()

# Select longitudinal visit-level data
pbcseq_visits <- survival::pbcseq |>
  dplyr::select(
    id,
    day,
    ascites,
    hepato,
    spiders,
    edema,
    bili,
    chol,
    albumin,
    alk.phos,
    ast,
    platelet,
    protime,
    stage
  )

# Save data
usethis::use_data(
  pbcseq_baseline,
  pbcseq_visits,
  overwrite = TRUE
)
