# Internal helpers shared by all compartmental models -----------------------

# Build the mutable working table: phenotypes + dynamics, tagged with the
# internal id (iid) so the original individual order can be restored later.
# copy() guards against modifying the dynamics slot in place via set()/`:=`.
.epi_working_table <- function(epop) {
  X <- cbind(as.data.table(epop@pheno), copy(epop@dynamics))
  X[, iid := epop@iid]
  X[]
}

# Fix the disease trajectory of the initial donors using a model-specific
# path generator (e.g. .generate_sir_path, .generate_sidr_path). Modifies X
# by reference.
.epi_init_donors <- function(X, path_fn, simParam) {
  purrr::walk(X[, .I[donor == 1L]], \(i) {
    data.table::set(X, i, c("status", simParam$timings),
                    path_fn(0.0, X, i, simParam))
  })
  invisible(X)
}

# Recombine step: restore the population's own row order, drop the borrowed
# phenotype columns, iid, and per-iteration bookkeeping, add the parasites
# flag, report, and write back into the population.
.finalise_epi_dynamics <- function(epop, Y) {
  orig_cols <- names(epop@dynamics)

  # Restore epop's individual order. iid is NOT guaranteed to be sorted (it can
  # be arbitrary after crossing/subsetting), so we reorder to match epop@iid
  # rather than sorting ascending.
  Y <- Y[match(epop@iid, iid)]

  Y[, parasites := !is.na(Tinf)]

  # Keep only the original dynamics columns (minus group_inf) plus parasites.
  # This automatically discards pheno cols, iid, and event_rate, and is robust
  # to however as.data.table() named the phenotype columns.
  keep <- union(setdiff(orig_cols, "group_inf"), "parasites")
  Y <- Y[, ..keep]

  final_t <- max(Y$Tdeath, na.rm = TRUE) |> signif(5)
  message(sprintf("- Final t = %f, values are:", final_t),
          paste(capture.output(table(Y$status)), collapse = ", "), "\n")

  epop@dynamics <- Y[]
  epop
}

.safe_sample <- function(x, ...) {
  if (length(x) > 1L) {
    sample(x, ...)
  } else if (length(x) <= 1L) {
    x
  }
}