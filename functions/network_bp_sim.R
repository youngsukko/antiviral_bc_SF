# ==============================================================================
# Integer encoding legend
# ------------------------------------------------------------------------------
# status            : 1L = S, 2L = I, 3L = R
# contact_type      : 0L = seed, 1L = household, 2L = community
# treat_reason      : 1L = self, 2L = household, 3L = community
# quarantine_reason : 1L = self, 2L = household, 3L = community
# asymptomatic      : 0L = symptomatic, 1L = asymptomatic
# ==============================================================================


# ==============================================================================
# precompute_neighbors()
# ==============================================================================
precompute_neighbors <- function(net_data) {

  n_people <- nrow(net_data$nodes)

  build_neighbor_list <- function(edges, n) {
    if (nrow(edges) == 0) {
      out   <- vector("list", n)
      out[] <- list(integer(0))
      return(out)
    }
    both <- c(edges$from, edges$to)
    nbrs <- c(edges$to,   edges$from)
    raw  <- split(nbrs, both)
    out  <- vector("list", n)
    pid  <- as.integer(names(raw))
    for (k in seq_along(pid)) out[[pid[k]]] <- as.integer(raw[[k]])
    out[sapply(out, is.null)] <- list(integer(0))
    out
  }

  hh_nbr   <- build_neighbor_list(net_data$hh_edges,  n_people)
  comm_nbr <- build_neighbor_list(net_data$comm_edges, n_people)

  lapply(seq_len(n_people), function(i)
    list(household = hh_nbr[[i]], community = comm_nbr[[i]])
  )
}


# ==============================================================================
# network_bp_sim()
#
# Arguments (new vs original):
#   seed_symptomatic : logical (default TRUE)
#     TRUE  — seed case is symptomatic  (asymptomatic = 0L)
#             triggers self-treatment and ring PEP on day antiviral_start
#     FALSE — seed case is asymptomatic (asymptomatic = 1L)
#             cannot self-treat (no symptom onset); ring is triggered from
#             time_infectious instead of time_onset
# ==============================================================================
network_bp_sim <- function(
    neighbors_by_type,
    node_df,
    R0,
    epidemic_prob_threshold    = 0.999,
    generation_time_fn,
    infection_to_onset_fn,
    infectious_before_onset_fn,
    prop_asymptomatic,
    p_inf_household,
    p_inf_community,
    antiviral_start,
    logistical_delay_fn,
    drug_eff_inf_fn,
    drug_eff_inf_hh_fn,
    drug_eff_trans_fn,
    eff_inf_e_max,
    eff_inf_hh_e_max,
    eff_trans_e_max,
    prob_treat_self,
    prob_treat_household,
    prob_treat_community,
    time_to_quarantine_fn,
    prob_quarantine_self,
    prob_quarantine_household,
    prob_quarantine_community,
    quarantine_efficacy,
    seeding_cases,
    seed_symptomatic = TRUE,    # NEW — controls seed case symptom status
    t0   = 0,
    seed = 42
) {
  set.seed(seed)

  n_people <- nrow(node_df)

  # --------------------------------------------------------------------------
  # Epidemic threshold (Whittle 1955, asymptomatic-only)
  # --------------------------------------------------------------------------
  if (R0 > 1) {
    epidemic_queue_threshold <- ceiling(
      log(1 - epidemic_prob_threshold) / log(1 / R0)
    )
  } else {
    epidemic_queue_threshold <- Inf
  }

  # --------------------------------------------------------------------------
  # Initialize tdf
  # --------------------------------------------------------------------------
  tdf <- data.frame(
    person_id             = node_df$person_id,
    household             = node_df$household,
    status                = rep(1L, n_people),
    generation            = NA_integer_,
    ancestor_id           = NA_integer_,
    contact_type          = NA_integer_,
    time_infection        = NA_real_,
    time_onset            = NA_real_,
    time_infectious       = NA_real_,
    asymptomatic          = NA_integer_,
    treated               = 0L,
    time_treated          = NA_real_,
    treat_reason          = NA_integer_,
    efficacy_at_infection = NA_real_,
    efficacy_transmission = NA_real_,
    quarantined           = 0L,
    time_quarantined      = NA_real_,
    quarantine_reason     = NA_integer_,
    offspring_generated   = FALSE,
    stringsAsFactors      = FALSE
  )

  # --------------------------------------------------------------------------
  # Seed case initialization
  # seed_symptomatic controls whether the index case has symptoms.
  #   TRUE  : asymptomatic = 0L (symptomatic)
  #           time_onset set normally; ring can be triggered from onset
  #   FALSE : asymptomatic = 1L (asymptomatic)
  #           time_onset = NA; ring triggered from time_infectious only
  # --------------------------------------------------------------------------
  seed_ids <- sample(seq_len(n_people), seeding_cases, replace = FALSE)

  for (sid in seed_ids) {
    tdf$status[sid]       <- 2L
    tdf$generation[sid]   <- 1L
    tdf$contact_type[sid] <- 0L
    tdf$time_infection[sid] <- t0

    if (seed_symptomatic) {
      # Symptomatic seed
      tdf$asymptomatic[sid]    <- 0L
      tdf$time_onset[sid]      <- t0 + infection_to_onset_fn(1)
      tdf$time_infectious[sid] <- t0 +
        (tdf$time_onset[sid] - t0) * infectious_before_onset_fn(1)
    } else {
      # Asymptomatic seed — no symptom onset, infectious from time_infectious
      tdf$asymptomatic[sid]    <- 1L
      tdf$time_onset[sid]      <- NA_real_
      tdf$time_infectious[sid] <- t0 +
        infectious_before_onset_fn(1) * infection_to_onset_fn(1)
    }
  }

  active_queue <- seed_ids
  is_epidemic  <- FALSE

  # --------------------------------------------------------------------------
  # Contact group template
  # --------------------------------------------------------------------------
  contact_groups_template <- list(
    list(ctype = 1L, sar = p_inf_household,
         p_treat      = prob_treat_household,
         p_quarantine = prob_quarantine_household,
         treat_code   = 2L, quar_code = 2L),
    list(ctype = 2L, sar = p_inf_community,
         p_treat      = prob_treat_community,
         p_quarantine = prob_quarantine_community,
         treat_code   = 3L, quar_code = 3L)
  )

  # --------------------------------------------------------------------------
  # Main loop
  # --------------------------------------------------------------------------
  while (length(active_queue) > 0) {

    # Epidemic check: asymptomatic active cases only
    n_asymp_active <- sum(tdf$asymptomatic[active_queue] == 1L, na.rm = TRUE)
    if (n_asymp_active >= epidemic_queue_threshold) {
      is_epidemic <- TRUE
      break
    }

    active_queue <- active_queue[order(tdf$time_infection[active_queue])]
    idx          <- active_queue[1]
    active_queue <- active_queue[-1]

    is_asymp         <- tdf$asymptomatic[idx]
    t_onset_idx      <- tdf$time_onset[idx]
    t_infectious_idx <- tdf$time_infectious[idx]

    # Reference time for generation time anchor:
    # symptomatic  -> onset time
    # asymptomatic -> time_infectious
    t_idx   <- if (is_asymp == 0L) t_onset_idx else t_infectious_idx
    gen_idx <- tdf$generation[idx]

    # Step 1: Logistical delay
    logistical_delay <- logistical_delay_fn(1)

    # Step 2: Self-treatment and self-quarantine (symptomatic only)
    if (is_asymp == 0L && !is.na(t_idx) && t_idx >= antiviral_start) {

      if (tdf$treated[idx] == 0L && rbinom(1, 1, prob_treat_self) == 1L) {
        tdf$treated[idx]      <- 1L
        tdf$time_treated[idx] <- t_onset_idx + logistical_delay
        tdf$treat_reason[idx] <- 1L
      }

      if (tdf$quarantined[idx] == 0L && rbinom(1, 1, prob_quarantine_self) == 1L) {
        tdf$quarantined[idx]       <- 1L
        tdf$time_quarantined[idx]  <- t_onset_idx + time_to_quarantine_fn(1)
        tdf$quarantine_reason[idx] <- 1L
      }
    }

    # Step 3: Source-side transmission reduction (symptomatic treated only)
    efficacy_trans_idx <- 0
    if (tdf$treated[idx] == 1L && !is.na(tdf$time_treated[idx]) && is_asymp == 0L) {
      t_onset_to_drug    <- tdf$time_treated[idx] - t_idx
      efficacy_trans_idx <- drug_eff_trans_fn(t_onset_to_drug)
      tdf$efficacy_transmission[idx] <- efficacy_trans_idx / eff_trans_e_max
    }

    # Step 4: Ring trigger condition
    # Symptomatic  -> triggered from onset (t_onset_idx) if t_idx >= antiviral_start
    # Asymptomatic -> cannot trigger a ring (undetected)
    can_trigger_ring <- (is_asymp == 0L) && !is.na(t_idx) && (t_idx >= antiviral_start)
    t_ring_absolute  <- if (can_trigger_ring) t_onset_idx + logistical_delay else Inf

    nbrs_idx <- neighbors_by_type[[idx]]

    for (grp in contact_groups_template) {

      nbr_ids <- if (grp$ctype == 1L) nbrs_idx$household else nbrs_idx$community
      if (length(nbr_ids) == 0L) next

      for (nbr_id in nbr_ids) {

        # (4a) Ring PEP / quarantine proposal
        if (can_trigger_ring) {
          if (tdf$treated[nbr_id] == 0L && rbinom(1, 1, grp$p_treat) == 1L) {
            tdf$treated[nbr_id]      <- 1L
            tdf$time_treated[nbr_id] <- t_ring_absolute
            tdf$treat_reason[nbr_id] <- grp$treat_code
          }
          if (tdf$quarantined[nbr_id] == 0L && rbinom(1, 1, grp$p_quarantine) == 1L) {
            tdf$quarantined[nbr_id]       <- 1L
            tdf$time_quarantined[nbr_id]  <- t_ring_absolute
            tdf$quarantine_reason[nbr_id] <- grp$quar_code
          }
        }

        # (4b) Transmission — susceptibles only
        if (tdf$status[nbr_id] != 1L) next

        # [Gate 1] SAR
        if (rbinom(1, 1, grp$sar) == 0L) next

        t_infection_nbr <- t_idx + generation_time_fn(1)

        # [Gate 2] Source-side transmission reduction
        if (rbinom(1, 1, efficacy_trans_idx) == 1L) next

        # [Gate 3] Quarantine block
        if (tdf$quarantined[idx] == 1L && !is.na(tdf$time_quarantined[idx])) {
          if (t_infection_nbr > tdf$time_quarantined[idx]) {
            if (rbinom(1, 1, quarantine_efficacy) == 1L) next
          }
        }

        # [Gate 4] Target-side infection prevention
        pep_prevented       <- FALSE
        efficacy_inf_record <- NA_real_
        if (tdf$treated[nbr_id] == 1L && !is.na(tdf$time_treated[nbr_id])) {
          t_since_pep <- tdf$time_treated[nbr_id] - t_idx
          if (grp$ctype == 1L) {
            efficacy_inf        <- drug_eff_inf_hh_fn(t_since_pep)
            efficacy_inf_record <- if (eff_inf_hh_e_max > 0)
              efficacy_inf / eff_inf_hh_e_max else 0
          } else {
            efficacy_inf        <- drug_eff_inf_fn(t_since_pep)
            efficacy_inf_record <- if (eff_inf_e_max > 0)
              efficacy_inf / eff_inf_e_max else 0
          }
          if (rbinom(1, 1, efficacy_inf) == 1L) pep_prevented <- TRUE
        }
        if (pep_prevented) next

        # Infection confirmed
        tdf$status[nbr_id]         <- 2L
        tdf$generation[nbr_id]     <- gen_idx + 1L
        tdf$ancestor_id[nbr_id]    <- idx
        tdf$contact_type[nbr_id]   <- grp$ctype
        tdf$time_infection[nbr_id] <- t_infection_nbr

        is_asymp_nbr             <- rbinom(1, 1, prop_asymptomatic)
        tdf$asymptomatic[nbr_id] <- is_asymp_nbr
        if (is_asymp_nbr == 0L) {
          tdf$time_onset[nbr_id]      <- t_infection_nbr + infection_to_onset_fn(1)
          tdf$time_infectious[nbr_id] <- t_infection_nbr +
            (tdf$time_onset[nbr_id] - t_infection_nbr) * infectious_before_onset_fn(1)
        } else {
          tdf$time_infectious[nbr_id] <- t_infection_nbr +
            infectious_before_onset_fn(1) * infection_to_onset_fn(1)
        }

        if (!is.na(efficacy_inf_record))
          tdf$efficacy_at_infection[nbr_id] <- efficacy_inf_record

        active_queue <- c(active_queue, nbr_id)

      } # neighbor loop
    } # contact group loop

    tdf$offspring_generated[idx] <- TRUE
    tdf$status[idx]              <- 3L

  } # main loop

  tdf$offspring_generated <- NULL

  infected_df <- tdf[!is.na(tdf$time_infection), ]
  infected_df <- infected_df[order(infected_df$time_infection,
                                   infected_df$person_id), ]

  list(
    full          = tdf,
    infected      = infected_df,
    is_epidemic   = is_epidemic,
    active_at_end = length(active_queue)
  )
}