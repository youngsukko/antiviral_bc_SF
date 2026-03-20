# ==============================================================================
# Integer encoding legend (used throughout both functions)
# ------------------------------------------------------------------------------
# status            : 1L = S (susceptible), 2L = I (infected), 3L = R (recovered)
# contact_type      : 0L = seed, 1L = household, 2L = community
# treat_reason      : 1L = self, 2L = household, 3L = community
# quarantine_reason : 1L = self, 2L = household, 3L = community
# asymptomatic      : 0L = symptomatic, 1L = asymptomatic
# ==============================================================================


# ==============================================================================
# Antiviral efficacy functions (global defaults; may be overridden inside
# run_scenario() by locally-scoped closures).
#
# Infection prevention — community contacts (logistic decay in time since PEP):
#   E(t) = e_max / (1 + exp(kappa * (t - tau)))
#
# Infection prevention — household contacts (separate logistic; household
#   prophylaxis trials consistently show higher efficacy than community trials):
#   E_hh(t) = e_max_hh / (1 + exp(kappa_hh * (t - tau_hh)))
#   Default parameters (Hayden 2004 household PEP data):
#     e_max = 0.63, kappa = 6.06, tau = 0.89
#
# Transmission reduction — source-side (anchored to symptom onset):
#   E(t) = e_max / (1 + exp(kappa * (t - tau)))
# ==============================================================================
drug_eff_inf_fn <- function(t)
  antiviral_eff_inf_e_max / (1 + exp(antiviral_eff_inf_kappa * (t - antiviral_eff_inf_tau)))

# Household-specific infection prevention efficacy function.
# Parameters reflect higher household PEP efficacy observed in dedicated
# household trials (e.g., Hayden 2004) relative to community-level estimates.
drug_eff_inf_hh_fn <- function(t)
  antiviral_eff_inf_hh_e_max / (1 + exp(antiviral_eff_inf_hh_kappa * (t - antiviral_eff_inf_hh_tau)))

drug_eff_trans_fn <- function(t)
  antiviral_eff_trans_e_max / (1 + exp(antiviral_eff_trans_kappa * (t - antiviral_eff_trans_tau)))


# ==============================================================================
# precompute_neighbors()
#
# Builds per-person household and community neighbor lists from a net_data
# object. Call this ONCE before running repeated simulations on the same
# network — the result can be passed directly into network_bp_sim().
#
# Arguments:
#   net_data : output of create_sf_network() — needs $hh_edges, $comm_edges,
#              and $nodes (for n_people)
#
# Returns:
#   A list of length n_people. Each element is a list with two integer vectors:
#     $household  : neighbor person_ids connected by a household edge
#     $community  : neighbor person_ids connected by a community edge
# ==============================================================================
precompute_neighbors <- function(net_data) {
  
  n_people <- nrow(net_data$nodes)
  
  # Helper: build a person -> neighbor list from an undirected edge data.frame.
  # Each undirected edge (from, to) is expanded in both directions so every
  # person finds their neighbors when indexed by person_id.
  build_neighbor_list <- function(edges, n) {
    if (nrow(edges) == 0) {
      out    <- vector("list", n)
      out[]  <- list(integer(0))
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
  
  hh_nbr   <- build_neighbor_list(net_data$hh_edges,   n_people)
  comm_nbr <- build_neighbor_list(net_data$comm_edges,  n_people)
  
  lapply(seq_len(n_people), function(i) {
    list(household = hh_nbr[[i]], community = comm_nbr[[i]])
  })
}


# ==============================================================================
# network_bp_sim()
#
# Stochastic branching process on an explicit contact network.
# Transmission occurs only along igraph edges.
# Household and community contacts are handled separately with their own SAR,
# treatment probability, and quarantine probability.
#
# Epidemic termination (Whittle 1955 — asymptomatic-only threshold):
#   Symptomatic cases are subject to strong intervention (self-quarantine,
#   self-treatment, ring PEP triggered by symptom onset). Their contribution
#   to onward transmission is therefore substantially suppressed, making them
#   poor proxies for the uncontrolled epidemic risk. The threshold is therefore
#   evaluated against the number of ASYMPTOMATIC active infections only:
#
#       i_asymp = |{j in active_queue : tdf$asymptomatic[j] == 1L}|
#
#   Under a branching process approximation, the probability that all i_asymp
#   lineages go extinct without causing a major epidemic is (1/R0)^i_asymp.
#   The simulation is declared an epidemic and terminated early when:
#
#       1 - (1/R0)^i_asymp  >=  epidemic_prob_threshold
#   <=> (1/R0)^i_asymp      <=  1 - epidemic_prob_threshold
#   <=> i_asymp             >=  log(1 - epidemic_prob_threshold) / log(1/R0)
#
#   Default threshold = 0.999 (i.e., 99.9% confidence of epidemic).
#   The minimum asymptomatic queue size that triggers this is precomputed once
#   as epidemic_queue_threshold and checked at the top of each loop iteration.
#
# Household vs. community infection prevention efficacy:
#   Pruning gate [4] (target-side infection prevention) uses two separate
#   logistic efficacy functions:
#     - Household contacts : drug_eff_inf_hh_fn  (higher e_max; household PEP trials)
#     - Community contacts : drug_eff_inf_fn      (community-level estimates)
#   The corresponding e_max denominators are eff_inf_hh_e_max and eff_inf_e_max.
#
# Intervention cascade for each infected person (idx):
#   Step 1: Draw logistical delay for this person's ring
#   Step 2: Self-treatment / self-quarantine (if symptomatic + AV active)
#   Step 3: Compute source-side transmission reduction efficacy
#   Step 4: For each neighbor:
#     (4a) Propose ring PEP / quarantine to the neighbor
#     (4b) Attempt transmission through 4 independent pruning gates:
#          [1] SAR Bernoulli trial
#          [2] Source-side transmission reduction (treated index)
#          [3] Quarantine block (only for exposures after quarantine start)
#          [4] Target-side infection prevention (neighbor received PEP)
#              — uses household or community efficacy function by contact type
#
# Key design:
#   tdf is pre-allocated for all N people (person_id == row index -> O(1) lookup).
#   neighbors_by_type is pre-built by precompute_neighbors() and passed in,
#   so repeated calls on the same network incur no redundant computation.
#
# Arguments:
#   neighbors_by_type           : output of precompute_neighbors()
#   node_df                     : net_data$nodes — needs $person_id, $household
#   R0                          : basic reproduction number (used for epidemic threshold)
#   epidemic_prob_threshold     : P(epidemic) threshold to declare early termination
#                                 (default 0.999); evaluated on asymptomatic queue only
#   generation_time_fn          : function(n) — draw n generation times
#   infection_to_onset_fn       : function(n) — draw n incubation periods
#   infectious_before_onset_fn  : function(n) — pulling latent period from incubation period (proportion)
#   prop_asymptomatic           : proportion of infections that are asymptomatic
#   p_inf_household             : per-contact SAR within household
#   p_inf_community             : per-contact SAR in community
#   antiviral_start             : day AV program begins
#   logistical_delay_fn         : function(n) — draw n logistical delays
#   drug_eff_inf_fn             : function(t) — community infection prevention efficacy at t
#   drug_eff_inf_hh_fn          : function(t) — household infection prevention efficacy at t
#   drug_eff_trans_fn           : function(t) — transmission reduction efficacy at t
#   eff_inf_e_max               : denominator for relative community infection prevention efficacy
#   eff_inf_hh_e_max            : denominator for relative household infection prevention efficacy
#   eff_trans_e_max             : denominator for relative transmission reduction efficacy
#   prob_treat_self/household/community    : PEP acceptance probabilities by reason
#   time_to_quarantine_fn       : function(n) — draw n quarantine delays
#   prob_quarantine_self/household/community : quarantine probabilities by reason
#   quarantine_efficacy         : P(block | transmission would occur after quarantine)
#   seeding_cases               : number of index cases at t0
#   t0                          : simulation start time (default 0)
#   seed                        : RNG seed
#
# Returns:
#   list(
#     full          = tdf for all N people,
#     infected      = tdf for infected only, sorted by time_infection,
#     is_epidemic   = logical — TRUE if epidemic threshold was reached,
#     active_at_end = number of active infections when loop exited
#   )
# ==============================================================================
network_bp_sim <- function(
    neighbors_by_type,
    node_df,
    R0,
    epidemic_prob_threshold = 0.999,
    generation_time_fn, infection_to_onset_fn, infectious_before_onset_fn, prop_asymptomatic,
    p_inf_household, p_inf_community,
    antiviral_start, logistical_delay_fn,
    drug_eff_inf_fn, drug_eff_inf_hh_fn, drug_eff_trans_fn,
    eff_inf_e_max, eff_inf_hh_e_max, eff_trans_e_max,
    prob_treat_self, prob_treat_household, prob_treat_community,
    time_to_quarantine_fn,
    prob_quarantine_self, prob_quarantine_household, prob_quarantine_community,
    quarantine_efficacy,
    seeding_cases,
    t0 = 0, seed = 42
) {
  set.seed(seed)
  
  n_people <- nrow(node_df)
  
  # --------------------------------------------------------------------------
  # Precompute the asymptomatic queue size that triggers epidemic declaration.
  #
  # Threshold is based solely on asymptomatic active infections because
  # symptomatic cases trigger strong intervention (self-isolation, ring PEP)
  # and are not reliable proxies for uncontrolled epidemic risk under the
  # intervention model.
  #
  # We want the smallest integer i such that:
  #   1 - (1/R0)^i >= epidemic_prob_threshold
  # Rearranging:
  #   i >= log(1 - epidemic_prob_threshold) / log(1/R0)
  #
  # ceiling() ensures we require strictly at least that many active cases.
  # R0 must be > 1 for a major epidemic to be possible; if R0 <= 1 the
  # threshold is never triggered (set to Inf so the loop always runs to
  # natural extinction).
  # --------------------------------------------------------------------------
  if (R0 > 1) {
    epidemic_queue_threshold <- ceiling(
      log(1 - epidemic_prob_threshold) / log(1 / R0)
    )
  } else {
    epidemic_queue_threshold <- Inf
  }
  
  # --------------------------------------------------------------------------
  # Initialize tdf for all N people.
  # person_id == row index (1-based) -> O(1) row access throughout.
  # All categorical columns use integer encoding (see legend at top of file).
  # --------------------------------------------------------------------------
  tdf <- data.frame(
    person_id             = node_df$person_id,
    household             = node_df$household,
    status                = rep(1L, n_people),   # 1L = S
    generation            = NA_integer_,
    ancestor_id           = NA_integer_,
    contact_type          = NA_integer_,          # 0=seed, 1=household, 2=community
    time_infection        = NA_real_,
    time_onset            = NA_real_,
    time_infectious       = NA_real_,
    asymptomatic          = NA_integer_,          # 0=symptomatic, 1=asymptomatic
    treated               = 0L,
    time_treated          = NA_real_,
    treat_reason          = NA_integer_,          # 1=self, 2=household, 3=community
    efficacy_at_infection = NA_real_,
    efficacy_transmission = NA_real_,
    quarantined           = 0L,
    time_quarantined      = NA_real_,
    quarantine_reason     = NA_integer_,          # 1=self, 2=household, 3=community
    offspring_generated   = FALSE,
    stringsAsFactors      = FALSE
  )
  
  # --------------------------------------------------------------------------
  # Seed case initialization.
  # Seed cases are forced symptomatic (asymptomatic = 0L) so they are
  # observable by the surveillance system and can trigger a ring.
  # --------------------------------------------------------------------------
  seed_ids <- sample(seq_len(n_people), seeding_cases, replace = FALSE)
  
  for (sid in seed_ids) {
    tdf$status[sid]          <- 2L   # I
    tdf$generation[sid]      <- 1L
    tdf$contact_type[sid]    <- 0L   # seed
    tdf$time_infection[sid]  <- t0
    tdf$asymptomatic[sid]    <- 0L
    tdf$time_onset[sid]      <- t0 + infection_to_onset_fn(1)
    tdf$time_infectious[sid] <- t0 + (tdf$time_onset[sid] - t0)*infectious_before_onset_fn(1)
  }
  
  active_queue <- seed_ids
  is_epidemic  <- FALSE
  
  # --------------------------------------------------------------------------
  # Contact group definitions.
  # Built once outside the main loop; contact_type codes match the legend.
  # --------------------------------------------------------------------------
  contact_groups_template <- list(
    list(ctype = 1L, sar = p_inf_household,
         p_treat      = prob_treat_household,
         p_quarantine = prob_quarantine_household,
         treat_code   = 2L,
         quar_code    = 2L),
    list(ctype = 2L, sar = p_inf_community,
         p_treat      = prob_treat_community,
         p_quarantine = prob_quarantine_community,
         treat_code   = 3L,
         quar_code    = 3L)
  )
  
  # --------------------------------------------------------------------------
  # Main simulation loop.
  # Exits when:
  #   (a) active_queue is empty -> outbreak contained (natural extinction)
  #   (b) count of ASYMPTOMATIC cases in active_queue >=
  #       epidemic_queue_threshold -> epidemic declared
  #       Rationale: symptomatic cases are intervention-suppressed and should
  #       not drive the R0-based epidemic probability calculation.
  #       (Whittle 1955: probability of contained outcome is now < 0.1%)
  # --------------------------------------------------------------------------
  while (length(active_queue) > 0) {
    
    # Epidemic check: count only asymptomatic active infections.
    # Symptomatic individuals trigger ring PEP and self-isolate; their
    # transmission potential is substantially curtailed by the intervention
    # model and they are therefore excluded from the Whittle threshold.
    n_asymp_active <- sum(tdf$asymptomatic[active_queue] == 1L, na.rm = TRUE)
    if (n_asymp_active >= epidemic_queue_threshold) {
      is_epidemic <- TRUE
      break
    }
    
    active_queue <- active_queue[order(tdf$time_infection[active_queue])]
    idx          <- active_queue[1]
    active_queue <- active_queue[-1]
    
    # t_idx            <- tdf$time_infection[idx]
    is_asymp         <- tdf$asymptomatic[idx]
    t_onset_idx      <- tdf$time_onset[idx]
    t_infectious_idx <- tdf$time_infectious[idx]
    t_idx            <- if (is_asymp == 0L) t_onset_idx else t_infectious_idx
    gen_idx          <- tdf$generation[idx]
    
    # ------------------------------------------------------------------
    # Step 1: Draw logistical delay for this person's ring.
    # ------------------------------------------------------------------
    logistical_delay <- logistical_delay_fn(1)
    
    # ------------------------------------------------------------------
    # Step 2: Self-treatment and self-quarantine.
    # ------------------------------------------------------------------
    if (is_asymp == 0L && t_idx >= antiviral_start) {
      
      if (tdf$treated[idx] == 0L && rbinom(1, 1, prob_treat_self) == 1L) {
        tdf$treated[idx]      <- 1L
        tdf$time_treated[idx] <- t_onset_idx + logistical_delay
        tdf$treat_reason[idx] <- 1L   # self
      }
      
      if (tdf$quarantined[idx] == 0L && rbinom(1, 1, prob_quarantine_self) == 1L) {
        tdf$quarantined[idx]       <- 1L
        tdf$time_quarantined[idx]  <- t_onset_idx + time_to_quarantine_fn(1)
        tdf$quarantine_reason[idx] <- 1L   # self
      }
    }
    
    # ------------------------------------------------------------------
    # Step 3: Source-side transmission reduction efficacy.
    # ------------------------------------------------------------------
    efficacy_trans_idx <- 0
    if (tdf$treated[idx] == 1L && !is.na(tdf$time_treated[idx]) && is_asymp == 0L) {
      # t_ref       <- if (is_asymp == 0L) t_onset_idx else t_infectious_idx
      t_onset_to_drug    <- tdf$time_treated[idx] - t_idx
      efficacy_trans_idx <- drug_eff_trans_fn(t_onset_to_drug)
      tdf$efficacy_transmission[idx] <- efficacy_trans_idx / eff_trans_e_max
    }
    
    # ------------------------------------------------------------------
    # Step 4: Process all neighbors.
    # ------------------------------------------------------------------
    can_trigger_ring <- (is_asymp == 0L) && (t_idx >= antiviral_start)
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
        
        # (4b) Transmission attempt: susceptibles only
        if (tdf$status[nbr_id] != 1L) next   # 1L = S
        
        # [Pruning 1] SAR-based transmission attempt
        if (rbinom(1, 1, grp$sar) == 0L) next
        
        t_infection_nbr <- t_idx + generation_time_fn(1)
        
        # [Pruning 2] Source-side transmission reduction
        if (rbinom(1, 1, efficacy_trans_idx) == 1L) next
        
        # [Pruning 3] Quarantine block
        if (tdf$quarantined[idx] == 1L && !is.na(tdf$time_quarantined[idx])) {
          if (t_infection_nbr > tdf$time_quarantined[idx]) {
            if (rbinom(1, 1, quarantine_efficacy) == 1L) next
          }
        }
        
        # [Pruning 4] Target-side infection prevention.
        # Household contacts receive drug_eff_inf_hh_fn (higher efficacy
        # consistent with household PEP trial data); community contacts
        # receive drug_eff_inf_fn (community-level estimates).
        pep_prevented       <- FALSE
        efficacy_inf_record <- NA_real_
        if (tdf$treated[nbr_id] == 1L && !is.na(tdf$time_treated[nbr_id])) {
          # t_ref       <- if (is_asymp == 0L) t_onset_idx else t_infectious_idx
          t_since_pep <- tdf$time_treated[nbr_id] - t_idx
          if (grp$ctype == 1L) {
            # Household contact: use household-specific efficacy function
            efficacy_inf        <- drug_eff_inf_hh_fn(t_since_pep)
            efficacy_inf_record <- efficacy_inf / eff_inf_hh_e_max
          } else {
            # Community contact: use community-level efficacy function
            efficacy_inf        <- drug_eff_inf_fn(t_since_pep)
            efficacy_inf_record <- efficacy_inf / eff_inf_e_max
          }
          if (rbinom(1, 1, efficacy_inf) == 1L) pep_prevented <- TRUE
        }
        if (pep_prevented) next
        
        # All pruning passed — infection confirmed
        tdf$status[nbr_id]         <- 2L   # I
        tdf$generation[nbr_id]     <- gen_idx + 1L
        tdf$ancestor_id[nbr_id]    <- idx
        tdf$contact_type[nbr_id]   <- grp$ctype
        tdf$time_infection[nbr_id] <- t_infection_nbr
        
        is_asymp_nbr             <- rbinom(1, 1, prop_asymptomatic)
        tdf$asymptomatic[nbr_id] <- is_asymp_nbr
        if (is_asymp_nbr == 0L) {
          tdf$time_onset[nbr_id]      <- t_infection_nbr + infection_to_onset_fn(1)
          tdf$time_infectious[nbr_id] <- t_infection_nbr + (tdf$time_onset[nbr_id] - t_infection_nbr)*infectious_before_onset_fn(1)
        } else {
          tdf$time_infectious[nbr_id] <- t_infection_nbr + infectious_before_onset_fn(1)*infection_to_onset_fn(1)
        }
        
        if (!is.na(efficacy_inf_record))
          tdf$efficacy_at_infection[nbr_id] <- efficacy_inf_record
        
        active_queue <- c(active_queue, nbr_id)
        
      } # end neighbor loop
    } # end contact group loop
    
    tdf$offspring_generated[idx] <- TRUE
    tdf$status[idx]              <- 3L   # R
    
  } # end main loop
  
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