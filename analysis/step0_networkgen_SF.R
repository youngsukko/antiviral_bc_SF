# ==============================================================================
# scp_build_network_sf.R
#
# Purpose:
#   Build a realistic contact network for San Francisco using empirical
#   household data (df_merged_sf) and the Prem age-structured contact matrix.
#
# Key difference from scp_build_network.R (RTI + Prem version):
#
#   [Old]  O(N^2) pairwise enumeration
#     Iterates over all i-j pairs directly.
#     N=250 -> ~31,000 pairs; N=10,000 -> ~50,000,000 pairs.
#     Infeasible for SF full population (~870,000 people).
#
#   [New]  O(N * mean_contacts) Poisson sampling
#     For each person i:
#       For each age group ag (1-16):
#         Draw k ~ Poisson(Prem[ag_i, ag])  <- number of contacts with this group
#         Randomly sample k people from age group ag
#         Remove self and household members
#     Scales linearly with N regardless of population size.
#
#   [Parallelization]
#     Edge generation for each person i is independent.
#     Person list is split into chunks and distributed across cores via
#     parallel::parLapply. Duplicate edges (i-j and j-i) are prevented
#     by keeping only pairs where j > i.
#
#   [Contact model note]
#     Household edges  : full clique (all members connected, edge_type = "home")
#     Community edges  : density-dependent Poisson sampling based on Prem matrix
#                        (Prem m_ij values used directly as Poisson lambda;
#                         does NOT divide by n_j, so this is density-dependent,
#                         not frequency-dependent)
#
# Required inputs:
#   df_merged_sf : data frame already loaded in R environment
#                  Required columns: hh_id (household ID), agep (age in years)
#   Prem Excel files in: data/Prem_contact/
#
# Outputs (output/ folder):
#   net_data_sf.rds                 : network object for simulation
#   plot_sf_network_age.png         : network visualization colored by age group
#   plot_sf_degree_distribution.png : degree distribution (household / community / total)
#   plot_sf_age_degree.png          : mean community degree by age group
# ==============================================================================

library(dplyr)
library(tidyr)
library(readxl)
library(igraph)
library(ggplot2)
library(parallel)
library(arrow)

# ==============================================================================
# [Configuration Block] — Modify only this section
# ==============================================================================

# --- Input data (name of object already loaded in R environment) --------------
df_merged_sf <- read_parquet("output/df_hh_sf.parquet")
input_df_name <- "df_merged_sf"   # must contain hh_id and agep columns

# --- Prem contact matrix file paths ------------------------------------------
prem_dir        <- "data/Prem_contact"
prem_country    <- "United States of America"
prem_components <- c("school", "work", "other")   # components to sum for community matrix

# --- Population size cap (NULL = use all households) -------------------------
# If running the full SF population is too slow or memory-intensive,
# set this to a smaller number (e.g., 5000) to sample that many households.
# NULL uses the entire df_merged_sf.
max_households <- NULL
# max_households <- 10000

# --- Parallelization settings ------------------------------------------------
n_cores    <- max(1L, detectCores() - 1L)   # reserve one core for the OS
chunk_size <- 500L   # number of people processed per core per batch
# reduce if memory is tight; increase if you have RAM to spare

# --- Random seed -------------------------------------------------------------
network_seed <- 42

# --- Output directory --------------------------------------------------------
output_dir <- "output"
figure_dir <- "figures"
if (!dir.exists(output_dir)) dir.create(output_dir)


# ==============================================================================
# Function definitions
# ==============================================================================

# ------------------------------------------------------------------------------
# load_prem_matrices()
#
# Load Prem contact matrices from Excel files and sum specified components
# into a single community contact matrix.
#
# Returns: list(community, school, work, other, home)
#   community : summed 16x16 matrix (e.g., school + work + other)
#   others    : individual component matrices (for diagnostics)
# ------------------------------------------------------------------------------
load_prem_matrices <- function(data_dir,
                               country    = "United States of America",
                               components = c("school", "work", "other")) {
  cat("Loading Prem contact matrices...\n")
  
  file_map <- list(
    school = "MUestimates_school_2.xlsx",
    work   = "MUestimates_work_2.xlsx",
    other  = "MUestimates_other_locations_2.xlsx",
    home   = "MUestimates_home_2.xlsx"
  )
  
  load_one <- function(comp) {
    fpath <- file.path(data_dir, file_map[[comp]])
    if (!file.exists(fpath)) {
      warning(sprintf("File not found: %s", fpath))
      return(NULL)
    }
    mat <- as.matrix(read_excel(fpath, sheet = country, col_names = FALSE))
    storage.mode(mat) <- "numeric"
    cat(sprintf("  Loaded: %s (%d x %d)\n", comp, nrow(mat), ncol(mat)))
    mat
  }
  
  mats <- lapply(
    setNames(c("school", "work", "other", "home"),
             c("school", "work", "other", "home")),
    load_one
  )
  
  # Sum specified components into the community matrix
  comm_mat <- matrix(0, 16, 16)
  for (comp in components)
    if (!is.null(mats[[comp]])) comm_mat <- comm_mat + mats[[comp]]
  
  cat(sprintf("  Community matrix: sum of [%s]\n", paste(components, collapse = " + ")))
  
  list(
    community = comm_mat,
    school    = mats$school,
    work      = mats$work,
    other     = mats$other,
    home      = mats$home
  )
}


# ------------------------------------------------------------------------------
# prepare_node_df()
#
# Extract and prepare the node data frame from df_merged_sf.
# Assigns sequential person_id (1 to N) that aligns with igraph vertex indices.
# Optionally samples max_hh households if the full population is too large.
#
# Returns: data.frame(person_id, hh_id, agep, age_group)
#   person_id  : sequential integer ID (1 to N), matches igraph vertex number
#   hh_id      : original household ID from input data
#   agep       : age in years
#   age_group  : 5-year age band index (1 = 0-4y, 2 = 5-9y, ..., 16 = 75+y)
# ------------------------------------------------------------------------------
prepare_node_df <- function(df, max_hh = NULL, seed = 42) {
  
  cat("Preparing node data frame...\n")
  cat(sprintf("  Input rows      : %d\n", nrow(df)))
  
  all_hh_ids <- unique(df$hh_id)
  cat(sprintf("  Total households: %d\n", length(all_hh_ids)))
  
  # Optionally subsample households
  if (!is.null(max_hh) && max_hh < length(all_hh_ids)) {
    set.seed(seed)
    sampled_hh <- sample(all_hh_ids, max_hh)
    df <- df[df$hh_id %in% sampled_hh, ]
    cat(sprintf("  Household cap applied: %d households\n", max_hh))
  }
  
  # Re-assign sequential person_id starting from 1
  # This ensures person_id == igraph vertex index (required for O(1) tdf lookup)
  node_df <- df %>%
    select(hh_id, agep) %>%
    mutate(
      person_id = row_number(),
      age_group = pmin(floor(agep / 5) + 1L, 16L)   # cap at group 16 (75+)
    ) %>%
    select(person_id, hh_id, agep, age_group) %>%
    as.data.frame()
  
  hh_sizes <- table(node_df$hh_id)
  cat(sprintf("  Final population : %d people\n", nrow(node_df)))
  cat(sprintf("  Final households : %d\n", length(hh_sizes)))
  cat(sprintf("  Mean hh size     : %.2f\n", mean(hh_sizes)))
  cat(sprintf("  HH size range    : %d - %d\n", min(hh_sizes), max(hh_sizes)))
  
  node_df
}


# ------------------------------------------------------------------------------
# create_household_edges()
#
# Connect all members within each household as a complete clique (edge_type = "home").
# Every household member pair gets an edge unconditionally — no probability involved.
# Runs in O(sum of hh_size^2) which is effectively O(N) for typical household sizes.
#
# Returns: data.frame(from, to, edge_type = "home")
# ------------------------------------------------------------------------------
create_household_edges <- function(node_df) {
  
  hh_list <- split(node_df$person_id, node_df$hh_id)
  
  edge_list <- lapply(hh_list, function(members) {
    if (length(members) < 2) return(NULL)   # single-person household: no edges
    pairs <- combn(members, 2)
    data.frame(
      from      = pairs[1, ],
      to        = pairs[2, ],
      edge_type = 1L,
      stringsAsFactors = FALSE
    )
  })
  
  result <- bind_rows(Filter(Negate(is.null), edge_list))
  cat(sprintf("  Household edges created: %d\n", nrow(result)))
  result
}


# ------------------------------------------------------------------------------
# create_community_edges_fast()
#
# Generate community (non-household) edges using Poisson sampling from the
# Prem contact matrix. Runs in O(N * mean_contacts) — suitable for large N.
#
# Algorithm (density-dependent contact model):
#   For each person i:
#     For each age group ag in 1:16:
#       lambda = Prem[age_group_i, ag]   <- expected contacts with this age group
#       k ~ Poisson(lambda)              <- actual number of contacts drawn
#       Sample k people from age group ag (excluding self and household members)
#       Keep only pairs where sampled_id > i  (prevents duplicate edges)
#
# Note on contact model:
#   This is DENSITY-DEPENDENT: Prem m_ij is used directly as Poisson lambda.
#   The number of contacts a person makes does NOT depend on group size n_j.
#   (Contrast: frequency-dependent would use lambda = m_ij / n_j * N)
#
# Parallelization:
#   People are split into chunks of size `chunk_size`.
#   Each chunk is processed independently on a separate core via parLapply.
#   Each chunk receives a unique random seed for reproducibility.
#
# Arguments:
#   node_df               : data.frame with person_id, hh_id, age_group
#   prem_community_matrix : 16x16 Prem contact matrix (community = school+work+other)
#   n_cores               : number of parallel cores to use
#   chunk_size            : number of people per processing chunk
#   seed                  : base random seed (chunk i uses seed + i)
#
# Returns: data.frame(from, to, edge_type = "community")
# ------------------------------------------------------------------------------
create_community_edges_fast <- function(node_df,
                                        prem_community_matrix,
                                        n_cores    = 1L,
                                        chunk_size = 500L,
                                        seed       = 42) {
  N         <- nrow(node_df)
  age_group <- node_df$age_group
  hh_id     <- node_df$hh_id
  
  # Pre-build per-age-group member index lists for O(1) candidate lookup
  ag_members <- lapply(1:16, function(ag) which(age_group == ag))
  names(ag_members) <- as.character(1:16)
  
  cat(sprintf("  Building community edges... (%d people | %d cores | chunk = %d)\n",
              N, n_cores, chunk_size))
  
  # Split all person indices into chunks
  all_ids  <- seq_len(N)
  chunks   <- split(all_ids, ceiling(all_ids / chunk_size))
  n_chunks <- length(chunks)
  cat(sprintf("  Total chunks to process: %d\n", n_chunks))
  
  # Inner function: process one chunk of people
  # Returns a data.frame of (from, to, edge_type) for this chunk
  process_chunk <- function(chunk_ids, chunk_seed) {
    set.seed(chunk_seed)
    from_list <- list()
    to_list   <- list()
    
    for (i in chunk_ids) {
      ag_i <- age_group[i]
      hh_i <- hh_id[i]
      
      for (ag_j in 1:16) {
        # Expected number of contacts with age group ag_j (Poisson lambda)
        lambda <- prem_community_matrix[ag_i, ag_j]
        if (lambda <= 0) next
        
        k <- rpois(1, lambda)   # actual number of contacts this person makes today
        if (k == 0) next
        
        # Candidate pool: members of age group ag_j, excluding self and household
        candidates <- ag_members[[ag_j]]
        candidates <- candidates[candidates != i & hh_id[candidates] != hh_i]
        if (length(candidates) == 0) next
        
        # Sample min(k, available) people without replacement
        k_actual <- min(k, length(candidates))
        sampled  <- sample(candidates, k_actual, replace = FALSE)
        
        # Keep only j > i to avoid creating duplicate edges (i-j and j-i)
        keep <- sampled[sampled > i]
        if (length(keep) > 0) {
          from_list[[length(from_list) + 1]] <- rep(i, length(keep))
          to_list[[length(to_list)   + 1]]   <- keep
        }
      }
    }
    
    if (length(from_list) == 0) return(NULL)
    
    data.frame(
      from      = unlist(from_list),
      to        = unlist(to_list),
      edge_type = 2L,
      stringsAsFactors = FALSE
    )
  }
  
  t0 <- proc.time()
  
  if (n_cores > 1) {
    # Parallel execution: spin up PSOCK cluster
    cl <- makeCluster(n_cores, type = "PSOCK")
    on.exit(stopCluster(cl), add = TRUE)
    
    # Export shared objects needed by each worker
    clusterExport(
      cl,
      varlist = c("age_group", "hh_id", "ag_members", "prem_community_matrix"),
      envir   = environment()
    )
    
    # Pack chunk_ids and chunk_seed together into a single list per chunk,
    # then use parLapply (parMapply does not exist in the parallel package)
    chunk_args <- lapply(seq_along(chunks), function(ci) {
      list(chunk_ids = chunks[[ci]], chunk_seed = seed + ci)
    })
    
    results <- parLapply(cl, chunk_args, function(args) {
      process_chunk(args$chunk_ids, args$chunk_seed)
    })
    
  } else {
    # Single-core execution with progress reporting
    results <- vector("list", n_chunks)
    for (ci in seq_along(chunks)) {
      if (ci %% 50 == 0)
        cat(sprintf("    Processing chunk %d / %d...\n", ci, n_chunks))
      results[[ci]] <- process_chunk(chunks[[ci]], seed + ci)
    }
  }
  
  elapsed <- (proc.time() - t0)[["elapsed"]]
  
  comm_edges <- bind_rows(Filter(Negate(is.null), results))
  
  if (nrow(comm_edges) == 0) {
    cat("  WARNING: No community edges generated.\n")
    return(data.frame(from = integer(0), to = integer(0),
                      edge_type = integer(0), stringsAsFactors = FALSE))
  }
  
  # Safety deduplication (i < j constraint above should prevent most duplicates)
  comm_edges <- comm_edges[!duplicated(paste(comm_edges$from, comm_edges$to)), ]
  
  cat(sprintf("  Community edges created: %d (%.1f sec)\n", nrow(comm_edges), elapsed))
  comm_edges
}


# ------------------------------------------------------------------------------
# create_sf_network()
#
# Master function: calls prepare_node_df, create_household_edges, and
# create_community_edges_fast in sequence, then assembles the igraph object.
#
# Returns: list compatible with network_bp_sim() as net_data argument
#   $graph        : igraph object (vertex attributes: person_id, hh_id, agep, age_group)
#   $nodes        : data.frame (person_id, hh_id, agep, age_group, household)
#   $hh_edges     : household edge data.frame (edge_type = "home")
#   $comm_edges   : community edge data.frame (edge_type = "community")
#   $all_edges    : combined edge data.frame
#   $n_households : number of households
#   $network_type : "SF-Prem"
#   $prem_matrices: loaded Prem matrices (for diagnostics)
# ------------------------------------------------------------------------------
create_sf_network <- function(df, prem_mats,
                              max_hh     = NULL,
                              n_cores    = 1L,
                              chunk_size = 500L,
                              seed       = 42) {
  cat(strrep("=", 62), "\n")
  cat(" Building SF-Prem Network\n")
  cat(strrep("=", 62), "\n\n")
  
  node_df <- prepare_node_df(df, max_hh = max_hh, seed = seed)
  N       <- nrow(node_df)
  
  cat("\nStep 1: Household edges (full clique per household)...\n")
  hh_edges <- create_household_edges(node_df)
  
  cat("\nStep 2: Community edges (Prem Poisson sampling)...\n")
  comm_edges <- create_community_edges_fast(
    node_df               = node_df,
    prem_community_matrix = prem_mats$community,
    n_cores               = n_cores,
    chunk_size            = chunk_size,
    seed                  = seed
  )
  
  all_edges <- bind_rows(hh_edges, comm_edges)
  
  cat("\nStep 3: Building igraph object...\n")
  g <- graph_from_data_frame(d = all_edges, vertices = node_df, directed = FALSE)
  
  # Add 'household' column for compatibility with network_bp_sim()
  node_df$household <- node_df$hh_id
  
  n_hh <- length(unique(node_df$hh_id))
  
  cat(sprintf("\n  Population       : %d\n", N))
  cat(sprintf("  Households       : %d\n", n_hh))
  cat(sprintf("  Household edges  : %d\n", nrow(hh_edges)))
  cat(sprintf("  Community edges  : %d\n", nrow(comm_edges)))
  cat(sprintf("  Total edges      : %d\n", ecount(g)))
  cat(sprintf("  Mean degree      : %.2f\n", mean(degree(g))))
  
  list(
    graph         = g,
    nodes         = node_df,
    hh_edges      = hh_edges,
    comm_edges    = comm_edges,
    all_edges     = all_edges,
    n_households  = n_hh,
    network_type  = "SF-Prem",
    prem_matrices = prem_mats
  )
}


# ==============================================================================
# [Section 1] Load Prem matrices and input data
# ==============================================================================

prem_mats <- load_prem_matrices(prem_dir, prem_country, prem_components)

if (!exists(input_df_name)) {
  stop(sprintf("Object '%s' not found in R environment. Please load the data first.",
               input_df_name))
}
df_input <- get(input_df_name)

cat(sprintf("\nInput data: %s  (%d rows x %d cols)\n",
            input_df_name, nrow(df_input), ncol(df_input)))
cat("Columns used: hh_id, agep\n\n")

# Report population scale before building
n_hh_total  <- length(unique(df_input$hh_id))
n_pop_total <- nrow(df_input)
cat(sprintf("Total households: %d | Total population: %d\n", n_hh_total, n_pop_total))

if (!is.null(max_households)) {
  est_pop <- round(n_pop_total * max_households / n_hh_total)
  cat(sprintf("After cap (~%d HH): estimated ~%d people\n", max_households, est_pop))
} else {
  cat(sprintf("Using full population: %d people\n", n_pop_total))
  if (n_pop_total > 50000) {
    cat("  ** WARNING: population > 50,000.\n")
    cat("     Consider setting max_households to a smaller value for testing.\n\n")
  }
}
cat(sprintf("Parallel: %d cores | chunk_size: %d\n\n", n_cores, chunk_size))


# ==============================================================================
# [Section 2] Build network
# ==============================================================================

t_start <- proc.time()

net_data <- create_sf_network(
  df         = df_input,
  prem_mats  = prem_mats,
  max_hh     = max_households,
  n_cores    = n_cores,
  chunk_size = chunk_size,
  seed       = network_seed
)

elapsed_total <- (proc.time() - t_start)[["elapsed"]]
cat(sprintf("\nTotal build time: %.1f sec (%.1f min)\n",
            elapsed_total, elapsed_total / 60))

rds_path <- file.path(output_dir, "net_data_sf.rds")
saveRDS(net_data, rds_path)
cat(sprintf("Saved: %s\n", rds_path))


# ==============================================================================
# [Section 3] Visualization
# ==============================================================================

library(ggplot2)
library(dplyr)
library(igraph)

# --------------------------------------------------------------------------
# 3-1. Household contact degree distribution
# --------------------------------------------------------------------------
hh_degree <- data.frame(
  person_id = c(net_data$hh_edges$from, net_data$hh_edges$to)
) %>%
  count(person_id, name = "degree") %>%
  right_join(data.frame(person_id = net_data$nodes$person_id), by = "person_id") %>%
  mutate(degree = ifelse(is.na(degree), 0L, degree))

p_hh_deg <- ggplot(hh_degree, aes(x = degree)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white", alpha = 0.85) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(
    title    = "Household Contact Degree Distribution",
    subtitle = sprintf("mean = %.2f | N = %d people", mean(hh_degree$degree), nrow(hh_degree)),
    x = "Number of household contacts",
    y = "Number of individuals"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

path_hh_deg <- file.path(figure_dir, "plot_sf_hh_degree.png")
ggsave(path_hh_deg, p_hh_deg, width = 7, height = 5, dpi = 150)
cat(sprintf("Saved: %s\n", path_hh_deg))


# --------------------------------------------------------------------------
# 3-2. Community contact degree distribution
# --------------------------------------------------------------------------
comm_degree <- data.frame(
  person_id = c(net_data$comm_edges$from, net_data$comm_edges$to)
) %>%
  count(person_id, name = "degree") %>%
  right_join(data.frame(person_id = net_data$nodes$person_id), by = "person_id") %>%
  mutate(degree = ifelse(is.na(degree), 0L, degree))

p_comm_deg <- ggplot(comm_degree, aes(x = degree)) +
  geom_histogram(binwidth = 1, fill = "tomato", color = "white", alpha = 0.85) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(
    title    = "Community Contact Degree Distribution",
    subtitle = sprintf("mean = %.2f | N = %d people", mean(comm_degree$degree), nrow(comm_degree)),
    x = "Number of community contacts",
    y = "Number of individuals"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

path_comm_deg <- file.path(figure_dir, "plot_sf_comm_degree.png")
ggsave(path_comm_deg, p_comm_deg, width = 7, height = 5, dpi = 150)
cat(sprintf("Saved: %s\n", path_comm_deg))


# --------------------------------------------------------------------------
# 3-3. Network visualization — 100 sampled households
# --------------------------------------------------------------------------
set.seed(123)
sampled_hh <- sample(unique(net_data$nodes$household), min(100, net_data$n_households))

sub_nodes <- net_data$nodes %>%
  filter(household %in% sampled_hh)
sub_ids <- sub_nodes$person_id  # integer vector

sub_hh_edges <- net_data$hh_edges %>%
  filter(from %in% sub_ids & to %in% sub_ids)

sub_comm_edges <- net_data$comm_edges %>%
  filter(from %in% sub_ids & to %in% sub_ids)

sub_edges <- bind_rows(sub_hh_edges, sub_comm_edges)

g_sub <- graph_from_data_frame(
  d        = sub_edges[, c("from", "to", "edge_type")],
  directed = FALSE,
  vertices = data.frame(name = sub_nodes$person_id,
                        household = sub_nodes$household)
)

set.seed(123)
layout_coords <- layout_with_fr(g_sub)

# node_plot: integer person_id 그대로 사용
node_plot <- sub_nodes %>%
  mutate(
    x  = layout_coords[, 1],
    y  = layout_coords[, 2]
  )

# graph_from_data_frame은 vertex name을 character로 바꾸므로
# edge_df의 from/to도 character로 맞춰서 join
edge_df <- as_data_frame(g_sub, what = "edges") %>%
  mutate(from_int = as.integer(from), to_int = as.integer(to))

edge_plot <- edge_df %>%
  left_join(node_plot %>% select(person_id, x, y), by = c("from_int" = "person_id")) %>%
  rename(x_start = x, y_start = y) %>%
  left_join(node_plot %>% select(person_id, x, y), by = c("to_int" = "person_id")) %>%
  rename(x_end = x, y_end = y)

p_net <- ggplot() +
  geom_segment(
    data  = edge_plot %>% filter(edge_type == 2L),
    aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
    color = "tomato", alpha = 0.2, linewidth = 0.3
  ) +
  geom_segment(
    data  = edge_plot %>% filter(edge_type == 1L),
    aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
    color = "steelblue", alpha = 0.8, linewidth = 0.7
  ) +
  geom_point(
    data  = node_plot,
    aes(x = x, y = y, color = factor(household)),
    size  = 2.0, show.legend = FALSE
  ) +
  labs(
    title    = "Network Visualization — 100 Sampled Households",
    subtitle = sprintf("Blue edges: household  |  Red edges: community  |  Color: household ID\n%d nodes, %d edges shown",
                       nrow(node_plot), nrow(edge_plot)),
    x = NULL, y = NULL
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40")
  )

path_net <- file.path(figure_dir, "plot_sf_network_sample.png")
ggsave(path_net, p_net, width = 10, height = 10, dpi = 150)
cat(sprintf("Saved: %s\n", path_net))

library(ggforce)

# --------------------------------------------------------------------------
# 6-3. Network visualization — synthetic example network
# --------------------------------------------------------------------------
set.seed(123)

# --- Parameters ---
n_hh_vis    <- 100
hh_size_min <- 2
hh_size_max <- 5
jitter_r    <- 0.28   # member scatter radius around household center

# --- Synthetic Prem-style contact matrix ---
prem_syn <- matrix(0.3, nrow = 16, ncol = 16)
diag(prem_syn) <- 3.0
prem_syn[3:4, 3:4] <- prem_syn[3:4, 3:4] + 2.5
prem_syn[5:9, 5:9] <- prem_syn[5:9, 5:9] + 1.5
prem_syn[3:4, 5:9] <- prem_syn[3:4, 5:9] + 0.8
prem_syn[5:9, 3:4] <- prem_syn[5:9, 3:4] + 0.8

# --- Node generation ---
hh_sizes_vis <- sample(hh_size_min:hh_size_max, n_hh_vis, replace = TRUE)
age_group_raw <- rgamma(sum(hh_sizes_vis), shape = 3, rate = 0.3)
age_groups_vis <- pmin(pmax(round(age_group_raw), 1L), 16L)

syn_nodes <- data.frame(
  person_id = seq_len(sum(hh_sizes_vis)),
  household = rep(seq_len(n_hh_vis), times = hh_sizes_vis),
  age_group = age_groups_vis
)
N_vis <- nrow(syn_nodes)

# --- Household grid centers (fixed spacing = 1.0) ---
grid_cols  <- ceiling(sqrt(n_hh_vis))
grid_rows  <- ceiling(n_hh_vis / grid_cols)
hh_centers <- data.frame(
  household = seq_len(n_hh_vis),
  cx = ((seq_len(n_hh_vis) - 1) %% grid_cols),
  cy = floor((seq_len(n_hh_vis) - 1) / grid_cols)
)

# Place members around household center with tight jitter
syn_nodes <- syn_nodes %>%
  left_join(hh_centers, by = "household") %>%
  mutate(
    angle = runif(N_vis, 0, 2 * pi),
    r     = runif(N_vis, 0, jitter_r),
    x     = cx + r * cos(angle),
    y     = cy + r * sin(angle)
  )

# --- Household edges ---
syn_hh_edges <- do.call(rbind, lapply(seq_len(n_hh_vis), function(hh) {
  members <- syn_nodes$person_id[syn_nodes$household == hh]
  if (length(members) < 2) return(NULL)
  pairs <- combn(members, 2)
  data.frame(from = pairs[1,], to = pairs[2,], stringsAsFactors = FALSE)
}))

# --- Community edges (Prem Poisson sampling) ---
ag_members_vis <- lapply(1:16, function(ag)
  syn_nodes$person_id[syn_nodes$age_group == ag]
)
comm_from <- integer(0); comm_to <- integer(0)
for (i in seq_len(N_vis)) {
  ag_i <- syn_nodes$age_group[i]; hh_i <- syn_nodes$household[i]
  for (ag_j in 1:16) {
    k <- rpois(1, prem_syn[ag_i, ag_j])
    if (k == 0) next
    candidates <- ag_members_vis[[ag_j]]
    candidates <- candidates[candidates != i &
                               syn_nodes$household[candidates] != hh_i &
                               candidates > i]
    if (length(candidates) == 0) next
    sampled   <- sample(candidates, min(k, length(candidates)), replace = FALSE)
    comm_from <- c(comm_from, rep(i, length(sampled)))
    comm_to   <- c(comm_to, sampled)
  }
}
syn_comm_edges <- data.frame(from = comm_from, to = comm_to,
                             stringsAsFactors = FALSE)

# --- Edge coordinate tables ---
node_xy <- syn_nodes %>% select(person_id, x, y)

hh_edge_plot <- syn_hh_edges %>%
  left_join(node_xy, by = c("from" = "person_id")) %>% rename(x_start = x, y_start = y) %>%
  left_join(node_xy, by = c("to"   = "person_id")) %>% rename(x_end   = x, y_end   = y)

comm_edge_plot <- syn_comm_edges %>%
  left_join(node_xy, by = c("from" = "person_id")) %>% rename(x_start = x, y_start = y) %>%
  left_join(node_xy, by = c("to"   = "person_id")) %>% rename(x_end   = x, y_end   = y)

# --- Age group color palette ---
age_pal        <- colorRampPalette(c("#4575b4","#91bfdb","#fee090","#fc8d59","#d73027"))(16)
age_labels_vis <- paste0(seq(0,75,5), "-", seq(4,79,5), "y")

# --- Plot ---
p_net <- ggplot() +
  # Community edges (behind everything)
  geom_segment(
    data = comm_edge_plot,
    aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
    color = "gray60", alpha = 0.12, linewidth = 0.2
  ) +
  # Household convex hull bubbles — works for 1, 2, or 3+ members
  geom_mark_ellipse(
    data        = syn_nodes,
    aes(x = x, y = y, group = factor(household)),
    fill        = "gray88",
    color       = "gray65",
    alpha       = 0.35,
    expand      = unit(4, "mm"),   # padding around points
    show.legend = FALSE
  ) +
  # Household edges
  geom_segment(
    data = hh_edge_plot,
    aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
    color = "gray30", alpha = 0.8, linewidth = 0.7
  ) +
  # Nodes colored by age group
  geom_point(
    data  = syn_nodes,
    aes(x = x, y = y, fill = factor(age_group)),
    shape = 21, size = 4.0, color = "white", stroke = 0.4,
    show.legend = TRUE
  ) +
  scale_fill_manual(
    values = age_pal,
    labels = age_labels_vis,
    name   = "Age group"
  ) +
  guides(fill = guide_legend(ncol = 2, override.aes = list(size = 3))) +
  coord_fixed() +   # keep grid square — prevents ellipse distortion
  labs(
    title    = "Example Contact Network Structure (Synthetic)",
    subtitle = sprintf(
      "%d households (%d–%d members)  |  Node color: age group  |  Gray lines: community contacts\nFor illustration only — same contact model as simulation",
      n_hh_vis, hh_size_min, hh_size_max
    ),
    x = NULL, y = NULL
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5, color = "gray40", size = 9),
    legend.position = "right",
    legend.text     = element_text(size = 8),
    legend.title    = element_text(size = 9, face = "bold")
  )

path_net <- file.path(figure_dir, "plot_sf_network_sample.png")
ggsave(path_net, p_net, width = 11, height = 10, dpi = 150)
cat(sprintf("Saved: %s\n", path_net))