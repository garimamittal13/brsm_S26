# ============================================================
# MST Event Boundary Analysis – BRSM 2026
# ============================================================

# SET THIS to the folder containing your brsm_data/ folder
setwd("~/Desktop")   # <-- change if needed

# ⚠️  WARNING: participant 00033 in Both_item_task/both_data has TWO task files
# and TWO test files (timestamps 11h18 and 12h22). Delete the earlier/incomplete
# session files before running, otherwise that participant's trials will be doubled.


# Install any missing packages silently before loading
required_pkgs <- c("tidyverse", "ggplot2", "patchwork", "knitr", "scales",
                   "lme4", "emmeans", "BayesFactor", "Matrix")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org", quiet = TRUE)
}

library(tidyverse)
library(ggplot2)
library(patchwork)
library(knitr)
library(scales)
library(lme4)        # GLMMs
library(emmeans)     # marginal means from models
library(BayesFactor) # Bayesian t-tests for null results

theme_set(theme_minimal(base_size = 12))

RESULTS_DIR <- "results"
dir.create(RESULTS_DIR, showWarnings = FALSE)

# --- Object bins (shared by item_only and Both; 192 objects) ---
obj_bins <- read_delim(
  "brsm_data/item_only/Set6 bins.txt",
  delim = "\t", col_names = c("img_num", "bin"), show_col_types = FALSE
) %>% mutate(stim_type = "Objects")

# --- Scene bins (shared by item_only and Both; 192 scenes) ---
scene_bins <- read_delim(
  "brsm_data/item_only/SetScC bins.txt",
  delim = "\t", col_names = c("img_num", "bin"), show_col_types = FALSE
) %>% mutate(stim_type = "Scenes")

# --- Task-only object bins (384 objects) ---
obj_bins_task <- read_delim(
  "brsm_data/task_only/Set6 bins_ob.txt",
  delim = "\t", col_names = c("img_num", "bin"), show_col_types = FALSE
) %>% mutate(stim_type = "Objects")

# Combine for item_only / Both
bins_item_both <- bind_rows(obj_bins, scene_bins)

# --- Scenes mapping (correct response for shoebox question) ---
scenes_map <- read_delim(
  "brsm_data/item_only/scenes_mapping.txt",
  delim = " ", col_names = c("scene_path", "correct_key"),
  show_col_types = FALSE
) %>%
  mutate(
    scene_filename = basename(scene_path),
    scene_img_num  = as.integer(str_extract(scene_filename, "\\d+"))
  ) %>%
  distinct(scene_img_num, .keep_all = TRUE)

cat("item_only/Both bins:", nrow(bins_item_both), "rows\n")
cat("task_only bins:", nrow(obj_bins_task), "rows\n")
cat("Scenes mapping:", nrow(scenes_map), "entries\n")

# Helper: read one CSV, attach participant_id and condition
# Force all columns to character to avoid type clashes across files
read_mst <- function(fpath, condition) {
  df <- read_csv(fpath, col_types = cols(.default = col_character()),
                 show_col_types = FALSE)
  df$participant_id <- str_extract(basename(fpath), "^\\d{5}")
  df$condition      <- condition
  df
}

# Gather file lists per condition
conditions <- list(
  list(label = "Item_Only",  folder = "brsm_data/item_only/item_only_data"),
  list(label = "Both",       folder = "brsm_data/Both_item_task/both_data"),
  list(label = "Task_Only",  folder = "brsm_data/task_only/task_only_data")
)

task_list <- list()
test_list <- list()

for (cond in conditions) {
  files <- list.files(cond$folder, full.names = TRUE, pattern = "\\.csv$")
  task_files <- files[str_detect(basename(files), "_task_")]
  test_files <- files[str_detect(basename(files), "_test_")]
  task_list <- c(task_list, lapply(task_files, read_mst, condition = cond$label))
  test_list <- c(test_list, lapply(test_files, read_mst, condition = cond$label))
}

task_raw <- bind_rows(task_list)
test_raw <- bind_rows(test_list)

cat("Task rows loaded:", nrow(task_raw), "\n")
cat("Test rows loaded:", nrow(test_raw), "\n")
cat("Unique participants (task):",
    n_distinct(paste(task_raw$participant_id, task_raw$condition)), "\n")
cat("Unique participants (test):",
    n_distinct(paste(test_raw$participant_id, test_raw$condition)), "\n")

task_trials <- task_raw %>%
  filter(str_detect(image_path, "Objects|Scenes"))

cat("Trial rows after filtering:", nrow(task_trials), "\n")

task_trials <- task_trials %>%
  mutate(
    rt_during = suppressWarnings(as.numeric(`key_resp_9.rt`)),
    rt_after  = suppressWarnings(as.numeric(`key_resp_8.rt`)),
    encoding_rt = case_when(
      !is.na(rt_during) ~ rt_during,
      !is.na(rt_after)  ~ 3 + rt_after,
      TRUE              ~ NA_real_
    ),
    # Capture the actual key pressed (from either response window)
    enc_key = case_when(
      !is.na(`key_resp_9.keys`) & `key_resp_9.keys` != "" &
        `key_resp_9.keys` != "None" ~ `key_resp_9.keys`,
      !is.na(`key_resp_8.keys`) & `key_resp_8.keys` != "" &
        `key_resp_8.keys` != "None" ~ `key_resp_8.keys`,
      TRUE ~ NA_character_
    )
  )

# Per professor's suggestion: Mid = positions 3, 4, 5 only.
# Positions 2 and 6 are boundary-adjacent and excluded from Mid
# to ensure clean separation from boundary effects.
task_trials <- task_trials %>%
  group_by(participant_id, condition) %>%
  mutate(
    trial_idx    = row_number(),
    pos_in_event = ((trial_idx - 1) %% 7) + 1,    # 1..7
    event_pos    = case_when(
      pos_in_event == 1             ~ "Post",
      pos_in_event == 7             ~ "Pre",
      pos_in_event %in% c(3, 4, 5) ~ "Mid",
      TRUE                          ~ NA_character_  # positions 2 & 6 excluded
    )
  ) %>%
  ungroup()

# Normalise image filename: extract just the basename
task_trials <- task_trials %>%
  mutate(
    filename  = str_replace_all(basename(image_path), "\\\\", "/"),
    img_num   = as.integer(str_extract(filename, "\\d+")),
    stim_type = if_else(str_detect(image_path, "Objects"), "Objects", "Scenes")
  )

cat("Trials per participant (first 5):\n")
task_trials %>% count(participant_id, condition) %>% head(5) %>% print()
cat("\nEvent position counts (positions 2 & 6 should be NA):\n")
task_trials %>% count(pos_in_event, event_pos) %>% print()

task_trials %>%
  filter(!is.na(event_pos)) %>%
  count(event_pos) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  kable(caption = "Event position distribution among included trials (expected: ~20\\% Pre, ~60\\% Mid, ~20\\% Post after excluding positions 2 & 6)",
        booktabs = TRUE)

cat("\nNote: positions 2 and 6 are excluded (boundary-adjacent).\n")
cat("Excluded trials:", sum(is.na(task_trials$event_pos)), "\n")
cat("Included trials:", sum(!is.na(task_trials$event_pos)), "\n")

# --- 1. Encoding response rate (all conditions) ---
enc_response_rate <- task_trials %>%
  group_by(participant_id, condition) %>%
  summarise(
    n_trials    = n(),
    n_responded = sum(!is.na(enc_key)),
    resp_rate   = n_responded / n_trials,
    .groups = "drop"
  )

# --- 2. Scene encoding accuracy (item_only and Both only, which have scenes) ---
scene_trials <- task_trials %>%
  filter(stim_type == "Scenes", condition %in% c("Item_Only", "Both")) %>%
  left_join(scenes_map, by = c("img_num" = "scene_img_num")) %>%
  filter(!is.na(correct_key), !is.na(enc_key)) %>%
  mutate(correct = (enc_key == correct_key))

scene_acc <- scene_trials %>%
  group_by(participant_id, condition) %>%
  summarise(
    n_scene_trials = n(),
    n_correct      = sum(correct),
    scene_accuracy = n_correct / n_scene_trials,
    .groups = "drop"
  )

cat("Participants with response rate:", nrow(enc_response_rate), "\n")
cat("Participants with scene accuracy:", nrow(scene_acc), "\n")

test_trials <- test_raw %>%
  filter(!is.na(image_path), image_path != "") %>%
  mutate(
    filename  = str_replace_all(basename(image_path), "\\\\", "/"),
    folder    = str_extract(image_path, "^[^/\\\\]+"),
    item_type = case_when(
      folder == "Foils"                ~ "Foil",
      str_detect(filename, "b\\.jpg$") ~ "Lure",
      str_detect(filename, "a\\.jpg$") ~ "Target",
      TRUE                             ~ "Unknown"
    ),
    response = `key_resp_3.keys`,
    test_rt  = suppressWarnings(as.numeric(`key_resp_3.rt`))
  )

test_trials %>% count(item_type) %>%
  kable(caption = "Item type counts in test data", booktabs = TRUE)

test_trials <- test_trials %>%
  mutate(
    img_num   = as.integer(str_extract(filename, "\\d+")),
    stim_type = case_when(
      str_detect(image_path, "Objects") ~ "Objects",
      str_detect(image_path, "Scenes")  ~ "Scenes",
      TRUE ~ "Foil"
    )
  )

test_item_both <- test_trials %>%
  filter(condition %in% c("Item_Only", "Both")) %>%
  left_join(bins_item_both, by = c("img_num", "stim_type"))

test_task_only <- test_trials %>%
  filter(condition == "Task_Only") %>%
  left_join(obj_bins_task %>% select(img_num, bin), by = "img_num")

test_trials <- bind_rows(test_item_both, test_task_only)

cat("Lures with bin assigned:",
    sum(!is.na(test_trials$bin) & test_trials$item_type == "Lure"), "\n")

enc_lookup <- task_trials %>%
  filter(!is.na(event_pos)) %>%
  select(participant_id, condition, img_num, stim_type, event_pos) %>%
  distinct()

test_trials <- test_trials %>%
  left_join(enc_lookup,
            by = c("participant_id", "condition", "img_num", "stim_type"),
            suffix = c("", "_enc"))

cat("Targets with event_pos:",
    sum(!is.na(test_trials$event_pos) & test_trials$item_type == "Target"), "\n")
cat("Lures with event_pos:",
    sum(!is.na(test_trials$event_pos) & test_trials$item_type == "Lure"), "\n")

test_trials <- test_trials %>%
  mutate(
    correct_class = case_when(
      item_type == "Target" & response == "o" ~ TRUE,
      item_type == "Lure"   & response == "s" ~ TRUE,
      item_type == "Foil"   & response == "n" ~ TRUE,
      item_type %in% c("Target", "Lure", "Foil") ~ FALSE,
      TRUE ~ NA
    )
  )

# Participant-level test accuracy by item type
test_acc_by_type <- test_trials %>%
  filter(item_type %in% c("Target", "Lure", "Foil"), !is.na(correct_class)) %>%
  group_by(participant_id, condition, item_type) %>%
  summarise(
    n_trials  = n(),
    n_correct = sum(correct_class),
    accuracy  = n_correct / n_trials,
    .groups   = "drop"
  )

# Overall test accuracy per participant
test_acc_overall <- test_trials %>%
  filter(item_type %in% c("Target", "Lure", "Foil"), !is.na(correct_class)) %>%
  group_by(participant_id, condition) %>%
  summarise(
    n_trials  = n(),
    n_correct = sum(correct_class),
    accuracy  = n_correct / n_trials,
    .groups   = "drop"
  )

cat("Mean overall test accuracy:", round(mean(test_acc_overall$accuracy), 3), "\n")

n_task <- task_raw %>%
  distinct(participant_id, condition) %>%
  mutate(phase = "Encoding (Task)")

n_test <- test_raw %>%
  distinct(participant_id, condition) %>%
  mutate(phase = "Recognition (Test)")

n_both_phases <- n_task %>%
  inner_join(n_test %>% select(participant_id, condition),
             by = c("participant_id", "condition")) %>%
  mutate(phase = "Both Phases")

n_counts <- bind_rows(n_task, n_test, n_both_phases) %>%
  count(condition, phase) %>%
  mutate(phase = factor(phase, levels = c("Encoding (Task)",
                                          "Recognition (Test)",
                                          "Both Phases")))

p_n <- ggplot(n_counts, aes(x = condition, y = n, fill = phase)) +
  geom_col(position = position_dodge(0.8), width = 0.6) +
  geom_text(aes(label = n), position = position_dodge(0.8), vjust = -0.3,
            size = 3.5) +
  scale_fill_brewer(palette = "Pastel2") +
  labs(
    title = "Participant Count by Condition and Phase",
    subtitle = "Verifies complete data availability",
    x = "Condition", y = "Number of Participants", fill = "Phase"
  ) +
  theme(legend.position = "bottom")

print(p_n)
ggsave(file.path(RESULTS_DIR, "01_participant_counts.png"), p_n,
       width = 8, height = 5, dpi = 300)

kable(n_counts %>% pivot_wider(names_from = phase, values_from = n),
      caption = "Participant Counts by Condition", booktabs = TRUE)

rr_summary <- enc_response_rate %>%
  group_by(condition) %>%
  summarise(
    N  = n(),
    M  = mean(resp_rate),
    SD = sd(resp_rate),
    SE = sd(resp_rate) / sqrt(n()),
    .groups = "drop"
  )

p_rr <- ggplot(rr_summary, aes(x = condition, y = M, fill = condition)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width = 0.15) +
  geom_jitter(data = enc_response_rate, aes(x = condition, y = resp_rate),
              inherit.aes = FALSE, width = 0.15, alpha = 0.3, size = 1.2) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "grey50") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1.05)) +
  labs(
    title = "Manipulation Check: Encoding Response Rate by Condition",
    subtitle = "Proportion of trials with a key press; dots = individual participants",
    x = "Condition", y = "Response Rate"
  ) +
  theme(legend.position = "none")

print(p_rr)
ggsave(file.path(RESULTS_DIR, "02_encoding_response_rate.png"), p_rr,
       width = 7, height = 5, dpi = 300)

kable(rr_summary, digits = 3,
      caption = "Encoding Response Rate by Condition", booktabs = TRUE)

scene_acc_summary <- scene_acc %>%
  group_by(condition) %>%
  summarise(
    N  = n(),
    M  = mean(scene_accuracy),
    SD = sd(scene_accuracy),
    SE = sd(scene_accuracy) / sqrt(n()),
    .groups = "drop"
  )

p_sacc <- ggplot(scene_acc_summary, aes(x = condition, y = M, fill = condition)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width = 0.15) +
  geom_jitter(data = scene_acc, aes(x = condition, y = scene_accuracy),
              inherit.aes = FALSE, width = 0.15, alpha = 0.3, size = 1.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1.05)) +
  labs(
    title = "Manipulation Check: Scene Encoding Accuracy by Condition",
    subtitle = "Red dashed line = chance (50\\%); dots = individual participants",
    x = "Condition", y = "Scene Encoding Accuracy"
  ) +
  theme(legend.position = "none")

print(p_sacc)
ggsave(file.path(RESULTS_DIR, "03_scene_encoding_accuracy.png"), p_sacc,
       width = 7, height = 5, dpi = 300)

kable(scene_acc_summary, digits = 3,
      caption = "Scene Encoding Accuracy by Condition", booktabs = TRUE)

test_acc_summary <- test_acc_by_type %>%
  group_by(condition, item_type) %>%
  summarise(
    N  = n(),
    M  = mean(accuracy),
    SE = sd(accuracy) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(item_type = factor(item_type, levels = c("Target", "Lure", "Foil")))

p_tacc <- ggplot(test_acc_summary,
                 aes(x = item_type, y = M, fill = condition)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_errorbar(aes(ymin = M - SE, ymax = M + SE),
                position = position_dodge(0.7), width = 0.15) +
  geom_hline(yintercept = 1/3, linetype = "dashed", color = "red") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1.05)) +
  labs(
    title = "Test-Phase Classification Accuracy by Item Type and Condition",
    subtitle = "Target->old; Lure->similar; Foil->new. Red line = chance (33\\%)",
    x = "Item Type", y = "Classification Accuracy", fill = "Condition"
  )

print(p_tacc)
ggsave(file.path(RESULTS_DIR, "04_test_accuracy_by_type.png"), p_tacc,
       width = 9, height = 5, dpi = 300)

# Overall test accuracy table
test_acc_overall_summary <- test_acc_overall %>%
  group_by(condition) %>%
  summarise(
    N  = n(),
    M  = round(mean(accuracy), 3),
    SD = round(sd(accuracy), 3),
    .groups = "drop"
  )
kable(test_acc_overall_summary,
      caption = "Overall Test Accuracy by Condition", booktabs = TRUE)

enc_rt_part <- task_trials %>%
  filter(!is.na(encoding_rt), !is.na(event_pos)) %>%
  group_by(participant_id, condition, event_pos) %>%
  summarise(mean_rt = mean(encoding_rt), .groups = "drop")

enc_rt_summary <- enc_rt_part %>%
  group_by(event_pos) %>%
  summarise(
    M  = mean(mean_rt),
    SE = sd(mean_rt) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(event_pos = factor(event_pos, levels = c("Pre", "Mid", "Post")))

p1 <- ggplot(enc_rt_summary, aes(x = event_pos, y = M, fill = event_pos)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width = 0.15) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "H1: Encoding RT by Event Position",
    x = "Event Position", y = "Mean Encoding RT (s)", fill = "Position"
  )

print(p1)
ggsave(file.path(RESULTS_DIR, "05_encoding_rt_overall.png"), p1,
       width = 7, height = 5, dpi = 300)

enc_rt_cond <- enc_rt_part %>%
  group_by(condition, event_pos) %>%
  summarise(
    M  = mean(mean_rt),
    SE = sd(mean_rt) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(event_pos = factor(event_pos, levels = c("Pre", "Mid", "Post")))

p2 <- ggplot(enc_rt_cond,
             aes(x = event_pos, y = M, fill = condition)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_errorbar(aes(ymin = M - SE, ymax = M + SE),
                position = position_dodge(0.7), width = 0.15) +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "H1 x H5: Encoding RT by Event Position and Condition",
    x = "Event Position", y = "Mean Encoding RT (s)", fill = "Condition"
  )

print(p2)
ggsave(file.path(RESULTS_DIR, "06_encoding_rt_by_condition.png"), p2,
       width = 9, height = 5, dpi = 300)

tc <- task_trials %>%
  filter(!is.na(encoding_rt)) %>%
  group_by(trial_idx) %>%
  summarise(M = mean(encoding_rt), .groups = "drop")

boundaries <- seq(1, 280, by = 7)

p_tc <- ggplot(tc, aes(x = trial_idx, y = M)) +
  geom_line(alpha = 0.6) +
  geom_vline(xintercept = boundaries, color = "red", alpha = 0.25,
             linetype = "dashed") +
  geom_smooth(method = "loess", span = 0.1, se = FALSE, color = "steelblue") +
  labs(
    title = "Encoding RT Time-Course (All Participants Averaged)",
    subtitle = "Red dashed lines = event boundaries (every 7 items)",
    x = "Trial Position (1-280)", y = "Mean RT (s)"
  )

print(p_tc)
ggsave(file.path(RESULTS_DIR, "07_encoding_rt_timecourse.png"), p_tc,
       width = 12, height = 5, dpi = 300)

enc_rt_part_filt <- enc_rt_part %>% filter(!is.na(event_pos))

enc_rt_anova <- aov(mean_rt ~ event_pos + Error(participant_id / event_pos),
                    data = enc_rt_part_filt)
summary(enc_rt_anova)

pairwise.t.test(enc_rt_part_filt$mean_rt, enc_rt_part_filt$event_pos,
                paired = TRUE, p.adjust.method = "holm")

enc_rt_twoway <- aov(mean_rt ~ event_pos * condition +
                       Error(participant_id / event_pos),
                     data = enc_rt_part_filt)
summary(enc_rt_twoway)

# Trial-level data with valid RT and event_pos
rt_trial <- task_trials %>%
  filter(!is.na(encoding_rt), !is.na(event_pos)) %>%
  mutate(event_pos = factor(event_pos, levels = c("Mid", "Post", "Pre")))
# Mid as reference level so Post vs Mid is the key contrast

rt_glmm <- glmer(
  encoding_rt ~ event_pos + (1 | participant_id),
  data   = rt_trial,
  family = Gamma(link = "log")
)

summary(rt_glmm)

# Estimated marginal means (back-transformed to seconds)
rt_emm <- emmeans(rt_glmm, ~ event_pos, type = "response")
rt_emm_df <- as.data.frame(rt_emm) %>%
  mutate(event_pos = factor(event_pos, levels = c("Pre", "Mid", "Post")))

cat("GLMM estimated marginal means (seconds):\n")
print(rt_emm_df)

# Pairwise contrasts from GLMM
pairs(rt_emm, adjust = "holm")

# GLMM with condition interaction (H5)
rt_glmm_cond <- glmer(
  encoding_rt ~ event_pos * condition + (1 | participant_id),
  data   = rt_trial,
  family = Gamma(link = "log")
)
summary(rt_glmm_cond)

# LRT: does adding condition interaction improve fit?
anova(rt_glmm, rt_glmm_cond)

p1_glmm <- ggplot(rt_emm_df,
                  aes(x = event_pos, y = response, fill = event_pos)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.15) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title    = "H1: Encoding RT by Event Position (Gamma GLMM)",
    subtitle = "Model-estimated marginal means; error bars = 95% CI",
    x = "Event Position", y = "Estimated Mean RT (s)", fill = "Position"
  )

print(p1_glmm)
ggsave(file.path(RESULTS_DIR, "05b_encoding_rt_glmm.png"), p1_glmm,
       width = 7, height = 5, dpi = 300)

resp_counts <- test_trials %>%
  filter(item_type != "Unknown", !is.na(response)) %>%
  count(item_type, response) %>%
  group_by(item_type) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    response  = factor(response, levels = c("o", "s", "n"),
                       labels = c("Old", "Similar", "New")),
    item_type = factor(item_type, levels = c("Target", "Lure", "Foil"))
  )

p3 <- ggplot(resp_counts, aes(x = item_type, y = prop, fill = response)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("Old" = "#1b9e77",
                               "Similar" = "#d95f02",
                               "New" = "#7570b3")) +
  labs(
    title = "Response Distribution by Item Type",
    x = "Item Type", y = "Proportion", fill = "Response"
  )

print(p3)
ggsave(file.path(RESULTS_DIR, "08_response_distribution.png"), p3,
       width = 8, height = 5, dpi = 300)

resp_by_pos <- test_trials %>%
  filter(item_type %in% c("Target", "Lure"),
         !is.na(event_pos), !is.na(response)) %>%
  group_by(participant_id, condition, item_type, event_pos, response) %>%
  summarise(n_resp = n(), .groups = "drop") %>%
  group_by(participant_id, condition, item_type, event_pos) %>%
  mutate(prop = n_resp / sum(n_resp)) %>%
  ungroup()

resp_pos_summary <- resp_by_pos %>%
  group_by(item_type, event_pos, response) %>%
  summarise(
    M  = mean(prop),
    SE = sd(prop) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    event_pos = factor(event_pos, levels = c("Pre", "Mid", "Post")),
    response  = factor(response, levels = c("o", "s", "n"),
                       labels = c("Old", "Similar", "New"))
  )

p3b <- ggplot(resp_pos_summary,
              aes(x = event_pos, y = M, fill = response)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = M - SE, ymax = M + SE),
                position = position_dodge(0.8), width = 0.15) +
  facet_wrap(~ item_type) +
  scale_fill_manual(values = c("Old" = "#1b9e77",
                               "Similar" = "#d95f02",
                               "New" = "#7570b3")) +
  labs(
    title = "Response Proportions by Event Position and Item Type",
    x = "Event Position", y = "Mean Proportion", fill = "Response"
  )

print(p3b)
ggsave(file.path(RESULTS_DIR, "09_response_proportions_by_position.png"), p3b,
       width = 10, height = 5, dpi = 300)

p_old_target <- test_trials %>%
  filter(item_type == "Target", !is.na(event_pos)) %>%
  group_by(participant_id, condition, event_pos) %>%
  summarise(p_old = mean(response == "o", na.rm = TRUE), .groups = "drop")

p_old_foil <- test_trials %>%
  filter(item_type == "Foil") %>%
  group_by(participant_id, condition) %>%
  summarise(p_old_foil = mean(response == "o", na.rm = TRUE), .groups = "drop")

rec_df <- p_old_target %>%
  left_join(p_old_foil, by = c("participant_id", "condition")) %>%
  mutate(REC = p_old - p_old_foil)

rec_summary <- rec_df %>%
  group_by(event_pos) %>%
  summarise(M = mean(REC, na.rm = TRUE),
            SE = sd(REC, na.rm = TRUE) / sqrt(n()),
            .groups = "drop") %>%
  mutate(event_pos = factor(event_pos, levels = c("Pre", "Mid", "Post")))

p4 <- ggplot(rec_summary, aes(x = event_pos, y = M, fill = event_pos)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "H2: Recognition Memory Index (REC) by Event Position",
    subtitle = "REC = P(old|Target) - P(old|Foil)",
    x = "Event Position", y = "REC"
  )

print(p4)
ggsave(file.path(RESULTS_DIR, "10_REC_by_event_position.png"), p4,
       width = 7, height = 5, dpi = 300)

rec_cond <- rec_df %>%
  group_by(condition, event_pos) %>%
  summarise(M = mean(REC, na.rm = TRUE),
            SE = sd(REC, na.rm = TRUE) / sqrt(n()),
            .groups = "drop") %>%
  mutate(event_pos = factor(event_pos, levels = c("Pre", "Mid", "Post")))

p5 <- ggplot(rec_cond,
             aes(x = event_pos, y = M, fill = condition)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_errorbar(aes(ymin = M - SE, ymax = M + SE),
                position = position_dodge(0.7), width = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "H2 x H5: REC by Event Position and Condition",
    x = "Event Position", y = "REC", fill = "Condition"
  )

print(p5)
ggsave(file.path(RESULTS_DIR, "11_REC_by_condition.png"), p5,
       width = 9, height = 5, dpi = 300)

rec_anova <- aov(REC ~ event_pos + Error(participant_id / event_pos),
                 data = rec_df)
summary(rec_anova)

pairwise.t.test(rec_df$REC, rec_df$event_pos,
                paired = TRUE, p.adjust.method = "holm")

rec_twoway <- aov(REC ~ event_pos * condition +
                    Error(participant_id / event_pos),
                  data = rec_df)
summary(rec_twoway)

# Trial-level binary: 1 = responded "Old", 0 = otherwise
target_trials_glmm <- test_trials %>%
  filter(item_type == "Target", !is.na(event_pos), !is.na(response)) %>%
  mutate(
    resp_old  = as.integer(response == "o"),
    event_pos = factor(event_pos, levels = c("Mid", "Post", "Pre"))
  )

foil_trials_glmm <- test_trials %>%
  filter(item_type == "Foil", !is.na(response)) %>%
  mutate(resp_old = as.integer(response == "o"))

# Foil baseline model (intercept only — no event_pos for foils)
foil_glmm <- glmer(
  resp_old ~ 1 + (1 | participant_id),
  data   = foil_trials_glmm,
  family = binomial
)
p_old_foil_est <- plogis(fixef(foil_glmm)[["(Intercept)"]])
cat("Estimated P(Old | Foil):", round(p_old_foil_est, 4), "\n")

# Target model: event_pos fixed effect
rec_glmm <- glmer(
  resp_old ~ event_pos + (1 | participant_id),
  data   = target_trials_glmm,
  family = binomial
)
summary(rec_glmm)

# Marginal means on probability scale, then subtract foil baseline → REC
rec_emm <- emmeans(rec_glmm, ~ event_pos, type = "response")
rec_emm_df <- as.data.frame(rec_emm) %>%
  mutate(
    REC      = prob - p_old_foil_est,
    REC_low  = asymp.LCL - p_old_foil_est,
    REC_high = asymp.UCL - p_old_foil_est,
    event_pos = factor(event_pos, levels = c("Pre", "Mid", "Post"))
  )

cat("GLMM-derived REC estimates:\n")
print(rec_emm_df %>% select(event_pos, prob, REC, REC_low, REC_high))

# Pairwise contrasts (log-odds scale, Holm-corrected)
pairs(rec_emm, adjust = "holm")

# GLMM with condition interaction (H5)
rec_glmm_cond <- glmer(
  resp_old ~ event_pos * condition + (1 | participant_id),
  data   = target_trials_glmm,
  family = binomial
)
summary(rec_glmm_cond)
anova(rec_glmm, rec_glmm_cond)

# Bayesian paired t-test for Post vs Mid REC (H2 null check)
rec_post <- rec_df %>% filter(event_pos == "Post") %>% pull(REC)
rec_mid  <- rec_df %>% filter(event_pos == "Mid")  %>% pull(REC)
rec_pre  <- rec_df %>% filter(event_pos == "Pre")  %>% pull(REC)

bf_rec_post_mid <- ttestBF(x = rec_post, y = rec_mid, paired = TRUE)
bf_rec_post_pre <- ttestBF(x = rec_post, y = rec_pre, paired = TRUE)

cat("Bayes Factor (Post vs Mid REC):", extractBF(bf_rec_post_mid)$bf, "\n")
cat("Bayes Factor (Post vs Pre REC):", extractBF(bf_rec_post_pre)$bf, "\n")
cat("BF > 3: evidence for H1 (effect exists)\n")
cat("BF < 1/3: evidence for H0 (null result)\n")

p4_glmm <- ggplot(rec_emm_df,
                  aes(x = event_pos, y = REC, fill = event_pos)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = REC_low, ymax = REC_high), width = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title    = "H2: REC by Event Position (Binomial GLMM)",
    subtitle = "GLMM marginal means − foil baseline; error bars = 95% CI",
    x = "Event Position", y = "REC"
  )

print(p4_glmm)
ggsave(file.path(RESULTS_DIR, "10b_REC_glmm.png"), p4_glmm,
       width = 7, height = 5, dpi = 300)

p_sim_lure <- test_trials %>%
  filter(item_type == "Lure", !is.na(event_pos)) %>%
  group_by(participant_id, condition, event_pos) %>%
  summarise(p_sim = mean(response == "s", na.rm = TRUE), .groups = "drop")

p_sim_foil <- test_trials %>%
  filter(item_type == "Foil") %>%
  group_by(participant_id, condition) %>%
  summarise(p_sim_foil = mean(response == "s", na.rm = TRUE), .groups = "drop")

ldi_df <- p_sim_lure %>%
  left_join(p_sim_foil, by = c("participant_id", "condition")) %>%
  mutate(LDI = p_sim - p_sim_foil)

ldi_summary <- ldi_df %>%
  group_by(event_pos) %>%
  summarise(M = mean(LDI, na.rm = TRUE),
            SE = sd(LDI, na.rm = TRUE) / sqrt(n()),
            .groups = "drop") %>%
  mutate(event_pos = factor(event_pos, levels = c("Pre", "Mid", "Post")))

p6 <- ggplot(ldi_summary, aes(x = event_pos, y = M, fill = event_pos)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "H3: Lure Discrimination Index (LDI) by Event Position",
    subtitle = "LDI = P(similar|Lure) - P(similar|Foil)",
    x = "Event Position", y = "LDI"
  )

print(p6)
ggsave(file.path(RESULTS_DIR, "12_LDI_by_event_position.png"), p6,
       width = 7, height = 5, dpi = 300)

ldi_cond <- ldi_df %>%
  group_by(condition, event_pos) %>%
  summarise(M = mean(LDI, na.rm = TRUE),
            SE = sd(LDI, na.rm = TRUE) / sqrt(n()),
            .groups = "drop") %>%
  mutate(event_pos = factor(event_pos, levels = c("Pre", "Mid", "Post")))

p7 <- ggplot(ldi_cond,
             aes(x = event_pos, y = M, fill = condition)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_errorbar(aes(ymin = M - SE, ymax = M + SE),
                position = position_dodge(0.7), width = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "H3 x H5: LDI by Event Position and Condition",
    x = "Event Position", y = "LDI", fill = "Condition"
  )

print(p7)
ggsave(file.path(RESULTS_DIR, "13_LDI_by_condition.png"), p7,
       width = 9, height = 5, dpi = 300)

ldi_anova <- aov(LDI ~ event_pos + Error(participant_id / event_pos),
                 data = ldi_df)
summary(ldi_anova)

pairwise.t.test(ldi_df$LDI, ldi_df$event_pos,
                paired = TRUE, p.adjust.method = "holm")

ldi_twoway <- aov(LDI ~ event_pos * condition +
                    Error(participant_id / event_pos),
                  data = ldi_df)
summary(ldi_twoway)

lure_trials_glmm <- test_trials %>%
  filter(item_type == "Lure", !is.na(event_pos), !is.na(response)) %>%
  mutate(
    resp_sim  = as.integer(response == "s"),
    event_pos = factor(event_pos, levels = c("Mid", "Post", "Pre"))
  )

# Foil baseline for "Similar" responses
foil_sim_glmm <- glmer(
  as.integer(response == "s") ~ 1 + (1 | participant_id),
  data   = test_trials %>% filter(item_type == "Foil", !is.na(response)),
  family = binomial
)
p_sim_foil_est <- plogis(fixef(foil_sim_glmm)[["(Intercept)"]])
cat("Estimated P(Similar | Foil):", round(p_sim_foil_est, 4), "\n")

# Lure model: event_pos fixed effect
ldi_glmm <- glmer(
  resp_sim ~ event_pos + (1 | participant_id),
  data   = lure_trials_glmm,
  family = binomial
)
summary(ldi_glmm)

ldi_emm <- emmeans(ldi_glmm, ~ event_pos, type = "response")
ldi_emm_df <- as.data.frame(ldi_emm) %>%
  mutate(
    LDI      = prob - p_sim_foil_est,
    LDI_low  = asymp.LCL - p_sim_foil_est,
    LDI_high = asymp.UCL - p_sim_foil_est,
    event_pos = factor(event_pos, levels = c("Pre", "Mid", "Post"))
  )

cat("GLMM-derived LDI estimates:\n")
print(ldi_emm_df %>% select(event_pos, prob, LDI, LDI_low, LDI_high))

pairs(ldi_emm, adjust = "holm")

ldi_glmm_cond <- glmer(
  resp_sim ~ event_pos * condition + (1 | participant_id),
  data   = lure_trials_glmm,
  family = binomial
)
summary(ldi_glmm_cond)
anova(ldi_glmm, ldi_glmm_cond)

# Bayesian paired t-tests for LDI — H3 predicted a null, so BF is crucial
ldi_post <- ldi_df %>% filter(event_pos == "Post") %>% pull(LDI)
ldi_mid  <- ldi_df %>% filter(event_pos == "Mid")  %>% pull(LDI)
ldi_pre  <- ldi_df %>% filter(event_pos == "Pre")  %>% pull(LDI)

bf_ldi_post_mid <- ttestBF(x = ldi_post, y = ldi_mid, paired = TRUE)
bf_ldi_post_pre <- ttestBF(x = ldi_post, y = ldi_pre, paired = TRUE)

cat("Bayes Factor (Post vs Mid LDI):", extractBF(bf_ldi_post_mid)$bf, "\n")
cat("Bayes Factor (Post vs Pre LDI):", extractBF(bf_ldi_post_pre)$bf, "\n")
cat("BF < 1/3: substantial evidence for null (no boundary effect on LDI)\n")
cat("1/3 < BF < 3: inconclusive\n")
cat("BF > 3: evidence for boundary effect\n")

p6_glmm <- ggplot(ldi_emm_df,
                  aes(x = event_pos, y = LDI, fill = event_pos)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = LDI_low, ymax = LDI_high), width = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title    = "H3: LDI by Event Position (Binomial GLMM)",
    subtitle = "GLMM marginal means − foil baseline; error bars = 95% CI",
    x = "Event Position", y = "LDI"
  )

print(p6_glmm)
ggsave(file.path(RESULTS_DIR, "12b_LDI_glmm.png"), p6_glmm,
       width = 7, height = 5, dpi = 300)

p_sim_lure_bin <- test_trials %>%
  filter(item_type == "Lure", !is.na(event_pos), !is.na(bin)) %>%
  group_by(participant_id, condition, event_pos, bin) %>%
  summarise(p_sim = mean(response == "s", na.rm = TRUE), .groups = "drop")

ldi_bin_df <- p_sim_lure_bin %>%
  left_join(p_sim_foil, by = c("participant_id", "condition")) %>%
  mutate(LDI = p_sim - p_sim_foil)

ldi_bin_summary <- ldi_bin_df %>%
  group_by(event_pos, bin) %>%
  summarise(M = mean(LDI, na.rm = TRUE),
            SE = sd(LDI, na.rm = TRUE) / sqrt(n()),
            .groups = "drop") %>%
  mutate(event_pos = factor(event_pos, levels = c("Pre", "Mid", "Post")),
         bin = as.factor(bin))

p8 <- ggplot(ldi_bin_summary,
             aes(x = bin, y = M, color = event_pos, group = event_pos)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "H3 + H4: LDI Across Lure Bins and Event Positions",
    subtitle = "Bin 1 = most similar (hardest), Bin 5 = least similar (easiest)",
    x = "Lure Similarity Bin", y = "LDI", color = "Event\nPosition"
  )

print(p8)
ggsave(file.path(RESULTS_DIR, "14_LDI_by_bin_and_event_position.png"), p8,
       width = 9, height = 5, dpi = 300)

ldi_bin_cond <- ldi_bin_df %>%
  group_by(condition, event_pos, bin) %>%
  summarise(M = mean(LDI, na.rm = TRUE),
            SE = sd(LDI, na.rm = TRUE) / sqrt(n()),
            .groups = "drop") %>%
  mutate(event_pos = factor(event_pos, levels = c("Pre", "Mid", "Post")),
         bin = as.factor(bin))

p9 <- ggplot(ldi_bin_cond,
             aes(x = bin, y = M, color = event_pos, group = event_pos)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width = 0.12) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ condition) +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "H3 + H4 + H5: LDI by Bin, Event Position, and Condition",
    x = "Lure Similarity Bin", y = "LDI", color = "Event\nPosition"
  )

print(p9)
ggsave(file.path(RESULTS_DIR, "15_LDI_by_bin_condition.png"), p9,
       width = 12, height = 5, dpi = 300)

ldi_bin_collapsed <- test_trials %>%
  filter(item_type == "Lure", !is.na(bin)) %>%
  group_by(participant_id, condition, bin) %>%
  summarise(p_sim = mean(response == "s", na.rm = TRUE), .groups = "drop") %>%
  left_join(p_sim_foil, by = c("participant_id", "condition")) %>%
  mutate(LDI = p_sim - p_sim_foil)

ldi_bin_col_summary <- ldi_bin_collapsed %>%
  group_by(bin) %>%
  summarise(
    M  = mean(LDI, na.rm = TRUE),
    SE = sd(LDI, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(bin = as.factor(bin))

p10 <- ggplot(ldi_bin_col_summary, aes(x = bin, y = M, group = 1)) +
  geom_line(linewidth = 1, color = "steelblue") +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width = 0.15,
                color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "H4: LDI Across Lure Similarity Bins (Collapsed)",
    subtitle = "Bin 1 = most similar (hardest), Bin 5 = least similar (easiest)",
    x = "Lure Similarity Bin", y = "LDI"
  )

print(p10)
ggsave(file.path(RESULTS_DIR, "16_LDI_by_bin_collapsed.png"), p10,
       width = 7, height = 5, dpi = 300)

ldi_bin_col_cond <- ldi_bin_collapsed %>%
  group_by(condition, bin) %>%
  summarise(
    M  = mean(LDI, na.rm = TRUE),
    SE = sd(LDI, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(bin = as.factor(bin))

p10b <- ggplot(ldi_bin_col_cond,
               aes(x = bin, y = M, color = condition, group = condition)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width = 0.12) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "H4 x H5: LDI Across Bins by Condition",
    x = "Lure Similarity Bin", y = "LDI", color = "Condition"
  )

print(p10b)
ggsave(file.path(RESULTS_DIR, "17_LDI_by_bin_condition_collapsed.png"), p10b,
       width = 9, height = 5, dpi = 300)

ldi_bin_anova <- aov(LDI ~ factor(bin) + Error(participant_id / factor(bin)),
                     data = ldi_bin_collapsed)
summary(ldi_bin_anova)

pairwise.t.test(ldi_bin_collapsed$LDI, ldi_bin_collapsed$bin,
                paired = TRUE, p.adjust.method = "holm")

ldi_bin_pos_anova <- aov(
  LDI ~ factor(bin) * event_pos +
    Error(participant_id / (factor(bin) * event_pos)),
  data = ldi_bin_df
)
summary(ldi_bin_pos_anova)

lure_bin_trials <- test_trials %>%
  filter(item_type == "Lure", !is.na(event_pos), !is.na(bin),
         !is.na(response)) %>%
  mutate(
    resp_sim  = as.integer(response == "s"),
    event_pos = factor(event_pos, levels = c("Mid", "Post", "Pre")),
    bin_num   = as.numeric(bin)   # numeric for slope estimation
  )

# Main effects model
ldi_bin_glmm <- glmer(
  resp_sim ~ bin_num + event_pos + (1 | participant_id),
  data   = lure_bin_trials,
  family = binomial
)
summary(ldi_bin_glmm)

# Interaction model
ldi_bin_glmm_int <- glmer(
  resp_sim ~ bin_num * event_pos + (1 | participant_id),
  data   = lure_bin_trials,
  family = binomial
)
summary(ldi_bin_glmm_int)

# LRT: does the interaction improve fit?
anova(ldi_bin_glmm, ldi_bin_glmm_int)

# Marginal means by bin and event_pos from the interaction model
ldi_bin_emm <- emmeans(ldi_bin_glmm_int,
                       ~ bin_num | event_pos,
                       at   = list(bin_num = 1:5),
                       type = "response")

ldi_bin_emm_df <- as.data.frame(ldi_bin_emm) %>%
  mutate(
    LDI       = prob - p_sim_foil_est,
    LDI_low   = asymp.LCL - p_sim_foil_est,
    LDI_high  = asymp.UCL - p_sim_foil_est,
    bin       = as.factor(bin_num),
    event_pos = factor(event_pos, levels = c("Pre", "Mid", "Post"))
  )

cat("GLMM-derived LDI by bin and event position:\n")
print(ldi_bin_emm_df %>% select(bin, event_pos, LDI, LDI_low, LDI_high))

p8_glmm <- ggplot(ldi_bin_emm_df,
                  aes(x = bin, y = LDI,
                      color = event_pos, group = event_pos)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = LDI_low, ymax = LDI_high), width = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_brewer(palette = "Set2") +
  labs(
    title    = "H3 + H4: LDI by Bin and Event Position (Binomial GLMM)",
    subtitle = "Model-estimated marginal means; Bin 1 = hardest, Bin 5 = easiest",
    x = "Lure Similarity Bin", y = "LDI", color = "Event\nPosition"
  )

print(p8_glmm)
ggsave(file.path(RESULTS_DIR, "14b_LDI_bin_glmm.png"), p8_glmm,
       width = 9, height = 5, dpi = 300)

test_rt_part <- test_trials %>%
  filter(item_type %in% c("Target", "Lure"),
         !is.na(event_pos), !is.na(test_rt)) %>%
  group_by(participant_id, condition, item_type, event_pos) %>%
  summarise(mean_rt = mean(test_rt), .groups = "drop")

test_rt_summary <- test_rt_part %>%
  group_by(item_type, event_pos) %>%
  summarise(M = mean(mean_rt),
            SE = sd(mean_rt) / sqrt(n()),
            .groups = "drop") %>%
  mutate(event_pos = factor(event_pos, levels = c("Pre", "Mid", "Post")))

p11 <- ggplot(test_rt_summary,
              aes(x = event_pos, y = M, fill = item_type)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_errorbar(aes(ymin = M - SE, ymax = M + SE),
                position = position_dodge(0.7), width = 0.15) +
  scale_fill_brewer(palette = "Pastel1") +
  labs(
    title = "Exploratory: Test RT by Event Position and Item Type",
    x = "Event Position", y = "Mean Test RT (s)", fill = "Item Type"
  )

print(p11)
ggsave(file.path(RESULTS_DIR, "18_test_rt_by_event_position.png"), p11,
       width = 9, height = 5, dpi = 300)

rec_overall <- rec_df %>%
  group_by(participant_id, condition) %>%
  summarise(REC = mean(REC, na.rm = TRUE), .groups = "drop")

ldi_overall <- ldi_df %>%
  group_by(participant_id, condition) %>%
  summarise(LDI = mean(LDI, na.rm = TRUE), .groups = "drop")

rec_ldi <- rec_overall %>%
  left_join(ldi_overall, by = c("participant_id", "condition"))

p12 <- ggplot(rec_ldi, aes(x = REC, y = LDI, color = condition)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Participant-Level REC vs LDI by Condition",
    x = "Recognition Memory (REC)", y = "Lure Discrimination (LDI)",
    color = "Condition"
  )

print(p12)
ggsave(file.path(RESULTS_DIR, "19_REC_vs_LDI_scatter.png"), p12,
       width = 8, height = 6, dpi = 300)

p13a <- ggplot(rec_overall, aes(x = condition, y = REC, fill = condition)) +
  geom_violin(alpha = 0.4) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "REC Distribution", x = "Condition", y = "REC") +
  theme(legend.position = "none")

p13b <- ggplot(ldi_overall, aes(x = condition, y = LDI, fill = condition)) +
  geom_violin(alpha = 0.4) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "LDI Distribution", x = "Condition", y = "LDI") +
  theme(legend.position = "none")

p13 <- p13a + p13b +
  plot_annotation(title = "Participant-Level Memory Metrics by Condition")

print(p13)
ggsave(file.path(RESULTS_DIR, "20_REC_LDI_violin.png"), p13,
       width = 12, height = 5, dpi = 300)

compute_diff_ci <- function(df, value_col, metric_label) {
  wide <- df %>%
    filter(event_pos %in% c("Post", "Mid")) %>%
    select(participant_id, condition, event_pos, value = all_of(value_col)) %>%
    pivot_wider(names_from = event_pos, values_from = value) %>%
    mutate(diff = Post - Mid)
  
  data.frame(
    metric = metric_label,
    M      = mean(wide$diff, na.rm = TRUE),
    SE     = sd(wide$diff, na.rm = TRUE) / sqrt(sum(!is.na(wide$diff))),
    N      = sum(!is.na(wide$diff))
  )
}

diff_rt  <- compute_diff_ci(enc_rt_part, "mean_rt", "Encoding RT (s)")
diff_rec <- compute_diff_ci(rec_df,      "REC",     "REC")
diff_ldi <- compute_diff_ci(ldi_df,      "LDI",     "LDI")

effect_df <- bind_rows(diff_rt, diff_rec, diff_ldi) %>%
  mutate(
    lower = M - qt(0.975, N - 1) * SE,
    upper = M + qt(0.975, N - 1) * SE,
    metric = factor(metric, levels = c("Encoding RT (s)", "REC", "LDI"))
  )

p14 <- ggplot(effect_df, aes(x = metric, y = M)) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  size = 0.8, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(
    title = "Event Boundary Effect: Post - Mid Difference (95% CI)",
    subtitle = "Positive = Post > Mid; Red dashed line = no difference",
    x = NULL, y = "Post - Mid Difference"
  )

print(p14)
ggsave(file.path(RESULTS_DIR, "21_effect_size_summary.png"), p14,
       width = 8, height = 4, dpi = 300)

kable(effect_df %>% select(metric, M, SE, lower, upper, N),
      digits = 4,
      caption = "Post - Mid Difference Summary (95\\% CI)", booktabs = TRUE)

# Descriptive stats for Encoding RT
enc_table <- enc_rt_part %>%
  group_by(condition, event_pos) %>%
  summarise(
    N  = n(),
    M  = round(mean(mean_rt), 4),
    SD = round(sd(mean_rt), 4),
    .groups = "drop"
  )
kable(enc_table, caption = "Encoding RT Descriptives", booktabs = TRUE)

# REC descriptive stats
rec_table <- rec_df %>%
  group_by(condition, event_pos) %>%
  summarise(
    N  = n(),
    M  = round(mean(REC, na.rm = TRUE), 4),
    SD = round(sd(REC, na.rm = TRUE), 4),
    .groups = "drop"
  )
kable(rec_table, caption = "REC Descriptives", booktabs = TRUE)

# LDI descriptive stats
ldi_table <- ldi_df %>%
  group_by(condition, event_pos) %>%
  summarise(
    N  = n(),
    M  = round(mean(LDI, na.rm = TRUE), 4),
    SD = round(sd(LDI, na.rm = TRUE), 4),
    .groups = "drop"
  )
kable(ldi_table, caption = "LDI Descriptives", booktabs = TRUE)

# LDI by bin descriptives
ldi_bin_table <- ldi_bin_collapsed %>%
  group_by(bin) %>%
  summarise(
    N  = n(),
    M  = round(mean(LDI, na.rm = TRUE), 4),
    SD = round(sd(LDI, na.rm = TRUE), 4),
    .groups = "drop"
  )
kable(ldi_bin_table, caption = "LDI by Bin Descriptives (Collapsed)",
      booktabs = TRUE)

sessionInfo()


sink("~/Desktop/MST_results_output.txt")

cat("=== DESCRIPTIVES ===\n")
cat("P(Old|Foil):", round(p_old_foil_est, 4), "\n")
cat("P(Similar|Foil):", round(p_sim_foil_est, 4), "\n")

cat("\n=== H1: RT GLMM ===\n")
print(summary(rt_glmm))
print(rt_emm_df)
print(pairs(rt_emm, adjust = "holm"))
print(anova(rt_glmm, rt_glmm_cond))

cat("\n=== H1: RT RM-ANOVA ===\n")
print(summary(enc_rt_anova))

cat("\n=== H2: REC GLMM ===\n")
print(summary(rec_glmm))
print(rec_emm_df)
print(pairs(rec_emm, adjust = "holm"))
cat("BF Post vs Mid REC:", extractBF(bf_rec_post_mid)$bf, "\n")
cat("BF Post vs Pre REC:", extractBF(bf_rec_post_pre)$bf, "\n")
print(anova(rec_glmm, rec_glmm_cond))

cat("\n=== H2: REC RM-ANOVA ===\n")
print(summary(rec_anova))

cat("\n=== H3: LDI GLMM ===\n")
print(summary(ldi_glmm))
print(ldi_emm_df)
print(pairs(ldi_emm, adjust = "holm"))
cat("BF Post vs Mid LDI:", extractBF(bf_ldi_post_mid)$bf, "\n")
cat("BF Post vs Pre LDI:", extractBF(bf_ldi_post_pre)$bf, "\n")
print(anova(ldi_glmm, ldi_glmm_cond))

cat("\n=== H3: LDI RM-ANOVA ===\n")
print(summary(ldi_anova))

cat("\n=== H4: BIN GLMM ===\n")
print(summary(ldi_bin_glmm))
print(anova(ldi_bin_glmm, ldi_bin_glmm_int))

cat("\n=== H4: BIN RM-ANOVA ===\n")
print(summary(ldi_bin_anova))

sink()
cat("Saved to ~/Desktop/MST_results_output.txt\n")