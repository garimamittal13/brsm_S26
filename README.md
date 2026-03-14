# Behavioral Data Analysis of Memory Fidelity (Mnemonic Similarity Task)

A data analysis project studying how event boundaries influence recognition memory and pattern separation using behavioral experiment logs from the Mnemonic Similarity Task (MST).

## Overview

This project analyzes behavioral data from a replication and extension of Morse et al. (2023), using the Mnemonic Similarity Task (MST) to study how event segmentation affects memory.

The core goals were to:
- clean and structure raw experimental logs
- engineer trial-level behavioral variables
- compute memory-related indices
- run statistical analyses to test hypotheses about encoding and recognition

## Dataset

The dataset consists of behavioral logs collected from **160 participants** across three experimental conditions:
- **Both** (n = 51)
- **Item_Only** (n = 56)
- **Task_Only** (n = 53)

Each participant produced:
- one **task/encoding CSV**
- one **test-phase CSV**

The experiment included:
- **40 events of 7 items each** during encoding
- **150 test images** per participant:
  - 60 targets
  - 60 lures
  - 30 foils

## Problem Statement

The project investigates whether **event boundaries** help or hurt memory.

In particular, it tests whether:
- post-boundary items are slower to encode
- post-boundary items have poorer recognition
- event boundaries affect lure discrimination / pattern separation
- lure similarity bins show a graded difficulty effect

## Data Processing Pipeline

The preprocessing workflow included:
- parsing raw participant CSV logs
- classifying stimuli as **target / lure / foil**
- deriving **event position labels**: pre-boundary, mid, post-boundary
- computing reaction times from multiple response fields
- assigning lure similarity bins
- validating dataset integrity through manipulation checks and summary plots

## Key Metrics

Two main behavioral metrics were analyzed:

- **Recognition Memory Index (REC)**  
  Measures target recognition corrected by foil false alarms

- **Lure Discrimination Index (LDI)**  
  Measures the ability to distinguish similar lures from previously seen items

## Methods

Statistical analyses included:
- **Generalized Linear Mixed Models (GLMMs)**
- **Repeated-measures ANOVA**
- **Holm-Bonferroni corrected comparisons**
- **Bayesian paired t-tests**

## Main Findings

- **Post-boundary items were slower to encode**
- **Recognition memory (REC) dropped for post-boundary items**
- **Lure discrimination (LDI) did not significantly differ across event positions**
- **LDI increased monotonically from similarity Bin 1 to Bin 5**
- Boundary effects differed across experimental conditions

## Repository Structure

- `brsm_data/` — raw and processed data files
- analysis scripts / notebooks — preprocessing, modeling, and visualization
- `README.md` — project overview

## Tech Stack

- **R**
- `lme4`
- `emmeans`
- `ggplot2`

## Why this project matters

This repository demonstrates an end-to-end behavioral data analysis workflow:
1. cleaning messy experimental logs
2. engineering trial-level variables
3. validating dataset quality
4. applying mixed-effects statistical models
5. interpreting human memory behavior through data
