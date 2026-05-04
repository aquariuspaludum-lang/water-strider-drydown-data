README
======================================================================
Dataset: Stage-specific effects of habitat desiccation on larval
development, survival, and wing morph determination in four Canadian
water strider species (Hemiptera: Gerridae)

Authors: Manabu Kishi, John R. Spence
Corresponding author: Manabu Kishi (kishi@hotmail.co.jp)
======================================================================


----------------------------------------------------------------------
1. STUDY OVERVIEW
----------------------------------------------------------------------
This dataset contains experimental data on the effects of habitat
desiccation on larval development, survival, flight muscle condition,
wing morph determination, and preoviposition period in four Canadian
water strider species. Nymphs were reared under eight desiccation
conditions in which damp paper was substituted for a water surface
during one specific larval instar, or throughout all instars.

Study species (ordered from more permanent to less permanent habitat):
  - Gerris pingreensis   Drake & Hottes  (large permanent lakes)
  - Gerris buenoi        Kirkaldy        (permanent to temporary)
  - Aquarius remigis     Say             (streams)
  - Limnoporus dissortis Drake & Harris  (temporary shallow ponds)

Collection sites:
  - G. pingreensis, G. buenoi, L. dissortis:
      George Lake Field Site, ~100 km north-west of Edmonton,
      central Alberta, Canada
  - A. remigis:
      Stream near Banff, west Alberta, Canada

Rearing conditions: 19L:5D photoperiod, 20 +/- 2 degrees C


----------------------------------------------------------------------
2. EXPERIMENTAL DESIGN — DESICCATION CONDITIONS
----------------------------------------------------------------------
Eight desiccation conditions were applied. In all conditions, nymphs
were held on a water surface EXCEPT during the instar(s) indicated,
when damp paper was substituted for a water surface.

  Condition 1  Water throughout all instars (Control A)
  Condition 2  Damp paper during 1st instar only
  Condition 3  Damp paper during 2nd instar only
  Condition 4  Damp paper during 3rd instar only
  Condition 5  Damp paper during 4th instar only
  Condition 6  Damp paper during early 5th instar only
  Condition 7  Damp paper during late 5th instar only
  Condition 8  Damp paper throughout all instars (Control D)

Wing morph data (Wing_types.csv) are summarised at the 4-group level:
  Group 1  Condition 1 only         (Control A)
  Group 2  Conditions 2, 3, 4       (damp during 1st-3rd instar)
  Group 3  Conditions 5, 6, 7       (damp during 4th-5th instar)
  Group 4  Condition 8 only         (Control D)


----------------------------------------------------------------------
3. STARTING SAMPLE SIZES
----------------------------------------------------------------------
The following starting N values (total nymphs entering the experiment
per species) are used by the R script to compute 1st-instar mortality.
These values represent n_per_replicate_case x n_replicates.

  G. pingreensis : 325  (65 per case x 5 replicates)
  G. buenoi      : 236  (59 per case x 4 replicates)
  A. remigis     : 207  (23 per case x 9 replicates)
  L. dissortis   : 165  (33 per case x 5 replicates)

Starting N is divided equally across the 8 conditions within each
species (i.e., same N per condition). First-instar deaths are not
recorded in the nymphal period CSV files; they are back-calculated
from starting N minus the number of individuals appearing at 2nd
instar in the nymphal period files.


----------------------------------------------------------------------
4. FILE DESCRIPTIONS
----------------------------------------------------------------------

4.1  Nymphal period data (4 files)
----------------------------------------------------------------------
Files:
  Nymphal_period_G_pingreensis.csv  (12,380 rows)
  Nymphal_period_G_buenoi.csv       ( 6,532 rows)
  Nymphal_period_A_remigis.csv      ( 6,709 rows)
  Nymphal_period_L_dissortis.csv    ( 4,162 rows)

Columns:
  Species        Species name (character)
  Replicate      Replicate case number (integer)
  Stage          Larval stage at which the individual was recorded;
                 values are "2", "3", "4", "5", "Adult" (character)
  Condition      Desiccation condition 1-8 (integer; see Section 2)
  Days_to_stage  Days from hatching to reaching this stage (numeric)

Structure:
  Each individual that survived to adult emergence contributes five
  rows (Stage = "2", "3", "4", "5", "Adult"), one per stage reached.
  Individuals that died before reaching a given stage do not appear
  at that stage or any subsequent stage.
  First-instar deaths are not included (see Section 3).

Note:
  Days_to_stage at Stage "Adult" equals the total nymphal period
  (used in Analysis A). Days_to_stage at intermediate stages equals
  the cumulative days from hatching to moulting into that stage.


4.2  Wing morph data (1 file)
----------------------------------------------------------------------
File: Wing_types.csv  (32 rows)

Columns:
  Species        Species name: "G_buenoi" or "A_remigis" (character)
  Sex            "female" or "male" (character)
  Group          Condition group 1-4 (integer; see Section 2)
  Wing_type      "long_wing" or "short_wing" (character)
  N_of_individuals  Count of individuals (integer)

Notes:
  - Wing morph data are available for G. buenoi and A. remigis only.
    G. pingreensis is apterous (no wing morph variation).
    L. dissortis is always long-winged.
  - Data are aggregated at the 4-group level (not 8 conditions).
    This is a limitation of the original data collection.
  - 2 species x 2 sexes x 4 groups x 2 wing types = 32 rows.


4.3  Flight muscle condition data (3 files)
----------------------------------------------------------------------
Files:
  flight_muscle_G_buenoi.csv    (8 rows)
  flight_muscle_A_remigis.csv   (8 rows)
  flight_muscle_L_dissortis.csv (8 rows)

Columns:
  Condition    Desiccation condition 1-8 (integer; see Section 2)
  Developed    Number of individuals with developed flight muscles
  Intermediate Number of individuals with intermediate flight muscles
  Histolysed   Number of individuals with histolysed flight muscles

Notes:
  - Flight muscle data are available for G. buenoi, A. remigis, and
    L. dissortis only. G. pingreensis is apterous and was excluded.
  - In G. buenoi, no individuals had Developed-type muscles
    (Developed = 0 across all conditions).
  - In L. dissortis, no individuals had Histolysed-type muscles
    (Histolysed = 0 across all conditions).
  - Each row represents aggregated counts for all individuals
    dissected within that condition.


4.4  Preoviposition period data (4 files)
----------------------------------------------------------------------
Files:
  preoviposition_period_G_buenoi_long_wing.csv   (49 rows)
  preoviposition_period_G_buenoi_short_wing.csv  (90 rows)
  preoviposition_period_A_remigis_long_wing.csv  (73 rows)
  preoviposition_period_A_remigis_short_wing.csv (59 rows)

Columns:
  Condition       Desiccation condition 1-8 (integer; see Section 2)
  Preoviposition  Days from adult emergence to first oviposition
                  (individual-level, numeric)

Notes:
  - Preoviposition period data are available for G. buenoi and
    A. remigis only. G. pingreensis and L. dissortis are univoltine
    and show reproductive dormancy; preoviposition data were not
    collected for these species.
  - Data are split by wing morph (long_wing / short_wing).
  - Each row represents one female individual.


4.5  R analysis script (1 file)
----------------------------------------------------------------------
File: analysis_water_strider.R

Statistical analyses used in the manuscript:
  Analysis A  Nymphal period — GLMM (Gamma distribution,
               Condition as fixed effect, Replicate as random effect)
  Analysis B  Stage-specific mortality — GLM (binomial);
               mortality is computed directly from the nymphal period
               CSV files combined with the starting N (see Section 3)
  Analysis C  Wing morph proportion — logistic regression
               (Group x Sex; G. buenoi and A. remigis only)
  Analysis D  Preoviposition period — GLM (Gamma distribution,
               Condition + Wing morph)
  Analysis E  Flight muscle condition (intra-species) — chi-square
  Analysis F  Flight muscle condition (inter-species) — chi-square
               or Fisher's exact test

Required R packages:
  lme4, lmerTest, emmeans, multcomp, tidyverse


----------------------------------------------------------------------
5. SPECIES-SPECIFIC NOTES
----------------------------------------------------------------------
  G. pingreensis  Apterous (no wing morphs); no flight muscle data;
                  univoltine; no preoviposition data.
  G. buenoi       Long-winged and short-winged morphs; multivoltine;
                  no Developed-type flight muscles observed.
  A. remigis      Long-winged and short-winged morphs; multivoltine.
  L. dissortis    Always long-winged; univoltine; no Histolysed-type
                  flight muscles observed; no preoviposition data.


----------------------------------------------------------------------
6. DATA AVAILABILITY
----------------------------------------------------------------------
Repository name : water-strider-drydown-data
URL             : https://github.com/aquariuspaludum-lang/water-strider-drydown-data



----------------------------------------------------------------------
7. RELATED PUBLICATION
----------------------------------------------------------------------
Kishi M. and Spence J.R. (in preparation). Stage-specific effects of
habitat desiccation on larval development, survival, and wing morph
determination in four Canadian water strider species
(Hemiptera: Gerridae).

Target journal: Ecological Entomology


----------------------------------------------------------------------
8. LICENCE
----------------------------------------------------------------------
These data are released under the Creative Commons Attribution 4.0
International licence (CC BY 4.0).
https://creativecommons.org/licenses/by/4.0/
