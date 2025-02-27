Title: 
Assessment of Case Influence in the Lasso with a Case-weight Adjusted Solution Path

Introduction:
This archive contains the code and datasets that can be used to reproduce all figures and tables in the paper.

Configuration:
The following configurations are used by the authors to run the code.
  
  1. R version 4.4.2 
  2. tidyverse 2.0.0
  3. glmnet 4.1-8
  4. lars 1.3
  5. MASS 7.3-60
  6. latex2exp 0.9.6
  7. mnormt 2.1.1
  8. ggpubr 0.6.0
  9. viridis 0.6.4
  10. hrbrthemes 0.8.0
  11. readxl 1.4.3
  12. reshape2 1.4.4
  13. scales 1.3.0

Structure and contents of the archive:

- caseweightlasso2024.R: for the core Cookdislasso function and some auxiliary functions that may be called by all other code.
- Code4Figures1-3.R: for Figures 1-3
- Code4Figures5-13.R: for Figures 5-13
- Code4Tables2,5-13.R: for Tables 2, 5-13
- Concrete_example(Figure4).R: for Figure 4 of the concrete compressive strength data example
- Diabetes_example(Table3).R: for Table 3 of the diabetes data example 
- Glioblastoma_example(Table4).R: for Table 4 of the glioblastoma gene expression data example

- The datasets used are stored in "datasets" subfolder.
- Tables are saved in "results" subfolder.
  
Usage:
To reproduce each table/figure, you need to change your working directory in R to this folder and run the scripts. 


Runtime estimate:
+--------------------------------+-----------------+------------------------+
|                                | MBP 32GB M2 MAX | Windows 16GB AMD 5600X |
+--------------------------------+-----------------+------------------------+
| Code4Figures1-3.R              |      < 1min     |         < 1min         |
+--------------------------------+-----------------+------------------------+
| Code4Figures5-13.R             |      < 1min     |         < 2min         |
+--------------------------------+-----------------+------------------------+
| Code4Tables2,5-13.R            |      25hrs      |          45hrs         |
+--------------------------------+-----------------+------------------------+
| Concrete_example(Figure4).R    |       8min      |          12min         |
+--------------------------------+-----------------+------------------------+
| Diabetes_example(Table3).R     |      2.5min     |          4min          |
+--------------------------------+-----------------+------------------------+
| Glioblastoma_example(Table4).R |       5min      |          8min          |
+--------------------------------+-----------------+------------------------+
