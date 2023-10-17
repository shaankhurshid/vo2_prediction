# vo2_prediction
Data analysis scripts for Deep ECG-VO2 study

# citation
Khurshid S, Churchill TW, Diamant N, Di Achille P, Reeder C, Singh P, Friedman SF, Wasfy MM, Alba GA, Maron BA, Systrom DM, Wertheim BM, Ellinor PT, Ho JE, Baggish AL, Batra P, Lubitz SA, Guseh JS. Eur J Prev Cardiol 2023 (https://doi.org/10.1093/eurjpc/zwad321)

# notes
- script 'vo2_ecg_pclr_v2_test.R' fits models and evaluates on the internal MGH holdout set
- script 'vo2_ecg_pclr_v2_validate.R' fits models and evaluates on the external BWH test set
- scripts assume the presence of PCLR-based ECG representations (model here: https://github.com/broadinstitute/ml4h/tree/master/model_zoo/PCLR)
- with PCLR representations created, model outputs can be regenerated using the coefficients and scale factors (see Supplementary Materials in publication above)
