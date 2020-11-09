##########
#
# Nacer et al.
# Pan-cancer application of a lung-adenocarcinoma-derived 
# gene-expression-based prognostic predictor
#
# Names for figure legends
#
##########


full_names <- c('THCA_TCGA' = 'Thyroid carcinoma (THCA)', 
                'KIP_TCGA' = 'Kidney renal papillary cell carcinoma (KIP)', 
                'PRAD_TCGA' = 'Prostate adenocarcinoma (PRAD)',
                'LUAD_TCGA' = 'Lung adenocarcinoma (LUAD)', 
                'KICh_TCGA' = 'Kidney chromophobe (KICh)', 
                'KICC_TCGA' = 'Kidney renal clear cell carcinoma (KICC)',
                # 'mixed_TCGA' = 'Kidney cancer mixed subtype (KImixed)',
                'PCPG_TCGA' = 'Pheochromocytoma and paraganglioma (PCPG)', 
                'LGG_TCGA' = 'Brain lower grade glioma (LGG)', 
                'ACC_TCGA' = 'Adrenocortical carcinoma (ACC)',
                'LIHC_TCGA' = 'Liver hepatocellular carcinoma (LIHC)', 
                'PAAD_TCGA' = 'Pancreatic adenocarcinoma (PAAD)', 
                'BRCA_TCGA' = 'Breast invasive carcinoma (BRCA)',
                'MESO_TCGA' = 'Mesothelioma (MESO)', 
                'CHOL_TCGA' = 'Cholangiocarcinoma (CHOL)', 
                'THYM_TCGA' = 'Thymoma (THYM)', 
                'STAD_TCGA' = 'Stomach adenocarcinoma (STAD)',
                'LUSC_TCGA' = 'Lung squamous cell carcinoma (LUSC)', 
                'UVM_TCGA' = 'Uveal melanoma (UVM)', 
                'SARC_TCGA' = 'Sarcoma (SARC)', 
                'BLCA_TCGA' = 'Bladder urothelial carcinoma (BLCA)',
                'OV_TCGA' = 'Ovarian serous cystadenocarcinoma (OV)', 
                'GBM_TCGA' = 'Glioblastoma multiforme (GBM)', 
                'HNSC_TCGA' = 'Head and neck squamous cell carcinoma (HNSC)',
                'ESCA_TCGA' = 'Esophageal carcinoma (ESCA)', 
                'COAD_TCGA' = 'Colon adenocarcinoma (COAD)', 
                'READ_TCGA' = 'Rectum adenocarcinoma (READ)',
                'UCEC_TCGA' = 'Uterine corpus endometrial carcinoma (UCEC)', 
                'UCS_TCGA' = 'Uterine carcinosarcoma (UCS)',
                'CESC_TCGA' = 'Cervical squamous cell carcinoma and endocervical adenocarcinoma (CESC)',
                'TGCT_TCGA' = 'Testicular germ cell tumors (TGCT)', 
                'SKCM_TCGA' = 'Skin cutaneous melanoma (SKCM)', 
                'DLBC_TCGA' = 'Lymphoid neoplasm diffuse large B-cell lymphoma (DLBC)',
                'BRCA_GOBO' = '[GOBO] Breast invasive carcinoma (BRCA)', 
                'BRCA_SCAN-B' = '[SCAN-B] Breast invasive carcinoma (BRCA)',
                'BRCA_Tabchy' = '[Tabchy] Breast cancer',
                'BRCA_Iwamoto' = '[Iwamoto] Breast cancer',
                'BRCA_Hess' = '[Hess] Breast cancer',
                'MEL_Gide' = '[Gide] Melanoma',
                'MEL_Liu' = '[Liu] Melanoma',
                'MEL_Riaz' = '[Riaz] Melanoma')

acronym_only <- c('THCA_TCGA' = 'THCA', 'KIP_TCGA' = 'KIP', 'PRAD_TCGA' = 'PRAD',
                  'LUAD_TCGA' = 'LUAD', 'KICh_TCGA' = 'KICh', 'KICC_TCGA' = 'KICC',
                  'PCPG_TCGA' = 'PCPG', 'LGG_TCGA' = 'LGG', 'ACC_TCGA' = 'ACC',
                  'LIHC_TCGA' = 'LIHC', 'PAAD_TCGA' = 'PAAD', 'BRCA_TCGA' = 'BRCA',
                  'MESO_TCGA' = 'MESO', 'CHOL_TCGA' = 'CHOL', 'THYM_TCGA' = 'THYM', 'STAD_TCGA' = 'STAD',
                  'LUSC_TCGA' = 'LUSC', 'UVM_TCGA' = 'UVM', 'SARC_TCGA' = 'SARC', 'BLCA_TCGA' = 'BLCA',
                  'OV_TCGA' = 'OV', 'GBM_TCGA' = 'GBM', 'HNSC_TCGA' = 'HNSC',
                  'ESCA_TCGA' = 'ESCA', 'COAD_TCGA' = 'COAD', 'READ_TCGA' = 'READ',
                  'UCEC_TCGA' = 'UCEC', 'UCS_TCGA' = 'UCS', 'CESC_TCGA' = 'CESC',
                  'TGCT_TCGA' = 'TGCT', 'SKCM_TCGA' = 'SKCM', 'DLBC_TCGA' = 'DLBC',
                  'BRCA_GOBO' = '[GOBO] BRCA', 'BRCA_SCAN-B' = '[SCAN-B] BRCA',
                  'BRCA_Tabchy' = '[Tabchy] BRCA',
                  'BRCA_Iwamoto' = '[Iwamoto] BRCA', 'BRCA_Hess' = '[Hess] BRCA',
                  'MEL_Gide' = '[Gide] Melanoma', 'MEL_Liu' = '[Liu] Melanoma',
                  'MEL_Riaz' = '[Riaz] Melanoma')

acronym_ast_os_clams <- c('THCA_TCGA' = 'THCA', 'KIP_TCGA' = '*KIP', 'PRAD_TCGA' = 'PRAD',
                        'LUAD_TCGA' = '*LUAD', 'KICh_TCGA' = 'KICh', 'KICC_TCGA' = 'KICC',
                        'PCPG_TCGA' = 'PCPG', 'LGG_TCGA' = '*LGG', 'ACC_TCGA' = 'ACC',
                        'LIHC_TCGA' = '*LIHC', 'PAAD_TCGA' = 'PAAD', 'BRCA_TCGA' = 'BRCA',
                        'MESO_TCGA' = 'MESO', 'CHOL_TCGA' = 'CHOL', 'THYM_TCGA' = 'THYM', 'STAD_TCGA' = 'STAD',
                        'LUSC_TCGA' = 'LUSC', 'UVM_TCGA' = 'UVM', 'SARC_TCGA' = 'SARC', 'BLCA_TCGA' = 'BLCA',
                        'OV_TCGA' = 'OV', 'GBM_TCGA' = 'GBM', 'HNSC_TCGA' = 'HNSC',
                        'ESCA_TCGA' = 'ESCA', 'COAD_TCGA' = 'COAD', 'READ_TCGA' = 'READ',
                        'UCEC_TCGA' = 'UCEC', 'UCS_TCGA' = 'UCS', 'CESC_TCGA' = 'CESC',
                        'TGCT_TCGA' = 'TGCT', 'SKCM_TCGA' = 'SKCM', 'DLBC_TCGA' = 'DLBC',
                        'BRCA_GOBO' = '[GOBO] *BRCA', 'BRCA_SCAN-B' = '[SCAN-B] *BRCA',
                        'BRCA_Tabchy' = '[Tabchy] Breast cancer',
                        'BRCA_Iwamoto' = '[Iwamoto] Breast cancer', 'BRCA_Hess' = '[Hess] Breast cancer',
                        'MEL_Gide' = '[Gide] Melanoma', 'MEL_Liu' = '[Liu] Melanoma',
                        'MEL_Riaz' = '[Riaz] Melanoma')