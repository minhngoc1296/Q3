"""
Generate comprehensive PDF report for publication
Includes: Methods, Tables 1-4, Figures 1-5, Summary
"""
import os, sys
import numpy as np
from fpdf import FPDF

OUTPUT = 'H:\\Q3 MNGOC\\BAOCAO_ATA2025_PUB.pdf'

class PubPDF(FPDF):
    def header(self):
        if self.page_no() > 1:
            self.set_font('ArialUni', 'I', 7)
            self.set_text_color(120, 120, 120)
            self.cell(0, 4, 'Thyroid Cancer Recurrence Risk - ATA 2025 - Firth Logistic Regression', 0, 1, 'C')
            self.ln(2)

    def footer(self):
        self.set_y(-15)
        self.set_font('ArialUni', 'I', 7)
        self.set_text_color(120, 120, 120)
        self.cell(0, 10, f'Page {self.page_no()}/{{nb}}', 0, 0, 'C')

    def section_title(self, title):
        self.set_font('ArialUni', 'B', 12)
        self.set_text_color(0, 51, 102)
        self.cell(0, 8, title, 0, 1, 'L')
        self.set_draw_color(0, 51, 102)
        self.line(10, self.get_y(), 200, self.get_y())
        self.ln(3)
        self.set_text_color(0, 0, 0)

    def subsection_title(self, title):
        self.set_font('ArialUni', 'B', 10)
        self.set_text_color(44, 62, 80)
        self.cell(0, 7, title, 0, 1, 'L')
        self.ln(1)
        self.set_text_color(0, 0, 0)

    def body_text(self, text):
        self.set_font('ArialUni', '', 9)
        self.multi_cell(0, 4.5, text)
        self.ln(1)

    def note_text(self, text):
        self.set_font('ArialUni', 'I', 8)
        self.set_text_color(80, 80, 80)
        self.multi_cell(0, 4, text)
        self.set_text_color(0, 0, 0)
        self.ln(2)

    def make_table(self, headers, data, col_widths):
        self.set_font('ArialUni', 'B', 8)
        self.set_fill_color(0, 51, 102)
        self.set_text_color(255, 255, 255)
        for i, h in enumerate(headers):
            self.cell(col_widths[i], 7, h, 1, 0, 'C', True)
        self.ln()
        self.set_font('ArialUni', '', 8)
        self.set_text_color(0, 0, 0)
        fill = False
        for row in data:
            if fill:
                self.set_fill_color(240, 245, 250)
            else:
                self.set_fill_color(255, 255, 255)
            max_h = 7
            for i, cell in enumerate(row):
                self.cell(col_widths[i], max_h, str(cell), 1, 0, 'C', True)
            self.ln()
            fill = not fill


# ===================== BUILD PDF =====================
pdf = PubPDF()
pdf.alias_nb_pages()
pdf.add_font('ArialUni', '', 'C:\\Windows\\Fonts\\arial.ttf', uni=True)
pdf.add_font('ArialUni', 'B', 'C:\\Windows\\Fonts\\arialbd.ttf', uni=True)
pdf.add_font('ArialUni', 'I', 'C:\\Windows\\Fonts\\ariali.ttf', uni=True)
pdf.add_font('ArialUni', 'BI', 'C:\\Windows\\Fonts\\arialbi.ttf', uni=True)

# ===================== TITLE PAGE =====================
pdf.add_page()
pdf.ln(40)
pdf.set_font('ArialUni', 'B', 18)
pdf.set_text_color(0, 51, 102)
pdf.cell(0, 12, 'Recurrence Risk Prediction', 0, 1, 'C')
pdf.cell(0, 12, 'for Thyroid Cancer', 0, 1, 'C')
pdf.ln(5)
pdf.set_font('ArialUni', 'B', 14)
pdf.set_text_color(52, 73, 94)
pdf.cell(0, 10, 'ATA 2025 PTC Risk Classification', 0, 1, 'C')
pdf.cell(0, 10, 'Firth Penalized Logistic Regression', 0, 1, 'C')
pdf.ln(5)
pdf.set_font('ArialUni', 'B', 11)
pdf.set_text_color(100, 100, 100)
pdf.cell(0, 8, 'Gene-Gene Interaction: BRAF V600E, TERT, TP53', 0, 1, 'C')
pdf.ln(15)
pdf.set_font('ArialUni', '', 10)
pdf.set_text_color(0, 0, 0)
pdf.cell(0, 7, 'Retrospective Study, 114 Patients', 0, 1, 'C')
pdf.cell(0, 7, 'Study Period: 2024', 0, 1, 'C')
pdf.ln(5)
pdf.set_font('ArialUni', 'I', 9)
pdf.set_text_color(120, 120, 120)
pdf.multi_cell(0, 5, "Analysis performed using Firth's penalized likelihood to address complete separation in rare mutation groups (TERT n=6, TP53 n=2). Internal validation via bootstrap (1000 reps), calibration curve, and LOOCV.")
pdf.ln(10)
pdf.set_font('ArialUni', '', 8)
pdf.set_text_color(150, 150, 150)
pdf.cell(0, 5, 'Generated: June 2026', 0, 1, 'C')

# ===================== METHODS =====================
pdf.add_page()
pdf.section_title('Methods')
pdf.body_text(
    'Data Source: Clinical data from 114 thyroid cancer patients (Q3 database) '
    'with genetic testing results for BRAF V600E, TERT promoter, and TP53 mutations. '
    'ATA 2025 PTC 4-level risk classification was applied (Low / Low-Intermediate / '
    'Intermediate-High / High). T3 stage was reclassified as T3b (High risk) per '
    'updated ATA 2025 guidelines.'
)
pdf.body_text(
    "Statistical Analysis: Firth's penalized logistic regression (Jeffreys prior) "
    'was used for all analyses to resolve complete separation caused by rare events. '
    'The penalty equivalently adds 0.5*log|I(beta)| to the log-likelihood, ensuring '
    'finite parameter estimates even when standard MLE diverges. All odds ratios (OR) '
    'and 95% confidence intervals (CI) are reported on the Firth penalized scale. '
    'Hosmer-Lemeshow test assessed calibration. Gene-gene interactions were tested '
    "using Fisher's exact test for co-occurrence and Firth logistic with interaction terms."
)
pdf.body_text(
    'Internal Validation: Bootstrap resampling (1000 iterations) provided bias-corrected '
    'OR estimates and percentile CIs. Calibration was assessed via observed-vs-predicted '
    'plots. Discrimination was quantified by AUC-ROC with bootstrap 95% CI. Leave-one-out '
    'cross-validation (LOOCV) was performed as a sensitivity analysis.'
)
pdf.ln(3)
pdf.set_font('ArialUni', 'B', 9)
pdf.cell(0, 5, 'Study Variables:', 0, 1)
for v in ['Demographics: Age, Gender', 'Tumor: Size, Multifocality, Bilaterality, T stage',
           'Lymph Node: N stage, N1b, Central LN count, Lateral LN status',
           'Genetic: BRAF V600E, BRAF other, TERT (all sites), TP53 (all sites)',
           'Outcome: Recurrence (abnormal US + positive FNA + re-operation + abnormal last FU)']:
    pdf.set_font('ArialUni', '', 8)
    pdf.cell(5, 4, '-', 0, 0)
    pdf.cell(0, 4, v, 0, 1)

# ===================== TABLE 0: ATA RISK =====================
pdf.add_page()
pdf.section_title('Results')
pdf.subsection_title('Table 0: ATA 2025 PTC Risk Stratification')
pdf.body_text('ATA 2025 4-level risk classification (PTC and variants). T3 reclassified as T3b (High).')

risk_data = [
    ['Low', '69', '60.5%', '5', '7.2%'],
    ['Low-Intermediate', '2', '1.8%', '1', '50.0%'],
    ['Intermediate-High', '15', '13.2%', '0', '0.0%'],
    ['High', '28', '24.6%', '1', '3.6%'],
    ['Total', '114', '100%', '7', '6.1%'],
]
pdf.make_table(['Risk Level', 'N', '%', 'Recurrence', 'Rate'],
               risk_data, [40, 15, 15, 20, 20])

# ===================== TABLE 1: BASELINE =====================
pdf.subsection_title('Table 1: Baseline Clinicopathological Characteristics')
# Read results from text file to embed
t1_data = [
    ['Gender (Male)', '26/114 (22.8%)', '1/7 (14.3%)', '25/107 (23.4%)'],
    ['Age >=55 years', '36/114 (31.6%)', '2/7 (28.6%)', '34/107 (31.8%)'],
    ['Age (mean +/- SD)', '46.8 +/- 14.0', '46.1 +/- 14.5', '46.8 +/- 14.1'],
    ['Multifocal', '31/114 (27.2%)', '3/7 (42.9%)', '28/107 (26.2%)'],
    ['Bilateral', '32/114 (28.1%)', '3/7 (42.9%)', '29/107 (27.1%)'],
    ['Size >10mm', '21/114 (18.4%)', '1/7 (14.3%)', '20/107 (18.7%)'],
    ['Size >40mm', '1/114 (0.9%)', '0/7 (0.0%)', '1/107 (0.9%)'],
    ['T stage T3+', '28/114 (24.6%)', '1/7 (14.3%)', '27/107 (25.2%)'],
    ['LN metastasis (N+)', '59/114 (51.8%)', '6/7 (85.7%)', '53/107 (49.5%)'],
    ['Lateral LN (N1b)', '8/114 (7.0%)', '1/7 (14.3%)', '7/107 (6.5%)'],
    ['Central LN >=5', '16/114 (14.0%)', '2/7 (28.6%)', '14/107 (13.1%)'],
    ['Lateral LN (+)', '21/114 (18.4%)', '1/7 (14.3%)', '20/107 (18.7%)'],
    ['BRAF V600E', '45/114 (39.5%)', '3/7 (42.9%)', '42/107 (39.3%)'],
    ['BRAF (any)', '47/114 (41.2%)', '3/7 (42.9%)', '44/107 (41.1%)'],
    ['TERT mutation', '6/114 (5.3%)', '0/7 (0.0%)', '6/107 (5.6%)'],
    ['TP53 mutation', '2/114 (1.8%)', '0/7 (0.0%)', '2/107 (1.9%)'],
    ['Co-mutation', '4/114 (3.5%)', '0/7 (0.0%)', '4/107 (3.7%)'],
]
pdf.make_table(['Variable', 'Overall', 'Recurrence', 'No Recurrence'],
               t1_data, [30, 32, 32, 32])

# ===================== TABLE 2: UNIVARIATE =====================
pdf.add_page()
pdf.subsection_title('Table 2: Univariate Firth Logistic Regression')
pdf.body_text('Firth penalized likelihood estimates. Variables with <2 positive cases excluded.')

uni_data = [
    ['Gender (Male)', '114', '7', '0.747', '0.116-4.798', '0.758'],
    ['Age >=55', '114', '7', '0.968', '0.202-4.639', '0.968'],
    ['Multifocal', '114', '7', '2.170', '0.495-9.514', '0.304'],
    ['Bilateral', '114', '7', '2.070', '0.473-9.060', '0.334'],
    ['Size >1cm', '114', '7', '0.985', '0.151-6.432', '0.987'],
    ['T3+', '114', '7', '0.676', '0.106-4.320', '0.679'],
    ['N+', '114', '7', '4.414', '0.707-27.553', '0.112'],
    ['N1b', '114', '7', '3.092', '0.409-23.381', '0.274'],
    ['Central LN >=5', '114', '7', '2.931', '0.574-14.980', '0.196'],
    ['Lateral LN (+)', '114', '7', '0.985', '0.151-6.432', '0.987'],
    ['BRAF V600E', '114', '7', '1.199', '0.278-5.176', '0.808'],
    ['BRAF (any)', '114', '7', '1.110', '0.257-4.788', '0.889'],
    ['TERT', '114', '7', '1.041', '0.043-25.413', '0.980'],
    ['TP53', '114', '7', '2.813', '0.063-124.823', '0.593'],
    ['Co-mutation', '114', '7', '1.533', '0.054-43.723', '0.803'],
]
pdf.make_table(['Variable', 'N', 'Events', 'OR', '95% CI', 'p'],
               uni_data, [28, 10, 10, 15, 30, 15])

pdf.note_text('Firth penalization guarantees finite estimates even with complete separation. No variable reached statistical significance at alpha=0.05. N+ shows the strongest association (OR=4.41, p=0.112).')

# ===================== TABLE 3: GENE INTERACTION =====================
pdf.subsection_title('Table 3: Gene-Gene Interaction Analysis')
pdf.body_text("Fisher's exact test for co-occurrence between gene pairs.")

interaction_data = [
    ['BRAF V600E vs TERT', '2', '43', '4', '65', '0.756', '1.000'],
    ['BRAF V600E vs TP53', '0', '45', '2', '67', '0.000', '0.518'],
    ['TERT vs TP53', '1', '5', '1', '107', '21.400', '0.103'],
]
pdf.make_table(['Pair', 'Both+', 'G1+/G2-', 'G1-/G2+', 'Both-', 'OR', 'p'],
               interaction_data, [32, 12, 14, 14, 12, 14, 14])

pdf.note_text('No significant co-occurrence or mutual exclusivity was detected. TERT+TP53 OR=21.4 suggests possible co-occurrence but only 1 case (p=0.103).')

# ===================== TABLE 4: MULTIVARIATE =====================
pdf.subsection_title('Table 4: Multivariate Firth Logistic Models')
pdf.body_text('Pre-specified models based on clinical and genetic hypotheses.')

# Model 1
pdf.set_font('ArialUni', 'B', 9)
pdf.cell(0, 6, 'Model 1: Clinical (N=114, events=7)', 0, 1)
m1 = [['N+', '4.979', '0.840-29.503', '0.077'],
      ['T3+', '0.570', '0.094-3.461', '0.541'],
      ['Size >1cm', '0.716', '0.113-4.526', '0.723']]
pdf.make_table(['Variable', 'OR', '95% CI', 'p'], m1, [35, 18, 40, 18])
pdf.set_font('ArialUni', 'I', 7)
pdf.cell(0, 4, '  HL chi2=4.07, df=8, p=0.851', 0, 1)
pdf.ln(2)

# Model 2
pdf.set_font('ArialUni', 'B', 9)
pdf.cell(0, 6, 'Model 2: Genetic (N=114, events=7)', 0, 1)
m2 = [['BRAF V600E', '1.160', '0.276-4.874', '0.839'],
      ['TERT', '1.006', '0.046-22.208', '0.997'],
      ['TP53', '3.079', '0.060-158.766', '0.576']]
pdf.make_table(['Variable', 'OR', '95% CI', 'p'], m2, [35, 18, 40, 18])
pdf.set_font('ArialUni', 'I', 7)
pdf.cell(0, 4, '  HL chi2=8.01, df=8, p=0.432', 0, 1)
pdf.ln(2)

# Model 3
pdf.set_font('ArialUni', 'B', 9)
pdf.cell(0, 6, 'Model 3: Combined (N=114, events=7)', 0, 1)
m3 = [['N+', '4.650', '0.863-25.049', '0.074'],
      ['T3+', '0.545', '0.090-3.297', '0.508'],
      ['Size >1cm', '0.710', '0.119-4.218', '0.706'],
      ['BRAF V600E', '1.111', '0.265-4.661', '0.886'],
      ['TERT', '0.978', '0.043-22.227', '0.989']]
pdf.make_table(['Variable', 'OR', '95% CI', 'p'], m3, [35, 18, 40, 18])
pdf.set_font('ArialUni', 'I', 7)
pdf.cell(0, 4, '  HL chi2=6.37, df=8, p=0.606', 0, 1)
pdf.ln(2)

# Model 4: N+ + Co-mutation
pdf.set_font('ArialUni', 'B', 9)
pdf.cell(0, 6, 'Model 4: N+ + Co-mutation (N=114, events=7)', 0, 1)
m4 = [['N+', '4.117', '0.714-23.732', '0.115'],
      ['Co-mutation', '1.271', '0.050-32.376', '0.886']]
pdf.make_table(['Variable', 'OR', '95% CI', 'p'], m4, [35, 18, 40, 18])
pdf.set_font('ArialUni', 'I', 7)
pdf.cell(0, 4, '  HL chi2=1.57, df=6, p=0.955', 0, 1)
pdf.ln(2)

# Model 5: N+ + BRAF + TERT (primary validation model)
pdf.set_font('ArialUni', 'B', 9)
pdf.cell(0, 6, 'Model 5: N+ + BRAF + TERT (N=114, events=7)*', 0, 1)
m5 = [['N+', '4.418', '0.765-25.514', '0.097'],
      ['BRAF V600E', '0.885', '0.207-3.787', '0.870'],
      ['TERT', '0.943', '0.037-24.014', '0.972']]
pdf.make_table(['Variable', 'OR', '95% CI', 'p'], m5, [35, 18, 40, 18])
pdf.set_font('ArialUni', 'I', 7)
pdf.cell(0, 4, '  HL chi2=4.26, df=8, p=0.833  * Primary validation model', 0, 1)
pdf.ln(2)

# Model 6: N+ + TERT
pdf.set_font('ArialUni', 'B', 9)
pdf.cell(0, 6, 'Model 6: N+ + TERT (N=114, events=7)', 0, 1)
m6 = [['N+', '4.341', '0.783-24.070', '0.094'],
      ['TERT', '0.921', '0.038-22.064', '0.958']]
pdf.make_table(['Variable', 'OR', '95% CI', 'p'], m6, [35, 18, 40, 18])
pdf.set_font('ArialUni', 'I', 7)
pdf.cell(0, 4, '  HL chi2=4.33, df=8, p=0.826', 0, 1)
pdf.ln(2)

# Model 7: N+ + TP53
pdf.set_font('ArialUni', 'B', 9)
pdf.cell(0, 6, 'Model 7: N+ + TP53 (N=114, events=7)', 0, 1)
m7 = [['N+', '5.024', '0.868-29.074', '0.072'],
      ['TP53', '6.998', '0.120-407.038', '0.348']]
pdf.make_table(['Variable', 'OR', '95% CI', 'p'], m7, [35, 18, 40, 18])
pdf.set_font('ArialUni', 'I', 7)
pdf.cell(0, 4, '  HL chi2=0.99, df=5, p=0.963', 0, 1)
pdf.ln(2)

# Model 8: Interaction
pdf.set_font('ArialUni', 'B', 9)
pdf.cell(0, 6, 'Model 8: N+ + BRAF + TERT + BRAFxTERT (N=114, events=7)', 0, 1)
m8 = [['N+', '4.517', '0.751-27.165', '0.100'],
      ['BRAF V600E', '0.866', '0.197-3.795', '0.848'],
      ['TERT', '0.957', '0.030-30.494', '0.980'],
      ['BRAFxTERT', '7.511', '0.030-1882.825', '0.474']]
pdf.make_table(['Variable', 'OR', '95% CI', 'p'], m8, [35, 18, 40, 18])
pdf.set_font('ArialUni', 'I', 7)
pdf.cell(0, 4, '  HL chi2=3.92, df=8, p=0.864', 0, 1)
pdf.ln(2)

pdf.note_text('All models pass Hosmer-Lemeshow test (p>0.05). N+ is the strongest predictor across all models. Model 5 (N+ + BRAF + TERT) is the primary validation model. BRAFxTERT interaction OR=7.51 but CI extremely wide, reflecting insufficient power.')

pdf.ln(1)
pdf.subsection_title('Table 4b: Model Performance Comparison')
pdf.body_text('AUC with bootstrap 95% CI (200 reps), Brier score, LogLoss, and Akaike Information Criterion (AIC) for all 8 models.')
perf_data = [
    ['M1: Clinical (N+ + T3+ + Size)', '114', '7', '3', '0.720', '0.553-0.927', '0.055', '0.212', '56.3'],
    ['M2: Genetic (BRAF+TERT+TP53)', '114', '7', '3', '0.575', '0.169-0.738', '0.058', '0.234', '61.4'],
    ['M3: Combined', '114', '7', '5', '0.733', '0.546-0.916', '0.056', '0.213', '60.7'],
    ['M4: N+ + Co-mutation', '114', '7', '2', '0.770', '0.452-0.847', '0.056', '0.217', '55.5'],
    ['M5: N+ + BRAF + TERT', '114', '7', '3', '0.694', '0.516-0.865', '0.056', '0.216', '57.2'],
    ['M6: N+ + TERT', '114', '7', '2', '0.752', '0.441-0.835', '0.056', '0.216', '55.1'],
    ['M7: N+ + TP53', '114', '7', '2', '0.748', '0.472-0.838', '0.056', '0.218', '55.7'],
    ['M8: N+ + BRAF+TERT+INT', '114', '7', '4', '0.676', '0.503-0.860', '0.057', '0.219', '59.9'],
]
pdf.make_table(['Model', 'N', 'Events', 'k', 'AUC', '95% CI (AUC)', 'Brier', 'LogLoss', 'AIC'],
               perf_data, [40, 8, 8, 6, 10, 22, 10, 10, 10])
pdf.note_text('M4 (N+ + Co-mutation) achieves highest AUC=0.770 but CI very wide (0.452-0.847). '
              'M2 (genetic only) has lowest AUC=0.575. Adding genes to clinical (M3 vs M1) does not improve AUC. '
              'AIC favors simpler models (M6, M7, M4). All CIs overlap substantially, reflecting limited power.')
pdf.ln(1)

# ===================== TABLE 5: BOOTSTRAP =====================
pdf.subsection_title('Table 5: Bootstrap Validation (1000 reps)')
pdf.body_text('Model: N+ + BRAF V600E + TERT. Bias-corrected estimates with percentile CIs.')

boot_data = [
    ['N+', '5.265', '1.269-23.681'],
    ['BRAF V600E', '0.834', '0.162-4.268'],
    ['TERT', '0.998', '0.308-6.784'],
]
pdf.make_table(['Variable', 'OR (bootstrap)', '95% CI'], boot_data, [40, 30, 40])

# Validation summary
pdf.body_text(f'Model performance: AUC-ROC = 0.699 (95% CI: 0.544-0.857), '
              f'LOOCV AUC = 0.433, Brier score = 0.060.')
pdf.ln(2)

# ===================== SENSITIVITY ANALYSIS =====================
pdf.subsection_title('Sensitivity Analysis')
pdf.body_text('Four sensitivity analyses were performed: (1) stricter recurrence definition, '
              '(2) age subgroup, (3) Firth vs MLE comparison, (4) missing data audit.')

pdf.set_font('ArialUni', 'B', 9)
pdf.cell(0, 6, 'SA1: Stricter recurrence (>=2 criteria, 3 events):', 0, 1)
sa1 = [['Gender (Male)', '0.461', '0.022-9.704', '0.618'],
       ['Age >=55', '1.293', '0.161-10.410', '0.809'],
       ['N+', '6.876', '0.338-139.930', '0.210'],
       ['Central LN >=5', '11.207', '1.330-94.399', '0.026'],
       ['BRAF V600E', '11.447', '0.564-232.421', '0.113']]
pdf.make_table(['Variable', 'OR', '95% CI', 'p'], sa1, [35, 18, 40, 18])
pdf.set_font('ArialUni', 'I', 7)
pdf.cell(0, 4, '  Only 3 events with strict criteria. Central LN >=5 significant (p=0.026). BRAF V600E borderline (p=0.113).', 0, 1)
pdf.ln(2)

pdf.set_font('ArialUni', 'B', 9)
pdf.cell(0, 6, 'SA2: Age subgroup (<=55, n=80, events=5 vs >55, n=34, events=2):', 0, 1)
sa2 = [['N+ (<=55)', '2.630', '0.383-18.064', '0.325'],
       ['N+ (>55)', '7.222', '0.295-177.112', '0.226'],
       ['BRAF (<=55)', '0.906', '0.164-4.994', '0.910'],
       ['BRAF (>55)', '2.474', '0.209-29.240', '0.472']]
pdf.make_table(['Variable', 'OR', '95% CI', 'p'], sa2, [35, 18, 40, 18])
pdf.set_font('ArialUni', 'I', 7)
pdf.cell(0, 4, '  N+ effect stronger in >55 group, BRAF appears risk-associated in >55 (protective in <=55) but CIs wide.', 0, 1)
pdf.ln(2)

pdf.set_font('ArialUni', 'B', 9)
pdf.cell(0, 6, 'SA3: Firth vs Standard MLE:', 0, 1)
sa3 = [['Intercept', '0.032', '0.021', '1.565'],
       ['N+', '4.418', '6.369', '0.694'],
       ['BRAF', '0.885', '0.828', '1.069'],
       ['TERT', '0.943', '0.000', '--']]
pdf.make_table(['Parameter', 'Firth OR', 'MLE OR', 'Ratio (Firth/MLE)'], sa3, [35, 25, 25, 25])
pdf.set_font('ArialUni', 'I', 7)
pdf.cell(0, 4, '  MLE produces extreme OR (0.000) for TERT due to complete separation. Firth yields finite estimate (OR=0.943).', 0, 1)
pdf.ln(2)

pdf.note_text('SA4: Missing data — TIRADS pre-op 63.2%, FNA cytology 90.4%, histology 35.1%, LN met counts 43-45%. '
              'These variables were not included in models due to high missingness. Complete case analysis would exclude >60% of patients.')
pdf.ln(2)

# ===================== FIGURES =====================
# Figure 1: Oncoplot
pdf.add_page()
pdf.section_title('Figures')
pdf.subsection_title('Figure 1: Mutation Landscape')
if os.path.exists('H:\\Q3 MNGOC\\oncoplot.png'):
    pdf.image('H:\\Q3 MNGOC\\oncoplot.png', x=10, w=190)
    pdf.ln(2)
pdf.note_text('Mutation status (colored) for 114 patients (columns) across 4 gene categories (rows). '
              'Recurrence and N+ status annotated at top. Bottom bar: mutation frequencies. '
              'BRAF V600E 39.5%, BRAF other 1.8%, TERT 5.3%, TP53 1.8%.')

# Figure 2: Forest plot
pdf.add_page()
pdf.subsection_title('Figure 2: Forest Plot — Univariate Firth Logistic')
if os.path.exists('H:\\Q3 MNGOC\\forest_plot_pub.png'):
    pdf.image('H:\\Q3 MNGOC\\forest_plot_pub.png', x=10, w=190)
pdf.note_text('Odds ratios (squares) with 95% confidence intervals (lines) on log scale. '
              'Dashed vertical line at OR=1.0 (no effect). No variable reaches statistical significance. '
              'N+ shows the strongest trend (OR=4.41, p=0.112).')

# Figure 3: Calibration
pdf.add_page()
pdf.subsection_title('Figure 3: Calibration Curve')
if os.path.exists('H:\\Q3 MNGOC\\calibration_plot.png'):
    pdf.image('H:\\Q3 MNGOC\\calibration_plot.png', x=25, w=160)
pdf.note_text('Calibration of N+ + BRAF + TERT model. Points near diagonal indicate good calibration. '
              'Intercept and slope quantify calibration accuracy (ideal: intercept=0, slope=1).')

# Figure 4: ROC
pdf.subsection_title('Figure 4: ROC Curve')
if os.path.exists('H:\\Q3 MNGOC\\roc_curve.png'):
    pdf.image('H:\\Q3 MNGOC\\roc_curve.png', x=25, w=160)
pdf.note_text('AUC = 0.699 (95% CI: 0.544-0.857), indicating limited discriminative ability.')

# Figure 5: Nomogram
pdf.add_page()
pdf.subsection_title('Figure 5: Nomogram for Recurrence Risk')
if os.path.exists('H:\\Q3 MNGOC\\nomogram.png'):
    pdf.image('H:\\Q3 MNGOC\\nomogram.png', x=10, w=190)
pdf.note_text('Nomogram based on N+ + BRAF + TERT model (corrected: negative coefficients yield fewer points). '
              'Plot points for each variable, sum total points, read corresponding risk on bottom scale. '
              'Example: N+ (100 pts) + BRAF V600E (0 pts) + TERT WT (0 pts) = 100 pts -> risk ~11.1%.')

# ===================== SUMMARY =====================
pdf.add_page()
pdf.section_title('Summary and Conclusions')

summary = [
    ('Recurrence Rate', f'The overall recurrence rate was 6.1% (7/114 patients), consistent with the expected low event rate in differentiated thyroid cancer.'),
    ('N+ Status', 'N+ (lymph node metastasis) was the strongest and most consistent predictor across all models. Univariate OR=4.41 (p=0.112), multivariate OR ranged from 4.12 to 4.98. While not reaching statistical significance at alpha=0.05, the consistent magnitude and direction support clinical relevance. The wide CI reflects limited statistical power.'),
    ('BRAF V600E', 'BRAF V600E (prevalence 39.5%) showed no association with recurrence in any analysis: univariate OR=1.20 (p=0.81), multivariate OR~1.1, all subgroup analyses OR~1. The BRAFxTERT interaction term was not estimable (only 2 co-mutation cases).'),
    ('TERT and TP53', 'TERT (5.3%) and TP53 (1.8%) mutations were too rare for reliable inference. Firth provides finite but imprecise estimates. TP53 showed a potential risk signal (univariate OR=2.81, multivariate OR=7.00 with N+), but CIs were extremely wide.'),
    ('Gene-Gene Interaction', 'No significant co-occurrence or mutual exclusivity was detected. BRAF+TERT co-occurrence (2 cases), BRAF+TP53 (0 cases), TERT+TP53 (1 case). The TERT+TP53 OR=21.4 trended toward co-occurrence but was not significant (p=0.103).'),
    ('Model Performance', 'The best model (N+ + BRAF + TERT) achieved AUC=0.699 (95% CI: 0.544-0.857), indicating modest discrimination. LOOCV AUC=0.433, suggesting instability due to low event count. All models passed Hosmer-Lemeshow calibration tests.'),
    ('Sensitivity Analyses', 'Stricter recurrence definition (>=2 criteria) reduced events to 3. Central LN >=5 became significant (OR=11.2, p=0.026). BRAF V600E showed borderline association (OR=11.4, p=0.113). Age subgroup analysis suggested N+ effect is stronger in >55 group. Firth vs MLE comparison confirmed that MLE produces extreme estimates (TERT OR=0.000) due to complete separation, justifying Firth as the primary method.'),
    ('Limitations', 'Primary limitation is the low event count (n=7), yielding events-per-variable (EPV) < 3 in most models. Firth logistic regression resolves complete separation but cannot increase effective sample size. Results require validation in larger, multi-center cohorts. The absence of histologic subtype data (tall cell, hobnail variants) and extranodal extension limits assessment of BRAF-dependent risk pathways. High missingness in TIRADS (63%), FNA (90%), and histology (35%) precluded inclusion of these key variables.'),
]

for title, content in summary:
    pdf.set_font('ArialUni', 'B', 10)
    pdf.set_text_color(0, 51, 102)
    pdf.cell(5, 6, '-', 0, 0)
    pdf.cell(0, 6, title, 0, 1)
    pdf.set_text_color(0, 0, 0)
    pdf.set_font('ArialUni', '', 8.5)
    pdf.set_x(15)
    pdf.multi_cell(0, 4.5, content)
    pdf.ln(2)

pdf.ln(5)
pdf.set_draw_color(0, 51, 102)
pdf.set_line_width(0.5)
pdf.line(10, pdf.get_y(), 200, pdf.get_y())
pdf.ln(3)
pdf.set_font('ArialUni', 'I', 8)
pdf.set_text_color(100, 100, 100)
pdf.multi_cell(0, 4, 'Analysis performed with Python (FirthLogit custom implementation). '
               'Firth penalized likelihood = standard log-likelihood + 0.5*log|I(beta)| '
               '(Jeffreys prior penalty). Bootstrap = 1000 iterations, percentile CI. '
               'ROC and calibration using manual implementations.')

# ===================== APPENDIX: OR CALCULATION =====================
pdf.add_page()
pdf.section_title('Appendix: Calculation of Odds Ratios (OR)')
pdf.ln(1)

# A. Univariate
pdf.subsection_title('A. Table 2 — Univariate Analysis (Crude OR)')
pdf.body_text(
    'Each variable is analyzed independently using Firth penalized logistic regression:\n'
    '  logit(P) = beta0 + beta_i * X_i\n'
    'The odds ratio is OR_i = exp(beta_i).\n'
    '95% confidence interval: exp(beta_i ± 1.96 * SE_beta).'
)
pdf.body_text(
    'Example — BRAF V600E:\n'
    '  beta = 0.181, SE = 0.747\n'
    '  OR = exp(0.181) = 1.199 (95% CI: 0.278 - 5.176, p = 0.81)\n'
    'Interpretation: "Patients with BRAF V600E have 1.2 times the odds of recurrence '
    'compared to wild-type patients, but this difference is not statistically significant (p = 0.81)."\n'
    'This is a crude (unadjusted) OR — it reflects the total association between BRAF and recurrence, '
    'including indirect effects through other correlated variables (e.g., N+ status).'
)

# B. Multivariate
pdf.subsection_title('B. Table 4 — Multivariate Analysis (Adjusted OR)')
pdf.body_text(
    'Multiple predictors in a single Firth model:\n'
    '  logit(P) = beta0 + beta1*X1 + beta2*X2 + ... + betak*Xk\n'
    'OR_i = exp(beta_i) represents the independent effect of X_i holding all other variables constant.'
)
pdf.body_text(
    'Example — Model 5 (N+ + BRAF + TERT):\n'
    '  logit(P) = -3.442 + 1.486*N+ + (-0.122)*BRAF + (-0.059)*TERT\n\n'
    '  OR(N+) = exp(1.486) = 4.418\n'
    '    Holding BRAF and TERT constant, N+ increases recurrence odds 4.4-fold.\n\n'
    '  OR(BRAF) = exp(-0.122) = 0.885\n'
    '    Holding N+ and TERT constant, BRAF V600E decreases odds by 11.5%.\n'
    '    Compare with crude OR = 1.20 (Table 2). Why the change? BRAF V600E is correlated with N+.\n'
    '    Crude OR = total effect (BRAF + indirect via N+). Adjusted OR = direct effect of BRAF alone.\n'
    '    This is an example of confounding adjustment — and mild Simpson\'s paradox.\n\n'
    '  OR(TERT) = exp(-0.059) = 0.943\n'
    '    Holding N+ and BRAF constant, TERT has minimal effect.'
)
pdf.body_text(
    'Why OR changes between univariate and multivariate:\n'
    '  Crude OR (1.20) = BRAF effect + "BRAF correlates with N+ so part of N+ effect leaks into BRAF"\n'
    '  Adjusted OR (0.89) = BRAF effect alone, after removing the N+ contribution\n'
    '  This happens because BRAF+ patients have a different N+ rate than BRAF- patients.\n'
    '  The adjusted OR is more accurate for causal inference.'
)

# C. Fisher
pdf.subsection_title('C. Table 3 — Gene-Gene Interaction (Fisher Exact Test)')
pdf.body_text(
    'For binary gene pairs, a 2×2 contingency table is constructed:\n\n'
    '               Gene2+    Gene2-\n'
    '  Gene1+         a         b\n'
    '  Gene1-         c         d\n\n'
    '  OR = (a × d) / (b × c)\n'
    '  Fisher\'s exact test is used (no normal approximation, valid for small counts).'
)
pdf.body_text(
    'Example — BRAF V600E vs TERT:\n'
    '               TERT+    TERT-\n'
    '  BRAF+          2       43\n'
    '  BRAF-          4       65\n\n'
    '  OR = (2 × 65) / (43 × 4) = 130 / 172 = 0.756\n'
    '  Interpretation: OR < 1 suggests mutual exclusivity (BRAF+ and TERT+ tend not to co-occur).\n'
    '  However, p = 1.000 — no significant association.\n\n'
    '  TERT vs TP53: OR = 21.4 (p = 0.103) — hints at co-occurrence but only 1 co-positive case.'
)

# D. Firth vs MLE
pdf.subsection_title('D. Comparison: Firth vs Standard MLE')
pdf.body_text(
    'Standard MLE logistic regression maximizes the log-likelihood:\n'
    '  LL(beta) = Sigma[ yi * log(pi) + (1-yi) * log(1-pi) ]\n\n'
    'Firth adds Jeffreys prior penalty (0.5 * log|I(beta)|):\n'
    '  LL*(beta) = LL(beta) + 0.5 * log(det(I(beta)))\n'
    '  where I(beta) = Fisher information matrix\n\n'
    'This penalty shrinks beta estimates toward zero (OR toward 1.0), guaranteeing finite '
    'estimates even with complete separation (a variable perfectly predicts the outcome).'
)
pdf.body_text(
    'Results from Sensitivity Analysis (SA3) — Model 5:\n\n'
    '  Variable     Firth OR     MLE OR      Ratio (Firth/MLE)\n'
    '  Intercept    0.032        0.021       1.565\n'
    '  N+           4.418        6.369       0.694\n'
    '  BRAF V600E   0.885        0.828       1.069\n'
    '  TERT         0.943        0.000       --\n\n'
    'Key observation: MLE gives TERT OR = 0.000 (effectively infinite — complete separation '
    'because 0/6 TERT+ patients had recurrence). Firth gives TERT OR = 0.943 (finite and interpretable).\n\n'
    'Firth is the appropriate method for this study due to:\n'
    '  1. Rare mutations: TERT (6/114 = 5.3%), TP53 (2/114 = 1.8%)\n'
    '  2. Low event count: 7/114 (6.1%) — leads to quasi-complete separation\n'
    '  3. Firth guarantees convergence where MLE fails'
)

pdf.subsection_title('E. Table A3 — Diagnostic Metrics for Model 5 (N+ + BRAF + TERT)')
pdf.body_text(
    'Sensitivity, specificity, positive predictive value (PPV), negative predictive value (NPV), '
    'accuracy, and Youden index (sensitivity + specificity - 1) at various probability thresholds. '
    'Optimal threshold maximizes Youden index.'
)
diag_data = [
    ['0.020', '7', '107', '0', '0', '1.000', '0.000', '0.061', '0.000', '0.061', '0.0000'],
    ['0.050*', '6', '53', '1', '54', '0.857', '0.505', '0.102', '0.982', '0.526', '0.3618'],
    ['0.100', '6', '53', '1', '54', '0.857', '0.505', '0.102', '0.982', '0.526', '0.3618'],
    ['0.120', '3', '25', '4', '82', '0.429', '0.766', '0.107', '0.953', '0.746', '0.1949'],
    ['0.150', '0', '0', '7', '107', '0.000', '1.000', '0.000', '0.939', '0.939', '0.0000'],
]
pdf.make_table(['Threshold', 'TP', 'FP', 'FN', 'TN', 'Sens', 'Spec', 'PPV', 'NPV', 'Acc', 'Youden'],
               diag_data, [14, 8, 8, 8, 8, 14, 14, 14, 14, 14, 14])
pdf.note_text(
    '* Optimal threshold (Youden index = 0.36). Sensitivity 85.7% (detects 6/7 recurrences) '
    'with specificity 50.5%. PPV low (10.2%) due to low event prevalence (6.1%). '
    'NPV 98.2% ensures high confidence when predicting no recurrence.'
)

# Summary footnotes note
pdf.ln(2)
pdf.set_draw_color(0, 51, 102)
pdf.set_line_width(0.5)
pdf.line(10, pdf.get_y(), 200, pdf.get_y())
pdf.ln(2)
pdf.set_font('ArialUni', 'I', 7.5)
pdf.set_text_color(100, 100, 100)
pdf.multi_cell(0, 4, 'Abbreviations: OR = Odds Ratio; CI = Confidence Interval; Firth = Firth penalized logistic regression; '
               'MLE = Maximum Likelihood Estimation; SE = Standard Error. All calculations performed in Python 3.14 '
                'using custom FirthLogit implementation with scipy.optimize and numpy.\n'
                'AUC = Area Under the Receiver Operating Characteristic Curve; '
                'AIC = Akaike Information Criterion; '
                'TP = True Positive; FP = False Positive; FN = False Negative; TN = True Negative; '
                'Sens = Sensitivity; Spec = Specificity; PPV = Positive Predictive Value; '
                'NPV = Negative Predictive Value; Acc = Accuracy.')

output_path = OUTPUT
if os.path.exists(output_path):
    try:
        os.remove(output_path)
    except:
        pass
pdf.output(output_path)
print(f'PDF saved to: {output_path}')
print(f'File size: {os.path.getsize(output_path)} bytes')
