import os
from fpdf import FPDF

class PDF(FPDF):
    def header(self):
        pass

    def footer(self):
        pass

    def section_title(self, title):
        self.set_font('ArialUni', 'B', 13)
        self.set_text_color(0, 51, 102)
        self.cell(0, 8, title, 0, 1, 'L')
        self.set_text_color(0, 0, 0)
        self.ln(1)

    def subsection_title(self, title):
        self.set_font('ArialUni', 'B', 11)
        self.set_text_color(0, 76, 153)
        self.cell(0, 7, title, 0, 1, 'L')
        self.set_text_color(0, 0, 0)
        self.ln(1)

    def make_table(self, headers, data, col_widths=None):
        if col_widths is None:
            w = (self.w - 20) / len(headers)
            col_widths = [w] * len(headers)
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
            self.set_fill_color(235, 245, 255) if fill else self.set_fill_color(255, 255, 255)
            for i, cell in enumerate(row):
                self.cell(col_widths[i], 7, str(cell), 1, 0, 'C', True)
            self.ln()
            fill = not fill

    def add_body_text(self, text):
        self.set_font('ArialUni', '', 9)
        self.multi_cell(0, 5, text)
        self.ln(1)

    def add_note_text(self, text):
        self.set_font('ArialUni', 'I', 8)
        self.set_text_color(100, 100, 100)
        self.multi_cell(0, 4.5, text)
        self.set_text_color(0, 0, 0)
        self.ln(1)

pdf = PDF()
pdf.alias_nb_pages()
pdf.add_font('ArialUni', '', 'C:\\Windows\\Fonts\\arial.ttf')
pdf.add_font('ArialUni', 'B', 'C:\\Windows\\Fonts\\arialbd.ttf')
pdf.add_font('ArialUni', 'I', 'C:\\Windows\\Fonts\\ariali.ttf')
pdf.add_font('ArialUni', 'BI', 'C:\\Windows\\Fonts\\arialbi.ttf')

# Title page
pdf.add_page()
pdf.ln(40)
pdf.set_font('ArialUni', 'B', 20)
pdf.set_text_color(0, 51, 102)
pdf.cell(0, 12, 'ANALYSIS REPORT', 0, 1, 'C')
pdf.ln(4)
pdf.cell(0, 12, 'RECURRENCE RISK PREDICTION FOR THYROID CANCER', 0, 1, 'C')
pdf.cell(0, 12, 'ACCORDING TO ATA 2025 GUIDELINES', 0, 1, 'C')
pdf.ln(20)
pdf.set_font('ArialUni', '', 12)
pdf.set_text_color(0, 0, 0)
pdf.cell(0, 8, 'Firth Penalized Logistic Regression', 0, 1, 'C')
pdf.cell(0, 8, 'Retrospective Study, 114 patients', 0, 1, 'C')
pdf.cell(0, 8, 'Study period: 2024', 0, 1, 'C')
pdf.ln(40)
pdf.set_font('ArialUni', 'I', 9)
pdf.set_text_color(100, 100, 100)
pdf.cell(0, 6, 'Firth method resolves complete separation in rare mutation groups', 0, 1, 'C')

# --- SECTION 1 ---
pdf.add_page()
pdf.section_title('1. Study Sample Characteristics')
pdf.add_body_text('Total patients: 114. Recurrence events: 7 (6.1%).')

pdf.subsection_title('1.1. ATA 2025 PTC Risk Stratification (4-level)')
pdf.set_font('ArialUni', 'I', 8)
pdf.set_text_color(100, 100, 100)
pdf.cell(0, 5, '(ATA 2025 PTC 4-level: Low / Low-Intermediate / Intermediate-High / High)', 0, 1, 'L')
pdf.set_text_color(0, 0, 0)
pdf.ln(1)
ata_data = [
    ['Low', '69', '60.5%', '5/69 (7.2%)'],
    ['Low-Intermediate', '2', '1.8%', '1/2 (50.0%)'],
    ['Intermediate-High', '15', '13.2%', '0/15 (0.0%)'],
    ['High', '28', '24.6%', '1/28 (3.6%)'],
    ['Total', '114', '100%', '7/114 (6.1%)'],
]
pdf.make_table(['Risk Level', 'N', '%', 'Recurrence'], ata_data, [45, 20, 20, 30])
pdf.ln(2)
pdf.add_note_text('Note: Recurrence rates do not increase monotonically across levels due to missing histopathological criteria (aggressive histology, margins, ETE, vascular invasion) unavailable in this dataset.')

pdf.subsection_title('1.2. Descriptive Statistics')
desc_data = [
    ['Multifocal (>1 nodule)', '31', '27.2%'],
    ['Bilateral involvement', '32', '28.1%'],
    ['Tumor size >1cm', '21', '18.4%'],
    ['Tumor size >4cm', '1', '0.9%'],
    ['T stage T3+', '28', '24.6%'],
    ['LN metastasis (N+)', '59', '51.8%'],
    ['Lateral neck meta. (N1b)', '8', '7.0%'],
    ['Central LN >=5', '16', '14.0%'],
    ['Lateral neck LN (+)', '21', '18.4%'],
    ['BRAF V600E', '45', '39.5%'],
    ['Other BRAF', '2', '1.8%'],
    ['TERT mutation', '6', '5.3%'],
    ['TP53 mutation', '2', '1.8%'],
    ['Any mutation', '51', '44.7%'],
    ['Co-mutation (>=2 genes)', '4', '3.5%'],
    ['BRAF + TERT', '3', '2.6%'],
    ['BRAF + TP53', '0', '0.0%'],
]
pdf.make_table(['Variable', 'N (n=114)', '%'], desc_data, [50, 30, 25])

# --- SECTION 2: UNIVARIATE (FIRTH) ---
pdf.add_page()
pdf.section_title('2. Univariate Firth Logistic Regression')
pdf.add_body_text('Firth penalized likelihood logistic regression resolves complete separation in mutation groups with zero events, producing finite odds ratios and confidence intervals for all variables. Variables with <2 positive cases (tumor size >4cm, BRAF+TP53) are excluded.')

uni_headers = ['Variable', 'N', 'Events', 'OR', '95% CI', 'p-value']
uni_data = [
    ['Multifocal (>1 nodule)', '114', '7', '2.170', '0.495 - 9.514', '0.304'],
    ['Bilateral involvement', '114', '7', '2.070', '0.473 - 9.060', '0.334'],
    ['Tumor size >1cm', '114', '7', '0.985', '0.151 - 6.432', '0.987'],
    ['Tumor size >4cm', '114', '7', '--', '--', 'excluded'],
    ['T stage T3+', '114', '7', '0.676', '0.106 - 4.320', '0.679'],
    ['LN metastasis (N+)', '114', '7', '4.414', '0.707 - 27.553', '0.112'],
    ['Lateral neck (N1b)', '114', '7', '3.092', '0.409 - 23.381', '0.274'],
    ['Central LN >=5', '114', '7', '2.931', '0.574 - 14.980', '0.196'],
    ['Lateral neck LN (+)', '114', '7', '0.985', '0.151 - 6.432', '0.987'],
    ['BRAF V600E', '114', '7', '1.199', '0.278 - 5.176', '0.808'],
    ['BRAF (any)', '114', '7', '1.110', '0.257 - 4.788', '0.889'],
    ['TERT mutation', '114', '7', '1.041', '0.043 - 25.413', '0.980'],
    ['TP53 mutation', '114', '7', '2.813', '0.063 - 124.823', '0.593'],
    ['Any mutation', '114', '7', '0.954', '0.222 - 4.110', '0.950'],
    ['Co-mutation (>=2)', '114', '7', '1.533', '0.054 - 43.723', '0.803'],
    ['BRAF + TERT', '114', '7', '1.990', '0.060 - 66.025', '0.700'],
    ['BRAF + TP53', '114', '7', '--', '--', 'excluded'],
]
pdf.make_table(uni_headers, uni_data, [38, 10, 10, 18, 35, 18])

pdf.ln(3)
pdf.add_note_text('Note: Firth penalization uses Jeffrey\'s prior penalty (0.5*log|I(beta)|), guaranteeing finite estimates even with complete separation. TERT (OR=1.04, CI=0.04-25.41) and TP53 (OR=2.81, CI=0.06-124.82) now have finite and interpretable estimates, though wide CIs reflect low event count.')

pdf.subsection_title('2.1. Key Observations')
notes = [
    'N+ is the strongest predictor (OR=4.41; p=0.112), approaching significance.',
    'Central LN >=5 shows moderate risk (OR=2.93), consistent with ATA guidelines.',
    'TP53 Firth OR=2.81 suggests possible risk elevation but CI is very wide.',
    'TERT OR=1.04 indicates no association with recurrence in this sample.',
    'BRAF V600E OR=1.20, not significant (p=0.81).',
]
for n in notes:
    pdf.cell(5, 4.5, '-', 0, 0)
    pdf.set_font('ArialUni', '', 8.5)
    pdf.multi_cell(0, 4.5, '  ' + n)
    pdf.ln(1)

# --- SECTION 3: MULTIVARIATE (FIRTH) ---
pdf.add_page()
pdf.section_title('3. Multivariate Firth Logistic Regression')
pdf.add_body_text('Four models fitted with Firth penalized logistic regression.')

m_headers = ['Variable', 'OR', '95% CI', 'p-value']

pdf.subsection_title('3.1. Model 1 - Clinical (N=114, events=7)')
m1_data = [
    ['Size >1cm', '0.716', '0.113 - 4.526', '0.723'],
    ['T3+', '0.570', '0.094 - 3.461', '0.541'],
    ['N+', '4.979', '0.840 - 29.503', '0.077'],
]
pdf.make_table(m_headers, m1_data, [30, 18, 38, 18])
pdf.ln(2)
pdf.set_font('ArialUni', 'I', 8)
pdf.cell(0, 4, 'Hosmer-Lemeshow: chi2=4.07, df=8, p=0.851', 0, 1, 'L')
pdf.set_font('ArialUni', '', 9)
pdf.ln(2)

pdf.subsection_title('3.2. Model 2 - Genetic (N=114, events=7)')
m2_data = [
    ['BRAF V600E', '1.180', '0.282 - 4.934', '0.820'],
    ['TERT', '1.040', '0.045 - 23.937', '0.980'],
]
pdf.make_table(m_headers, m2_data, [30, 18, 38, 18])
pdf.ln(2)
pdf.set_font('ArialUni', 'I', 8)
pdf.cell(0, 4, 'Hosmer-Lemeshow: chi2=7.86, df=8, p=0.447', 0, 1, 'L')
pdf.set_font('ArialUni', '', 9)
pdf.ln(2)

pdf.subsection_title('3.3. Model 3 - Combined (N=114, events=7)')
m3_data = [
    ['Size >1cm', '0.710', '0.119 - 4.218', '0.706'],
    ['T3+', '0.545', '0.090 - 3.297', '0.508'],
    ['N+', '4.650', '0.863 - 25.049', '0.074'],
    ['BRAF V600E', '1.111', '0.265 - 4.661', '0.886'],
    ['TERT', '0.978', '0.043 - 22.227', '0.989'],
]
pdf.make_table(m_headers, m3_data, [30, 18, 38, 18])
pdf.ln(2)
pdf.set_font('ArialUni', 'I', 8)
pdf.cell(0, 4, 'Hosmer-Lemeshow: chi2=6.37, df=8, p=0.606', 0, 1, 'L')
pdf.set_font('ArialUni', '', 9)
pdf.ln(2)

pdf.subsection_title('3.4. Model 4 - N+ and Co-mutation (N=114, events=7)')
m4_data = [
    ['N+', '4.118', '0.710 - 23.891', '0.115'],
    ['Co-mutation', '2.163', '0.065 - 71.702', '0.666'],
]
pdf.make_table(m_headers, m4_data, [30, 18, 38, 18])
pdf.ln(2)
pdf.set_font('ArialUni', 'I', 8)
pdf.cell(0, 4, 'Hosmer-Lemeshow: chi2=4.11, df=8, p=0.848', 0, 1, 'L')
pdf.set_font('ArialUni', '', 9)

pdf.subsection_title('3.5. Model 5 - TERT + N+ (N=114, events=7)')
m5_data = [
    ['TERT', '1.008', '0.040 - 25.601', '0.996'],
    ['N+', '4.343', '0.734 - 25.703', '0.106'],
]
pdf.make_table(m_headers, m5_data, [30, 18, 38, 18])
pdf.ln(2)
pdf.set_font('ArialUni', 'I', 8)
pdf.cell(0, 4, 'Hosmer-Lemeshow: chi2=7.08, df=8, p=0.528', 0, 1, 'L')
pdf.set_font('ArialUni', '', 9)
pdf.ln(2)

pdf.subsection_title('3.6. Model 6 - TP53 + N+ (N=114, events=7)')
m6_data = [
    ['TP53', '7.000', '0.120 - 407.327', '0.348'],
    ['N+', '4.252', '0.680 - 26.575', '0.122'],
]
pdf.make_table(m_headers, m6_data, [30, 18, 38, 18])
pdf.ln(2)
pdf.set_font('ArialUni', 'I', 8)
pdf.cell(0, 4, 'Hosmer-Lemeshow: chi2=4.62, df=8, p=0.797', 0, 1, 'L')
pdf.set_font('ArialUni', '', 9)
pdf.ln(2)

pdf.subsection_title('3.7. Model 7 - TERT + TP53 Mutual (N=114, events=7)')
m7_data = [
    ['TERT', '1.005', '0.043 - 23.572', '0.997'],
    ['TP53', '3.080', '0.062 - 152.409', '0.572'],
]
pdf.make_table(m_headers, m7_data, [30, 18, 38, 18])
pdf.ln(2)
pdf.set_font('ArialUni', 'I', 8)
pdf.cell(0, 4, 'Hosmer-Lemeshow: chi2=7.94, df=8, p=0.439', 0, 1, 'L')
pdf.set_font('ArialUni', '', 9)

# --- SECTION 3.5 (old) → now 3.8: BRAF SUBGROUP ANALYSIS ---
pdf.add_page()
pdf.section_title('3.8. Detailed Analysis of BRAF V600E')

pdf.add_body_text('BRAF V600E prevalence is 39.5% (45/114). Although known in the literature as a risk factor, BRAF V600E showed no significant association with recurrence in this cohort. Subgroup analyses were performed to explore potential effect modification.')

# Forest plot
pdf.subsection_title('3.8.1. Forest Plot: BRAF OR by Subgroup')
pdf.image('H:\\Q3 MNGOC\\braf_forest_plot.png', x=20, w=170)
pdf.ln(2)
pdf.add_note_text('Firth OR for BRAF V600E across clinical subgroups. All subgroups show OR close to 1.00, indicating no significant risk elevation regardless of N status, multifocality, or T stage. N0 and T3+ subgroups had only 1 event each, precluding estimation.')

# Stacked bar
pdf.subsection_title('3.8.2. Recurrence by BRAF and N Status')
pdf.image('H:\\Q3 MNGOC\\braf_stacked_bar.png', x=15, w=180)
pdf.ln(2)
pdf.add_note_text('Left: Recurrence rates are driven by N status, not by BRAF. In N+ patients, BRAF+ (10.7%) vs BRAF- (9.7%) are nearly identical. In N0 patients, BRAF+ (0%) vs BRAF- (2.6%). Right: BRAF distribution across N0 and N+ is balanced (17 vs 28).')

# Recurrence bar
pdf.subsection_title('3.8.3. Overall Recurrence by BRAF Status')
pdf.image('H:\\Q3 MNGOC\\braf_recurrence_bar.png', x=35, w=140)
pdf.ln(2)
pdf.add_note_text('Overall: BRAF- 4/69 (5.8%) vs BRAF+ 3/45 (6.7%). Firth OR=1.20, 95% CI 0.28-5.18, p=0.81. No statistically significant difference.')

# Heatmap
pdf.subsection_title('3.8.4. Mutation Co-occurrence Matrix')
pdf.image('H:\\Q3 MNGOC\\mutation_heatmap.png', x=35, w=140)
pdf.ln(2)
pdf.add_note_text('Co-mutation counts among 114 patients. BRAF+TERT: 2 cases (1.8%). BRAF+TP53: 0 cases. TERT+TP53: 1 case (0.9%). Very low co-occurrence rates limit the analysis of synergistic effects.')

pdf.subsection_title('3.8.5. Interpretation')
pdf.add_body_text('BRAF V600E does not independently predict recurrence in this dataset. Subgroup analyses confirm that N+ status is the dominant risk factor. Possible explanations: (1) inadequate power (only 7 events), (2) BRAF effect may depend on histologic subtype (tall cell, hobnail) unavailable in this dataset, (3) BRAF requires TERT co-mutation for risk elevation (only 2 patients with both). Larger studies incorporating histopathology are needed.')

# --- SECTION 4: SUMMARY ---
pdf.add_page()
pdf.section_title('4. Key Findings and Conclusions')

summary_points = [
    ('Firth Method', 'Firth penalized logistic regression was used to address complete separation caused by rare mutation groups (TERT n=6, TP53 n=2, co-mutation n=4). All estimates now have finite odds ratios and confidence intervals.'),
    ('Recurrence Rate', 'The overall recurrence rate was 6.1% (7/114 patients).'),
    ('N+ as Predictor', 'N+ is the strongest and most consistent predictor: univariate OR=4.41 (p=0.112); Model 1 OR=4.98 (p=0.077). The consistent trend across all models supports its clinical relevance.'),
    ('TERT and TP53', 'Firth provides finite OR for TERT (1.04) and TP53 (2.81), but wide CIs indicate insufficient data. TP53 OR=2.81 hints at risk elevation, while TERT shows no association. Models adjusting for N+ (Model 5: TERT OR=1.01, p=0.996; Model 6: TP53 OR=7.00, p=0.348) confirm these mutations do not add predictive value beyond N+ status.'),
    ('BRAF V600E', 'BRAF V600E (39.5% prevalence) is not associated with recurrence in this cohort (OR=1.20; p=0.81). Subgroup analysis confirmed no effect modification by N status, multifocality, or T stage. This may reflect limited power (7 events) or missing histopathologic data (tall cell variant, ETE) rather than true absence of effect. Only 2 patients had BRAF+TERT co-mutation, precluding analysis of this synergistic pair.'),
    ('Model Fit', 'All 7 multivariate models pass Hosmer-Lemeshow goodness-of-fit test (p>0.05), indicating adequate model specification.'),
    ('Limitations', 'The limited number of recurrence events (n=7) substantially reduces statistical power. Firth resolves separation but cannot increase effective sample size. Confidence intervals around rare mutation ORs remain very wide. Larger studies with Cox regression are recommended.'),
]

for title, content in summary_points:
    pdf.set_font('ArialUni', 'B', 9)
    pdf.cell(5, 5, '-', 0, 0)
    pdf.cell(0, 5, ' ' + title, 0, 1)
    pdf.set_font('ArialUni', '', 8.5)
    pdf.set_x(15)
    pdf.multi_cell(0, 4.5, '  ' + content)
    pdf.ln(2)

pdf.ln(5)
pdf.set_draw_color(0, 51, 102)
pdf.set_line_width(0.5)
pdf.line(10, pdf.get_y(), 200, pdf.get_y())
pdf.ln(3)
pdf.set_font('ArialUni', 'I', 8)
pdf.set_text_color(100, 100, 100)
pdf.multi_cell(0, 4.5, 'Firth logistic regression (Jeffreys prior). Data encoded from Q3 source file and genetic results file. Encoded data file: DATA_MAHOA_ATA2025.xlsx | Analysis file: RESULTS_ATA2025_FIRTH.txt')

output_path = 'H:\\Q3 MNGOC\\BAOCAO_ATA2025_EN_FIRTH.pdf'
if os.path.exists(output_path):
    output_path = 'H:\\Q3 MNGOC\\BAOCAO_ATA2025_EN_FIRTH_v2.pdf'
pdf.output(output_path)
print(f'PDF saved to: {output_path}')
print(f'File size: {os.path.getsize(output_path)} bytes')
