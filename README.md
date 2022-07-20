# Risk_Score

Calculate the Genetic Risk Score (GRS) for a list of subjects based on predefined risk weights.

Reads a (vcf) file with population/sample genotype information and outputs a genetic risk score based on predefined weights. Several different algorithms are available including a generic aggregate score and several specific, published algorithms.

Predefined weights are included based on the previously published algorithms, but the user can override these if desired.

Should be possible to install with:
'''
pip install git+https://github.com/carsten0202/risk_score
'''
