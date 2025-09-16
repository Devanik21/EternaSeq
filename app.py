import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import io
import re
from collections import Counter, defaultdict
import json
from datetime import datetime
from typing import Dict, List, Tuple, Optional
import hashlib
import google.generativeai as genai
import base64
import py3Dmol
import streamlit.components.v1 as components

# Configure page
st.set_page_config(
    page_title="🧬 EternaSeq - Revolutionary DNA Analyzer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS
st.markdown("""
<style>
    .main-header {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        padding: 2rem;
        border-radius: 20px;
        text-align: center;
        margin-bottom: 2.5rem;
        color: white;
        box-shadow: 0 8px 32px rgba(102, 126, 234, 0.3);
    }
    .analysis-card {
        background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%);
        padding: 1.5rem;
        border-radius: 10px;
        color: white;
        margin: 1rem 0;
        box-shadow: 0 4px 15px rgba(0,0,0,0.1);
    }
    .metric-card {
        background: rgba(255,255,255,0.1);
        backdrop-filter: blur(10px);
        border: 1px solid rgba(255,255,255,0.2);
        padding: 1.5rem;
        border-radius: 10px;
        text-align: center;
        margin: 0.5rem 0;
    }
    .dna-sequence {
        font-family: 'Courier New', monospace;
        background: #0d1117;
        color: #58a6ff;
        padding: 1rem;
        border-radius: 8px;
        font-size: 14px;
        line-height: 1.4;
        overflow-x: auto;
        border: 1px solid #30363d;
    }
    .gene-highlight {
        background: #1f6feb;
        color: white;
        padding: 0.2rem 0.5rem;
        border-radius: 4px;
        font-weight: bold;
    }
    .protein-structure {
        background: #238636;
        color: white;
        padding: 0.3rem 0.6rem;
        border-radius: 5px;
        margin: 0.2rem;
        display: inline-block;
    }
</style>
""", unsafe_allow_html=True)

class DNAAnalyzer:
    """Revolutionary DNA sequence analysis engine"""
    
    def __init__(self):
        self.codon_table = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }
        
        self.amino_acids = {
            'A': 'Alanine', 'R': 'Arginine', 'N': 'Asparagine', 'D': 'Aspartic acid',
            'C': 'Cysteine', 'E': 'Glutamic acid', 'Q': 'Glutamine', 'G': 'Glycine',
            'H': 'Histidine', 'I': 'Isoleucine', 'L': 'Leucine', 'K': 'Lysine',
            'M': 'Methionine', 'F': 'Phenylalanine', 'P': 'Proline', 'S': 'Serine',
            'T': 'Threonine', 'W': 'Tryptophan', 'Y': 'Tyrosine', 'V': 'Valine',
            '*': 'Stop codon'
        }
        
        # Known gene patterns and functions
        self.gene_database = {
            'CYTB': {'name': 'Cytochrome b', 'function': 'Electron transport chain', 'type': 'Mitochondrial'},
            'COX1': {'name': 'Cytochrome c oxidase I', 'function': 'Cellular respiration', 'type': 'Mitochondrial'},
            'ACTB': {'name': 'Beta-actin', 'function': 'Cytoskeleton structure', 'type': 'Nuclear'},
            'GAPDH': {'name': 'Glyceraldehyde-3-phosphate dehydrogenase', 'function': 'Glycolysis', 'type': 'Nuclear'},
            'TP53': {'name': 'Tumor protein p53', 'function': 'Cell cycle regulation', 'type': 'Nuclear'},
            'BRCA1': {'name': 'Breast cancer 1', 'function': 'DNA repair', 'type': 'Nuclear'},
            'APOE': {'name': 'Apolipoprotein E', 'function': 'Lipid transport', 'type': 'Nuclear'},
            'CFTR': {'name': 'Cystic fibrosis transmembrane conductance regulator', 'function': 'Ion transport', 'type': 'Nuclear'}
        }
        
        # Species-specific markers
        self.species_markers = {
            'human': ['ACTB', 'GAPDH', 'TP53', 'BRCA1'],
            'chimpanzee': ['CYTB', 'COX1', 'ACTB'],
            'mouse': ['ACTB', 'GAPDH'],
            'fruit_fly': ['ACTB'],
            'yeast': ['GAPDH']
        }

    def parse_sequence_file(self, file_content: str) -> Dict:
        """Parse uploaded sequence file"""
        lines = file_content.strip().split('\n')
        sequences = []
        
        for line in lines:
            if line.strip():
                parts = line.split('\t')
                if len(parts) >= 2:
                    sequence = parts[0].strip().upper()
                    class_label = parts[1].strip() if len(parts) > 1 else 'Unknown'
                    
                    # Clean sequence - keep only ATCG
                    clean_sequence = re.sub(r'[^ATCG]', '', sequence)
                    
                    if len(clean_sequence) > 0:
                        sequences.append({
                            'sequence': clean_sequence,
                            'class': class_label,
                            'length': len(clean_sequence)
                        })
        
        return {'sequences': sequences, 'total_sequences': len(sequences)}

    def analyze_composition(self, sequence: str) -> Dict:
        """Analyze nucleotide composition and patterns"""
        composition = Counter(sequence)
        length = len(sequence)
        
        # Basic composition
        gc_content = (composition['G'] + composition['C']) / length * 100 if length > 0 else 0
        at_content = (composition['A'] + composition['T']) / length * 100 if length > 0 else 0
        
        # Dinucleotide analysis
        dinucleotides = {}
        for i in range(len(sequence) - 1):
            dinuc = sequence[i:i+2]
            dinucleotides[dinuc] = dinucleotides.get(dinuc, 0) + 1
        
        # CpG islands detection
        cpg_count = dinucleotides.get('CG', 0)
        cpg_ratio = cpg_count / (length / 2) if length > 0 else 0
        
        return {
            'composition': dict(composition),
            'length': length,
            'gc_content': gc_content,
            'at_content': at_content,
            'dinucleotides': dinucleotides,
            'cpg_count': cpg_count,
            'cpg_ratio': cpg_ratio
        }

    def find_orfs(self, sequence: str) -> List[Dict]:
        """Find Open Reading Frames"""
        orfs = []
        
        for frame in range(3):
            for strand in [sequence, self.reverse_complement(sequence)]:
                strand_name = "+" if strand == sequence else "-"
                
                for i in range(frame, len(strand) - 2, 3):
                    codon = strand[i:i+3]
                    if len(codon) == 3 and codon == 'ATG':  # Start codon
                        # Look for stop codon
                        for j in range(i + 3, len(strand) - 2, 3):
                            stop_codon = strand[j:j+3]
                            if len(stop_codon) == 3 and stop_codon in ['TAA', 'TAG', 'TGA']:
                                orf_seq = strand[i:j+3]
                                if len(orf_seq) >= 30:  # Minimum ORF length
                                    protein_seq = self.translate_dna(orf_seq)
                                    orfs.append({
                                        'start': i + 1,
                                        'end': j + 3,
                                        'strand': strand_name,
                                        'frame': frame + 1,
                                        'length': len(orf_seq),
                                        'dna_sequence': orf_seq,
                                        'protein_sequence': protein_seq,
                                        'protein_length': len(protein_seq) - 1  # Exclude stop codon
                                    })
                                break
        
        # Sort by length (longest first)
        orfs.sort(key=lambda x: x['length'], reverse=True)
        return orfs

    def reverse_complement(self, sequence: str) -> str:
        """Generate reverse complement of DNA sequence"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return ''.join(complement.get(base, base) for base in reversed(sequence))

    def translate_dna(self, sequence: str) -> str:
        """Translate DNA sequence to protein"""
        protein = []
        for i in range(0, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if len(codon) == 3:
                amino_acid = self.codon_table.get(codon, 'X')
                protein.append(amino_acid)
        return ''.join(protein)

    def identify_genes(self, sequence: str) -> List[Dict]:
        """Identify potential genes using pattern matching"""
        genes_found = []
        
        # Simple gene identification based on known patterns
        for gene_id, gene_info in self.gene_database.items():
            # Look for gene-like patterns (this is simplified)
            pattern_score = 0
            
            # Check for start codons
            start_codons = sequence.count('ATG')
            if start_codons > 0:
                pattern_score += start_codons * 10
            
            # Check for stop codons
            stop_codons = sequence.count('TAA') + sequence.count('TAG') + sequence.count('TGA')
            if stop_codons > 0:
                pattern_score += stop_codons * 5
            
            # Check GC content (genes typically have specific GC patterns)
            gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
            if 40 <= gc_content <= 60:
                pattern_score += 20
            
            if pattern_score > 30:  # Threshold for gene-like sequence
                genes_found.append({
                    'gene_id': gene_id,
                    'name': gene_info['name'],
                    'function': gene_info['function'],
                    'type': gene_info['type'],
                    'confidence': min(95, pattern_score),
                    'predicted_location': f"1-{len(sequence)}"
                })
        
        return genes_found

    def species_classification(self, sequence: str) -> Dict:
        """Classify species based on sequence patterns"""
        species_scores = {}
        
        # Analyze sequence characteristics
        gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
        
        # Species-specific patterns
        if 45 <= gc_content <= 50:
            species_scores['human'] = 0.8
            species_scores['chimpanzee'] = 0.7
        elif 38 <= gc_content <= 45:
            species_scores['mouse'] = 0.6
        elif gc_content > 50:
            species_scores['fruit_fly'] = 0.5
        
        # Check for mitochondrial markers
        if 'ATGCCC' in sequence or 'CYTB' in sequence:
            species_scores['chimpanzee'] = species_scores.get('chimpanzee', 0) + 0.2
        
        # Normalize scores
        if species_scores:
            max_score = max(species_scores.values())
            for species in species_scores:
                species_scores[species] = species_scores[species] / max_score
        
        return species_scores

    def predict_protein_structure(self, protein_sequence: str) -> Dict:
        """Predict protein structure and properties (simulated)"""
        # Remove stop codon if present
        protein_sequence = protein_sequence.replace('*', '')

        if not protein_sequence:
            return {'error': 'Invalid or empty protein sequence'}
        
        # Amino acid composition
        aa_composition = Counter(protein_sequence)
        
        # Predict secondary structure (simplified simulation)
        structure_seq = []
        for aa in protein_sequence:
            # Assign structure based on amino acid properties (simplified)
            if aa in ['A', 'E', 'L', 'M', 'Q', 'K', 'R', 'H']: # Helix formers
                structure_seq.append('H')
            elif aa in ['F', 'I', 'Y', 'V', 'W', 'T', 'C']: # Sheet formers
                structure_seq.append('S')
            else: # Coil/turn formers
                structure_seq.append('C')
        
        # Smooth the structure to create longer segments
        smoothed_structure = list(structure_seq)
        for i in range(1, len(smoothed_structure) - 1):
            if smoothed_structure[i-1] == smoothed_structure[i+1]:
                smoothed_structure[i] = smoothed_structure[i-1]
        
        # Predict domains based on structure
        domains = []
        if 'HHHHH' in "".join(smoothed_structure):
            domains.append('Alpha-Helix Rich')
        if 'SSSSS' in "".join(smoothed_structure):
            domains.append('Beta-Sheet Rich')
        if len(protein_sequence) > 100 and 'C'*10 in "".join(smoothed_structure):
            domains.append('Disordered Region')
        if not domains:
            domains.append('Mixed Alpha/Beta')
        
        # Calculate properties
        molecular_weight = len(protein_sequence) * 110  # Average MW per residue
        
        return {
            'length': len(protein_sequence),
            'composition': dict(aa_composition),
            'molecular_weight': molecular_weight,
            'predicted_domains': domains,
            'secondary_structure_sequence': "".join(smoothed_structure),
            'secondary_structure_summary': {
                'helix_percent': smoothed_structure.count('H') / len(smoothed_structure) * 100 if smoothed_structure else 0,
                'sheet_percent': smoothed_structure.count('S') / len(smoothed_structure) * 100 if smoothed_structure else 0,
                'coil_percent': smoothed_structure.count('C') / len(smoothed_structure) * 100 if smoothed_structure else 0,
            }
        }

    def _score_guide_rna(self, guide_seq: str, full_sequence: str) -> Dict:
        """Scores a guide RNA based on on-target and off-target metrics."""
        score = 100.0
        
        # 1. GC content (ideal 40-60%)
        gc_content = (guide_seq.count('G') + guide_seq.count('C')) / len(guide_seq) * 100
        if not 40 <= gc_content <= 60:
            score -= abs(gc_content - 50)  # Penalize deviation from 50%
        
        # 2. Poly-T termination signal (bad)
        if 'TTTT' in guide_seq:
            score -= 40
            
        # 3. On-target score (simple heuristic: GC clamp at the end is good)
        if guide_seq.endswith('G'):
            score += 10
            
        # 4. Off-target score (simple version: count exact matches elsewhere)
        # A real tool would use BLAST or allow mismatches.
        off_target_hits = full_sequence.count(guide_seq) - 1
        score -= off_target_hits * 10 # Penalize each off-target hit
        
        return {
            'gc_content': gc_content,
            'on_target_score': max(0, min(100, score)),
            'off_target_hits': off_target_hits
        }

    def design_crispr_guides(self, sequence: str) -> List[Dict]:
        """Design and score guide RNAs for CRISPR/Cas9."""
        guides = []
        pam = 'GG' # We look for NGG
        
        # Search forward strand
        for match in re.finditer(pam, sequence):
            pam_start = match.start()
            if pam_start >= 21: # Need 20nt guide + 1nt N
                guide_start = pam_start - 21
                guide_seq = sequence[guide_start:guide_start+20]
                
                scores = self._score_guide_rna(guide_seq, sequence)
                overall_score = scores['on_target_score']
                
                if overall_score > 50: # Filter for decent guides
                    guides.append({
                        'sequence': guide_seq,
                        'location': guide_start,
                        'pam_location': pam_start,
                        'strand': '+',
                        'gc_content': scores['gc_content'],
                        'on_target_score': scores['on_target_score'],
                        'off_target_hits': scores['off_target_hits'],
                        'overall_score': overall_score
                    })

        # Search reverse strand
        rev_comp = self.reverse_complement(sequence)
        for match in re.finditer(pam, rev_comp):
            pam_start_rev = match.start()
            if pam_start_rev >= 21:
                guide_start_rev = pam_start_rev - 21
                guide_seq_rev = rev_comp[guide_start_rev:guide_start_rev+20]
                
                # Guide sequence should be reported as it appears on the forward strand
                guide_seq_fwd = self.reverse_complement(guide_seq_rev)
                
                scores = self._score_guide_rna(guide_seq_rev, rev_comp)
                overall_score = scores['on_target_score']
                
                if overall_score > 50:
                    # Location on forward strand
                    fwd_location = len(sequence) - (guide_start_rev + 20)
                    guides.append({
                        'sequence': guide_seq_fwd,
                        'location': fwd_location,
                        'pam_location': len(sequence) - pam_start_rev,
                        'strand': '-',
                        'gc_content': scores['gc_content'],
                        'on_target_score': scores['on_target_score'],
                        'off_target_hits': scores['off_target_hits'],
                        'overall_score': overall_score
                    })
        
        # Sort by overall score
        guides.sort(key=lambda x: x['overall_score'], reverse=True)
        return guides

    def evolutionary_analysis(self, sequence: str) -> Dict:
        """Perform evolutionary analysis"""
        # Calculate mutation hotspots
        mutation_hotspots = []
        
        # Look for repetitive elements
        for i in range(len(sequence) - 6):
            motif = sequence[i:i+6]
            count = sequence.count(motif)
            if count > 2:
                mutation_hotspots.append({
                    'motif': motif,
                    'frequency': count,
                    'positions': [j for j in range(len(sequence)-5) if sequence[j:j+6] == motif]
                })
        
        # Calculate evolutionary conservation score (simplified)
        conservation_score = 85 + np.random.normal(0, 10)  # Simulated
        conservation_score = max(0, min(100, conservation_score))
        
        return {
            'conservation_score': conservation_score,
            'mutation_hotspots': sorted(mutation_hotspots, key=lambda x: x['frequency'], reverse=True)[:5],
            'evolutionary_rate': 'Slow' if conservation_score > 80 else 'Moderate' if conservation_score > 60 else 'Fast'
        }

    def pharmacogenomics_analysis(self, sequence: str) -> Dict:
        """Analyze pharmacogenomic implications"""
        drug_targets = []
        
        # Look for known drug target patterns
        if 'ATGCCC' in sequence:
            drug_targets.append({
                'target': 'Cytochrome P450',
                'drugs': ['Warfarin', 'Clopidogrel', 'Codeine'],
                'significance': 'High',
                'recommendation': 'Dose adjustment may be required'
            })
        
        if sequence.count('CG') > len(sequence) * 0.05:
            drug_targets.append({
                'target': 'Methylation sites',
                'drugs': ['5-Azacytidine', 'Decitabine'],
                'significance': 'Medium',
                'recommendation': 'Monitor methylation status'
            })
        
        return {
            'drug_targets': drug_targets,
            'pharmacogenomic_score': len(drug_targets) * 25
        }

    def epigenetic_analysis(self, sequence: str) -> Dict:
        """Analyzes epigenetic features like methylation potential."""
        cpg_islands = []
        
        # A simple CpG island definition: region > 200bp, GC > 50%, Obs/Exp CpG > 0.6
        window_size = 200
        step_size = 50
        
        if len(sequence) < window_size:
            return {
                'cpg_islands_found': 0,
                'cpg_islands': [],
                'g_quadruplexes_found': 0,
                'g_quadruplexes': [],
                'overall_methylation_potential': 'Low'
            }

        for i in range(0, len(sequence) - window_size, step_size):
            window = sequence[i:i+window_size]
            g_count = window.count('G')
            c_count = window.count('C')
            gc_content = (g_count + c_count) / window_size * 100
            
            cpg_count = window.count('CG')
            # Observed/Expected ratio
            if g_count > 0 and c_count > 0:
                expected_cpg = (c_count * g_count) / window_size
                obs_exp_ratio = cpg_count / expected_cpg if expected_cpg > 0 else 0
            else:
                obs_exp_ratio = 0
                
            if gc_content > 50 and obs_exp_ratio > 0.6:
                cpg_islands.append({
                    'start': i,
                    'end': i + window_size,
                    'gc_content': gc_content,
                    'obs_exp_ratio': obs_exp_ratio
                })
                
        # G-quadruplex prediction (simple pattern)
        g_quadruplex_pattern = r'G{3,5}\w{1,7}G{3,5}\w{1,7}G{3,5}\w{1,7}G{3,5}'
        g_quadruplexes = []
        for match in re.finditer(g_quadruplex_pattern, sequence):
            g_quadruplexes.append({
                'start': match.start(),
                'end': match.end(),
                'sequence': match.group(0)
            })
            
        return {
            'cpg_islands_found': len(cpg_islands),
            'cpg_islands': cpg_islands,
            'g_quadruplexes_found': len(g_quadruplexes),
            'g_quadruplexes': g_quadruplexes,
            'overall_methylation_potential': 'High' if len(cpg_islands) > 2 else 'Moderate' if len(cpg_islands) > 0 else 'Low'
        }

    def predict_non_coding_rna(self, sequence: str) -> Dict:
        """Predicts non-coding RNAs (ncRNAs) like miRNA and lncRNA."""
        mirnas = []
        lncrnas = []
        
        # 1. miRNA prediction (simplified: look for hairpin potential)
        # Look for short inverted repeats which could form hairpins.
        window_size = 60
        min_hairpin_len = 15
        if len(sequence) > window_size:
            for i in range(len(sequence) - window_size):
                window = sequence[i:i+window_size]
                # Find a short inverted repeat within the window
                for l in range(min_hairpin_len, window_size // 2):
                    sub = window[:l]
                    rev_comp_sub = self.reverse_complement(sub)
                    if rev_comp_sub in window[l:]:
                        mirnas.append({
                            'start': i,
                            'end': i + window_size,
                            'potential_hairpin_seq': window,
                            'confidence': 70 + len(sub) # Simple confidence score
                        })
                        break # Move to next window
        
        # 2. lncRNA prediction (simplified: long regions with no long ORFs)
        orfs = self.find_orfs(sequence)
        # Create a map of where ORFs are
        orf_map = [0] * len(sequence)
        for orf in orfs:
            if orf['length'] > 150: # A somewhat significant ORF
                for i in range(orf['start'] - 1, orf['end']):
                    if i < len(orf_map):
                        orf_map[i] = 1
        
        # Find long stretches of non-coding regions
        in_lncrna = False
        lnc_start = 0
        for i in range(len(orf_map)):
            if orf_map[i] == 0 and not in_lncrna:
                in_lncrna = True
                lnc_start = i
            elif orf_map[i] == 1 and in_lncrna:
                in_lncrna = False
                length = i - lnc_start
                if length > 200: # lncRNA minimum length
                    lncrnas.append({
                        'start': lnc_start,
                        'end': i,
                        'length': length,
                        'gc_content': (sequence[lnc_start:i].count('G') + sequence[lnc_start:i].count('C')) / length * 100 if length > 0 else 0
                    })
        if in_lncrna:
            length = len(sequence) - lnc_start
            if length > 200:
                lncrnas.append({
                    'start': lnc_start,
                    'end': len(sequence),
                    'length': length,
                    'gc_content': (sequence[lnc_start:].count('G') + sequence[lnc_start:].count('C')) / length * 100 if length > 0 else 0
                })

        return {
            'mirnas_found': len(mirnas),
            'mirnas': sorted(mirnas, key=lambda x: x['confidence'], reverse=True)[:5],
            'lncrnas_found': len(lncrnas),
            'lncrnas': sorted(lncrnas, key=lambda x: x['length'], reverse=True)[:5]
        }

    def detect_viral_integration(self, sequence: str) -> Dict:
        """Detects potential viral integration sites."""
        viral_signatures = {
            'HPV16': 'ACCGACAGCTCAGAGGAGGAG', # Human Papillomavirus 16
            'HBV': 'AATTCCACCAAGCCTCCAA',     # Hepatitis B Virus
            'HIV-1': 'GGTCTCTCTGGTTAGACCAGA', # Human Immunodeficiency Virus 1
            'EBV': 'AGGAACCAAGTCGATCTTCC',    # Epstein-Barr Virus
        }
        integrations = []
        for virus, signature in viral_signatures.items():
            for match in re.finditer(signature, sequence):
                integrations.append({
                    'virus': virus,
                    'start': match.start(),
                    'end': match.end(),
                    'signature_sequence': signature,
                    'confidence': 95.0 # High confidence if exact signature is found
                })
        return {
            'integrations_found': len(integrations),
            'integrations': sorted(integrations, key=lambda x: x['confidence'], reverse=True)
        }

    def microbiome_composition_analysis(self, sequence: str) -> Dict:
        """Analyzes sequence for microbiome composition markers (simulated)."""
        # Simplified markers for major bacterial phyla (16S rRNA variable regions)
        microbiome_markers = {
            'Firmicutes': 'TGGAGCATGTGGTTTAATTC', # e.g., Lactobacillus, Clostridium
            'Bacteroidetes': 'GGTCTGAGAGGATGATCAGT', # e.g., Bacteroides
            'Actinobacteria': 'GCTGGCACGTAGGTAGCC',   # e.g., Bifidobacterium
            'Proteobacteria': 'CCAGACTCCTACGGGAGGC',  # e.g., E. coli, Salmonella
            'Verrucomicrobia': 'GTCGTCAGCTCGTGTCGTGA' # e.g., Akkermansia
        }
        
        composition = {}
        total_hits = 0
        
        for phylum, marker in microbiome_markers.items():
            hits = len(re.findall(marker, sequence))
            if hits > 0:
                composition[phylum] = hits
                total_hits += hits
        
        # Normalize to percentages
        if total_hits > 0:
            for phylum in composition:
                composition[phylum] = (composition[phylum] / total_hits) * 100
        
        # Add a 'dominant_phylum' field
        dominant_phylum = "Unknown"
        if composition:
            dominant_phylum = max(composition, key=composition.get)
            
        return {
            'composition': composition,
            'total_marker_hits': total_hits,
            'dominant_phylum': dominant_phylum,
            'diversity_score': len(composition) / len(microbiome_markers) * 100 if microbiome_markers else 0
        }

    def detect_structural_variants(self, sequence: str) -> Dict:
        """Detects potential structural variants like duplications and inversions."""
        variants = []
        
        # 1. Tandem Duplication detection (simplified)
        # Look for long repeated segments next to each other
        window = 50
        if len(sequence) > window * 2:
            for i in range(len(sequence) - window * 2):
                segment1 = sequence[i:i+window]
                segment2 = sequence[i+window:i+window*2]
                if segment1 == segment2:
                    variants.append({
                        'type': 'Tandem Duplication',
                        'start': i,
                        'end': i + window * 2,
                        'length': window,
                        'sequence': segment1,
                        'confidence': 80.0
                    })
        
        # 2. Inversion detection (simplified)
        # Look for a segment followed by its reverse complement
        window = 40
        if len(sequence) > window * 2:
            for i in range(len(sequence) - window * 2):
                segment = sequence[i:i+window]
                rev_comp_segment = self.reverse_complement(segment)
                # Look for the reverse complement nearby
                next_segment = sequence[i+window:i+window*2]
                if rev_comp_segment == next_segment:
                    variants.append({
                        'type': 'Inversion',
                        'start': i,
                        'end': i + window * 2,
                        'length': window * 2,
                        'sequence': segment + '...' + next_segment,
                        'confidence': 75.0
                    })

        return {
            'variants_found': len(variants),
            'variants': sorted(variants, key=lambda x: x['start'])
        }

    def comparative_genomics(self, sequence: str) -> Dict:
        """Performs a simplified comparative genomics analysis."""
        # Simulated database of orthologous genes from different species
        ortholog_db = {
            'Human_TP53': 'ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATGGATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCTGTCCCTTCCCAGAAAACCTACCAGGGCAGCTACGGTTTCCGTCTGGGCTTCTTGCATTCTGGGACAGCCAAGTCTGTGACTTGCACGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACCTGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATGGTCTGGCCCCTCCTCAGCATCTTATCCGAGTGGAAGGAAATTTGCGTGTGGAGTATTTGGATGACAGAAACACTTTTCGACATAGTGTGGTGGTGCCCTATGAGCCGCCTGAGGTTGGCTCTGACTGTACCACCATCCACTACAACTACATGTGTAACAGTTCCTGCATGGGCGGCATGAACCGGAGGCCCATCCTCACCATCATCACACTGGAAGACTCCAGTGGTAATCTACTGGGACGGAACAGCTTTGAGGTGCGTGTTTGTGCCTGTCCTGGGAGAGACCGGCGCACAGAGGAAGAGAATCTCCGCAAGAAAGGGGAGCCTCACCACGAGCTGCCCCCAGGGAGCACTAAGCGAGCACTGCCCAACAACACCAGCTCCTCTCCCCAGCCAAAGAAGAAACCACTGGATGGAGAATATTTCACCCTTCAGATCCGTGGGCGTGAGCGCTTCGAGATGTTCCGAGAGCTGAATGAGGCCTTAGAGTACCCGGACAGGCCT',
            'Mouse_Trp53': 'ATGGAGGATTCGCAGTCAGATCCTAGCATTCGAGCCCCCTCTGAGTCAGGAAGCATTTTCAGGCCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATGGATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCTGTCCCTTCCCAGAAAACCTACCAGGGCAGCTACGGTTTCCGTCTGGGCTTCTTGCATTCTGGGACAGCCAAGTCTGTGACTTGCACGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACCTGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATGGTCTGGCCCCTCCTCAGCATCTTATCCGAGTGGAAGGAAATTTGCGTGTGGAGTATTTGGATGACAGAAACACTTTTCGACATAGTGTGGTGGTGCCCTATGAGCCGCCTGAGGTTGGCTCTGACTGTACCACCATCCACTACAACTACATGTGTAACAGTTCCTGCATGGGCGGCATGAACCGGAGGCCCATCCTCACCATCATCACACTGGAAGACTCCAGTGGTAATCTACTGGGACGGAACAGCTTTGAGGTGCGTGTTTGTGCCTGTCCTGGGAGAGACCGGCGCACAGAGGAAGAGAATCTCCGCAAGAAAGGGGAGCCTCACCACGAGCTGCCCCCAGGGAGCACTAAGCGAGCACTGCCCAACAACACCAGCTCCTCTCCCCAGCCAAAGAAGAAACCACTGGATGGAGAATATTTCACCCTTCAGATCCGTGGGCGTGAGCGCTTCGAGATGTTCCGAGAGCTGAATGAGGCCTTAGAGTACCCGGACAGGCCT',
            'Zebrafish_tp53': 'ATGGAGGACTCTCAATCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATGGATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCTGTCCCTTCCCAGAAAACCTACCAGGGCAGCTACGGTTTCCGTCTGGGCTTCTTGCATTCTGGGACAGCCAAGTCTGTGACTTGCACGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACCTGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATGGTCTGGCCCCTCCTCAGCATCTTATCCGAGTGGAAGGAAATTTGCGTGTGGAGTATTTGGATGACAGAAACACTTTTCGACATAGTGTGGTGGTGCCCTATGAGCCGCCTGAGGTTGGCTCTGACTGTACCACCATCCACTACAACTACATGTGTAACAGTTCCTGCATGGGCGGCATGAACCGGAGGCCCATCCTCACCATCATCACACTGGAAGACTCCAGTGGTAATCTACTGGGACGGAACAGCTTTGAGGTGCGTGTTTGTGCCTGTCCTGGGAGAGACCGGCGCACAGAGGAAGAGAATCTCCGCAAGAAAGGGGAGCCTCACCACGAGCTGCCCCCAGGGAGCACTAAGCGAGCACTGCCCAACAACACCAGCTCCTCTCCCCAGCCAAAGAAGAAACCACTGGATGGAGAATATTTCACCCTTCAGATCCGTGGGCGTGAGCGCTTCGAGATGTTCCGAGAGCTGAATGAGGCCTTAGAGTACCCGGACAGGCCT'
        }
        
        homology_results = []
        
        # Simplified alignment: find longest common substring using dynamic programming
        for species_gene, ortholog_seq in ortholog_db.items():
            s1 = sequence.upper()
            s2 = ortholog_seq.upper()
            m = [[0] * (1 + len(s2)) for _ in range(1 + len(s1))]
            longest, x_longest = 0, 0
            for x in range(1, 1 + len(s1)):
                for y in range(1, 1 + len(s2)):
                    if s1[x - 1] == s2[y - 1]:
                        m[x][y] = m[x - 1][y - 1] + 1
                        if m[x][y] > longest:
                            longest = m[x][y]
                            x_longest = x
                    else:
                        m[x][y] = 0
            
            match_length = longest
            
            if match_length > 20: # Minimum significant match
                identity = (match_length / min(len(s1), len(s2))) * 100 if min(len(s1), len(s2)) > 0 else 0
                homology_results.append({
                    'ortholog': species_gene,
                    'identity': identity,
                    'match_length': match_length,
                    'e_value': 10 ** (-identity/10) # Simulated E-value
                })
                
        return {
            'homology_hits': sorted(homology_results, key=lambda x: x['identity'], reverse=True)
        }

    def predict_3d_folding(self, sequence: str) -> Dict:
        """Predicts 3D genome folding patterns like TADs and loops (simulated)."""
        # Simplified insulator motifs (like CTCF binding sites)
        insulator_motifs = ['CCCTCCTCCCC', 'GCGCGGGGGC']
        promoter_like = 'TATAAT'
        enhancer_like = 'TGTGTC'
        
        tads = []
        loops = []
        
        # 1. Find TAD boundaries
        boundaries = [0]
        for motif in insulator_motifs:
            for match in re.finditer(motif, sequence):
                boundaries.append(match.start())
        boundaries.append(len(sequence))
        boundaries = sorted(list(set(boundaries)))
        
        for i in range(len(boundaries) - 1):
            start, end = boundaries[i], boundaries[i+1]
            if end - start > 1000: # Minimum TAD size
                tads.append({
                    'start': start,
                    'end': end,
                    'length': end - start,
                    'type': 'Predicted TAD'
                })

        # 2. Find loops within TADs
        promoters = [m.start() for m in re.finditer(promoter_like, sequence)]
        enhancers = [m.start() for m in re.finditer(enhancer_like, sequence)]
        
        for tad in tads:
            tad_promoters = [p for p in promoters if tad['start'] <= p < tad['end']]
            tad_enhancers = [e for e in enhancers if tad['start'] <= e < tad['end']]
            
            # Connect each promoter to the nearest enhancer in the same TAD
            for p_pos in tad_promoters:
                closest_e = -1
                min_dist = float('inf')
                for e_pos in tad_enhancers:
                    dist = abs(p_pos - e_pos)
                    if dist < min_dist and dist > 100: # Must be a long-range interaction
                        min_dist = dist
                        closest_e = e_pos
                
                if closest_e != -1:
                    loops.append({
                        'promoter_pos': p_pos,
                        'enhancer_pos': closest_e,
                        'distance': min_dist,
                        'tad_start': tad['start'],
                        'tad_end': tad['end'],
                        'strength': np.random.uniform(0.5, 1.0) # Simulated interaction strength
                    })

        return {
            'tads_found': len(tads),
            'tads': sorted(tads, key=lambda x: x['length'], reverse=True),
            'loops_found': len(loops),
            'loops': sorted(loops, key=lambda x: x['strength'], reverse=True)
        }

    def infer_regulatory_network(self, sequence: str, orfs: List[Dict]) -> Dict:
        """Infers a gene regulatory network from TFBS and gene locations (simulated)."""
        # Simplified Transcription Factor Binding Site (TFBS) motifs
        tfbs_db = {
            'SP1': 'GGGCGG',
            'AP-1': 'TGACTCA',
            'NF-kB': 'GGGACTTTCC',
            'CREB': 'TGACGTCA',
            'MYC': 'CACGTG'
        }
        
        regulations = []
        
        # Find ORFs to act as target genes
        target_genes = [{'id': f"ORF_{i+1}", 'start': orf['start']} for i, orf in enumerate(orfs)]
        
        if not target_genes:
            return {'regulations': []}
            
        for tf, motif in tfbs_db.items():
            for match in re.finditer(motif, sequence):
                tfbs_pos = match.start()
                
                # Find the closest downstream gene
                closest_gene = None
                min_dist = float('inf')
                
                for gene in target_genes:
                    dist = gene['start'] - tfbs_pos
                    if 0 < dist < 5000: # Look in a 5kb upstream window
                        if dist < min_dist:
                            min_dist = dist
                            closest_gene = gene['id']
                
                if closest_gene:
                    regulations.append({
                        'tf': tf,
                        'target_gene': closest_gene,
                        'tfbs_location': tfbs_pos,
                        'distance_to_gene': min_dist,
                        'regulation_type': 'Activation' if np.random.rand() > 0.3 else 'Repression' # Simulated
                    })
                    
        return {'regulations': regulations}

    def design_synthetic_circuits(self, sequence: str, orfs: List[Dict]) -> Dict:
        """Proposes synthetic biology circuits based on the sequence (simulated)."""
        # Library of synthetic biology parts
        parts_db = {
            'promoters': [{'id': 'pTet', 'sequence': 'TCCCTATCAGTGATAGAGATTGACATCCCTATCAGTGATAGAGATACTGAGCAC', 'strength': 95}],
            'rbs': [{'id': 'RBS_strong', 'sequence': 'AGGAGG', 'strength': 90}],
            'terminators': [{'id': 'Term_double', 'sequence': 'CCAGGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCGGTGAACGCTCT', 'efficiency': 98}],
            'reporters': [{'id': 'GFP', 'sequence': 'ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAATTGGAATACAACTATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATCAAAGTTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA'}]
        }
        designs = []
        # Design 1: Reporter construct
        p, r, g, t = parts_db['promoters'][0], parts_db['rbs'][0], parts_db['reporters'][0], parts_db['terminators'][0]
        design1_parts = [
            {'part': 'Promoter', 'id': p['id'], 'length': len(p['sequence'])},
            {'part': 'RBS', 'id': r['id'], 'length': len(r['sequence'])},
            {'part': 'Reporter Gene', 'id': g['id'], 'length': len(g['sequence'])},
            {'part': 'Terminator', 'id': t['id'], 'length': len(t['sequence'])}
        ]
        designs.append({
            'name': 'High-Expression GFP Reporter', 'description': 'A standard construct for testing gene expression levels.',
            'parts': design1_parts,
            'total_length': sum(part['length'] for part in design1_parts)
        })
        # Design 2: Regulate a native gene
        if orfs:
            design2_parts = [
                {'part': 'Promoter', 'id': p['id'], 'length': len(p['sequence'])},
                {'part': 'Native Gene', 'id': 'Largest ORF', 'length': orfs[0]['length']},
                {'part': 'Terminator', 'id': t['id'], 'length': len(t['sequence'])}
            ]
            designs.append({
                'name': 'Inducible Control of Native Gene', 'description': f"Places the largest identified ORF (length {orfs[0]['length']} bp) under the control of the inducible pTet promoter.",
                'parts': design2_parts,
                'total_length': sum(part['length'] for part in design2_parts)
            })
        return {'designs': designs}

    def generate_pdb_from_structure(self, protein_sequence: str, structure_sequence: str) -> str:
        """Generates a simulated PDB file string from a protein and its secondary structure."""
        pdb_lines = []
        atom_index = 1
        residue_index = 1
        x, y, z = 0.0, 0.0, 0.0
        
        # Helix parameters
        helix_radius = 2.5
        helix_pitch = 1.5 # Rise per residue
        
        # Sheet parameters
        sheet_step = 3.5
        
        for aa, struct in zip(protein_sequence, structure_sequence):
            if aa == '*': continue # Skip stop codons
            
            res_name = self.amino_acids.get(aa, "UNK").upper()[:3]
            
            pdb_line = f"ATOM  {atom_index:5d}  CA  {res_name} A{residue_index:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C  "
            pdb_lines.append(pdb_line)
            
            if struct == 'H': # Helix
                angle = (residue_index % 4) * (np.pi / 2) # simplified turn
                x += helix_pitch
                y = helix_radius * np.cos(angle)
                z = helix_radius * np.sin(angle)
            elif struct == 'S': # Sheet
                x += sheet_step * 0.8
                y += (residue_index % 2 - 0.5) * 2 * sheet_step * 0.2 # Zig-zag
            else: # Coil
                # Random walk
                x += np.random.uniform(-1.5, 1.5)
                y += np.random.uniform(-1.5, 1.5)
                z += np.random.uniform(-1.5, 1.5)
            
            atom_index += 1
            residue_index += 1
            
        return "\n".join(pdb_lines)

def create_sequence_visualization(sequence: str, max_length: int = 200) -> str:
    """Create formatted sequence visualization"""
    if len(sequence) <= max_length:
        display_seq = sequence
    else:
        display_seq = sequence[:max_length] + f"... ({len(sequence)-max_length} more nucleotides)"
    
    # Add spacing every 10 nucleotides
    formatted = ""
    for i, nucleotide in enumerate(display_seq):
        if i > 0 and i % 10 == 0:
            formatted += " "
        if i > 0 and i % 50 == 0:
            formatted += "\n"
        formatted += nucleotide
    
    return formatted

def visualize_protein_structure(structure_sequence: str) -> str:
    """Create an HTML visualization of the protein secondary structure."""
    html = "<div style='font-family: monospace; line-height: 20px; word-wrap: break-word;'>"
    color_map = {
        'H': '#FF6B6B', # Helix - red
        'S': '#4ECDC4', # Sheet - teal
        'C': '#45B7D1'  # Coil - blue
    }
    for element in structure_sequence:
        html += f"<span style='display: inline-block; width: 12px; height: 12px; background-color: {color_map.get(element, '#ccc')}; border-radius: 2px; margin: 1px;' title='{element}'></span>"
    html += "</div>"
    
    # Add a legend
    html += """
    <div style='margin-top: 10px; font-size: 12px;'>
        <span style='display: inline-block; width: 10px; height: 10px; background-color: #FF6B6B; margin-right: 5px;'></span> Helix (H)
        <span style='display: inline-block; width: 10px; height: 10px; background-color: #4ECDC4; margin-left: 10px; margin-right: 5px;'></span> Sheet (S)
        <span style='display: inline-block; width: 10px; height: 10px; background-color: #45B7D1; margin-left: 10px; margin-right: 5px;'></span> Coil (C)
    </div>
    """
    return html

def render_3d_protein_structure(pdb_string: str, width: int = 500, height: int = 400):
    """Renders a 3D protein structure using py3Dmol."""
    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_string, 'pdb')
    
    # Style the protein
    view.setStyle({'cartoon': {'color': 'spectrum'}})
    
    # Zoom to fit the structure
    view.zoomTo()
    
    # Generate the HTML for the viewer
    html_output = view.to_html()
    
    if html_output:
        # FIX: Explicitly convert the output to a string to avoid the TypeError
        components.html(str(html_output), height=int(height), width=int(width))
    else:
        st.warning("Could not generate a 3D visualization for this structure.")

def main():
    # App header
    st.markdown("""
    <div class="main-header">
        <h1>🧬 EternaSeq</h1>
        <h2>The Nobel Prize-Winning Genomic Discovery Engine</h2>
        <p>Unlock groundbreaking biological discoveries from any DNA sequence. Your next Nobel Prize starts here.</p>
    </div>
    """, unsafe_allow_html=True)

    with st.sidebar:
        st.title("🔬 Beta Features")
        st.info("Explore cutting-edge experimental analyses. These features are under active development.")
        
        beta_features = [
            "Select a Beta Feature...",
            "1. Epigenetic Analysis",
            "2. Non-coding RNA Prediction",
            "3. Viral Integration Site Detection",
            "4. Microbiome Composition Analysis",
            "5. Structural Variant Detection",
            "6. Comparative Genomics",
            "7. 3D Genome Folding Prediction",
            "8. Gene Regulatory Network Inference",
            "9. Synthetic Biology Circuit Design",
            "10. AI-Powered Drug Discovery"
        ]
        
        st.selectbox(
            "Choose a feature to run:",
            options=beta_features,
            key="selected_beta_feature"
        )


    # File upload section
    st.header("📁 DNA Sequence Upload")
    
    uploaded_file = st.file_uploader(
        "Upload DNA sequence file (e.g., chimpanzee.txt)",
        type=['txt', 'fasta', 'fa', 'seq'],
        help="File format: sequence[tab]class (e.g., ATGCCC...	4)"
    )

    # Sample data for demo
    col1, col2 = st.columns([3, 1])
    with col1:
        if not uploaded_file:
            st.info("💡 Upload a DNA sequence file or try the demo below")
    
    with col2:
        if st.button("🚀 Try Sample Data"):
            # Create sample chimpanzee sequence
            sample_content = "ATGCCCCAACTAAATACCGCCGTATGACCCACCATAATTACCCCCATACTCCTGACACTATTTCTCGTCACCCAACTAAAAATATTAAATTCAAATTACCATCTACCCCCCTCACCAAAACCCATAAAAATAAAAAACTACAATAAACCCTGAGAACCAAAATGAACGAAAATCTATTCGCTTCATTCGCTGCCCCCACAATCCTAG\t4"
            
            # Store in session state
            st.session_state.sample_data = sample_content
            st.success("Sample chimpanzee sequence loaded!")

    # Process uploaded file or sample data
    file_content = None
    if uploaded_file:
        file_content = str(uploaded_file.read(), 'utf-8')
    elif st.session_state.get('sample_data'):
        file_content = st.session_state.sample_data

    if file_content:
        # Initialize analyzer
        analyzer = DNAAnalyzer()
        
        # Parse file
        with st.spinner("🧬 Parsing DNA sequences..."):
            parsed_data = analyzer.parse_sequence_file(file_content)
        
        if parsed_data['total_sequences'] == 0:
            st.error("No valid DNA sequences found. Please check file format.")
            return

        # Display sequence info
        st.success(f"✅ Successfully loaded {parsed_data['total_sequences']} DNA sequence(s)")
        
        # Process each sequence
        for seq_idx, seq_data in enumerate(parsed_data['sequences']):
            sequence = seq_data['sequence']
            class_label = seq_data['class']
            
            # =================================================================
            # == Perform all analyses once before rendering tabs for efficiency ==
            # =================================================================
            with st.spinner(f"🔬 Analyzing Sequence {seq_idx + 1}... This may take a moment."):
                composition = analyzer.analyze_composition(sequence)
                gc_content = composition['gc_content']
                orfs = analyzer.find_orfs(sequence)
                genes_found = analyzer.identify_genes(sequence)
                species_scores = analyzer.species_classification(sequence)
                evo_analysis = analyzer.evolutionary_analysis(sequence)
                pharma_analysis = analyzer.pharmacogenomics_analysis(sequence)
                crispr_guides = analyzer.design_crispr_guides(sequence)
                
                # Run beta analysis if selected
                epigenetic_results = None
                ncrna_results = None
                viral_results = None
                microbiome_results = None
                sv_results = None
                comp_gen_results = None
                folding_results = None
                reg_network_results = None
                synth_bio_results = None
                selected_feature = st.session_state.get('selected_beta_feature', "Select a Beta Feature...")
                if selected_feature.startswith("1."):
                    epigenetic_results = analyzer.epigenetic_analysis(sequence)
                elif selected_feature.startswith("2."):
                    ncrna_results = analyzer.predict_non_coding_rna(sequence)
                elif selected_feature.startswith("3."):
                    viral_results = analyzer.detect_viral_integration(sequence)
                elif selected_feature.startswith("4."):
                    microbiome_results = analyzer.microbiome_composition_analysis(sequence)
                elif selected_feature.startswith("5."):
                    sv_results = analyzer.detect_structural_variants(sequence)
                elif selected_feature.startswith("6."):
                    comp_gen_results = analyzer.comparative_genomics(sequence)
                elif selected_feature.startswith("7."):
                    folding_results = analyzer.predict_3d_folding(sequence)
                elif selected_feature.startswith("8."):
                    reg_network_results = analyzer.infer_regulatory_network(sequence, orfs)
                elif selected_feature.startswith("9."):
                    synth_bio_results = analyzer.design_synthetic_circuits(sequence, orfs)
                # Metrics for advanced tab
                entropy = -sum(p * np.log2(p) for p in [sequence.count(n)/len(sequence) for n in 'ATCG'] if p > 0)
                words_3 = [sequence[i:i+3] for i in range(len(sequence)-2)]
                unique_3mers = len(set(words_3))
                
                max_repeat = 0
                if sequence:
                    current_repeat = 1
                    max_repeat = 1
                    for i in range(1, len(sequence)):
                        if sequence[i] == sequence[i-1]:
                            current_repeat += 1
                        else:
                            max_repeat = max(max_repeat, current_repeat)
                            current_repeat = 1
                    max_repeat = max(max_repeat, current_repeat)

            st.markdown(f"---")
            st.header(f"🔬 Analysis Results - Sequence {seq_idx + 1}")
            
            # Basic sequence info
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.metric("Sequence Length", f"{composition['length']:,} bp")
            with col2:
                st.metric("Class Label", class_label)
            with col3:
                st.metric("GC Content", f"{gc_content:.1f}%")
            with col4:
                st.metric("Complexity", "High" if len(set(sequence)) == 4 else "Medium")

            # Create tabs for different analyses
            tab1, tab2, tab3, tab4, tab5, tab6, tab7 = st.tabs([
                "🧬 Sequence Analysis", "🔍 Gene Identification", "🧪 Protein Prediction", 
                "🌿 Species Classification", "⚗️ Drug Targets", "🔬 Advanced Analysis", "✂️ Gene Editing (CRISPR)"
            ])

            with tab1:
                st.subheader("📊 Nucleotide Composition Analysis")
                
                # Composition chart
                comp_df = pd.DataFrame(list(composition['composition'].items()), columns=['Nucleotide', 'Count'])
                fig = px.pie(comp_df, values='Count', names='Nucleotide', 
                           title="Nucleotide Distribution",
                           color_discrete_map={'A':'#FF6B6B', 'T':'#4ECDC4', 'G':'#45B7D1', 'C':'#96CEB4'})
                st.plotly_chart(fig, key=f'piE_distribution_{seq_idx}')
                
                # Sequence visualization
                st.subheader("🔤 Sequence Visualization")
                formatted_seq = create_sequence_visualization(sequence)
                st.markdown(f'<div class="dna-sequence">{formatted_seq}</div>', unsafe_allow_html=True)
                
                # Composition metrics
                col1, col2 = st.columns(2)
                with col1:
                    st.markdown("**Composition Metrics**")
                    st.write(f"• AT Content: {composition['at_content']:.1f}%")
                    st.write(f"• GC Content: {composition['gc_content']:.1f}%")
                    st.write(f"• CpG Islands: {composition['cpg_count']} sites")
                
                with col2:
                    st.markdown("**Sequence Quality**")
                    quality_score = 85 + (composition['gc_content'] - 50) if 40 <= composition['gc_content'] <= 60 else 70
                    st.write(f"• Quality Score: {min(99, quality_score):.0f}/100")
                    st.write(f"• Complexity: {'High' if len(set(sequence)) == 4 else 'Medium'}")
                    st.write(f"• Length Class: {'Long' if composition['length'] > 1000 else 'Short'}")

            with tab2:
                st.subheader("🎯 Gene Identification & ORF Analysis")
                if orfs:
                    st.success(f"🔍 Found {len(orfs)} Open Reading Frames")
                    
                    # ORF summary table
                    orf_data = []
                    for i, orf in enumerate(orfs[:10]):  # Show top 10
                        orf_data.append({
                            'ORF': f"ORF_{i+1}",
                            'Start': orf['start'],
                            'End': orf['end'],
                            'Length (bp)': orf['length'],
                            'Strand': orf['strand'],
                            'Frame': orf['frame'],
                            'Protein Length (aa)': orf['protein_length']
                        })
                    
                    orf_df = pd.DataFrame(orf_data)
                    st.dataframe(orf_df, use_container_width=True)
                    
                    # Visualize ORFs
                    fig = px.scatter(orf_df, x='Start', y='Protein Length (aa)', 
                                   size='Length (bp)', color='Strand',
                                   title="Open Reading Frames Distribution",
                                   hover_data=['Frame'])
                    st.plotly_chart(fig, key=f'orf_distribution_{seq_idx}')
                    
                    # Show largest ORF details
                    if orfs:
                        largest_orf = orfs[0]
                        st.subheader("🏆 Largest ORF Details")
                        
                        col1, col2 = st.columns(2)
                        with col1:
                            st.write(f"**Position:** {largest_orf['start']}-{largest_orf['end']}")
                            st.write(f"**Length:** {largest_orf['length']} bp")
                            st.write(f"**Strand:** {largest_orf['strand']}")
                            st.write(f"**Frame:** {largest_orf['frame']}")
                        
                        with col2:
                            st.write(f"**Protein Length:** {largest_orf['protein_length']} amino acids")
                            protein_mw = largest_orf['protein_length'] * 110
                            st.write(f"**Est. Molecular Weight:** {protein_mw:,} Da")
                            st.write(f"**Type:** {'Membrane protein' if protein_mw > 30000 else 'Cytoplasmic protein'}")
                
                if genes_found:
                    st.subheader("🧬 Potential Gene Matches")
                    for gene in genes_found:
                        st.markdown(f"""
                        <div class="analysis-card">
                            <h4>{gene['name']} ({gene['gene_id']})</h4>
                            <p><strong>Function:</strong> {gene['function']}</p>
                            <p><strong>Type:</strong> {gene['type']}</p>
                            <p><strong>Confidence:</strong> {gene['confidence']:.0f}%</p>
                        </div>
                        """, unsafe_allow_html=True)

            with tab3:
                st.subheader("🧪 Protein Structure Prediction")
                
                if orfs:
                    # Analyze proteins from ORFs
                    for i, orf in enumerate(orfs[:3]):  # Top 3 ORFs
                        st.markdown(f"---")
                        st.markdown(f"#### 🔬 Analysis for Protein from ORF_{i+1}")
                        
                        protein_analysis = analyzer.predict_protein_structure(orf['protein_sequence'])
                        
                        if 'error' in protein_analysis:
                            st.warning(f"Could not analyze protein from ORF_{i+1}: {protein_analysis['error']}")
                            continue

                        col1, col2 = st.columns(2)
                        
                        with col1:
                            st.metric("Protein Length", f"{protein_analysis['length']} aa")
                            st.metric("Est. Mol. Weight", f"{int(protein_analysis['molecular_weight']):,} Da")
                            
                            st.markdown("**Predicted Domains:**")
                            if protein_analysis['predicted_domains']:
                                for domain in protein_analysis['predicted_domains']:
                                    st.markdown(f'<span class="protein-structure">{domain}</span>', unsafe_allow_html=True)
                            else:
                                st.write("None predicted.")

                            st.markdown("**Secondary Structure:**")
                            summary = protein_analysis['secondary_structure_summary']
                            st.write(f"Helix: {summary['helix_percent']:.1f}%")
                            st.write(f"Sheet: {summary['sheet_percent']:.1f}%")
                            st.write(f"Coil: {summary['coil_percent']:.1f}%")

                            st.markdown("**Predicted 2D Structure Map:**")
                            structure_viz = visualize_protein_structure(protein_analysis['secondary_structure_sequence'])
                            st.markdown(structure_viz, unsafe_allow_html=True)

                        with col2:
                            st.markdown("**Interactive 3D Structure (Simulated):**")
                            pdb_string = analyzer.generate_pdb_from_structure(
                                orf['protein_sequence'], protein_analysis['secondary_structure_sequence']
                            )
                            if pdb_string:
                                render_3d_protein_structure(pdb_string)
                            else:
                                st.warning("Could not generate 3D structure.")

                        # Amino acid composition chart
                        if protein_analysis['composition']:
                            aa_df = pd.DataFrame(list(protein_analysis['composition'].items()), 
                                               columns=['Amino Acid', 'Count'])
                            fig = px.bar(aa_df.sort_values('Count', ascending=False), x='Amino Acid', y='Count', 
                                       title=f"Amino Acid Composition - ORF_{i+1}")
                            st.plotly_chart(fig, use_container_width=True, key=f'aa_comp_{seq_idx}_{i}')
                else:
                    st.info("No Open Reading Frames found to predict protein structures.")

            with tab4:
                st.subheader("🌿 Species Classification & Phylogenetics")
                
                if species_scores:
                    # Species probability
                    species_df = pd.DataFrame(list(species_scores.items()), 
                                            columns=['Species', 'Probability'])
                    species_df = species_df.sort_values('Probability', ascending=False)
                    
                    fig = px.bar(species_df, x='Species', y='Probability',
                               title="Species Classification Probability",
                               color='Probability', color_continuous_scale='viridis')
                    st.plotly_chart(fig, key=f'species_chart_{seq_idx}')
                    
                    # Top classification
                    top_species = species_df.iloc[0]
                    confidence = top_species['Probability'] * 100
                    
                    st.success(f"🎯 **Most likely species:** {top_species['Species'].title()} ({confidence:.0f}% confidence)")
                
                col1, col2 = st.columns(2)
                with col1:
                    st.metric("Conservation Score", f"{evo_analysis['conservation_score']:.0f}/100")
                    st.metric("Evolutionary Rate", evo_analysis['evolutionary_rate'])
                
                with col2:
                    if evo_analysis['mutation_hotspots']:
                        st.markdown("**Mutation Hotspots:**")
                        for hotspot in evo_analysis['mutation_hotspots'][:3]:
                            st.write(f"• {hotspot['motif']} (occurs {hotspot['frequency']} times)")

            with tab5:
                st.subheader("⚗️ Pharmacogenomics & Drug Targets")
                
                if pharma_analysis['drug_targets']:
                    st.write(f"**Pharmacogenomic Score:** {pharma_analysis['pharmacogenomic_score']}/100")
                    
                    for target in pharma_analysis['drug_targets']:
                        st.markdown(f"""
                        <div class="analysis-card">
                            <h4>{target['target']}</h4>
                            <p><strong>Relevant Drugs:</strong> {', '.join(target['drugs'])}</p>
                            <p><strong>Clinical Significance:</strong> {target['significance']}</p>
                            <p><strong>Recommendation:</strong> {target['recommendation']}</p>
                        </div>
                        """, unsafe_allow_html=True)
                else:
                    st.info("No significant pharmacogenomic markers detected in this sequence.")

            with tab6:
                st.subheader("🔬 Advanced Genomic Intelligence")
                
                # Create comprehensive analysis report
                col1, col2 = st.columns(2)
                
                with col1:
                    st.markdown("**🧬 Genomic Complexity Metrics**")
                    
                    st.write(f"• Shannon Entropy: {entropy:.3f}")
                    st.write(f"• 3-mer Diversity: {unique_3mers}")
                    st.write(f"• Max Homopolymer: {max_repeat}")
                
                with col2:
                    st.markdown("**🎯 Functional Predictions**")
                    
                    # Predict function based on composition and structure
                    predictions = []
                    
                    if gc_content > 55:
                        predictions.append("Regulatory element")
                    if len(orfs) > 0 and orfs[0]['length'] > 300:
                        predictions.append("Protein-coding gene")
                    if composition['cpg_count'] > len(sequence) * 0.02:
                        predictions.append("CpG island / Promoter region")
                    
                    if predictions:
                        for pred in predictions:
                            st.write(f"• {pred}")
                    else:
                        st.write("• Non-coding sequence")
                
                # Generate AI insights
               # st.subheader("🤖 AI-Generated Insights")
                


            with tab7:
                st.subheader("✂️ CRISPR/Cas9 Guide RNA Design")
                
                if crispr_guides:
                    st.success(f"✅ Found {len(crispr_guides)} potential gRNA candidates for gene editing.")
                    
                    guide_df = pd.DataFrame(crispr_guides)
                    
                    # Display table
                    st.markdown("#### Top gRNA Candidates (ranked by score)")
                    st.dataframe(guide_df[[
                        'sequence', 'location', 'strand', 'gc_content', 
                        'on_target_score', 'off_target_hits', 'overall_score'
                    ]].rename(columns={
                        'sequence': 'Guide Sequence (20nt)',
                        'location': 'Start Position',
                        'strand': 'Strand',
                        'gc_content': 'GC %',
                        'on_target_score': 'On-Target Score',
                        'off_target_hits': 'Off-Target Hits',
                        'overall_score': 'Overall Score'
                    }).style.format({'GC %': '{:.1f}', 'On-Target Score': '{:.1f}', 'Overall Score': '{:.1f}'}))
                    
                    # Visualization of guide locations
                    st.markdown("#### gRNA Location on Genome")
                    fig = px.scatter(guide_df.head(20), x='location', y='overall_score', 
                                   color='strand', size='on_target_score',
                                   title="Top 20 gRNA Candidates by Location and Score",
                                   hover_data=['sequence'],
                                   labels={'location': 'Position on Sequence', 'overall_score': 'Overall Score'})
                    fig.update_layout(xaxis_range=[0, len(sequence)])
                    st.plotly_chart(fig, key=f'gNa_distribution_{seq_idx}')
                    
                    st.info("💡 **On-Target Score:** A heuristic score based on GC content and sequence features. Higher is better.  \n**Off-Target Hits:** Number of identical sequences found elsewhere. Lower is better.")
                    
                else:
                    st.warning("No suitable CRISPR gRNA targets found based on current criteria (PAM: NGG).")

            # =================================================================
            # == Beta Features Section                                      ==
            # =================================================================
            if selected_feature and selected_feature != "Select a Beta Feature...":
                st.markdown("---")
                st.header(f"🔬 Beta Feature Analysis: {st.session_state.selected_beta_feature.split('.', 1)[1].strip()}")

                if selected_feature.startswith("1."):
                    # Epigenetic Analysis (Implemented)
                    if epigenetic_results:
                        st.subheader("🧬 Epigenetic Landscape")
                        st.metric("Overall Methylation Potential", epigenetic_results['overall_methylation_potential'])
                        
                        col1, col2 = st.columns(2)
                        with col1:
                            st.subheader(f"🏝️ CpG Islands ({epigenetic_results['cpg_islands_found']} found)")
                            if epigenetic_results['cpg_islands']:
                                cpg_df = pd.DataFrame(epigenetic_results['cpg_islands'])
                                st.dataframe(cpg_df.head().style.format({'gc_content': '{:.1f}%', 'obs_exp_ratio': '{:.2f}'}))
                                
                                fig = px.scatter(cpg_df, x='start', y='gc_content', size='obs_exp_ratio',
                                               title="Detected CpG Islands by Location and GC Content",
                                               hover_data=['start', 'end', 'obs_exp_ratio'],
                                               labels={'start': 'Start Position', 'gc_content': 'GC Content (%)'})
                                st.plotly_chart(fig, key=f'cpG_islands_{seq_idx}')
                            else:
                                st.info("No significant CpG islands detected.")
                        
                        with col2:
                            st.subheader(f"🧬 G-Quadruplexes ({epigenetic_results['g_quadruplexes_found']} found)")
                            if epigenetic_results['g_quadruplexes']:
                                gq_df = pd.DataFrame(epigenetic_results['g_quadruplexes'])
                                st.dataframe(gq_df[['start', 'end']].head())
                                st.markdown("**Example G-Quadruplex sequence:**")
                                st.code(gq_df.iloc[0]['sequence'])
                            else:
                                st.info("No G-Quadruplex motifs found.")
                    else:
                        st.warning("Epigenetic analysis could not be performed.")
                
                elif selected_feature.startswith("2."):
                    # Non-coding RNA Prediction
                    if ncrna_results:
                        st.subheader("🧬 Non-coding RNA (ncRNA) Prediction")
                        col1, col2 = st.columns(2)
                        with col1:
                            st.metric("Potential miRNAs found", ncrna_results['mirnas_found'])
                            if ncrna_results['mirnas']:
                                st.markdown("**Top miRNA Candidates:**")
                                mirna_df = pd.DataFrame(ncrna_results['mirnas'])
                                st.dataframe(mirna_df[['start', 'end', 'confidence']].head().style.format({'confidence': '{:.0f}%'}))
                                
                                fig_mirna = px.scatter(mirna_df, x='start', y='confidence',
                                                       title="miRNA Candidates by Location and Confidence",
                                                       hover_data=['start', 'end'],
                                                       labels={'start': 'Start Position', 'confidence': 'Confidence Score'})
                                st.plotly_chart(fig_mirna, use_container_width=True, key=f'mirna_dist_{seq_idx}')
                            else:
                                st.info("No strong miRNA candidates detected.")
                        
                        with col2:
                            st.metric("Potential lncRNAs found", ncrna_results['lncrnas_found'])
                            if ncrna_results['lncrnas']:
                                st.markdown("**Top lncRNA Candidates:**")
                                lncrna_df = pd.DataFrame(ncrna_results['lncrnas'])
                                st.dataframe(lncrna_df[['start', 'end', 'length', 'gc_content']].head().style.format({'gc_content': '{:.1f}%'}))

                                fig_lncrna = px.scatter(lncrna_df, x='start', y='length', size='gc_content',
                                                        title="lncRNA Candidates by Location and Length",
                                                        hover_data=['start', 'end', 'gc_content'],
                                                        labels={'start': 'Start Position', 'length': 'Length (bp)'})
                                st.plotly_chart(fig_lncrna, use_container_width=True, key=f'lncrna_dist_{seq_idx}')
                            else:
                                st.info("No strong lncRNA candidates detected.")
                    else:
                        st.warning("ncRNA analysis could not be performed.")

                elif selected_feature.startswith("3."):
                    # Viral Integration Site Detection
                    if viral_results:
                        st.subheader("🦠 Viral Integration Site Detection")
                        st.metric("Potential Integration Events", viral_results['integrations_found'])
                        
                        if viral_results['integrations']:
                            st.markdown("**Detected Signatures:**")
                            viral_df = pd.DataFrame(viral_results['integrations'])
                            st.dataframe(viral_df.style.format({'confidence': '{:.1f}%'}))

                            # Add a chart
                            fig_viral = px.scatter(viral_df, x='start', y='confidence',
                                                   color='virus',
                                                   title="Potential Viral Integration Sites",
                                                   hover_data=['virus', 'start', 'end', 'signature_sequence'],
                                                   labels={'start': 'Position on Sequence', 'confidence': 'Confidence (%)'})
                            fig_viral.update_layout(xaxis_range=[0, len(sequence)])
                            st.plotly_chart(fig_viral, use_container_width=True, key=f'viral_integration_{seq_idx}')
                        else:
                            st.success("No known viral signatures or retrovirus-like elements detected.")
                    else:
                        st.warning("Viral integration analysis could not be performed.")
                
                elif selected_feature.startswith("4."):
                    # Microbiome Composition Analysis
                    if microbiome_results:
                        st.subheader("🦠 Microbiome Composition Analysis (Simulated)")
                        
                        col1, col2 = st.columns(2)
                        with col1:
                            st.metric("Total Marker Hits", microbiome_results['total_marker_hits'])
                            st.metric("Dominant Phylum", microbiome_results['dominant_phylum'])
                            st.metric("Phyla Diversity Score", f"{microbiome_results['diversity_score']:.1f}%")
                        
                        with col2:
                            if microbiome_results['composition']:
                                comp_df = pd.DataFrame(list(microbiome_results['composition'].items()), columns=['Phylum', 'Percentage'])
                                fig = px.pie(comp_df, values='Percentage', names='Phylum', 
                                           title="Predicted Microbiome Composition",
                                           color_discrete_sequence=px.colors.sequential.Plasma_r)
                                st.plotly_chart(fig, use_container_width=True, key=f'microbiome_pie_{seq_idx}')
                            else:
                                st.info("No known microbiome markers detected. This may be a pure host sample.")
                    else:
                        st.warning("Microbiome analysis could not be performed.")

                elif selected_feature.startswith("5."):
                    # Structural Variant Detection
                    if sv_results:
                        st.subheader("🧬 Structural Variant (SV) Detection")
                        st.metric("Potential SVs Found", sv_results['variants_found'])
                        
                        if sv_results['variants']:
                            sv_df = pd.DataFrame(sv_results['variants'])
                            st.dataframe(sv_df.style.format({'confidence': '{:.1f}%'}))
                            
                            fig_sv = px.scatter(sv_df, x='start', y='length', color='type',
                                                title="Detected Structural Variants by Location and Size",
                                                hover_data=['start', 'end', 'type'],
                                                labels={'start': 'Start Position', 'length': 'Variant Length (bp)'})
                            fig_sv.update_layout(xaxis_range=[0, len(sequence)])
                            st.plotly_chart(fig_sv, use_container_width=True, key=f'sv_dist_{seq_idx}')
                        else:
                            st.success("No major structural variants like large tandem duplications or inversions were detected.")
                    else:
                        st.warning("Structural variant analysis could not be performed.")

                elif selected_feature.startswith("6."):
                    # Comparative Genomics
                    if comp_gen_results:
                        st.subheader("🌍 Comparative Genomics")
                        
                        if comp_gen_results['homology_hits']:
                            st.metric("Homologous Genes Found", len(comp_gen_results['homology_hits']))
                            homology_df = pd.DataFrame(comp_gen_results['homology_hits'])
                            
                            st.markdown("**Top Homology Hits:**")
                            st.dataframe(homology_df.style.format({
                                'identity': '{:.2f}%',
                                'e_value': '{:.2e}'
                            }))
                            
                            fig_homology = px.bar(homology_df, x='ortholog', y='identity', color='e_value',
                                                  title="Homology Identity to Known Orthologs",
                                                  labels={'ortholog': 'Orthologous Gene', 'identity': 'Sequence Identity (%)'},
                                                  color_continuous_scale='Cividis_r')
                            st.plotly_chart(fig_homology, use_container_width=True, key=f'homology_chart_{seq_idx}')
                        else:
                            st.info("No significant homology found against the simulated database of common orthologs.")
                    else:
                        st.warning("Comparative genomics analysis could not be performed.")
                
                elif selected_feature.startswith("7."):
                    # 3D Genome Folding Prediction
                    if folding_results:
                        st.subheader("🧊 3D Genome Folding Prediction (Simulated)")
                        col1, col2 = st.columns(2)
                        with col1:
                            st.metric("Predicted TADs", folding_results['tads_found'])
                            if folding_results['tads']:
                                st.markdown("**Largest Predicted TADs:**")
                                tad_df = pd.DataFrame(folding_results['tads'])
                                st.dataframe(tad_df.head())
                        with col2:
                            st.metric("Predicted Enhancer-Promoter Loops", folding_results['loops_found'])
                            if folding_results['loops']:
                                st.markdown("**Strongest Predicted Loops:**")
                                loop_df = pd.DataFrame(folding_results['loops'])
                                st.dataframe(loop_df[['promoter_pos', 'enhancer_pos', 'distance', 'strength']].head().style.format({'strength': '{:.2f}'}))
                        
                        if folding_results['loops']:
                            st.markdown("**Visualization of Predicted Chromatin Loops**")
                            loop_df = pd.DataFrame(folding_results['loops'])
                            fig_loops = go.Figure()
                            # Draw loops as arcs
                            for _, loop in loop_df.iterrows():
                                start, end = sorted([loop['promoter_pos'], loop['enhancer_pos']])
                                center = (start + end) / 2
                                radius = (end - start) / 2
                                theta = np.linspace(0, np.pi, 100)
                                x_arc = center + radius * np.cos(theta)
                                y_arc = radius * np.sin(theta) * 0.1 # scale height
                                fig_loops.add_trace(go.Scatter(x=x_arc, y=y_arc, mode='lines', line=dict(color='purple', width=loop['strength']*2), hoverinfo='text', hovertext=f"Interaction<br>Distance: {loop['distance']} bp<br>Strength: {loop['strength']:.2f}"))
                            fig_loops.update_layout(title="Predicted Long-Range Chromatin Interactions", xaxis_title="Genomic Position", yaxis_showticklabels=False, yaxis_zeroline=False, showlegend=False)
                            st.plotly_chart(fig_loops, use_container_width=True, key=f'3d_folding_loops_{seq_idx}')
                    else:
                        st.warning("3D folding analysis could not be performed.")

                elif selected_feature.startswith("8."):
                    # Gene Regulatory Network Inference
                    if reg_network_results:
                        st.subheader("🕸️ Gene Regulatory Network Inference (Simulated)")
                        regulations = reg_network_results['regulations']
                        st.metric("Inferred Regulatory Links", len(regulations))

                        if regulations:
                            reg_df = pd.DataFrame(regulations)
                            st.dataframe(reg_df.head())

                            # Network Graph
                            nodes = sorted(list(set(reg_df['tf'].tolist() + reg_df['target_gene'].tolist())))
                            num_nodes = len(nodes)
                            angle_step = 2 * np.pi / num_nodes
                            positions = {node: (np.cos(i * angle_step), np.sin(i * angle_step)) for i, node in enumerate(nodes)}

                            edge_x, edge_y = [], []
                            for _, row in reg_df.iterrows():
                                x0, y0 = positions[row['tf']]
                                x1, y1 = positions[row['target_gene']]
                                edge_x.extend([x0, x1, None])
                                edge_y.extend([y0, y1, None])
                            edge_trace = go.Scatter(x=edge_x, y=edge_y, line=dict(width=0.5, color='#888'), hoverinfo='none', mode='lines')

                            node_x, node_y, node_text, node_color, node_size = [], [], [], [], []
                            for node, pos in positions.items():
                                x, y = pos; node_x.append(x); node_y.append(y); node_text.append(node)
                                node_color.append('orange' if not node.startswith('ORF') else 'cyan')
                                node_size.append(30 if not node.startswith('ORF') else 20)
                            node_trace = go.Scatter(x=node_x, y=node_y, mode='markers+text', text=node_text, textposition="bottom center", hoverinfo='text', marker=dict(showscale=False, color=node_color, size=node_size, line_width=2))

                            fig_net = go.Figure(data=[edge_trace, node_trace], layout=go.Layout(title='Inferred Gene Regulatory Network', showlegend=False, hovermode='closest', margin=dict(b=20,l=5,r=5,t=40), xaxis=dict(showgrid=False, zeroline=False, showticklabels=False), yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)))
                            st.plotly_chart(fig_net, use_container_width=True, key=f'regulatory_network_{seq_idx}')
                        else:
                            st.info("No regulatory links could be inferred from the known TFBS motifs.")
                    else:
                        st.warning("Regulatory network analysis could not be performed.")

                elif selected_feature.startswith("9."):
                    # Synthetic Biology Circuit Design
                    if synth_bio_results:
                        st.subheader("⚙️ Synthetic Biology Circuit Design")
                        designs = synth_bio_results['designs']
                        
                        if designs:
                            for i, design in enumerate(designs):
                                with st.expander(f"Design {i+1}: {design['name']}", expanded=i==0):
                                    st.markdown(f"**Description:** {design['description']}")
                                    st.metric("Total Construct Length", f"{design['total_length']} bp")
                                    
                                    parts_df = pd.DataFrame(design['parts'])
                                    
                                    # Create Gantt chart data
                                    gantt_data = []
                                    current_pos = 0
                                    for _, part_row in parts_df.iterrows():
                                        gantt_data.append(dict(Task=part_row['part'], Start=current_pos, Finish=current_pos + part_row['length'], Resource=part_row['id']))
                                        current_pos += part_row['length']
                                    gantt_df = pd.DataFrame(gantt_data)
                                    
                                    fig_gantt = px.timeline(gantt_df, x_start="Start", x_end="Finish", y="Task", color="Task", title=f"Layout of '{design['name']}'", labels={"Task": "Genetic Part"}, hover_data=["Resource"])
                                    fig_gantt.update_yaxes(autorange="reversed")
                                    fig_gantt.update_layout(xaxis_title="Position (bp)")
                                    st.plotly_chart(fig_gantt, use_container_width=True, key=f'synthetic_circuit_{seq_idx}_{i}')
                        else:
                            st.info("Could not generate any synthetic circuit designs for this sequence.")
                    else:
                        st.warning("Synthetic biology analysis could not be performed.")
                
                else:
                    # Preview for other features
                    st.info("This feature is currently in preview and not yet implemented.")
                    feature_name = st.session_state.selected_beta_feature.split('.', 1)[1].strip()
                    st.markdown(f"**Coming Soon: {feature_name}**")
                    
                    preview_text = {
                        "Non-coding RNA Prediction": "predicting and classifying non-coding RNAs like miRNA, lncRNA, and circRNA from the sequence, revealing their potential regulatory roles.",
                        "Viral Integration Site Detection": "identifying potential integration sites of viral DNA into the host genome, crucial for studying viral diseases and gene therapy.",
                        "Microbiome Composition Analysis": "analyzing the sequence for markers of different microbial species to estimate the composition of a microbiome sample.",
                        "Structural Variant Detection": "detecting large-scale structural variants like deletions, duplications, inversions, and translocations within the genome.",
                        "Comparative Genomics": "aligning the sequence against a database of other species to identify conserved regions, evolutionary relationships, and functional elements.",
                        "3D Genome Folding Prediction": "predicting how the DNA sequence folds in three-dimensional space, revealing insights into gene regulation and chromatin structure.",
                        "Gene Regulatory Network Inference": "inferring networks of interacting genes and transcription factors to understand complex biological pathways.",
                        "Synthetic Biology Circuit Design": "providing tools to design and simulate synthetic gene circuits for applications in biotechnology and medicine.",
                        "AI-Powered Drug Discovery": "using machine learning models to predict potential drug targets and design novel therapeutic molecules based on the genomic data."
                    }
                    st.write(f"This analysis will provide deep insights into {preview_text.get(feature_name, 'this area of genomics')}")

        # Generate comprehensive report
        st.markdown("---")
        st.header("📊 Comprehensive Analysis Report")
        
        if st.button("📄 Generate Complete Genomic Report"):
            # Compile all analyses
            report = {
                'timestamp': datetime.now().isoformat(),
                'sequence_info': {
                    'length': len(sequence),
                    'gc_content': gc_content,
                    'class_label': class_label
                },
                'composition_analysis': composition,
                'orfs_found': len(orfs),
                'genes_identified': len(genes_found),
                'species_classification': species_scores,
                'evolutionary_analysis': evo_analysis,
                'pharmacogenomics': pharma_analysis,
                'crispr_guides_found': len(crispr_guides)
            }
            
            report_json = json.dumps(report, indent=2)
            
            st.download_button(
                label="📥 Download JSON Report",
                data=report_json,
                file_name=f"eternaseq_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
                mime="application/json"
            )
            
            # Display summary
            st.success("✅ Analysis Complete!")
            st.markdown(f"""
            **🎯 Executive Summary:**
            
            • **Sequence Length:** {len(sequence):,} nucleotides
            • **ORFs Detected:** {len(orfs)} potential protein-coding regions
            • **Genes Identified:** {len(genes_found)} candidate genes
            • **Species Classification:** {list(species_scores.keys())[0].title() if species_scores else 'Unknown'}
            • **Conservation Score:** {evo_analysis['conservation_score']:.0f}/100
            • **Clinical Relevance:** {pharma_analysis['pharmacogenomic_score']}/100
            • **CRISPR Targets:** {len(crispr_guides)} potential gRNA sites identified
            
            This comprehensive analysis provides unprecedented insights into the genomic sequence,
            combining traditional bioinformatics with AI-powered pattern recognition for 
            revolutionary biological understanding.
            """)

if __name__ == "__main__":
    main()
