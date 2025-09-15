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
    
    # FIX: Ensure the generated HTML is a valid string before rendering.
    # This prevents a TypeError if py3Dmol returns None or an empty object.
    if html_output:
        components.html(html_output, height=int(height), width=int(width))
    else:
        st.warning("Could not generate a 3D visualization for this structure.")

def get_gemini_insights(api_key: str, analysis_summary: Dict) -> List[str]:
    """Generate insights using Gemini AI"""
    if not api_key:
        return ["Please enter a Google AI API Key in the sidebar to generate AI insights."]
    
    try:
        genai.configure(api_key=api_key)
        model = genai.GenerativeModel('gemini-1.5-flash-latest')
        
        # Create a detailed prompt
        prompt = f"""
        You are a world-class geneticist and bioinformatician. Analyze the following DNA sequence summary and provide 4-5 deep, insightful conclusions for a genomics researcher. Be concise, use scientific language, and format your output as a bulleted list (using '*' for bullets).

        **Sequence Summary:**
        - **Length:** {analysis_summary['length']} bp
        - **GC Content:** {analysis_summary['gc_content']:.1f}%
        - **Identified ORFs:** {analysis_summary['orfs_count']}
        - **Largest ORF:** {analysis_summary['largest_orf_bp']} bp (protein: {analysis_summary['largest_orf_aa']} aa)
        - **Potential Genes Found:** {analysis_summary['genes_found']}
        - **Predicted Species:** {analysis_summary['predicted_species']} ({analysis_summary['species_confidence']:.0f}% confidence)
        - **Evolutionary Conservation Score:** {analysis_summary['conservation_score']:.0f}/100 (Rate: {analysis_summary['evolutionary_rate']})
        - **Pharmacogenomic Relevance:** {analysis_summary['pharma_score']}/100, Targets: {analysis_summary['pharma_targets']}
        - **Shannon Entropy:** {analysis_summary['entropy']:.3f}

        **Your Insights (as a bulleted list):**
        """
        
        response = model.generate_content(prompt)
        
        # Clean up the response
        insights = response.text.strip().split('\n')
        cleaned_insights = [re.sub(r'^\s*[\*\-]\s*', '', i).strip() for i in insights if i.strip()]
        return cleaned_insights

    except Exception as e:
        st.error(f"Gemini API Error: {e}")
        return ["Could not generate AI insights. Please check your API key and network connection."]

def main():
    # App header
    st.markdown("""
    <div class="main-header">
        <h1>🧬 EternaSeq</h1>
        <h2>The Nobel Prize-Winning Genomic Discovery Engine</h2>
        <p>Unlock groundbreaking biological discoveries from any DNA sequence. Your next Nobel Prize starts here.</p>
    </div>
    """, unsafe_allow_html=True)

    # Sidebar for API Key
    st.sidebar.header("💎 Gemini AI Integration")
    api_key = st.sidebar.text_input("Enter your Google AI API Key", type="password", help="Get your key from https://aistudio.google.com/app/apikey")

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
                st.subheader("🤖 AI-Generated Insights")
                
                with st.spinner("🧠 Consulting Gemini for deep insights..."):
                    # Prepare summary for Gemini
                    top_species_name = "Unknown"
                    top_species_prob = 0
                    if species_scores:
                        top_species_df = pd.DataFrame(list(species_scores.items()), columns=['Species', 'Probability']).sort_values('Probability', ascending=False)
                        if not top_species_df.empty:
                            top_species_name = top_species_df.iloc[0]['Species'].title()
                            top_species_prob = top_species_df.iloc[0]['Probability'] * 100

                    analysis_summary = {
                        'length': len(sequence),
                        'gc_content': gc_content,
                        'orfs_count': len(orfs),
                        'largest_orf_bp': orfs[0]['length'] if orfs else 0,
                        'largest_orf_aa': orfs[0]['protein_length'] if orfs else 0,
                        'genes_found': ', '.join([g['name'] for g in genes_found]) if genes_found else 'None',
                        'predicted_species': top_species_name,
                        'species_confidence': top_species_prob,
                        'conservation_score': evo_analysis['conservation_score'],
                        'evolutionary_rate': evo_analysis['evolutionary_rate'],
                        'pharma_score': pharma_analysis['pharmacogenomic_score'],
                        'pharma_targets': ', '.join([t['target'] for t in pharma_analysis['drug_targets']]) if pharma_analysis['drug_targets'] else 'None',
                        'entropy': entropy,
                    }
                    
                    insights = get_gemini_insights(api_key, analysis_summary)
                    for insight in insights:
                        st.info(f"💡 {insight}")

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
                    st.plotly_chart(fig, use_container_width=True)
                    
                    st.info("💡 **On-Target Score:** A heuristic score based on GC content and sequence features. Higher is better.  \n**Off-Target Hits:** Number of identical sequences found elsewhere. Lower is better.")
                    
                else:
                    st.warning("No suitable CRISPR gRNA targets found based on current criteria (PAM: NGG).")

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
