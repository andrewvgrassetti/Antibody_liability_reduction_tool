# Antibody Liability Reduction Tool

A Python pipeline for systematically identifying and reducing surface-exposed liabilities (hydrophobic patches, positive charge clusters) in antibody VH sequences while maintaining humanness and stability.

## Overview

This tool takes an antibody VH sequence with high polyreactivity and:

1. Numbers residues (ANARCI/IMGT) and classifies surface exposure
2. Identifies surface-exposed hydrophobic and positively charged residues
3. Generates humanness-constrained mutations using Abysis frequency data
4. Evaluates candidates through TAP (surface patches), DeepSP (stability), and OASis (humanness)
5. Expands validated point mutations into double/triple combinations via Bayesian optimization
6. Re-evaluates and ranks all combinatorial candidates

## Status

🚧 Under development
