# Sequence Alignment Using Edit Distance Algorithms
This repository contains implementations of two of the most standard sequence alignment algorithms used in bioinformatics:
- **Needleman-Wunsch** (Global Alignment)
- **Smith-Waterman** (Local Alignment) 

## Objectives:
The `objective` of this repo is to:
- Understand and implement **global and local alignment techniques** using dynamic programming.
- Explore the principles of **minimum edit distances** in the context of DNA sequence alignment.
- **Analyze the impact** of scoring parameters (match, mismatch, gap penalties) on the alignment results.

## Features
- **Customizable scoring**: Set your own match score, mismatch penalty, and gap penalty.  
- **Visual output**: Heatmaps of the scoring matrix with the optimal alignment path plotted.  
- **Execution time**: Measures the time taken to compute the alignment. 

## Usage .py
After installation of the required python packages [described on the `requirements.txt`], run the program **alignment.py**.

You will be asked to:
- Choose the algorithm: **Needleman-Wunsch** [0] or **Smith-Waterman** [1].
- Input the **sequences** to align/compare.
- Input the **scoring parameters**: match score, mismatch penalty, and gap penalty

The program will return the `aligned sequences`, `score`, `computation time` and `heatmap`.

## Usage jupyter-notebook
For personal use, we recommend running the .py file. This notebook was created mainly to visualize the 40+ graphs easily, without saving all images to a folder.

    ><(((('>
