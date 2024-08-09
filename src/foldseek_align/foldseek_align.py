#!/usr/bin/env python3
"""
This module provides functionality to perform Foldseek analysis in two modes:
all-vs-all and pdb-vs-pdb.

It includes a FoldseekAnalyzer class that handles the entire process of running Foldseek,
creating pairwise FASTA files, and generating fident matrices or scores based on the mode.
"""

import os
import subprocess
import csv
import numpy as np
from typing import List, Dict, Tuple
import logging
from multiprocessing import Pool, cpu_count
import argparse
import sys
import shutil

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


class FoldseekAnalyzer:
    """
    A class to perform Foldseek analysis in all-vs-all or pdb-vs-pdb mode.

    Attributes:
        input_dir (str): Directory containing input PDB files (for all-vs-all mode).
        output_dir (str): Output directory for all results.
        cpu_cores (int): Number of CPU cores to use for analysis.
        foldseek_path (str): Path to the Foldseek executable.
    """

    def __init__(self, input_dir: str, output_dir: str, cpu_percentage: float = 25):
        """
        Initialize the FoldseekAnalyzer.

        Args:
            input_dir (str): Directory containing input PDB files (for all-vs-all mode).
            output_dir (str): Output directory for all results.
            cpu_percentage (float): Percentage of CPU cores to use (default: 25).
        """
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.cpu_cores = max(1, int(cpu_count() * (cpu_percentage / 100)))
        os.makedirs(self.output_dir, exist_ok=True)

        self.foldseek_path = shutil.which("foldseek")
        if not self.foldseek_path:
            raise FileNotFoundError("Foldseek not found in system path. Please ensure it's installed via conda.")

    def run_all_vs_all_analysis(self, create_matrix: bool = False, generate_pairwise: bool = True) -> None:
        """
        Run the complete Foldseek all-vs-all analysis pipeline.

        Args:
            create_matrix (bool): Whether to create fident matrix (default: False).
            generate_pairwise (bool): Whether to generate pairwise alignments (default: True).
        """
        results_file = os.path.join(self.output_dir, "foldseek_results.tsv")
        try:
            self._run_foldseek_all_vs_all(results_file)
            if create_matrix:
                self._create_fident_matrix(results_file)
            if generate_pairwise:
                pairwise_dir = os.path.join(self.output_dir, 'pairwise_fastas')
                os.makedirs(pairwise_dir, exist_ok=True)
                self._create_pairwise_fastas(results_file, pairwise_dir)
        except Exception as e:
            logging.error(f"An error occurred during all-vs-all analysis: {str(e)}")
            raise

    def run_pdb_vs_pdb_analysis(self, input_pdb1: str, input_pdb2: str, output_fasta: str,
                                output_fident: bool = False, generate_pairwise: bool = True) -> None:
        """
        Run the Foldseek pdb-vs-pdb analysis pipeline.

        Args:
            input_pdb1 (str): Path to the first input PDB file.
            input_pdb2 (str): Path to the second input PDB file.
            output_fasta (str): Path to the output FASTA file.
            output_fident (bool): Whether to output fident score (default: False).
            generate_pairwise (bool): Whether to generate pairwise alignment (default: True).
        """
        try:
            self._run_foldseek_pdb_vs_pdb(input_pdb1, input_pdb2, output_fasta)
            if output_fident:
                self._save_fident_score(input_pdb1, input_pdb2)
            if generate_pairwise:
                pairwise_dir = os.path.dirname(output_fasta)
                self._create_pairwise_fastas_pdb_vs_pdb(input_pdb1, input_pdb2, pairwise_dir)
        except Exception as e:
            logging.error(f"An error occurred during pdb-vs-pdb analysis: {str(e)}")
            raise

    def _run_foldseek_all_vs_all(self, results_file: str) -> None:
        """
        Run Foldseek all-vs-all search.

        Args:
            results_file (str): Path to save the Foldseek results.
        """
        cmd = [
            self.foldseek_path,
            "easy-search",
            self.input_dir,
            self.input_dir,
            results_file,
            f"{self.output_dir}/tmp",
            "--format-output",
            "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,taln",
            "-a",
            "--threads",
            str(self.cpu_cores)
        ]
        try:
            subprocess.run(cmd, check=True, stderr=subprocess.PIPE, text=True)
            logging.info(f"Foldseek all-vs-all search completed. Results saved in {results_file}")
        except subprocess.CalledProcessError as e:
            logging.error(f"Error running Foldseek: {e.stderr}")
            raise

    def _run_foldseek_pdb_vs_pdb(self, input_pdb1: str, input_pdb2: str, output_fasta: str) -> None:
        """
        Run Foldseek pdb-vs-pdb comparison.

        Args:
            input_pdb1 (str): Path to the first input PDB file.
            input_pdb2 (str): Path to the second input PDB file.
            output_fasta (str): Path to save the output FASTA file.
        """
        cmd = [
            self.foldseek_path,
            "easy-search",
            input_pdb1,
            input_pdb2,
            output_fasta,
            f"{self.output_dir}/tmp",
            "--format-output",
            "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,taln",
            "-a",
            "--threads",
            str(self.cpu_cores)
        ]
        try:
            subprocess.run(cmd, check=True, stderr=subprocess.PIPE, text=True)
            logging.info(f"Foldseek pdb-vs-pdb comparison completed. Results saved in {output_fasta}")
        except subprocess.CalledProcessError as e:
            logging.error(f"Error running Foldseek: {e.stderr}")
            raise

    def _create_fident_matrix(self, results_file: str) -> None:
        """
        Create a fident matrix from Foldseek results.

        Args:
            results_file (str): Path to the Foldseek results file.
        """
        fident_scores = self._calculate_fident_scores(results_file)
        protein_names = sorted(list(set(k[0] for k in fident_scores.keys())))
        proteins = [os.path.splitext(os.path.basename(name))[0] for name in protein_names]
        matrix_size = len(proteins)

        fident_matrix = np.ones((matrix_size, matrix_size))

        for i, p1 in enumerate(proteins):
            for j, p2 in enumerate(proteins):
                if i != j:
                    fident_matrix[i, j] = fident_scores.get((protein_names[i], protein_names[j]), 0.0)

        fident_matrix_file = os.path.join(self.output_dir, 'fident_matrix.csv')
        self._write_matrix_to_csv(fident_matrix_file, fident_matrix, proteins)
        logging.info(f"Fident score matrix created: {fident_matrix_file}")

    def _write_matrix_to_csv(self, file_path: str, matrix: np.ndarray, proteins: List[str]) -> None:
        """
        Write a matrix to a CSV file.

        Args:
            file_path (str): Path to the output CSV file.
            matrix (np.ndarray): The matrix to write.
            proteins (List[str]): List of protein names.
        """
        with open(file_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([''] + proteins)
            for i, row in enumerate(matrix):
                writer.writerow([proteins[i]] + list(row))

    def _calculate_fident_scores(self, results_file: str) -> Dict[Tuple[str, str], float]:
        """
        Calculate fident scores from Foldseek results.

        Args:
            results_file (str): Path to the Foldseek results file.

        Returns:
            Dict[Tuple[str, str], float]: Dictionary of pairwise fident scores.
        """
        fident_scores = {}
        with open(results_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                query, target, fident = row[:3]
                fident_score = float(fident)
                fident_scores[(query, target)] = fident_score
                fident_scores[(target, query)] = fident_score
        return fident_scores

    def _create_pairwise_fastas(self, results_file: str, pairwise_dir: str) -> None:
        """
        Create pairwise FASTA files from Foldseek results.

        Args:
            results_file (str): Path to the Foldseek results file.
            pairwise_dir (str): Directory to store pairwise FASTA files.
        """
        with open(results_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            with Pool(processes=self.cpu_cores) as pool:
                pool.starmap(self._create_single_fasta,
                             [(row[0], row[1], row[12], row[13], pairwise_dir) for row in reader])
        logging.info(f"Pairwise FASTA files created in {pairwise_dir}")

    def _create_single_fasta(self, query: str, target: str, qaln: str, taln: str, pairwise_dir: str) -> None:
        """
        Create a single pairwise FASTA file.

        Args:
            query (str): Query sequence name.
            target (str): Target sequence name.
            qaln (str): Query aligned sequence.
            taln (str): Target aligned sequence.
            pairwise_dir (str): Directory to store the FASTA file.
        """
        query_base = os.path.splitext(os.path.basename(query))[0]
        target_base = os.path.splitext(os.path.basename(target))[0]

        fasta_filename = f"{query_base}_{target_base}.fasta"
        fasta_path = os.path.join(pairwise_dir, fasta_filename)

        with open(fasta_path, 'w') as fasta_file:
            fasta_file.write(f">{query_base}\n{qaln}\n")
            fasta_file.write(f">{target_base}\n{taln}\n")

    def _create_pairwise_fastas_pdb_vs_pdb(self, input_pdb1: str, input_pdb2: str, pairwise_dir: str) -> None:
        """
        Create pairwise FASTA file for pdb-vs-pdb comparison.

        Args:
            input_pdb1 (str): Path to the first input PDB file.
            input_pdb2 (str): Path to the second input PDB file.
            pairwise_dir (str): Directory to store the pairwise FASTA file.
        """
        # Implementation similar to _create_single_fasta, but for pdb-vs-pdb mode
        # You'll need to extract the aligned sequences from the Foldseek output
        pass

    def _save_fident_score(self, input_pdb1: str, input_pdb2: str) -> None:
        """
        Save fident score for pdb-vs-pdb comparison.

        Args:
            input_pdb1 (str): Path to the first input PDB file.
            input_pdb2 (str): Path to the second input PDB file.
        """
        # Implementation to extract and save the fident score
        # You'll need to parse the Foldseek output to get the fident score
        pass


def main():
    parser = argparse.ArgumentParser(
        description="""
        Perform Foldseek analysis in two modes: all_vs_all or pdb_vs_pdb.

        This script runs Foldseek on PDB files, either performing all-vs-all alignment
        or comparing two specific PDB files. It can generate fident matrices, pairwise
        FASTA files, and fident scores based on the chosen mode and options.

        Foldseek must be installed and available in the system path.
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("--mode",
                        choices=["all_vs_all", "pdb_vs_pdb"],
                        required=True,
                        help="Mode of operation: 'all_vs_all' or 'pdb_vs_pdb'")

    parser.add_argument("--cpu_percentage",
                        type=float,
                        default=25,
                        help="Percentage of CPU cores to use for the analysis (default: 25). "
                             "Increase for faster processing, decrease if you're running out of memory.")

    # All vs All mode arguments
    all_vs_all_group = parser.add_argument_group("All vs All mode arguments")
    all_vs_all_group.add_argument("--input_dir",
                                  help="Directory containing input PDB files for all_vs_all mode.")
    all_vs_all_group.add_argument("--output_dir",
                                  help="Output directory for all_vs_all mode results.")
    all_vs_all_group.add_argument("--create_matrix",
                                  action="store_true",
                                  help="Create fident matrix in all_vs_all mode.")

    # PDB vs PDB mode arguments
    pdb_vs_pdb_group = parser.add_argument_group("PDB vs PDB mode arguments")
    pdb_vs_pdb_group.add_argument("--input_pdb1",
                                  help="First input PDB file for pdb_vs_pdb mode.")
    pdb_vs_pdb_group.add_argument("--input_pdb2",
                                  help="Second input PDB file for pdb_vs_pdb mode.")
    pdb_vs_pdb_group.add_argument("--output_fasta",
                                  help="Output FASTA file for pdb_vs_pdb mode.")
    pdb_vs_pdb_group.add_argument("--output_fident",
                                  action="store_true",
                                  help="Save fident score as a txt file in pdb_vs_pdb mode.")

    # Common arguments
    parser.add_argument("--pairwise",
                        action="store_true",
                        default=True,
                        help="Generate pairwise alignment FASTA files (default: True).")

    args = parser.parse_args()

    if args.mode == "all_vs_all":
        if not args.input_dir or not args.output_dir:
            parser.error("all_vs_all mode requires --input_dir and --output_dir")

        analyzer = FoldseekAnalyzer(args.input_dir, args.output_dir, args.cpu_percentage)
        analyzer.run_all_vs_all_analysis(create_matrix=args.create_matrix, generate_pairwise=args.pairwise)

    elif args.mode == "pdb_vs_pdb":
        if not args.input_pdb1 or not args.input_pdb2 or not args.output_fasta:
            parser.error("pdb_vs_pdb mode requires --input_pdb1, --input_pdb2, and --output_fasta")

        analyzer = FoldseekAnalyzer(None, os.path.dirname(args.output_fasta), args.cpu_percentage)
        analyzer.run_pdb_vs_pdb_analysis(args.input_pdb1, args.input_pdb2, args.output_fasta,
                                         output_fident=args.output_fident, generate_pairwise=args.pairwise)


if __name__ == "__main__":
    main()