#!/usr/bin/env python3
"""
Wrapper for MAGeCK RRA and MLE analysis.

This wrapper interfaces with the official MAGeCK tool (Li et al., Genome Biology 2014)
to perform:
- RRA (Robust Rank Aggregation) - pairwise comparisons
- MLE (Maximum Likelihood Estimation) - multi-condition comparisons

MLE is particularly useful for comparing multiple conditions simultaneously
against a control group.

MAGeCK invocation is resolved robustly (see ``resolve_mageck_invocation``):
1. An explicit binary path (``mageck_binary``), if given and present.
2. ``mageck`` found on PATH.
3. ``conda run -n <env> mageck`` as a fallback, if a conda env is given.
"""

import shutil
import subprocess
import pandas as pd
from pathlib import Path
from typing import List, Optional, Dict


class MAGeCKNotFoundError(RuntimeError):
    """Raised when no usable MAGeCK installation can be located."""


def resolve_mageck_invocation(conda_env: Optional[str] = None,
                              mageck_binary: Optional[str] = None) -> List[str]:
    """Return the command prefix used to invoke MAGeCK, failing early if absent.

    Resolution order: explicit binary -> mageck on PATH -> conda run -n env.
    Raises MAGeCKNotFoundError with a clear message if none work.
    """
    # 1. Explicit binary path
    if mageck_binary:
        if Path(mageck_binary).exists() or shutil.which(mageck_binary):
            return [mageck_binary]
        raise MAGeCKNotFoundError(f"Specified MAGeCK binary not found: {mageck_binary}")

    # 2. mageck on PATH
    found = shutil.which("mageck")
    if found:
        return [found]

    # 3. conda env fallback
    if conda_env:
        if shutil.which("conda"):
            return ["conda", "run", "-n", conda_env, "mageck"]
        raise MAGeCKNotFoundError(
            "MAGeCK is not on PATH and 'conda' was not found to use env "
            f"'{conda_env}'. Install MAGeCK or point --mageck-binary at it."
        )

    raise MAGeCKNotFoundError(
        "MAGeCK not found on PATH. Install it (e.g. `conda create -n mageck "
        "-c bioconda mageck`) and pass --conda-env mageck, or point "
        "--mageck-binary at the executable."
    )


def get_mageck_version(invocation: List[str]) -> str:
    """Return the MAGeCK version string for a resolved invocation, or 'unknown'."""
    try:
        result = subprocess.run(
            invocation + ["--version"],
            capture_output=True, text=True, timeout=60,
        )
        out = (result.stdout or result.stderr or "").strip()
        return out.splitlines()[0] if out else "unknown"
    except Exception:
        return "unknown"


class MAGeCKWrapper:
    """
    Wrapper for running MAGeCK RRA / MLE analysis.

    Uses the official MAGeCK implementation, located via
    ``resolve_mageck_invocation``.
    """

    def __init__(self, conda_env: Optional[str] = None,
                 mageck_binary: Optional[str] = None,
                 verify: bool = True):
        """
        Args:
            conda_env: Name of conda environment with MAGeCK installed (fallback).
            mageck_binary: Explicit path to the mageck executable (highest priority).
            verify: If True, resolve and version-check MAGeCK now and fail early.
        """
        self.conda_env = conda_env
        self.mageck_binary = mageck_binary
        self._invocation = None
        if verify:
            self._invocation = resolve_mageck_invocation(conda_env, mageck_binary)
            version = get_mageck_version(self._invocation)
            print(f"  Using MAGeCK: {' '.join(self._invocation)} (version: {version})")

    @property
    def invocation(self) -> List[str]:
        if self._invocation is None:
            self._invocation = resolve_mageck_invocation(self.conda_env, self.mageck_binary)
        return self._invocation

    def run_test(
        self,
        count_matrix: str,
        control_samples: List[str],
        treatment_samples: List[str],
        output_prefix: str,
        normalization: str = "median",
        control_sgrna: Optional[str] = None,
        fdr_threshold: float = 0.25
    ) -> dict:
        """
        Run MAGeCK test (RRA analysis).

        Args:
            count_matrix: Path to count matrix file (tab-separated)
            control_samples: List of control sample column names
            treatment_samples: List of treatment sample column names
            output_prefix: Output file prefix
            normalization: Normalization method ('median', 'total', 'control')
            control_sgrna: Optional path to control sgRNA list for normalization
            fdr_threshold: FDR threshold for calling significance (default: 0.25)

        Returns:
            Dictionary with paths to output files
        """
        cmd = self.invocation + [
            "test",
            "-k", count_matrix,
            "-t", ",".join(treatment_samples),
            "-c", ",".join(control_samples),
            "-n", output_prefix,
            "--norm-method", normalization,
            "--adjust-method", "fdr",
        ]

        if control_sgrna:
            cmd.extend(["--control-sgrna", control_sgrna])

        cmd.extend(["--gene-test-fdr-threshold", str(fdr_threshold)])

        print(f"\nRunning MAGeCK test...")
        print(f"  Command: {' '.join(cmd)}")

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            print(result.stdout)
            return {
                'gene_summary': f"{output_prefix}.gene_summary.txt",
                'sgrna_summary': f"{output_prefix}.sgrna_summary.txt",
                'log': f"{output_prefix}.log",
            }
        except subprocess.CalledProcessError as e:
            print(f"ERROR running MAGeCK:")
            print(f"  Return code: {e.returncode}")
            print(f"  STDOUT: {e.stdout}")
            print(f"  STDERR: {e.stderr}")
            raise

    def parse_gene_results(
        self,
        gene_summary_file: str,
        direction: str = "both",
        fdr_threshold: float = 0.05
    ) -> pd.DataFrame:
        """
        Parse MAGeCK gene summary file.

        Args:
            gene_summary_file: Path to .gene_summary.txt file
            direction: Which hits to return ('neg', 'pos', or 'both')
            fdr_threshold: FDR cutoff used to flag the Significance column

        Returns:
            DataFrame with gene results
        """
        gene_df = pd.read_csv(gene_summary_file, sep='\t')

        gene_df = gene_df.rename(columns={
            'id': 'Gene',
            'num': 'Num_sgRNAs',
            'neg|score': 'RRA_Score_Neg',
            'neg|p-value': 'PValue_Neg',
            'neg|fdr': 'FDR_Neg',
            'neg|rank': 'Rank_Neg',
            'pos|score': 'RRA_Score_Pos',
            'pos|p-value': 'PValue_Pos',
            'pos|fdr': 'FDR_Pos',
            'pos|rank': 'Rank_Pos',
            'neg|lfc': 'LFC_Neg',
            'pos|lfc': 'LFC_Pos'
        })

        if direction == "neg":
            gene_df['RRA_Score'] = gene_df['RRA_Score_Neg']
            gene_df['FDR'] = gene_df['FDR_Neg']
            gene_df['logFC'] = gene_df['LFC_Neg']
            gene_df['Direction'] = 'depleted'
        elif direction == "pos":
            gene_df['RRA_Score'] = gene_df['RRA_Score_Pos']
            gene_df['FDR'] = gene_df['FDR_Pos']
            gene_df['logFC'] = gene_df['LFC_Pos']
            gene_df['Direction'] = 'enriched'
        else:  # both
            gene_df['FDR'] = gene_df[['FDR_Neg', 'FDR_Pos']].min(axis=1)
            gene_df['Direction'] = gene_df.apply(
                lambda row: 'depleted' if row['FDR_Neg'] < row['FDR_Pos'] else 'enriched',
                axis=1
            )
            gene_df['logFC'] = gene_df.apply(
                lambda row: row['LFC_Neg'] if row['Direction'] == 'depleted' else row['LFC_Pos'],
                axis=1
            )
            gene_df['RRA_Score'] = gene_df.apply(
                lambda row: row['RRA_Score_Neg'] if row['Direction'] == 'depleted' else row['RRA_Score_Pos'],
                axis=1
            )

        gene_df['Significance'] = gene_df['FDR'].apply(
            lambda x: 'significant' if x <= fdr_threshold else 'not_significant'
        )

        return gene_df

    def parse_sgrna_results(self, sgrna_summary_file: str) -> pd.DataFrame:
        """Parse MAGeCK sgRNA summary file."""
        sgrna_df = pd.read_csv(sgrna_summary_file, sep='\t')
        sgrna_df = sgrna_df.rename(columns={
            'sgrna': 'sgRNA',
            'Gene': 'Gene',
            'control_count': 'Control_Mean',
            'treatment_count': 'Treatment_Mean',
            'LFC': 'logFC',
            'control_var': 'Control_Var',
            'treatment_var': 'Treatment_Var',
            'p.low': 'PValue_Low',
            'p.high': 'PValue_High',
            'p.twosided': 'PValue',
            'FDR': 'FDR'
        })
        return sgrna_df

    def generate_design_matrix(
        self,
        control_samples: List[str],
        condition_groups: Dict[str, List[str]],
        output_file: str
    ) -> str:
        """
        Generate design matrix for MAGeCK MLE.

        Args:
            control_samples: List of control sample names
            condition_groups: Dict mapping condition names to sample lists
            output_file: Path to save design matrix

        Returns:
            Path to design matrix file
        """
        all_samples = control_samples.copy()
        for condition_samples in condition_groups.values():
            all_samples.extend(condition_samples)

        design_data = []
        for sample in all_samples:
            row = {'Samples': sample, 'baseline': 1}
            for condition_name, condition_samples in condition_groups.items():
                row[condition_name] = 1 if sample in condition_samples else 0
            design_data.append(row)

        design_df = pd.DataFrame(design_data)
        column_order = ['Samples', 'baseline'] + list(condition_groups.keys())
        design_df = design_df[column_order]
        design_df.to_csv(output_file, sep='\t', index=False)

        print(f"\n  Design matrix saved: {output_file}")
        print(f"  Samples: {len(all_samples)}")
        print(f"  Conditions: {len(condition_groups)} (+ baseline)")

        return output_file

    def run_mle(
        self,
        count_matrix: str,
        design_matrix: str,
        output_prefix: str,
        normalization: str = "median",
        control_sgrna: Optional[str] = None
    ) -> dict:
        """
        Run MAGeCK MLE (Maximum Likelihood Estimation).

        Args:
            count_matrix: Path to count matrix file (tab-separated)
            design_matrix: Path to design matrix file
            output_prefix: Output file prefix
            normalization: Normalization method ('median', 'total', 'control')
            control_sgrna: Optional path to control sgRNA list for normalization

        Returns:
            Dictionary with paths to output files
        """
        cmd = self.invocation + [
            "mle",
            "-k", count_matrix,
            "-d", design_matrix,
            "-n", output_prefix,
            "--norm-method", normalization,
        ]

        if control_sgrna:
            cmd.extend(["--control-sgrna", control_sgrna])

        print(f"\nRunning MAGeCK MLE...")
        print(f"  Command: {' '.join(cmd)}")

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            print(result.stdout)
            return {
                'gene_summary': f"{output_prefix}.gene_summary.txt",
                'sgrna_summary': f"{output_prefix}.sgrna_summary.txt",
                'log': f"{output_prefix}.log",
            }
        except subprocess.CalledProcessError as e:
            print(f"ERROR running MAGeCK MLE:")
            print(f"  Return code: {e.returncode}")
            print(f"  STDOUT: {e.stdout}")
            print(f"  STDERR: {e.stderr}")
            raise

    def parse_mle_gene_results(
        self,
        gene_summary_file: str,
        fdr_threshold: float = 0.05
    ) -> pd.DataFrame:
        """
        Parse MAGeCK MLE gene summary file.

        Args:
            gene_summary_file: Path to .gene_summary.txt file
            fdr_threshold: FDR threshold for significance
        """
        gene_df = pd.read_csv(gene_summary_file, sep='\t')

        if 'gene' in gene_df.columns:
            gene_df = gene_df.rename(columns={'gene': 'Gene'})

        fdr_cols = [col for col in gene_df.columns if col.endswith('|fdr')]
        for fdr_col in fdr_cols:
            condition = fdr_col.replace('|fdr', '')
            sig_col = f"{condition}|significant"
            gene_df[sig_col] = gene_df[fdr_col].apply(
                lambda x: 'significant' if x <= fdr_threshold else 'not_significant'
            )

        return gene_df

    def parse_mle_sgrna_results(self, sgrna_summary_file: str) -> pd.DataFrame:
        """Parse MAGeCK MLE sgRNA summary file."""
        return pd.read_csv(sgrna_summary_file, sep='\t')


def run_celltype_comparison(
    count_matrix_file: str,
    control_samples: List[str],
    treatment_samples: List[str],
    output_prefix: str,
    conda_env: Optional[str] = None,
    mageck_binary: Optional[str] = None,
    normalization: str = "median",
    control_sgrna: Optional[str] = None,
    fdr_threshold: float = 0.05
) -> dict:
    """
    Convenience function to run MAGeCK test for a single comparison.

    Args:
        count_matrix_file: Path to count matrix
        control_samples: List of control sample names
        treatment_samples: List of treatment sample names
        output_prefix: Output file prefix
        conda_env: Conda environment name (fallback for locating MAGeCK)
        mageck_binary: Explicit path to mageck executable
        normalization: Normalization method
        control_sgrna: Optional control sgRNA list for control normalization
        fdr_threshold: FDR threshold used both for MAGeCK and for flagging hits

    Returns:
        Dictionary with parsed results
    """
    wrapper = MAGeCKWrapper(conda_env=conda_env, mageck_binary=mageck_binary)

    output_files = wrapper.run_test(
        count_matrix=count_matrix_file,
        control_samples=control_samples,
        treatment_samples=treatment_samples,
        output_prefix=output_prefix,
        normalization=normalization,
        control_sgrna=control_sgrna,
        fdr_threshold=fdr_threshold,
    )

    gene_results = wrapper.parse_gene_results(
        output_files['gene_summary'], fdr_threshold=fdr_threshold)
    sgrna_results = wrapper.parse_sgrna_results(output_files['sgrna_summary'])

    return {
        'gene_results': gene_results,
        'sgrna_results': sgrna_results,
        'output_files': output_files
    }


def run_mle_multi_condition(
    count_matrix_file: str,
    control_samples: List[str],
    condition_groups: Dict[str, List[str]],
    output_prefix: str,
    conda_env: Optional[str] = None,
    mageck_binary: Optional[str] = None,
    normalization: str = "median",
    control_sgrna: Optional[str] = None,
    fdr_threshold: float = 0.05
) -> dict:
    """
    Run MAGeCK MLE for multiple conditions with auto-generated design matrix.

    Args:
        count_matrix_file: Path to count matrix
        control_samples: List of control sample names
        condition_groups: Dict mapping condition names to sample lists
        output_prefix: Output file prefix
        conda_env: Conda environment name (fallback for locating MAGeCK)
        mageck_binary: Explicit path to mageck executable
        normalization: Normalization method
        control_sgrna: Optional control sgRNA list for control normalization
        fdr_threshold: FDR threshold used for flagging hits

    Returns:
        Dictionary with parsed results and design matrix
    """
    wrapper = MAGeCKWrapper(conda_env=conda_env, mageck_binary=mageck_binary)

    design_matrix_file = f"{output_prefix}_design_matrix.txt"
    wrapper.generate_design_matrix(
        control_samples=control_samples,
        condition_groups=condition_groups,
        output_file=design_matrix_file
    )

    output_files = wrapper.run_mle(
        count_matrix=count_matrix_file,
        design_matrix=design_matrix_file,
        output_prefix=output_prefix,
        normalization=normalization,
        control_sgrna=control_sgrna,
    )

    gene_results = wrapper.parse_mle_gene_results(
        output_files['gene_summary'], fdr_threshold=fdr_threshold)
    sgrna_results = wrapper.parse_mle_sgrna_results(output_files['sgrna_summary'])

    return {
        'gene_results': gene_results,
        'sgrna_results': sgrna_results,
        'output_files': output_files,
        'design_matrix': design_matrix_file
    }
