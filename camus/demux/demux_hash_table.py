#!/usr/bin/env python3
"""
CAMUS cell-type demultiplexer - hash-table implementation.

Demultiplexes pooled CRISPR-screen FASTQs by a cell-type barcode and counts
sgRNAs per cell type. Uses dictionary lookup for O(1) gRNA matching instead of
iterating all gRNAs (50-100x faster than per-read Hamming iteration).

Correctness features:
- Ambiguous 1-mismatch gRNA variants are DROPPED, not silently overwritten.
- Barcode ties (a read equidistant from two cell-type barcodes) are rejected
  as ambiguous instead of being assigned arbitrarily.
- Optional base-quality filtering on the barcode and gRNA windows.
- Optional anchor/scaffold-motif search for the gRNA, for staggered reads where
  a fixed offset does not work.
- Warns when cell-type barcodes are too close to error-correct safely.
"""

import csv
import gzip
import os
from collections import defaultdict
import time


class FastDemultiplexer:
    def __init__(self, barcode_csv, barcode_start=22, barcode_length=8,
                 grna_start=12, grna_length=20, max_barcode_mismatches=1,
                 allow_grna_mismatch=True, min_qual=0,
                 grna_anchor=None, anchor_max_offset=40):
        """
        Args:
            barcode_csv: CSV with columns Celltype,Barcode.
            barcode_start/length: 0-based window for the cell-type barcode in R2.
            grna_start/length: 0-based window for the gRNA in R1 (fixed-offset mode).
            max_barcode_mismatches: max Hamming distance for barcode error correction.
            allow_grna_mismatch: enable precomputed 1-mismatch gRNA lookup.
            min_qual: minimum per-base Phred quality required across the barcode and
                gRNA windows (0 disables quality filtering). Assumes Phred+33.
            grna_anchor: optional scaffold/anchor sequence preceding the gRNA. When
                set, the gRNA is taken from the bases immediately following the first
                occurrence of this anchor within the first ``anchor_max_offset`` bases
                of R1 (overrides the fixed ``grna_start``).
            anchor_max_offset: how far into R1 to search for the anchor.
        """
        self.celltype_barcodes = self._load_barcodes(barcode_csv)
        self.barcode_start = barcode_start
        self.barcode_length = barcode_length
        self.grna_start = grna_start
        self.grna_length = grna_length
        self.max_barcode_mismatches = max_barcode_mismatches
        self.allow_grna_mismatch = allow_grna_mismatch
        self.min_qual = min_qual
        self.grna_anchor = grna_anchor.upper() if grna_anchor else None
        self.anchor_max_offset = anchor_max_offset

        # Track detailed statistics
        self.barcode_stats = defaultdict(lambda: defaultdict(int))
        self.celltype_counts = defaultdict(int)
        self.n_ambiguous_variants = 0

        print(f"Loaded {len(self.celltype_barcodes)} cell types: {list(self.celltype_barcodes.keys())}")
        if self.grna_anchor:
            print(f"gRNA extraction: anchor '{self.grna_anchor}' + {grna_length} bp "
                  f"(searched within first {anchor_max_offset} bp of R1)")
        else:
            print(f"gRNA extraction: positions {grna_start} to {grna_start + grna_length}")
        print(f"Barcode extraction: positions {barcode_start} to {barcode_start + barcode_length}")
        print(f"gRNA mismatch tolerance: {'1 mismatch' if allow_grna_mismatch else 'exact match only'}")
        if self.min_qual > 0:
            print(f"Base-quality filter: min Phred {self.min_qual} across barcode + gRNA windows")

        # Sanity-check that barcodes are far enough apart to error-correct safely.
        self._check_barcode_separation()

    def _load_barcodes(self, csv_file):
        """Load cell type barcodes."""
        barcodes = {}
        with open(csv_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                celltype = row['Celltype'].strip()
                barcode = row['Barcode'].strip().upper()
                barcodes[celltype] = barcode
        return barcodes

    def _check_barcode_separation(self):
        """Warn if any two cell-type barcodes are closer than 2*mismatches+1 apart.

        With error correction of up to ``max_barcode_mismatches``, two barcodes
        must differ by at least 2*max_barcode_mismatches+1 to be unambiguously
        correctable. Closer than that and some observed barcodes are equidistant
        from both references (and will be dropped as ties).
        """
        min_safe = 2 * self.max_barcode_mismatches + 1
        items = list(self.celltype_barcodes.items())
        warned = False
        for i in range(len(items)):
            for j in range(i + 1, len(items)):
                (ct1, bc1), (ct2, bc2) = items[i], items[j]
                d = self.hamming_distance(bc1, bc2)
                if d < min_safe:
                    print(f"  WARNING: barcodes for '{ct1}' and '{ct2}' differ by only "
                          f"{d} bp (need >= {min_safe} for {self.max_barcode_mismatches}-mismatch "
                          f"correction). Reads in between will be dropped as ambiguous.")
                    warned = True
        if warned:
            print("  -> Consider reducing --max-barcode-mismatches or redesigning barcodes.")

    def _load_library_as_dict(self, library_file):
        """
        Load gRNA library into hash tables for O(1) lookup.

        Returns:
            exact_match_dict: {gRNA_seq: gene} for exact matches.
            one_mismatch_dict: {gRNA_variant: (original_gRNA, gene)} for UNIQUE
                1-mismatch variants only.

        Ambiguous 1-mismatch variants (a variant reachable from two or more
        distinct sgRNAs, or colliding with an exact sgRNA) are dropped so the
        read is treated as no-match rather than mis-assigned.
        """
        print("\nBuilding gRNA hash tables...")
        start_time = time.time()

        exact_match_dict = {}
        with open(library_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    grna = parts[0].strip().upper()
                    gene = parts[1].strip()
                    if grna in exact_match_dict and exact_match_dict[grna] != gene:
                        print(f"  WARNING: sgRNA {grna} listed with two genes "
                              f"({exact_match_dict[grna]} / {gene}); keeping first.")
                        continue
                    exact_match_dict[grna] = gene

        one_mismatch_dict = {}
        n_ambiguous = 0
        if self.allow_grna_mismatch:
            # First pass: count how many distinct sgRNAs generate each variant.
            variant_sources = defaultdict(set)
            for grna in exact_match_dict:
                for pos in range(len(grna)):
                    original_base = grna[pos]
                    for new_base in ('A', 'C', 'G', 'T'):
                        if new_base != original_base:
                            variant = grna[:pos] + new_base + grna[pos + 1:]
                            variant_sources[variant].add(grna)

            # Second pass: keep a variant only if exactly one sgRNA produces it
            # and it does not collide with any exact sgRNA sequence.
            for variant, sources in variant_sources.items():
                if variant in exact_match_dict:
                    continue
                if len(sources) == 1:
                    grna = next(iter(sources))
                    one_mismatch_dict[variant] = (grna, exact_match_dict[grna])
                else:
                    n_ambiguous += 1

        elapsed = time.time() - start_time
        print(f"  Built hash tables in {elapsed:.2f} seconds")
        print(f"  Exact match entries: {len(exact_match_dict):,}")
        if self.allow_grna_mismatch:
            print(f"  Unique 1-mismatch variant entries: {len(one_mismatch_dict):,}")
            print(f"  Ambiguous variants dropped: {n_ambiguous:,}")
        print(f"  Total lookup entries: {len(exact_match_dict) + len(one_mismatch_dict):,}")

        self.n_ambiguous_variants = n_ambiguous
        return exact_match_dict, one_mismatch_dict

    def hamming_distance(self, s1, s2):
        """Calculate Hamming distance between two equal-length strings."""
        if len(s1) != len(s2):
            return float('inf')
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))

    @staticmethod
    def _min_phred(qual_window):
        """Minimum Phred+33 quality across a quality string window."""
        if not qual_window:
            return 0
        return min(ord(c) - 33 for c in qual_window)

    def extract_barcode(self, read2_seq):
        """Extract barcode from read2."""
        if len(read2_seq) < self.barcode_start + self.barcode_length:
            return None
        return read2_seq[self.barcode_start:self.barcode_start + self.barcode_length]

    def extract_grna(self, read1_seq):
        """Extract gRNA from read1.

        Uses anchor/scaffold search if ``grna_anchor`` is set, otherwise a fixed
        offset. Returns (grna_seq, start_index) or (None, None).
        """
        if self.grna_anchor:
            search_region = read1_seq[:self.anchor_max_offset + len(self.grna_anchor)].upper()
            idx = search_region.find(self.grna_anchor)
            if idx == -1:
                return None, None
            start = idx + len(self.grna_anchor)
            if len(read1_seq) < start + self.grna_length:
                return None, None
            return read1_seq[start:start + self.grna_length].upper(), start
        else:
            if len(read1_seq) < self.grna_start + self.grna_length:
                return None, None
            start = self.grna_start
            return read1_seq[start:start + self.grna_length].upper(), start

    def match_barcode(self, observed_bc):
        """Match barcode to cell type with error correction and tie rejection.

        Returns (celltype, distance). If two cell types are equidistant at the
        minimum distance, the read is ambiguous and (None, distance) is returned.
        """
        best_celltype = None
        best_distance = float('inf')
        second_best_distance = float('inf')

        for celltype, expected_bc in self.celltype_barcodes.items():
            distance = self.hamming_distance(observed_bc, expected_bc)
            if distance < best_distance:
                second_best_distance = best_distance
                best_distance = distance
                best_celltype = celltype
            elif distance < second_best_distance:
                second_best_distance = distance

        if best_distance > self.max_barcode_mismatches:
            return None, best_distance
        # Tie: the observed barcode is equally close to >= 2 cell types.
        if second_best_distance == best_distance:
            return None, best_distance
        return best_celltype, best_distance

    def match_grna_fast(self, observed_grna, exact_match_dict, one_mismatch_dict):
        """
        Fast gRNA matching using hash table lookup.

        Returns: (is_valid, matched_grna, gene)
        """
        if observed_grna in exact_match_dict:
            return True, observed_grna, exact_match_dict[observed_grna]

        if self.allow_grna_mismatch and observed_grna in one_mismatch_dict:
            original_grna, gene = one_mismatch_dict[observed_grna]
            return True, original_grna, gene

        return False, None, None

    def detect_control_sgrnas(self, exact_match_dict,
                              control_genes=("control", "non-targeting", "nontargeting",
                                             "ntc", "neg_control", "negctrl")):
        """Return the set of sgRNA sequences whose gene label marks them as
        non-targeting controls. Matching is case-insensitive and substring-based
        on the gene name."""
        keys = tuple(k.lower() for k in control_genes)
        controls = set()
        for grna, gene in exact_match_dict.items():
            g = str(gene).lower()
            if any(k in g for k in keys):
                controls.add(grna)
        return controls

    def process(self, read1_file, read2_file, library_file, output_prefix):
        """Process FASTQ files with fast hash table lookup."""
        print("\n" + "=" * 60)
        print("PROCESSING FASTQ FILES (HASH TABLE MODE)")
        print("=" * 60)

        exact_match_dict, one_mismatch_dict = self._load_library_as_dict(library_file)

        counts = {celltype: defaultdict(int) for celltype in self.celltype_barcodes.keys()}
        stats = {
            'total': 0,
            'assigned': defaultdict(int),
            'unassigned': 0,
            'ambiguous_barcode': 0,
            'no_barcode': 0,
            'no_grna': 0,
            'invalid_grna': 0,
            'low_quality': 0,
            'exact_match': 0,
            'one_mismatch': 0,
        }

        r1_opener = gzip.open if read1_file.endswith('.gz') else open
        r2_opener = gzip.open if read2_file.endswith('.gz') else open

        print("\nProcessing reads...")
        start_time = time.time()

        with r1_opener(read1_file, 'rt') as r1, r2_opener(read2_file, 'rt') as r2:
            while True:
                r1_lines = [r1.readline().rstrip('\n') for _ in range(4)]
                r2_lines = [r2.readline().rstrip('\n') for _ in range(4)]

                # End of file: header line empty
                if not r1_lines[0] or not r2_lines[0]:
                    break

                stats['total'] += 1

                if stats['total'] % 100000 == 0:
                    elapsed = time.time() - start_time
                    rate = stats['total'] / elapsed if elapsed > 0 else 0
                    print(f"  Processed {stats['total']:,} reads... ({rate:,.0f} reads/sec)", end='\r')

                r1_seq, r1_qual = r1_lines[1], r1_lines[3]
                r2_seq, r2_qual = r2_lines[1], r2_lines[3]

                # Extract barcode from read2
                barcode = self.extract_barcode(r2_seq)
                if not barcode:
                    stats['no_barcode'] += 1
                    continue

                # Extract gRNA from read1
                grna, grna_start = self.extract_grna(r1_seq)
                if not grna:
                    stats['no_grna'] += 1
                    continue

                # Base-quality filtering on barcode + gRNA windows
                if self.min_qual > 0:
                    bc_qual = r2_qual[self.barcode_start:self.barcode_start + self.barcode_length]
                    grna_qual = r1_qual[grna_start:grna_start + self.grna_length]
                    if (self._min_phred(bc_qual) < self.min_qual or
                            self._min_phred(grna_qual) < self.min_qual):
                        stats['low_quality'] += 1
                        continue

                # Fast hash table lookup for gRNA
                is_valid, matched_grna, gene = self.match_grna_fast(
                    grna, exact_match_dict, one_mismatch_dict
                )

                if not is_valid:
                    stats['invalid_grna'] += 1
                    continue

                if matched_grna == grna:
                    stats['exact_match'] += 1
                else:
                    stats['one_mismatch'] += 1

                # Match barcode to cell type (ties rejected)
                celltype, distance = self.match_barcode(barcode)

                self.barcode_stats[barcode]['total'] += 1

                if celltype:
                    counts[celltype][matched_grna] += 1
                    stats['assigned'][celltype] += 1
                    self.celltype_counts[celltype] += 1
                    self.barcode_stats[barcode]['assigned'] += 1
                    self.barcode_stats[barcode]['celltype'] = celltype
                    if distance == 0:
                        self.barcode_stats[barcode]['exact_match'] += 1
                    else:
                        self.barcode_stats[barcode]['corrected'] += 1
                else:
                    # distance <= max_mismatches but tied => ambiguous;
                    # distance > max_mismatches => simply unassigned.
                    if distance <= self.max_barcode_mismatches:
                        stats['ambiguous_barcode'] += 1
                        self.celltype_counts['ambiguous'] += 1
                        self.barcode_stats[barcode]['ambiguous'] += 1
                    else:
                        stats['unassigned'] += 1
                        self.celltype_counts['unassigned'] += 1
                        self.barcode_stats[barcode]['unassigned'] += 1

        elapsed = time.time() - start_time
        print(f"\n  Processed {stats['total']:,} total reads in {elapsed:.2f} seconds")
        if elapsed > 0:
            print(f"  Average speed: {stats['total']/elapsed:,.0f} reads/second")

        self._write_outputs(counts, exact_match_dict, output_prefix, stats)

        # Write a control-sgRNA list for MAGeCK control normalization, if any.
        controls = self.detect_control_sgrnas(exact_match_dict)
        if controls:
            control_file = f"{output_prefix}_control_sgrnas.txt"
            with open(control_file, 'w') as f:
                for grna in sorted(controls):
                    f.write(grna + "\n")
            print(f"✓ Control sgRNA list ({len(controls):,}): {control_file}")

        return stats

    def _write_outputs(self, counts, library_dict, output_prefix, stats):
        """Write count table and summary."""
        count_file = f"{output_prefix}_count.txt"

        with open(count_file, 'w') as f:
            celltypes = sorted(self.celltype_barcodes.keys())
            f.write("sgRNA\tGene\t" + "\t".join(celltypes) + "\n")
            for grna in sorted(library_dict.keys()):
                gene = library_dict[grna]
                celltype_counts = [str(counts[ct].get(grna, 0)) for ct in celltypes]
                f.write(f"{grna}\t{gene}\t" + "\t".join(celltype_counts) + "\n")

        print(f"\n✓ Count file: {count_file}")

        self._write_celltype_proportions(output_prefix, stats)
        self._write_barcode_performance(output_prefix, stats)

        summary_file = f"{output_prefix}_summary.txt"
        total = max(stats['total'], 1)
        with open(summary_file, 'w') as f:
            f.write("DEMULTIPLEXING SUMMARY (HASH TABLE MODE)\n")
            f.write("=" * 60 + "\n\n")
            f.write(f"Total reads: {stats['total']:,}\n")
            f.write(f"Assigned reads: {sum(stats['assigned'].values()):,}\n")
            f.write(f"Unassigned reads (barcode too far): {stats['unassigned']:,}\n")
            f.write(f"Ambiguous barcode (tie, dropped): {stats['ambiguous_barcode']:,}\n")
            f.write(f"No barcode: {stats['no_barcode']:,}\n")
            f.write(f"No gRNA: {stats['no_grna']:,}\n")
            f.write(f"Invalid gRNA: {stats['invalid_grna']:,}\n")
            f.write(f"Low quality (filtered): {stats['low_quality']:,}\n\n")

            f.write("gRNA MATCHING STATISTICS:\n")
            f.write(f"  Exact matches: {stats['exact_match']:,} "
                    f"({stats['exact_match']/total*100:.1f}%)\n")
            if self.allow_grna_mismatch:
                f.write(f"  1-mismatch matches: {stats['one_mismatch']:,} "
                        f"({stats['one_mismatch']/total*100:.1f}%)\n")
            f.write(f"  Total valid: {stats['exact_match'] + stats['one_mismatch']:,}\n\n")

            f.write("CELL TYPE DISTRIBUTION:\n")
            for celltype, count in sorted(stats['assigned'].items()):
                pct = (count / total * 100)
                f.write(f"  {celltype}: {count:,} ({pct:.1f}%)\n")

        print(f"✓ Summary: {summary_file}")

        print("\nSTATS:")
        print(f"  Total reads: {stats['total']:,}")
        print(f"  Valid gRNA matches: {stats['exact_match'] + stats['one_mismatch']:,}")
        print(f"    Exact: {stats['exact_match']:,} ({stats['exact_match']/total*100:.1f}%)")
        if self.allow_grna_mismatch:
            print(f"    1-mismatch: {stats['one_mismatch']:,} ({stats['one_mismatch']/total*100:.1f}%)")
        print(f"  Assigned to cell types: {sum(stats['assigned'].values()):,}")
        for celltype, count in sorted(stats['assigned'].items()):
            print(f"    {celltype}: {count:,}")
        print(f"  Ambiguous barcode (dropped): {stats['ambiguous_barcode']:,}")
        print(f"  Unassigned: {stats['unassigned']:,}")
        print(f"  Invalid gRNA: {stats['invalid_grna']:,}")
        if self.min_qual > 0:
            print(f"  Low quality (filtered): {stats['low_quality']:,}")

    def _write_celltype_proportions(self, output_prefix, stats):
        """Write cell type proportion report."""
        prop_file = f"{output_prefix}_celltype_proportions.txt"
        total_assigned = sum(stats['assigned'].values())
        total = max(stats['total'], 1)

        with open(prop_file, 'w') as f:
            f.write("CellType\tRead_Count\tPercentage\tExpected_Barcode\n")
            for celltype in sorted(self.celltype_barcodes.keys()):
                count = stats['assigned'].get(celltype, 0)
                pct = (count / total_assigned * 100) if total_assigned > 0 else 0
                expected_bc = self.celltype_barcodes[celltype]
                f.write(f"{celltype}\t{count}\t{pct:.2f}\t{expected_bc}\n")

            amb = stats['ambiguous_barcode']
            amb_pct = (amb / total * 100)
            f.write(f"ambiguous\t{amb}\t{amb_pct:.2f}\tN/A\n")

            unassigned_count = stats['unassigned']
            unassigned_pct = (unassigned_count / total * 100)
            f.write(f"unassigned\t{unassigned_count}\t{unassigned_pct:.2f}\tN/A\n")

        print(f"✓ Cell type proportions: {prop_file}")

    def _write_barcode_performance(self, output_prefix, stats):
        """Write barcode performance report."""
        perf_file = f"{output_prefix}_barcode_performance.txt"

        with open(perf_file, 'w') as f:
            f.write("Barcode\tCellType\tTotal_Reads\tAssigned_Reads\t"
                    "Assignment_Rate\tExact_Matches\tCorrected\n")

            top_barcodes = sorted(
                self.barcode_stats.items(),
                key=lambda x: x[1]['total'],
                reverse=True
            )[:20]

            for barcode, bc_stats in top_barcodes:
                total = bc_stats['total']
                assigned = bc_stats.get('assigned', 0)
                exact = bc_stats.get('exact_match', 0)
                corrected = bc_stats.get('corrected', 0)
                celltype = bc_stats.get('celltype', 'unassigned')
                assignment_rate = (assigned / total * 100) if total > 0 else 0
                f.write(f"{barcode}\t{celltype}\t{total}\t{assigned}\t"
                        f"{assignment_rate:.1f}%\t{exact}\t{corrected}\n")

        print(f"✓ Barcode performance: {perf_file}")


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description='CAMUS cell-type demultiplexer (hash-table implementation)'
    )
    parser.add_argument('--read1', required=True, help='Read 1 FASTQ (.fastq or .fastq.gz)')
    parser.add_argument('--read2', required=True, help='Read 2 FASTQ (.fastq or .fastq.gz)')
    parser.add_argument('--library', required=True, help='gRNA library (TSV: sgRNA<TAB>Gene)')
    parser.add_argument('--celltype-barcodes', required=True, help='Barcode CSV (Celltype,Barcode)')
    parser.add_argument('--output-prefix', required=True, help='Output prefix')
    parser.add_argument('--barcode-start', type=int, default=22, help='Barcode start position')
    parser.add_argument('--barcode-length', type=int, default=8, help='Barcode length')
    parser.add_argument('--grna-start', type=int, default=12, help='gRNA start position (fixed mode)')
    parser.add_argument('--grna-length', type=int, default=20, help='gRNA length')
    parser.add_argument('--max-barcode-mismatches', type=int, default=1,
                        help='Max barcode mismatches')
    parser.add_argument('--no-grna-mismatch', action='store_true',
                        help='Disable 1-mismatch tolerance for gRNA (exact match only)')
    parser.add_argument('--min-qual', type=int, default=0,
                        help='Minimum Phred quality across barcode + gRNA windows (0 = off)')
    parser.add_argument('--grna-anchor', type=str, default=None,
                        help='Scaffold/anchor sequence preceding the gRNA (enables anchor search)')
    parser.add_argument('--anchor-max-offset', type=int, default=40,
                        help='How far into R1 to search for the anchor')

    args = parser.parse_args()

    for f in [args.read1, args.read2, args.library, args.celltype_barcodes]:
        if not os.path.exists(f):
            print(f"ERROR: File not found: {f}")
            return

    output_dir = os.path.dirname(args.output_prefix)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    demux = FastDemultiplexer(
        barcode_csv=args.celltype_barcodes,
        barcode_start=args.barcode_start,
        barcode_length=args.barcode_length,
        grna_start=args.grna_start,
        grna_length=args.grna_length,
        max_barcode_mismatches=args.max_barcode_mismatches,
        allow_grna_mismatch=not args.no_grna_mismatch,
        min_qual=args.min_qual,
        grna_anchor=args.grna_anchor,
        anchor_max_offset=args.anchor_max_offset,
    )

    demux.process(
        read1_file=args.read1,
        read2_file=args.read2,
        library_file=args.library,
        output_prefix=args.output_prefix
    )

    print("\n" + "=" * 60)
    print("DEMULTIPLEXING COMPLETE!")
    print("=" * 60)


if __name__ == "__main__":
    main()
