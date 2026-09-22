"""
Integration tests for amino acid (protein) query input support in COMMEC.

Covers the protein screening path end-to-end using ScreenTesterFactory:
  - biorisk hmmscan runs directly on protein input (no 6-frame translation)
  - blastp replaces blastx for taxonomy search
  - blastn and RNA low-concern steps are skipped for protein-only queries
"""

from commec.config.result import ScreenStatus, ScreenStep
from commec.tests.screen_factory import ScreenTesterFactory


def test_protein_query_no_hits_passes(tmp_path):
    """Protein query with no biorisk or taxonomy hits should report PASS."""
    factory = ScreenTesterFactory("prot_no_hits", tmp_path)
    # Step 6: use is_protein=True so the factory writes protein chars and routes to blastp
    factory.add_query("prot1", 60, is_protein=True)
    result = factory.run()
    assert result.queries["prot1"].status.screen_status == ScreenStatus.PASS


def test_protein_query_biorisk_flag(tmp_path):
    """Protein query with a regulated biorisk hit should be FLAG."""
    factory = ScreenTesterFactory("prot_biorisk_flag", tmp_path)
    factory.add_query("prot1", 60, is_protein=True)
    # Step 6: biorisk hits for protein queries use query_name without frame suffix
    factory.add_hit(
        ScreenStep.BIORISK,
        "prot1",
        start=1,
        stop=50,
        title="dangerous_toxin",
        accession="TOX001",
        regulated=True,
    )
    result = factory.run()
    assert result.queries["prot1"].status.screen_status == ScreenStatus.FLAG


def test_protein_query_biorisk_warn(tmp_path):
    """Protein query with an unregulated biorisk hit should be WARNING."""
    factory = ScreenTesterFactory("prot_biorisk_warn", tmp_path)
    factory.add_query("prot1", 60, is_protein=True)
    # Step 6: regulated=False means the hit is a virulence factor, not a toxin
    factory.add_hit(
        ScreenStep.BIORISK,
        "prot1",
        start=1,
        stop=50,
        title="virulence_factor",
        accession="VIR001",
        regulated=False,
    )
    result = factory.run()
    assert result.queries["prot1"].status.screen_status == ScreenStatus.WARN


def test_protein_query_with_numeric_suffix_in_name(tmp_path):
    """Protein query named ORF_1 must not cause a KeyError in append_nt_querylength_info.

    The _base_name helper previously stripped the trailing _1 from ORF_1 because it
    looks like an NT frame suffix, producing ORF which is not a key in the queries dict.
    The fix checks the queries dict directly before attempting suffix stripping.
    """
    factory = ScreenTesterFactory("prot_numeric_suffix", tmp_path)
    factory.add_query("ORF_1", 60, is_protein=True)
    factory.add_hit(
        ScreenStep.BIORISK,
        "ORF_1",
        start=1,
        stop=50,
        title="dangerous_toxin",
        accession="TOX001",
        regulated=True,
    )
    result = factory.run()
    assert result.queries["ORF_1"].status.screen_status == ScreenStatus.FLAG


def test_protein_query_low_concern_hmmer_clears_taxonomy_flag(tmp_path):
    """Low-concern HMMER hit must clear a taxonomy flag on protein input.

    Previously, reset_query_statuses for LOW_CONCERN_DNA and LOW_CONCERN_RNA both
    wrote ScreenStatus.SKIP to the shared low_concern field before parse_low_concern_hits
    ran, so the HMMER low-concern step (which does work on protein input) never applied
    its results and taxonomy flags were never cleared.

    Note: biorisk hits are never clearable by the low-concern step (by design);
    only taxonomy hits are eligible for clearance.
    """
    factory = ScreenTesterFactory("prot_lc_clears_tax", tmp_path)
    factory.add_query("prot1", 60, is_protein=True)
    factory.add_hit(
        ScreenStep.TAXONOMY_AA,
        "prot1",
        start=1,
        stop=55,
        title="regulated pathogen protein",
        accession="REG001",
        taxid=12345,
        species="Dangerous species",
        genus="Dangerous",
        superkingdom="Bacteria",
        regulated=True,
    )
    factory.add_hit(
        ScreenStep.LOW_CONCERN_PROTEIN,
        "prot1",
        start=1,
        stop=55,
        title="housekeeping_protein",
        accession="HOUSE001",
    )
    result = factory.run()
    assert result.queries["prot1"].status.screen_status == ScreenStatus.CLEARED_FLAG


def test_protein_query_regulated_taxonomy_flags(tmp_path):
    """Protein query whose blastp best match is a regulated pathogen should be FLAG."""
    factory = ScreenTesterFactory("prot_tax_flag", tmp_path)
    factory.add_query("prot1", 60, is_protein=True)
    # Step 6: TAXONOMY_AA hit for a protein query is written to the blastp output file
    factory.add_hit(
        ScreenStep.TAXONOMY_AA,
        "prot1",
        start=1,
        stop=55,
        title="regulated pathogen protein",
        accession="REG001",
        taxid=12345,
        species="Dangerous species",
        genus="Dangerous",
        superkingdom="Bacteria",
        regulated=True,
    )
    result = factory.run()
    assert result.queries["prot1"].status.screen_status == ScreenStatus.FLAG
