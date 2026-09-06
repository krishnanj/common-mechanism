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
