"""Tests for the motif checker module."""

from __future__ import annotations

from antibody_liability_tool.mutations.motif_checker import MotifHit, check_motifs


class TestCheckMotifs:
    """Tests for check_motifs()."""

    def test_detect_n_glycosylation(self) -> None:
        """NAS triggers N-glycosylation (N[^P][ST])."""
        hits = check_motifs("XNASX", check_free_cys=False, check_oxidation_met=False)
        glyco = [h for h in hits if h.motif_name == "N-glycosylation"]
        assert len(glyco) >= 1
        assert glyco[0].matched_text == "NAS"

    def test_no_glycosylation_with_proline(self) -> None:
        """NPS should NOT trigger N-glycosylation (proline blocks it)."""
        hits = check_motifs("XNPSX", check_free_cys=False, check_oxidation_met=False)
        glyco = [h for h in hits if h.motif_name == "N-glycosylation"]
        assert len(glyco) == 0

    def test_detect_deamidation_ng(self) -> None:
        hits = check_motifs("XNGX", check_free_cys=False, check_oxidation_met=False)
        deamid = [h for h in hits if h.motif_name == "deamidation"]
        assert len(deamid) >= 1
        assert any(h.matched_text == "NG" for h in deamid)

    def test_detect_deamidation_ns(self) -> None:
        hits = check_motifs("XNSX", check_free_cys=False, check_oxidation_met=False)
        deamid = [h for h in hits if h.motif_name == "deamidation"]
        assert len(deamid) >= 1
        assert any(h.matched_text == "NS" for h in deamid)

    def test_detect_isomerization_dg(self) -> None:
        hits = check_motifs("XDGX", check_free_cys=False, check_oxidation_met=False)
        iso = [h for h in hits if h.motif_name == "isomerization"]
        assert len(iso) >= 1
        assert any(h.matched_text == "DG" for h in iso)

    def test_detect_isomerization_ds(self) -> None:
        hits = check_motifs("XDSX", check_free_cys=False, check_oxidation_met=False)
        iso = [h for h in hits if h.motif_name == "isomerization"]
        assert len(iso) >= 1
        assert any(h.matched_text == "DS" for h in iso)

    def test_detect_free_cysteine(self) -> None:
        hits = check_motifs("XACAX", check_free_cys=True, check_oxidation_met=False)
        cys = [h for h in hits if h.motif_name == "free_cysteine"]
        assert len(cys) == 1
        assert cys[0].matched_text == "C"

    def test_detect_surface_methionine(self) -> None:
        hits = check_motifs("XAMAX", check_free_cys=False, check_oxidation_met=True)
        met = [h for h in hits if h.motif_name == "oxidation_met"]
        assert len(met) == 1
        assert met[0].matched_text == "M"

    def test_surface_methionine_with_positions(self) -> None:
        """Only flag M at specified surface positions."""
        hits = check_motifs(
            "AMAG",
            check_free_cys=False,
            check_oxidation_met=True,
            surface_positions={0},
        )
        met = [h for h in hits if h.motif_name == "oxidation_met"]
        # M is at index 1, which is not in surface_positions={0}
        assert len(met) == 0

    def test_clean_sequence_no_motifs(self) -> None:
        """A clean sequence should trigger no motifs."""
        hits = check_motifs("AAAGGG", check_free_cys=True, check_oxidation_met=True)
        assert len(hits) == 0

    def test_motif_hit_dataclass(self) -> None:
        hit = MotifHit(
            motif_name="test",
            pattern="T",
            start_position=0,
            matched_text="T",
        )
        assert hit.motif_name == "test"
        assert hit.start_position == 0
