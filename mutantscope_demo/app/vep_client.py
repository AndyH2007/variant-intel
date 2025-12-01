# mutantscope_demo/app/vep_client.py
import json
from typing import Any, Dict, Optional, List
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, quote
from urllib.request import Request, urlopen

VEP_BASE_URL = "https://rest.ensembl.org/vep"

PLUGIN_PARAMETERS: Dict[str, Any] = {
    "AlphaMissense": 1,
    "AncestralAllele": 1,
    "Blosum62": 1,
    "CADD": 1,
    "ClinPred": 1,
    "Conservation": 1,
    "DosageSensitivity": 1,
    "EVE": 1,
    "Enformer": 1,
    "GO": 1,
    "GeneSplicer": 1,
    "Geno2MP": 1,
    "IntAct": 1,
    "LOEUF": 1,
    "LoF": 1,
    "MaveDB": 1,
    "MaxEntScan": 1,
    "NMD": 1,
    "OpenTargets": 1,
    "Phenotypes": 1,
    "REVEL": 1,
    "RiboseqORFs": 1,
    "SpliceAI": 2,
    "UTRAnnotator": 1,

    # core flags / fields
    "ambiguous_hgvs": 1,
    "appris": 1,
    "canonical": 1,
    "ccds": 1,
    "dbNSFP": "ALL",
    "dbscSNV": 1,
    "domains": 1,
    "failed": 1,
    "flag_pick": 1,
    "flag_pick_allele": 1,
    "flag_pick_allele_gene": 1,
    "ga4gh_vrs": 1,
    "gencode_basic": 0,
    "gencode_primary": 1,
    "hgvs": 1,
    "mane": 1,
    "merged": 1,
    "minimal": 1,
    "mirna": 1,
    "mutfunc": 1,
    "numbers": 1,
    "per_gene": 1,
    "pick": 1,
    "pick_allele": 1,
    "pick_allele_gene": 1,
    "protein": 1,
    "shift_3prime": 1,
    "shift_genomic": 1,
    "transcript_version": 1,
    "tsl": 1,
    "uniprot": 1,
    "variant_class": 1,
    "vcf_string": 1,
    "xref_refseq": 1,
}

def _build_query_params(extra: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    params = dict(PLUGIN_PARAMETERS)
    if extra:
        params.update(extra)
    return params

def fetch_vep_annotation(
    hgvs_notation: str,
    species: str = "human",
    extra_params: Optional[Dict[str, Any]] = None,
) -> List[Dict[str, Any]]:
    if not hgvs_notation:
        raise ValueError("HGVS notation is required.")

    encoded_hgvs = quote(hgvs_notation, safe="")
    params = _build_query_params(extra_params)
    qs = urlencode(params, doseq=True)
    url = f"{VEP_BASE_URL}/{quote(species, safe='')}/hgvs/{encoded_hgvs}?{qs}?refsq=1"

    req = Request(
        url,
        headers={
            "Content-Type": "application/json",
            "User-Agent": "precision-med-client/1.0",
        },
    )
    try:
        with urlopen(req) as resp:
            return json.load(resp)
    except HTTPError as e:
        body = e.read().decode("utf-8", errors="ignore")
        raise RuntimeError(f"VEP API returned HTTP {e.code}: {body or e.reason}") from e
    except URLError as e:
        raise RuntimeError(f"Failed to reach VEP API: {e.reason}") from e
