#!/usr/bin/env python3

from pyteomics import mzid

from typing import Dict
from typing import BinaryIO

# read identifications
def read_identifications(filename: str | BinaryIO, name: str) -> Dict[str, Dict]:
    """
    Returns a dictionary that proteins/peptides to scan numbers:
    Dict["name": str,
         "proteins": Dict[str, List[int]],
         "peptides": Dict[str, List[int]]
    """

    proteins_to_scannr = dict()
    peptides_to_scannr = dict()

    with mzid.read(filename) as reader:
        for s in reader:
            scan_nr = int(s["name"])
            for psm in s["SpectrumIdentificationItem"]:
                peptide = psm["PeptideSequence"]
                if peptide in peptides_to_scannr:
                    peptides_to_scannr[peptide].append(scan_nr)
                else:
                    peptides_to_scannr[peptide] = [scan_nr]
                for p in psm["PeptideEvidenceRef"]:
                    protein = p["accession"]
                    if protein in proteins_to_scannr:
                        proteins_to_scannr[protein].append(scan_nr)
                    else:
                        proteins_to_scannr[protein] = [scan_nr]
        reader.close()

    return {"name": name, "proteins": proteins_to_scannr, "peptides": peptides_to_scannr}