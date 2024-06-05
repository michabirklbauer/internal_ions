#!/usr/bin/env python3

import os
import shutil
import random
from datetime import datetime

from pyteomics import mzid
from psm_utils import io as psm_io

from typing import Dict
from typing import BinaryIO

# read identifications
def read_mzid(filename: str | BinaryIO, name: str) -> Dict[str, Dict]:
    """
    Returns a dictionary that proteins/peptides to scan numbers:
    Dict["name": str,
         "proteins": Dict[str, Set[int]],
         "peptides": Dict[str, Set[int]]
    """

    proteins_to_scannr = dict()
    peptides_to_scannr = dict()

    print("Read identifications in total:")

    with mzid.read(filename) as reader:
        nr_psms = 0
        for s in reader:
            scan_nr = int(s["name"])
            for psm in s["SpectrumIdentificationItem"]:
                peptide = psm["PeptideSequence"]
                if peptide in peptides_to_scannr:
                    peptides_to_scannr[peptide].add(scan_nr)
                else:
                    peptides_to_scannr[peptide] = {scan_nr}
                for p in psm["PeptideEvidenceRef"]:
                    protein = p["accession"]
                    if protein in proteins_to_scannr:
                        proteins_to_scannr[protein].add(scan_nr)
                    else:
                        proteins_to_scannr[protein] = {scan_nr}
                nr_psms += 1
                if nr_psms % 1000 == 0:
                    print(f"\t{nr_psms}")
        reader.close()

    print(f"\nFinished reading {nr_psms} identifications!")

    return {"name": name, "proteins": proteins_to_scannr, "peptides": peptides_to_scannr}

def read_other(filename: str | BinaryIO, name: str, verbose: bool = False) -> Dict[str, Dict]:
    """
    Returns a dictionary that proteins/peptides to scan numbers:
    Dict["name": str,
         "proteins": Dict[str, Set[int]],
         "peptides": Dict[str, Set[int]]
    """

    psms = list()
    if type(filename) == str:
        psms = psm_io.read_file(filename)
    else:
        tmp_dir_name = "tmp_fragannot_files_471635739"
        if os.path.exists(tmp_dir_name) and os.path.isdir(tmp_dir_name):
            shutil.rmtree(tmp_dir_name)
        os.makedirs(tmp_dir_name)
        output_name_prefix = tmp_dir_name + "/" + datetime.now().strftime("%b-%d-%Y_%H-%M-%S") + "_" + str(random.randint(10000, 99999))
        with open(output_name_prefix + filename.name, "wb") as f:
            f.write(filename.getbuffer())
            f.close()
        psms = psm_io.read_file(output_name_prefix + filename.name)
        try:
            os.remove(output_name_prefix + filename.name)
        except Exception as e:
            if verbose:
                print("Could not remove file: " + output_name_prefix + filename.name)
    if len(psms) == 0:
        print("Error: Couldn't read identifications file!")
        return {"name": name, "proteins": dict(), "peptides": dict()}

    proteins_to_scannr = dict()
    peptides_to_scannr = dict()

    print("Read identifications in total:")

    nr_psms = 0
    for psm in psms:
        scan_nr = int(psm["spectrum_id"])
        peptide = psm.peptidoform.sequence
        if peptide in peptides_to_scannr:
            peptides_to_scannr[peptide].add(scan_nr)
        else:
            peptides_to_scannr[peptide] = {scan_nr}
        for protein in psm["protein_list"]:
            if protein in proteins_to_scannr:
                proteins_to_scannr[protein].add(scan_nr)
            else:
                proteins_to_scannr[protein] = {scan_nr}
        nr_psms += 1
        if nr_psms % 1000 == 0:
            print(f"\t{nr_psms}")

    print(f"\nFinished reading {nr_psms} identifications!")

    return {"name": name, "proteins": proteins_to_scannr, "peptides": peptides_to_scannr}

def read_identifications(filename: str | BinaryIO, name: str, verbose: bool = False) -> Dict[str, Dict]:
    """
    Returns a dictionary that proteins/peptides to scan numbers:
    Dict["name": str,
         "proteins": Dict[str, Set[int]],
         "peptides": Dict[str, Set[int]]
    """

    if name.split(".")[-1].strip() == "mzid":
        return read_mzid(filename, name)
    return read_other(filename, name, verbose)
