import os
import shutil
import random
from datetime import datetime
from psm_utils import io as psm_io

from typing import Dict, BinaryIO
from fragannot.spectrumfile import SpectrumFile


def read_identifications(filename: str | BinaryIO,
                         filetype: str,
                         name: str,
                         spec_file: str | BinaryIO,
                         verbose: bool = False) -> Dict[str, Dict]:
    """
    Returns a dictionary that proteins/peptides to scan numbers:
    Dict["name": str,
         "proteins": Dict[str, Set[int]],
         "peptides": Dict[str, Set[int]]
    """

    if isinstance(filename, str):
        psms = psm_io.read_file(filename, filetype=filetype)
    else:
        tmp_dir_name = "tmp_fragannot_files_471635739"
        if os.path.exists(tmp_dir_name) and os.path.isdir(tmp_dir_name):
            shutil.rmtree(tmp_dir_name)
        os.makedirs(tmp_dir_name)
        output_name_prefix = tmp_dir_name + "/" + datetime.now().strftime("%b-%d-%Y_%H-%M-%S") + "_" + str(random.randint(10000, 99999))
        with open(output_name_prefix + filename.name, "wb") as f:
            f.write(filename.getbuffer())
        psms = psm_io.read_file(output_name_prefix + filename.name, filetype=filetype)
        try:
            os.remove(output_name_prefix + filename.name)
        except Exception:
            if verbose:
                print("Could not remove file: " + output_name_prefix + filename.name)
    if len(psms) == 0:
        print("Error: Couldn't read identifications file!")
        return {"name": name, "proteins": dict(), "peptides": dict()}

    proteins_to_scannr = dict()
    peptides_to_scannr = dict()
    peptide_to_peptidoforms = dict()
    proteins_to_peptides = dict()
    spec_file = SpectrumFile(spec_file)
    scan_numbers = {spec_id: i for i, spec_id in enumerate(spec_file.index)}

    print("Read identifications in total:")

    nr_psms = 0
    for psm in psms:
        peptide = psm.peptidoform.sequence
        scan_nr = scan_numbers[psm.spectrum_id]
        # begin parse necessary information
        # peptides_to_scannr
        if peptide in peptides_to_scannr:
            peptides_to_scannr[peptide].add(scan_nr)
        else:
            peptides_to_scannr[peptide] = {scan_nr}
        # proteins_to_scannr
        if psm["protein_list"] is not None:
            for protein in psm["protein_list"]:
                if protein in proteins_to_scannr:
                    proteins_to_scannr[protein].add(scan_nr)
                else:
                    proteins_to_scannr[protein] = {scan_nr}

        if peptide in peptide_to_peptidoforms:
            peptide_to_peptidoforms[peptide].add(psm.peptidoform.proforma)
        else:
            peptide_to_peptidoforms[peptide] = {psm.peptidoform.proforma}
        # proteins_to_peptides
        if psm["protein_list"] is not None:
            for protein in psm["protein_list"]:
                if protein in proteins_to_peptides:
                    proteins_to_peptides[protein].add(peptide)
                else:
                    proteins_to_peptides[protein] = {peptide}
        # end parse
        nr_psms += 1
        if nr_psms % 1000 == 0:
            print(f"\t{nr_psms}")

    print(f"\nFinished reading {nr_psms} identifications!")

    return {"name": name,
            "proteins_to_scannr": proteins_to_scannr,
            "peptides_to_scannr": peptides_to_scannr,
            # "scannr_to_peptidoforms": scannr_to_peptidoforms,
            "peptide_to_peptidoforms": peptide_to_peptidoforms,
            "proteins_to_peptides": proteins_to_peptides}
