#!/usr/bin/env python3

import os
import shutil
import random
from datetime import datetime

from psm_utils import io as psm_io

from typing import Tuple
from typing import Dict
from typing import BinaryIO

# parse scan number from the psm utils spectrum id
def parse_scannr(spectrum_id: str | int, i: int, pattern: str = "\\.\\d+\\.") -> Tuple[int, int]:
    """Parses the scan number from the params dictionary of the pyteomics mgf
    spectrum.

    Parameters
    ----------
    spectrum_id : str or int
        The spectrum id of a PSM as given by psm_utils.

    i : int
        The scan number to be returned in case of failure.

    pattern : str
        Regex pattern to use for parsing the scan number from the title if it
        can't be infered otherwise.

    Returns
    -------
    (exit_code, scan_nr) : Tuple
        A tuple with the exit code (0 if successful, 1 if parsing failed) at the
        first position [0] and the scan number at the second position [1].
    """

    if type(spectrum_id) == int:
        return spectrum_id

    # if there is a scan token in the title, try parse scan_nr
    if "scan" in spectrum_id:
        try:
            return (0, int(spectrum_id.split("scan=")[1].strip("\"")))
        except:
            pass

    # else try to parse by pattern
    try:
        scan_nr = re.findall(pattern, spectrum_id)[0]
        scan_nr = re.sub(r"[^0-9]", "", scan_nr)
        if len(scan_nr) > 0:
            return (0, int(scan_nr))
    except:
        pass

    # else try parse whole title
    try:
        return (0, int(float(spectrum_id)))
    except:
        pass

    # return unsuccessful parse
    return (1, i)

def read_identifications(filename: str | BinaryIO,
                         name: str,
                         pattern: str = "\\.\\d+\\.",
                         verbose: bool = False) -> Dict[str, Dict]:
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
    scannr_to_peptidoforms = dict()
    peptide_to_peptidoforms = dict()
    proteins_to_peptides = dict()

    print("Read identifications in total:")

    nr_psms = 0
    for psm in psms:
        # spectrum identfier, can also be str
        # see https://psm-utils.readthedocs.io/en/v1.2.0/api/psm_utils/#psm_utils.PSM
        # don't know how to handle tbh
        parsed_scan_nr = parse_scannr(psm["spectrum_id"], 0, pattern)
        scan_nr = parsed_scan_nr[1]
        if parsed_scan_nr[0] != 0:
            raise RuntimeError(f"Could not parse scan nr from spectrum id {psm['spectrum_id']}.")
        # this should return the unmodified peptide sequence
        # according to https://psm-utils.readthedocs.io/en/v1.2.0/api/psm_utils/#psm_utils.Peptidoform
        peptide = psm.peptidoform.sequence
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
        # scannr_to_peptidoforms
        # using psm.peptidoform.proforma
        # see https://psm-utils.readthedocs.io/en/v1.2.0/api/psm_utils/#psm_utils.Peptidoform
        if scan_nr in scannr_to_peptidoforms:
            scannr_to_peptidoforms[scan_nr].add(psm.peptidoform.proforma)
        else:
            scannr_to_peptidoforms[scan_nr] = {psm.peptidoform.proforma}
        # peptide_to_peptidoforms
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
            "scannr_to_peptidoforms": scannr_to_peptidoforms,
            "peptide_to_peptidoforms": peptide_to_peptidoforms,
            "proteins_to_peptides": proteins_to_peptides}
