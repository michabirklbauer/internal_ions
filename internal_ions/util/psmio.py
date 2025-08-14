import os
import streamlit as st
from psm_utils.io import read_file
from psm_utils.psm_list import PSMList
from tempfile import NamedTemporaryFile
from .spectrumio import SpectrumFile
import logging


@st.cache_data(hash_funcs={PSMList: lambda x: tuple(x.runs), SpectrumFile: lambda x: x.name})
def read_identifications(psms: PSMList,
                         name: str,
                         spec_file: SpectrumFile,
                         verbose: bool = False) -> dict[str, dict]:
    """
    Returns a dictionary that proteins/peptides to scan numbers:
    dict["name": str,
         "proteins": dict[str, Set[int]],
         "peptides": dict[str, Set[int]]
    """
    if len(psms) == 0:
        print("Error: Couldn't read identifications file!")
        return {"name": name, "proteins": {}, "peptides": {}}

    # proteins_to_scannr = dict()
    peptides_to_scannr = {}
    peptide_to_peptidoforms = {}
    peptidoforms_to_scannr = {}
    # proteins_to_peptides = dict()
    scan_numbers = {spec_id: i for i, spec_id in enumerate(spec_file.index)}
    for i in scan_numbers:
        print("Example ID", i)
        logging.debug("Example ID: %s", i)
        break

    print("Read identifications in total:")

    nr_psms = 0
    for psm in psms:
        peptide = psm.peptidoform.sequence
        scan_nr = scan_numbers[psm.spectrum_id]
        peptidoform = psm.peptidoform.proforma
        # begin parse necessary information
        # peptides_to_scannr
        if peptide in peptides_to_scannr:
            peptides_to_scannr[peptide].add(scan_nr)
        else:
            peptides_to_scannr[peptide] = {scan_nr}
        if peptidoform in peptidoforms_to_scannr:
            peptidoforms_to_scannr[peptidoform].add(scan_nr)
        else:
            peptidoforms_to_scannr[peptidoform] = {scan_nr}
        # proteins_to_scannr
        # if psm["protein_list"] is not None:
        #     for protein in psm["protein_list"]:
        #         if protein in proteins_to_scannr:
        #             proteins_to_scannr[protein].add(scan_nr)
        #         else:
        #             proteins_to_scannr[protein] = {scan_nr}

        if peptide in peptide_to_peptidoforms:
            peptide_to_peptidoforms[peptide].add(psm.peptidoform.proforma)
        else:
            peptide_to_peptidoforms[peptide] = {psm.peptidoform.proforma}
        # proteins_to_peptides
        # if psm["protein_list"] is not None:
        #     for protein in psm["protein_list"]:
        #         if protein in proteins_to_peptides:
        #             proteins_to_peptides[protein].add(peptide)
        #         else:
        #             proteins_to_peptides[protein] = {peptide}
        # end parse
        nr_psms += 1
        if nr_psms % 1000 == 0:
            print(f"\t{nr_psms}")

    print(f"\nFinished reading {nr_psms} identifications!")

    return {"name": name,
            # "proteins_to_scannr": proteins_to_scannr,
            "peptides_to_scannr": peptides_to_scannr,
            "peptidoforms_to_scannr": peptidoforms_to_scannr,
            # "scannr_to_peptidoforms": scannr_to_peptidoforms,
            "peptide_to_peptidoforms": peptide_to_peptidoforms,
            # "proteins_to_peptides": proteins_to_peptides
            }


@st.cache_data
def read_id_file(ident_file, file_format):
    # close the temp file before reading it with psm_utils to avoid permission error on Windows
    with NamedTemporaryFile(suffix=os.path.splitext(ident_file.name)[1], delete_on_close=False) as f:
        f.write(ident_file.getbuffer())
        f.close()
        return read_file(f.name, filetype=file_format)
