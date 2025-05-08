import streamlit as st

from fraggraph.combine_spectra import combine_spectra

from util.redirect import st_stdout
from util.spectrumio import filter_spectra
from util.tab3.plots import plot_consensus_spectrum
from util.tab3.plots import plot_spectra_chromatogram
from util.tab3.fraggraph import main as fraggraph_main

from util.constants import DIV_COLOR

from typing import List
from typing import Dict
from typing import Any


def get_minmz(selection: Dict[str, Any]) -> float:
    if len(selection.selection.box) == 0:
        return 0.0

    box = selection.selection.box[-1]
    return min(box["y"])


def get_maxmz(selection: Dict[str, Any]) -> float:
    if len(selection.selection.box) == 0:
        return 10000.0

    box = selection.selection.box[-1]
    return max(box["y"])


def get_minrt(selection: Dict[str, Any]) -> float:
    if len(selection.selection.box) == 0:
        return min([float(s["rt"]) for s in st.session_state["spectra"]["spectra"].values()])

    box = selection.selection.box[-1]
    return min(box["x"])


def get_maxrt(selection: Dict[str, Any]) -> float:
    if len(selection.selection.box) == 0:
        return max([float(s["rt"]) for s in st.session_state["spectra"]["spectra"].values()])

    box = selection.selection.box[-1]
    return max(box["x"])


def get_selected_scan_numbers(scan_number_str: str) -> List[int]:
    scan_nrs = set()
    for scan_nr in scan_number_str.split(","):
        if scan_nr.strip() != "":
            scan_nrs.add(int(scan_nr.strip()))
    return list(scan_nrs)


def main(argv=None) -> None:

    ############################################################################
    st.subheader("TITLE", divider=DIV_COLOR)
    st.markdown("Reannotation of fragments on a given spectrum (or merged spectra) using fragmentation graph.")

    ############################################################################
    st.subheader("Data Import", divider=DIV_COLOR)

    if "spectra" in st.session_state:
        st.success(f"Spectra from file \"{st.session_state['spectra']['name']}\" were successfully loaded!")

    ############################################################################
        st.subheader("Spectrum Selector", divider=DIV_COLOR)
        st.markdown("Using the following fields you can filter your spectra for further analyis.")

        # plot chromatogram
        spectra_chromatogram = st.plotly_chart(plot_spectra_chromatogram(st.session_state["spectra"]["spectra"]),
                                               use_container_width=True,
                                               theme="streamlit",
                                               key="chromatogram",
                                               on_select="rerun",
                                               selection_mode="box")

        spec_sel_col1, spec_sel_col2 = st.columns(2)

        with spec_sel_col1:
            first_scan = st.selectbox("Select the first scan number to analyze:",
                                      sorted(st.session_state["spectra"]["spectra"].keys()),
                                      index=0,
                                      help="Select the scan number where analysis should start.",
                                      key="first_scan_filter")

            min_mz = st.number_input("Select the minimum m/z for a spectrum:",
                                     min_value=0.0,
                                     max_value=10000.0,
                                     value=get_minmz(spectra_chromatogram),
                                     step=0.01,
                                     help="The minimum m/z.",
                                     key="min_mz_filter")

            min_rt = st.number_input("Select the minimum retention time for a spectrum:",
                                     min_value=min([float(s["rt"]) for s in st.session_state["spectra"]["spectra"].values()]),
                                     max_value=max([float(s["rt"]) for s in st.session_state["spectra"]["spectra"].values()]),
                                     value=get_minrt(spectra_chromatogram),
                                     step=0.01,
                                     help="The minimum rt.",
                                     key="min_rt_filter")

            # max_charge = st.number_input("Select the maximum charge for a spectrum:",
            #                              min_value = 1,
            #                              max_value = 20,
            #                              value = 4,
            #                              step = 1,
            #                              help = "The maximum charge.",
            #                              key="max_charge_filter")

        with spec_sel_col2:
            last_scan = st.selectbox("Select the last scan number to analyze:",
                                     sorted(st.session_state["spectra"]["spectra"].keys()),
                                     index=len(st.session_state["spectra"]["spectra"].keys()) - 1,
                                     help="Select the scan number where analysis should end.",
                                     key="last_scan_filter")

            max_mz = st.number_input("Select the maximum m/z for a spectrum:",
                                     min_value=0.0,
                                     max_value=10000.0,
                                     value=get_maxmz(spectra_chromatogram),
                                     step=0.01,
                                     help="The maximum m/z.",
                                     key="max_mz_filter")

            max_rt = st.number_input("Select the maximum retention time for a spectrum:",
                                     min_value=min([float(s["rt"]) for s in st.session_state["spectra"]["spectra"].values()]),
                                     max_value=max([float(s["rt"]) for s in st.session_state["spectra"]["spectra"].values()]),
                                     value=get_maxrt(spectra_chromatogram),
                                     step=0.01,
                                     help="The maximum rt.",
                                     key="max_rt_filter")

            # max_isotope = st.number_input("Select the maximum isotope for a spectrum:",
            #                               min_value = 1,
            #                               max_value = 20,
            #                               value = 4,
            #                               step = 1,
            #                               help = "The maximum isotope.",
            #                               key="max_isotope_filter")

        if "identifications" in st.session_state:
            st.success(f"Identifications from file \"{st.session_state['identifications']['name']}\" were successfully loaded!")

            protein_col, peptide_col = st.columns(2)

            with protein_col:
                st.selectbox("Select scans from a specific protein:",
                             st.session_state["identifications"]["proteins_to_scannr"].keys(),
                             key="selected_protein_scans",
                             index=None,
                             help="Select a protein of interest.")

            with peptide_col:
                st.selectbox("Select scans from a specific peptide:",
                             st.session_state["identifications"]["peptides_to_scannr"].keys(),
                             key="selected_peptide_scans",
                             index=None,
                             help="Select a peptide of interest.")

            cond1 = "selected_protein_scans" in st.session_state and st.session_state["selected_protein_scans"] is not None
            cond2 = "selected_peptide_scans" in st.session_state and st.session_state["selected_peptide_scans"] is not None
            if cond1 or cond2:
                selected_scans = set()
                if "selected_protein_scans" in st.session_state:
                    if st.session_state["selected_protein_scans"] is not None:
                        scans_from_protein_val = st.session_state["selected_protein_scans"]
                        selected_scans = selected_scans.union(st.session_state["identifications"]["proteins_to_scannr"][scans_from_protein_val])
                if "selected_peptide_scans" in st.session_state:
                    if st.session_state["selected_peptide_scans"] is not None:
                        scans_from_peptide_val = st.session_state["selected_peptide_scans"]
                        selected_scans = selected_scans.union(st.session_state["identifications"]["peptides_to_scannr"][scans_from_peptide_val])
                st.text_area("Currently selected mass spectra by means of scan numbers:",
                             value=",".join([str(s) for s in sorted(selected_scans)]),
                             key="selected_scan_numbers",
                             help="Currently selected mass spectra by means of scan numbers. Add or delete scan numbers here to modify.")

        else:
            st.info("No identifications file was provided! Filtering based on proteins/peptides not available unless a identifications file is uploaded in the \"Annotation\" tab!")

    ############################################################################
        st.subheader("Filter Spectra", divider=DIV_COLOR)

        filter_spectra_col1, filter_spectra_col2 = st.columns(2)

        with filter_spectra_col1:
            fg_mzd = st.number_input("Select mzd for combining spectra:",
                                     min_value=0.0,
                                     max_value=1.0,
                                     value=0.01,
                                     step=0.001,
                                     format="%0.3f",
                                     help="max m/z distance to merge peaks")

        with filter_spectra_col2:
            fg_cov = st.number_input("Select coverage for peak filtering:",
                                     min_value=0.0,
                                     max_value=1.0,
                                     value=0.25,
                                     step=0.001,
                                     format="%0.3f",
                                     help="minimum \"coverage\" to keep peaks in the consensus spectrum")

        l1, center_button, r1 = st.columns(3)

        with center_button:
            run_filter = st.button("Filter spectra and create consensus spectrum",
                                   type="primary",
                                   use_container_width=True)

        if run_filter:
            selected_scans_list = [i for i in range(int(first_scan), int(last_scan) + 1)]
            if "spectra" in st.session_state and st.session_state["spectra"] is not None:
                with st.status("Filtering spectra...") as filter_status:
                    if "selected_scan_numbers" in st.session_state:
                        if st.session_state["selected_scan_numbers"] is not None:
                            if st.session_state["selected_scan_numbers"].strip() != "":
                                selected_scans_list = get_selected_scan_numbers(st.session_state["selected_scan_numbers"])

                    filter_params = {"first_scan": first_scan,
                                     "last_scan": last_scan,
                                     "min_mz": min_mz,
                                     "max_mz": max_mz,
                                     "min_rt": min_rt,
                                     "max_rt": max_rt,
                                     # "max_charge": max_charge,
                                     # "max_isotope": max_isotope,
                                     "scans": selected_scans_list,
                                     }
                    with st_stdout("info"):
                        st.session_state["filtered_spectra"] = filter_spectra(st.session_state["spectra"]["spectra"],
                                                                              filter_params,
                                                                              st.session_state["spectra"]["name"])

                    filter_status.update(label="Successfully finished filtering spectra.", state="complete")
            else:
                st.error("Error reading spectra file! Please re-upload your file in the \"Annotation\" tab!", icon="🚨")

            # build consensus spectrum
            if "filtered_spectra" in st.session_state:
                with st.status("Combining spectra to consensus spectrum...") as consensus_status:
                    with st_stdout("info"):
                        # merge spectra into a single consensus spectrum
                        consensus_spectrum = combine_spectra(spectra_list=st.session_state["filtered_spectra"]["spectra"], mzd=fg_mzd)
                        print("Number of peaks in consensus spectrum before filtering: ", len(consensus_spectrum), "\n")
                        # filter out peaks that are found in less than cov of the spectra
                        consensus_spectrum = consensus_spectrum[consensus_spectrum["cov_spectra"] >= fg_cov]
                        consensus_spectrum = consensus_spectrum.reset_index(drop=True)
                        print("Number of peaks in consensus spectrum after filtering: ", len(consensus_spectrum))
                        st.session_state["consensus_spectrum"] = consensus_spectrum
                    st.success("Successfully created consensus spectrum!")
                    consensus_status.update(label="Successfully created consensus spectrum from "
                                                  f"filtered spectra from file {st.session_state['filtered_spectra']['name']}!",
                                            state="complete")

    else:
        st.error("No spectra file uploaded! Please upload a file in the \"Annotation\" tab!", icon="🚨")

    if "filtered_spectra" in st.session_state:
        st.subheader("Used Filter Parameters", divider=DIV_COLOR)
        st.markdown("Used filter parameters for current spectrum selection:")
        with st.expander("Click to show all filter parameters."):
            st.text(str(st.session_state["filtered_spectra"]["filter_params"]))

    if "consensus_spectrum" in st.session_state:

        ############################################################################
        st.subheader("Spectrum Viewer", divider=DIV_COLOR)
        # st.markdown("desc text")

        spectrum_plot = plot_consensus_spectrum(st.session_state["consensus_spectrum"]["its_mean"],
                                                st.session_state["consensus_spectrum"]["mz_mean"],
                                                st.session_state["consensus_spectrum"]["cov_spectra"]
                                                )
        st.plotly_chart(spectrum_plot, use_container_width=True)
        st.markdown("**Figure 1:** Displaying consensus spetrum.")

    if "consensus_spectrum" in st.session_state:

        ############################################################################
        st.subheader("Fraggraph Parameters", divider=DIV_COLOR)
        st.markdown("Please specify the parameters used for running Fraggraph.")

        if "selected_peptide_scans" in st.session_state and st.session_state["selected_peptide_scans"] is not None:
            selected_peptide = st.session_state["selected_peptide_scans"]
            possible_selection_values = [selected_peptide] + list(st.session_state["identifications"]["peptide_to_peptidoforms"][selected_peptide])
            fg_peptidoform1 = st.selectbox("Specify a peptidoform to consider:",
                                           possible_selection_values,
                                           index=None,
                                           help="Specify a peptidoform to consider for generating the fragment graph, e.g. "
                                                "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALRE.")

            fg_peptidoform2 = st.selectbox("Optionally, specify a peptidoform to compare to:",
                                           possible_selection_values,
                                           index=None,
                                           help="Specify a peptidoform to compare the fragment graph of the first peptidoform to, e.g. "
                                                "ARTKQTARKSTGGKAPRKQLATKAARKSAPAT[-79.966331]GGV[+79.966331]KKPHRYRPGTVALRE.")
        elif "selected_protein_scans" in st.session_state and st.session_state["selected_protein_scans"] is not None:
            selected_protein = st.session_state["selected_protein_scans"]
            possible_peptides = [peptidoform for peptidoform in st.session_state["identifications"]["proteins_to_peptides"][selected_protein]]
            fg_selected_peptide1 = st.selectbox("Select a peptide from the selected protein:",
                                                possible_peptides,
                                                index=None,
                                                help="Select a peptide to consider for generating the fragment graph, e.g. "
                                                     "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALRE.")

            if fg_selected_peptide1 is not None:
                possible_peptidoforms = [fg_selected_peptide1] + list(st.session_state["identifications"]["peptide_to_peptidoforms"][fg_selected_peptide1])
                fg_peptidoform1 = st.selectbox("Specify a peptidoform to consider:",
                                               possible_peptidoforms,
                                               index=None,
                                               help="Specify a peptidoform to consider for generating the fragment graph, e.g. "
                                                    "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALRE.")
                fg_peptidoform2 = st.selectbox("Optionally, specify a peptidoform to compare to:",
                                               possible_peptidoforms,
                                               index=None,
                                               help="Specify a peptidoform to compare the fragment graph of the first peptidoform to, e.g. "
                                                    "ARTKQTARKSTGGKAPRKQLATKAARKSAPAT[-79.966331]GGV[+79.966331]KKPHRYRPGTVALRE.")
        else:
            fg_peptidoform1 = st.text_input("Specify a peptidoform to consider:",
                                            value=None,
                                            help="Specify a peptidoform to consider for generating the fragment graph, e.g. "
                                                 "ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALRE.",
                                            placeholder="ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALRE")

            fg_peptidoform2 = st.text_input("Optionally, specify a peptidoform to compare to:",
                                            value=None,
                                            help="Specify a peptidoform to compare the fragment graph of the first peptidoform to, e.g. "
                                                 "ARTKQTARKSTGGKAPRKQLATKAARKSAPAT[-79.966331]GGV[+79.966331]KKPHRYRPGTVALRE.",
                                            placeholder="ARTKQTARKSTGGKAPRKQLATKAARKSAPAT[-79.966331]GGV[+79.966331]KKPHRYRPGTVALRE")

        fg_run_l, fg_run_center, fg_run_r = st.columns(3)

        with fg_run_center:
            run_fraggraph = st.button("Run Fraggraph!", type="primary", use_container_width=True)

        if run_fraggraph:
            fraggraph_main({"mzd": fg_mzd,
                            "cov": fg_cov,
                            "pep1": fg_peptidoform1,
                            "pep2": fg_peptidoform2})

    else:
        st.error("No consensus spectrum available. Please check spectra selection and create a consensus spectrum first", icon="🚨")
