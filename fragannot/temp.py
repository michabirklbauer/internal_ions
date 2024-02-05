##
import re
import numpy as np
from psm_utils import Peptidoform
from pyteomics import mass

ion_formulas = {
    "a": {"a": mass.Composition(formula="H-2O-1" + "C-1O-1")},
    "b": {"b": mass.Composition(formula="H-2O-1")},
    "x": {"x": mass.Composition(formula="H-2O-1" + "CO2")},
    "y": {"y": mass.Composition(formula="")},
    "cdot": {"cdot": mass.Composition(formula="H-2O-1" + "NH3")},
    "c": {"c": mass.Composition(formula="H-2O-1" + "NH4")},
    "c-1": {"c-1": mass.Composition(formula="H-2O-1" + "NH2")},
    "c+1": {"c+1": mass.Composition(formula="H-2O-1" + "NH5")},
    "zdot": {"zdot": mass.Composition(formula="H-2O-1" + "N-1" + "OH")},
    "z+1": {"z+1": mass.Composition(formula="H-2O-1" + "N-1" + "OH2")},  # z+M(H)
    "z+2": {"z+2": mass.Composition(formula="H-2O-1" + "N-1" + "OH3")},  # z+M(2H)
    "z+3": {"z+3": mass.Composition(formula="H-2O-1" + "N-1" + "OH4")},
    "c-zdot": {
        "c-zdot": mass.Composition(formula="H-2O-1" + "H-2O-1" + "OH5")
    },  # -O   = c + z (from msnbase issue 82)
    "c-z+1": {"c-z+1": mass.Composition(formula="H-2O-1" + "H-2O-1" + "OH6")},
    "cdot-zdot": {"cdot-zdot": mass.Composition(formula="H-2O-1" + "H-2O-1" + "OH4")},
    "cdot-z+1": {"cdot-z+1": mass.Composition(formula="H-2O-1" + "H-2O-1" + "OH5")},
    "n-n": {"n-n": mass.Composition(formula="P-1")},
    "b-y": {"b-y": mass.Composition(formula="H-2O-1" + "H-2O-1" + "")},
    "a-x": {"a-x": mass.Composition(formula="H-2O-1" + "H-2O-1" + "C-1O-1" + "CO2")},
}


def compute_theoretical_fragments():
    """Returns and set a list of m/z of fragment ions  and informations on the type/position of each fragments for a given peptidoform/proteoform"""

    # store the fragment types looked for:
    ionTypes = ["zdot", "c", "z+1"]

    sequence = "ARTKQTARKSTGGKAPRKQLATKAARKSAPSTGGVKKPHRYRPGTVALRE"
    modifications = {}
    frag_masses = {}

    # fragments masses:
    for ion_type in ionTypes:
        frag_masses_iontype = {}

        if "-" in ion_type:  # Internal Fragment
            # sum of all modification masses present in the internal fragment
            sum_mods = lambda modifications, i, j: sum(
                [mod["monoisotopicMassDelta"] for mod in modifications if i <= mod["location"] <= j]
            )  # sum mods delta for internal fragm ions
            # get all sub string of the peptide sequence
            sub_sequences = [
                (
                    sequence[i - 1 : j],
                    i,
                    j,
                    ion_type,
                    [mod["location"] for mod in modifications if i <= mod["location"] <= j],
                )
                for i in range(1, len(sequence))
                for j in range(i + 1, len(sequence))
            ]
            # compute internal frag masses
            frag_masses_iontype.update(
                {
                    ",".join(str(s) for s in seq[1:4]): round(
                        mass.fast_mass(
                            sequence=seq[0],
                            ion_type=ion_type,
                            ion_comp=ion_formulas[ion_type],
                        )
                        + sum_mods(modifications, seq[1], seq[2]),
                        4,
                    )
                    for seq in sub_sequences
                }
            )

        else:  # Terminal Fragment
            if any(i in ion_type for i in ["a", "b", "c"]):  # Nterm
                sum_mods = lambda modifications, i, j: sum(
                    [mod["monoisotopicMassDelta"] for mod in modifications if mod["location"] <= j]
                )
                sub_sequences = [
                    (
                        sequence[:j],
                        1,
                        j,
                        ion_type,
                        [mod["location"] for mod in modifications if mod["location"] <= j],
                    )
                    for j in range(2, len(sequence))
                ]
                frag_masses_iontype.update(
                    {
                        ",".join(str(s) for s in seq[1:4]): round(
                            mass.fast_mass(
                                sequence=seq[0],
                                ion_type=ion_type,
                                ion_comp=ion_formulas[ion_type],
                            )
                            + sum_mods(modifications, seq[1], seq[2]),
                            4,
                        )
                        for seq in sub_sequences
                    }
                )

            else:  # Cterm
                sum_mods = lambda modifications, i, j: sum(
                    [mod["monoisotopicMassDelta"] for mod in modifications if i <= mod["location"]]
                )
                sub_sequences = [
                    (
                        sequence[i - 1 :],
                        i,
                        len(sequence),
                        ion_type,
                        [mod["location"] for mod in modifications if i <= mod["location"]],
                    )
                    for i in range(1, len(sequence) + 1)
                ]
                frag_masses_iontype.update(
                    {
                        ",".join(str(s) for s in seq[1:4]): round(
                            mass.fast_mass(
                                sequence=seq[0],
                                ion_type=ion_type,
                                ion_comp=ion_formulas[ion_type],
                            )
                            + sum_mods(modifications, seq[1], seq[2]),
                            4,
                        )
                        for seq in sub_sequences
                    }
                )

        frag_masses[ion_type] = frag_masses_iontype

    return frag_masses


print(compute_theoretical_fragments())
# TODO remove addded H
ion_cap_formula = {
    "a": "H-2O-1" + "C-1O-1",
    "b": "H-2O-1",
    "x": "H-2O-1" + "CO2",
    "y": "",
    "c-1": "H-2O-1" + "NH3" + "H-1",
    "c": "H-2O-1" + "NH3",
    "cdot": "H-2O-1" + "NH3" + "H1",
    "c+1": "H-2O-1" + "NH5",
    "zdot": "H-2O-1" + "N-1" + "O",
    "z+1": "H-2O-1" + "N-1" + "OH",
    "z+2": "H-2O-1" + "N-1" + "OH2",
    "z+3": "H-2O-1" + "N-1" + "OH3",
    "t": "",
}

ion_cap_delta_mass = {
    name: mass.calculate_mass(formula, absolute=False) for (name, formula) in ion_cap_formula.items()
}


ion_direction = {
    "a": "n-term",
    "b": "n-term",
    "x": "c-term",
    "y": "c-term",
    "cdot": "n-term",
    "c": "n-term",
    "c-1": "n-term",
    "c+1": "n-term",
    "zdot": "c-term",
    "z+1": "c-term",
    "z+2": "c-term",
    "z+3": "c-term",
}


def parse_fragment_code(fragment_code: str):

    # test if fragment code format is valid*
    fragment_code_pattern = re.compile(".+(:).+(@)[0-9]+(:)[0-9]+(\()(\+|\-)[0-9](\))(\[(.*?)\])?")
    if bool(fragment_code_pattern.match(fragment_code)) == False:
        raise RuntimeError("Incorrect fragment code format: {0}".format(fragment_code))

    ## Parse fragment code

    start, end = [
        int(i) for i in re.search("(?<=\@)(.*?)(?=\()", fragment_code).group(1).split(":")
    ]  # Get start and end amino acid indexes
    ion_cap_start, ion_cap_end = [
        str(i) for i in re.search("^(.*?)(?=\@)", fragment_code).group(1).split(":")
    ]  # Get start and end ion caps name
    charge = int(re.search("(?<=\()(.*?)(?=\))", fragment_code).group(1))  # get charge state
    formula = re.search("(?<=\[)(.*?)(?=\])", fragment_code)
    if formula == None:
        formula = ""
    else:
        formula = str(re.search("(?<=\[)(.*?)(?=\])", fragment_code).group(1))

    return start, end, ion_cap_start, ion_cap_end, charge, formula


from pyteomics import parser


def theoretical_mass_to_charge(fragment_code, peptidoform):

    start, end, ion_cap_start, ion_cap_end, charge, formula = parse_fragment_code(fragment_code)

    # peptide and modification mass
    sequence = []
    mods = []
    for aa, mod in peptidoform.parsed_sequence[start - 1 : end]:
        sequence.append(aa)
        if not mod is None:
            mods.extend([m.mass for m in mod])

    # mass AA sequence
    ps = parser.parse("".join(sequence), show_unmodified_termini=True)
    P = mass.calculate_mass(parsed_sequence=ps)
    # mass modifications
    M = sum(mods)
    # mass start ion cap
    SI = ion_cap_delta_mass[ion_cap_start]
    # mass end ion cap
    EI = ion_cap_delta_mass[ion_cap_end]
    # hydrogen mass
    H = 1.00784
    # loss mass
    L = mass.calculate_mass(formula, absolute=True)

    print(P, M, SI, EI, H, charge, L)
    # Calculate fragment mass
    # TODO mass of hydrogen is considered in the ion cap formula ???
    fragment_mass = (P + M + SI + EI + (H * charge) - L) / np.abs(charge)

    return fragment_mass


print("\n\n")
print(ion_cap_delta_mass)
print("\n\n")
print(theoretical_mass_to_charge("t:cdot@1:3(+1)", peptidoform=Peptidoform("EGHERWGASRPGGDPSASRWR")))
print(theoretical_mass_to_charge("t:c@1:3(+1)", peptidoform=Peptidoform("EGHERWGASRPGGDPSASRWR")))

print(theoretical_mass_to_charge("x:t@17:21(+1)", peptidoform=Peptidoform("EGHERWGASRPGGDPSASRWR")))
print(theoretical_mass_to_charge("zdot:t@18:21(+1)", peptidoform=Peptidoform("EGHERWGASRPGGDPSASRWR")))
