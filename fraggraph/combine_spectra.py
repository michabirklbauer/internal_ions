import numpy as np
import pandas as pd
import pyteomics.mgf
import pyteomics.mzml
import math

#for testing
import plotly.express as px

def group_mz_values(x, mzd):
    mzdiff = np.diff(x)
    groups = np.cumsum(np.concatenate(([0], mzdiff >= mzd))) + 1
    return groups


def estimate_scattering(file_path = None, spectra_list = None, scattering_type = "mz"):
    """ Plot the distribution of the mz deviation or ppm deviation from each peaks across spectra in the list """

    # Load the spectra

    if file_path is not None:
        spectra = list(load_spectra(file_path))
    elif spectra_list is not None:
        spectra = spectra_list
    else:
        raise ValueError("Please provide either file_path or spectra_list.")

    # Deviation list
    deviations = []

    # For each spectrum:
    for sn in range(len(spectra)):
        spectrum = spectra[sn]

        # For each peak in the spectrum:
        for peak in spectrum['m/z array']:

            # Calculate the deviation to the closest peak in the other spectra
            for sn2 in range(sn, len(spectra)):

                spectrum2 = spectra[sn2]
                if scattering_type == "mz":
                    deviation = np.min(np.abs(peak - spectrum2['m/z array']))
                elif scattering_type == "ppm":
                    #calculate ppm deviation to the closest peak in the other spectra
                    deviation = np.min(np.abs(peak - spectrum2['m/z array']) / peak * 1e6)
                else:
                    raise ValueError("scattering_type must be either 'mz' or 'ppm'.")

                deviations.append(deviation)



            #print("Number of deviations:", len(deviations))

    threshold = 1 if scattering_type == "mz" else 100  # Adjust threshold for mz or ppm
    deviations = [x for x in deviations if x < threshold]

    # Plot the distribution
    bin_width = 0.001 if scattering_type == "mz" else 1  # Adjust bin width for mz or ppm
    nbins = math.ceil((max(deviations) - min(deviations)) / bin_width)
    x_label = 'm/z deviation (Da)' if scattering_type == "mz" else 'ppm deviation'
    histogram_dev = px.histogram(x=deviations, nbins=nbins, labels={'x': x_label})
    #improve to minimal aesthetic, black and white bars, no blue background
    histogram_dev.update_layout(plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)', font_color="black")
    #set bar color to black
    histogram_dev.update_traces(marker_color='black')



    #add a vertical line to the plot at 10ppm or 0.01Da
    if scattering_type == "mz":
        threshold = 0.01
    elif scattering_type == "ppm":
        threshold = 10
    histogram_dev.add_vline(x=threshold, line_width=3, line_dash="dash", line_color="red")


    #add annotation to the vertical line with the percentage of deviation below the threshold
    #calculate the percentage of deviation below the threshold
    percentage_below_threshold = len([x for x in deviations if x < threshold]) / len(deviations) * 100
    histogram_dev.add_annotation(x=threshold, y=0, text=f"{percentage_below_threshold:.2f}%", showarrow=True, arrowhead=1, yshift=10)


    histogram_dev.show()



def load_spectra(file, is_streamlit = False):

    if is_streamlit:
        file_path = file.name
    else:
        file_path = file

     # Read spectra from an .mgf or .mzML file
    if file_path.endswith('.mgf'):
        spectra = pyteomics.mgf.read(file)
    elif file_path.endswith('.mzML'):
        spectra = pyteomics.mzml.read(file)
    else:
        raise ValueError("Unsupported file format. Please provide .mgf or .mzML file.")

    return spectra

def combine_spectra(file_path = None,
                    spectra_list = None,
                    mzd = 0.01,
                    min_prop = 0.5,
                    weighted = False,
                    is_streamlit = False):



    # Load spectra from file
    if file_path is not None:

        spectra = list(load_spectra(file_path, is_streamlit))

        if len(spectra) == 1:
            return spectra[0]

        # Check if all spectra have the same MS level
        if len(set([s['ms level'] for s in spectra])) != 1:
            # If not, keep only MS2 spectra:
            spectra = [s for s in spectra if s['ms level'] == 2]
            print("Warning: not all spectra have the same MS level. Keeping only MS2 spectra.")

        #TODO first combine peak in each spectrum, then combine spectra
        for s in spectra:
            mz_groups = group_mz_values(s['m/z array'], mzd=mzd)
            s['m/z array'] = np.array([np.mean(s['m/z array'][mz_groups == i]) for i in range(1, max(mz_groups) + 1)])
            s['intensity array'] = np.array([np.sum(s['intensity array'][mz_groups == i]) for i in range(1, max(mz_groups) + 1)])

    elif spectra_list is not None:
        spectra = spectra_list
    else:
        raise ValueError("Please provide either file_path or spectra_list.")


    mz_data = np.concatenate([s['m/z array'] for s in spectra])
    intensity_data = np.concatenate([s['intensity array'] for s in spectra])

    keep = intensity_data > 0
    mz_data = mz_data[keep]
    intensity_data = intensity_data[keep]

    mz_order = np.argsort(mz_data)
    mz_data = mz_data[mz_order]
    intensity_data = intensity_data[mz_order]

    mz_groups = group_mz_values(mz_data, mzd=mzd)

    main_spectrum = {}

    main_spectrum['its_mean'] = [np.mean(intensity_data[mz_groups == i]) for i in range(1, max(mz_groups) + 1)]
    main_spectrum['its_median'] = [np.median(intensity_data[mz_groups == i]) for i in range(1, max(mz_groups) + 1)]
    main_spectrum['mz_mean'] = [np.mean(mz_data[mz_groups == i]) for i in range(1, max(mz_groups) + 1)]


    main_spectrum['mz_median'] = [np.median(mz_data[mz_groups == i]) for i in range(1, max(mz_groups) + 1)]
    main_spectrum['its_stdev'] = [np.std(intensity_data[mz_groups == i]) for i in range(1, max(mz_groups) + 1)]
    main_spectrum['n_peaks'] = [np.sum(mz_groups == i) for i in range(1, max(mz_groups) + 1)]
    #TODo check whether mz mean include zero value whe peak is not found in some spectra
    main_spectrum['cov_spectra'] = [np.sum(mz_groups == i) / len(spectra) for i in range(1, max(mz_groups) + 1)]


    #filter out rows by min_prop
    main_spectrum = pd.DataFrame(main_spectrum)
    main_spectrum = main_spectrum[main_spectrum['cov_spectra'] >= min_prop]

    #reset row numbers
    main_spectrum = main_spectrum.reset_index(drop=True)


    return main_spectrum

# Example usage
if __name__ == "__main__":
    input_file = "data/H3_tail_unmod/24unmod_Agc1e5Z10ETD30SA20Uscan3INJT150Res30000.mzML"

    #estimate_scattering(input_file, scattering_type="mz")
    #estimate_scattering(input_file, scattering_type="ppm")


    combined_spectrum = combine_spectra(input_file, mzd=0.15)
    print("mz_mean:", combined_spectrum['mz_mean'])
    print("its_mean:", combined_spectrum['its_mean'])
    print("its_median:", combined_spectrum['its_median'])
    print("n_peaks:", combined_spectrum['n_peaks'])


    #FOR TESTING#
    #plot the spectrum (use plotly)

    #plot as a barplot, with number of peak as annotation on top of the bar and error bars from the standard deviation
    fig = px.bar(x=combined_spectrum['mz_mean'], y=combined_spectrum['its_mean'], error_y=combined_spectrum['its_stdev'], labels={'x': 'm/z', 'y': 'Intensity'})
    #add the number of peaks as annotation
    fig.update_layout(annotations=[dict(x=combined_spectrum['mz_mean'][i], y=combined_spectrum['its_mean'][i], text=str(combined_spectrum['n_peaks'][i]), showarrow=False) for i in range(len(combined_spectrum['mz_mean']))])
    #add the number of spectra in which the peak is present as annotation
    #fig.update_layout(annotations=[dict(x=combined_spectrum['mz_mean'][i], y=combined_spectrum['its_mean'][i], text=str(combined_spectrum['cov_spectra'][i]), showarrow=False) for i in range(len(combined_spectrum['mz_mean']))])
    fig.show()
