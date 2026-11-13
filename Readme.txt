Organisation of Online Flow Cytometer Dataset Analysis
FSC Data

FSC data contain the FSC-A/H, SSC-A/H, and FITC-A/H measurements for every event recorded during a measurement.

The mbiomas toolbox provides functions to read, clean, gate, analyse, and plot experiments consisting of multiple flow cytometry measurements.
Important Data Structures
FCF – FCFiles

A Pandas DataFrame containing metadata about the experiment.

It includes information such as:

    Measurement time and number of events.
    Dilutions.
    Mean and median signal intensities.
    Sampling and regulation information.
    Fractions of gated populations.

FCD – FCData

A dictionary containing one entry for each measurement in the experiment.

Each entry contains the corresponding FSC data, including all measured events remaining after the cleaning procedures.
STALL – Statistics-All Cells

Contains statistical analyses of the datasets.

Each row represents one flow cytometry analysis and contains statistical information such as:

    Percentiles.
    Standard deviations.
    Entropy.
    Other statistical parameters.

Further Documentation

More detailed information about the data structures and how to use them can be found in the documentation of the mbiomas toolbox:

 https://gitlab.uliege.be/mipi/published-software/mbiomas-core/-/tree/92e072ff8032a7048bebea9ffe1738c548b5327a/mbiomas