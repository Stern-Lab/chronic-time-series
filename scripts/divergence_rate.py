############################################## 1. Import Libraries

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statistics import linear_regression
from scipy.stats import linregress
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
import os

############################################## 3. Data Processing Function (Main Function)

# def get_final_df(path, report_df, non_coding_regions, mutation_percentage_non_synonymous, mutation_percentage_synonymous):
#     """
#     Get the final DataFrame with divergence rates for non-synonymous and synonymous mutations.
#     :param path:
#     :param non_coding_regions:
#     :param mutation_percentage_non_synonymous:
#     :param mutation_percentage_synonymous:
#     :return:
#     """
#
#     # List
#     final_df = []
#
#     # A list of all patient_id and their timepoints
#     patient_timepoints = []
#
#     for patient_id_folder in os.listdir(base_directory):
#         patient_id_path = os.path.join(base_directory, patient_id_folder)
#         print(patient_id_path)
#         # Check if it's a directory
#         if os.path.isdir(patient_id_path):
#             # Iterate over timepoint folders
#             for timepoint_folder in os.listdir(patient_id_path):
#                 timepoint_path = os.path.join(patient_id_path, timepoint_folder)
#                 print(timepoint_path)
#                 # Check if it's a directory
#                 if os.path.isdir(timepoint_path):
#                     # Construct the file path
#                     file_path = os.path.join(timepoint_path, f'{patient_id_folder}_T{timepoint_folder}_merged.csv')
#                     print(file_path)
#                     # Read the CSV into a DataFrame
#                     df = pd.read_csv(file_path)
#
#                     # Ensure final_freq column contains floats
#                     df['final_freq'] = df['final_freq'].astype(float)
#
#                     # Extract patient_id and timepoint from the folder names
#                     patient_id = patient_id_folder
#                     timepoint = int(timepoint_folder)  # Convert to integer
#
#                     # Add patient_id and timepoint as new columns and convert timepoint to integer
#                     df['patient_id'] = patient_id
#                     df['timepoint'] = timepoint
#
#                     # Create a patinet_id and timepoint list
#                     patient_timepoints.append((patient_id, timepoint))
#
#                     # Drop rows with final_freq = -1
#                     df = df[df['final_freq'] != -1]
#
#                     # Label non-coding regions
#                     df = label_non_coding_regions(df, non_coding_regions)
#
#                     # Drop non-coding regions and other non-coding mutations
#                     df = df[df['type'] != 'non-coding']
#
#                     # Label synonymous, non-synonymous, and reversion mutations
#                     df = label_syn_non_syn_reversion(df)
#
#                     # Exclude indels (insertions + deletions) from df
#                     df = df[~df['mutation'].str.contains(r'\+|-')]
#
#                     # Count the number of non-synonymous and synonymous mutations (sanity check)
#                     n_non_synonymous, n_synonymous, n_reversion = count_syn_non_syn(df)
#
#                     # Filter out suspected problematic mutations
#                     df, filtered_out_df = TRY_filter_out_suspected_problematic_by_mutation(df, report_df)
#
#                     non_synonymous_divergence = calculate_divergence(df, 'non-synonymous', mutation_percentage_non_synonymous, genome_length_without_non_coding)
#                     synonymous_divergence = calculate_divergence(df, 'synonymous', mutation_percentage_synonymous, genome_length_without_non_coding)
#
#                     print(f"For {patient_id} the Non-synonymous divergence is: {non_synonymous_divergence}")
#                     print(f"For {patient_id} the Synonymous divergence is: {synonymous_divergence}")
#
#                     # Strip mutation column to only include reference and alternate nucleotides
#                     df[['original_nucleotide', 'mutated_nucleotide']] = df['mutation'].str.extract(r'([ACGT])\d+([ACGT])')
#
#                     # Create a new column for mutation_nucleotide
#                     df['mutation_nucleotide'] = df['original_nucleotide'] + df['mutated_nucleotide']
#
#
#                     # Create abriviated protein column for easier analysis
#                     protein_abbreviations = {
#                         'envelope protein': 'E',
#                         'membrane glycoprotein': 'M',
#                         'nucleocapsid phosphoprotein': 'N',
#                         'ORF10 protein': 'ORF10',
#                         'orf1ab polyprotein': 'ORF1ab',
#                         'ORF3a protein': 'ORF3a',
#                         'ORF6 protein': 'ORF6',
#                         'ORF7a protein': 'ORF7a',
#                         'ORF8 protein': 'ORF8',
#                         'surface glycoprotein': 'S'}
#
#                     df['protein'] = df['protein'].map(protein_abbreviations)
#
#                     # Create a new column for protein categories
#                     protein_categories = {
#                         'E': 'structural',
#                         'M': 'structural',
#                         'N': 'structural',
#                         'S': 'structural',
#                         'ORF3a': 'accessory',
#                         'ORF6': 'accessory',
#                         'ORF7a': 'accessory',
#                         'ORF8': 'accessory',
#                         'ORF10': 'accessory',
#                         'ORF1ab': 'replicase'
#                     }
#
#                     df['protein_category'] = df['protein'].map(protein_categories)
#
#                     # Append the final DataFrame to the list
#                     final_df.append(df)
#
#
#
#     # Concatenate the DataFrames in the list
#     final_df = pd.concat(final_df)
#     # Add missing timepoints for each patient
#     final_df = add_missing_timepoints(final_df, patient_timepoints)
#
#     # Save the final DataFrame only with columns of interest to a CSV file
#     final_df = final_df[['patient_id', 'timepoint', 'POS_x', 'mutation','final_freq', 'mutation_type', 'type', 'mutation_nucleotide', 'protein', 'protein_category']]
#     # Save the final DataFrame to a CSV file
#     final_df.to_csv('final_df_with_problematic_filter.csv', index=False)
#     return final_df
def get_final_df(path, report_df, non_coding_regions, mutation_percentage_non_synonymous, mutation_percentage_synonymous):
    """
    Get the final DataFrame with divergence rates for non-synonymous and synonymous mutations.
    :param path:
    :param non_coding_regions:
    :param mutation_percentage_non_synonymous:
    :param mutation_percentage_synonymous:
    :return:
    """

    # List to store final dataframes
    final_df_list = []

    # Collect all timepoints
    all_timepoints = []

    for patient_id_folder in os.listdir(base_directory):
        patient_id_path = os.path.join(base_directory, patient_id_folder)
        print(patient_id_path)
        # Check if it's a directory
        if os.path.isdir(patient_id_path):
            # Patient-specific data
            patient_df_list = []

            for timepoint_folder in os.listdir(patient_id_path):
                timepoint_path = os.path.join(patient_id_path, timepoint_folder)
                print(timepoint_path)
                # Check if it's a directory
                if os.path.isdir(timepoint_path):
                    # Construct the file path
                    file_path = os.path.join(timepoint_path, f'{patient_id_folder}_T{timepoint_folder}_merged.csv')
                    print(file_path)
                    # Read the CSV into a DataFrame
                    df = pd.read_csv(file_path)

                    # Ensure final_freq column contains floats
                    df['final_freq'] = df['final_freq'].astype(float)

                    # Extract patient_id and timepoint from the folder names
                    patient_id = patient_id_folder
                    timepoint = int(timepoint_folder)  # Convert to integer
                    patient_timepoint_list = [(patient_id, timepoint)]
                    # Add patient_id and timepoint as new columns
                    df['patient_id'] = patient_id
                    df['timepoint'] = timepoint

                    # Drop rows with final_freq = -1
                    df = df[df['final_freq'] != -1]

                    # Label non-coding regions
                    df = label_non_coding_regions(df, non_coding_regions)

                    # Drop non-coding regions and other non-coding mutations
                    df = df[df['type'] != 'non-coding']

                    # Label synonymous, non-synonymous, and reversion mutations
                    df = label_syn_non_syn_reversion(df)

                    # Exclude indels (insertions + deletions) from df
                    df = df[~df['mutation'].str.contains(r'\+|-')]

                    # Count the number of non-synonymous and synonymous mutations (sanity check)
                    n_non_synonymous, n_synonymous, n_reversion = count_syn_non_syn(df)

                    # Filter out suspected problematic mutations
                    df, filtered_out_df = TRY_filter_out_suspected_problematic_by_mutation(df, report_df)

                    non_synonymous_divergence = calculate_divergence(df, 'non-synonymous', mutation_percentage_non_synonymous, genome_length_without_non_coding)
                    synonymous_divergence = calculate_divergence(df, 'synonymous', mutation_percentage_synonymous, genome_length_without_non_coding)

                    print(f"For {patient_id} the Non-synonymous divergence is: {non_synonymous_divergence}")
                    print(f"For {patient_id} the Synonymous divergence is: {synonymous_divergence}")

                    # Strip mutation column to only include reference and alternate nucleotides
                    df[['original_nucleotide', 'mutated_nucleotide']] = df['mutation'].str.extract(r'([ACGT])\d+([ACGT])')

                    # Create a new column for mutation_nucleotide
                    df['mutation_nucleotide'] = df['original_nucleotide'] + df['mutated_nucleotide']

                    # Create abriviated protein column for easier analysis
                    protein_abbreviations = {
                        'envelope protein': 'E',
                        'membrane glycoprotein': 'M',
                        'nucleocapsid phosphoprotein': 'N',
                        'ORF10 protein': 'ORF10',
                        'orf1ab polyprotein': 'ORF1ab',
                        'ORF3a protein': 'ORF3a',
                        'ORF6 protein': 'ORF6',
                        'ORF7a protein': 'ORF7a',
                        'ORF8 protein': 'ORF8',
                        'surface glycoprotein': 'S'}

                    df['protein'] = df['protein'].map(protein_abbreviations)

                    # Create a new column for protein categories
                    protein_categories = {
                        'E': 'structural',
                        'M': 'structural',
                        'N': 'structural',
                        'S': 'structural',
                        'ORF3a': 'accessory',
                        'ORF6': 'accessory',
                        'ORF7a': 'accessory',
                        'ORF8': 'accessory',
                        'ORF10': 'accessory',
                        'ORF1ab': 'replicase'
                    }

                    df['protein_category'] = df['protein'].map(protein_categories)

                    # Append the patient's DataFrame to the list
                    patient_df_list.append(df)
                    # Collect the timepoints
                    all_timepoints += patient_timepoint_list

            # Concatenate the patient's DataFrames
            patient_df = pd.concat(patient_df_list)
            # Add missing timepoints
            patient_df = add_missing_timepoints(patient_df, all_timepoints)
            # Append the final patient DataFrame to the list
            final_df_list.append(patient_df)

    # Concatenate the DataFrames in the list
    final_df = pd.concat(final_df_list)
    # Ensure the 'timepoint' column is of integer type
    final_df['timepoint'] = final_df['timepoint'].astype(int)
    # Save the final DataFrame only with columns of interest to a CSV file
    final_df = final_df[['patient_id', 'timepoint', 'POS_x', 'mutation','final_freq', 'mutation_type', 'type', 'mutation_nucleotide', 'protein', 'protein_category']]
    # Save the final DataFrame to a CSV file
    final_df.to_csv(r"C:\Users\natal\OneDrive - Open University of Israel\שולחן העבודה\N1_OLD_reps_iVar_results_(0.01_100_50)_COMBINED_w_TechnionResequencing_W_INDELS\final_df_with_problematic_filter.csv", index=False)
    return final_df




##### Helper Functions #####

def label_non_coding_regions(data, non_coding_regions):
    """
    Labels the 'type' column in the DataFrame as 'non-coding' if the mutation position is within non-coding regions.

    Parameters:
    data (pd.DataFrame): The input DataFrame with mutation positions.
    non_coding_regions (list of tuples): List of tuples where each tuple defines the start and end of a non-coding region.

    Returns:
    pd.DataFrame: The DataFrame with an updated 'type' column.
    """
    # Create a list of non-coding positions
    non_coding = []
    for start, end in non_coding_regions:
        non_coding.extend(list(range(start, end)))

    # Check if the mutation is non-coding
    data['type'] = np.where(data['POS_x'].isin(non_coding), 'non-coding', 'NA')

    return data


def label_syn_non_syn_reversion(data):
    """
    Labels the 'type' column in the DataFrame as 'non-synonymous', 'synonymous', or 'reversion' based on the values of the 'protein_type' and 'mutation_type' columns.

    Parameters:
    data (pd.DataFrame): The input DataFrame with 'protein_type' and 'mutation_type' columns.

    Returns:
    pd.DataFrame: The DataFrame with an updated 'type' column.
    """
    types = []

    for idx, row in data.iterrows():
        protein = row['protein']
        mutation_type = row['mutation_type']

        if pd.notna(protein) and pd.notna(mutation_type):
            types.append('non-synonymous')
        elif pd.notna(protein) and pd.isna(mutation_type):
            types.append('synonymous')
        elif pd.isna(protein) and pd.isna(mutation_type):
            types.append('reversion')
        else:
            types.append('unknown')  # Added to handle any unexpected cases

    data['type'] = types

    return data

def count_syn_non_syn(data):
    """
    Count the number of non-synonymous and synonymous mutations in the data.
    :param data: Input DataFrame with 'type' column.
    :return: Returns the number of non-synonymous and synonymous mutations.
    """

    n_non_synonymous = data[data['type'] == 'non-synonymous'].shape[0]
    n_synonymous = data[data['type'] == 'synonymous'].shape[0]
    n_reversion = data[data['type'] == 'reversion'].shape[0]

    return n_non_synonymous, n_synonymous, n_reversion


def calculate_divergence(data, mutation_type, mutation_percentage, genome_length):
    """
    Calculate the divergence values for a given mutation type.
    :param genome_length: SARS-CoV-2 full genome length without non-coding regions
    :param mutation_percentage: The percentage of the genome where the mutation type can occur.
    :param data: Input DataFrame with 'type', 'final_freq' and 'timepoint' columns.
    :param mutation_type: The mutation type ('non-synonymous' or 'synonymous').
    :param n: The total number of possible mutations for the given type.
    :return:
    """
    # Calculate the normalization factor
    n = genome_length * mutation_percentage

    subset = data[data['type'] == mutation_type]
    summed_frequencies = subset.groupby('timepoint')['final_freq'].sum()
    divergence_values_normalized = summed_frequencies / n

    return divergence_values_normalized

def calculate_divergence_by_gene(data, mutation_type):


    # Filter the data based on the mutation type
    subset = data[data['type'] == mutation_type]

    # Group the data by gene and timepoint
    grouped = subset.groupby(['protein', 'timepoint'])['final_freq'].sum().reset_index()

    return grouped


def add_missing_timepoints(data, timepoints):
    """
    Add missing timepoints to the DataFrame.
    :param data: Input DataFrame with 'timepoint' column.
    :param timepoints: List of all timepoints.
    :return: Returns the DataFrame with missing timepoints added.
    """


   # Iterate over each patient
    for patient, patient_data in data.groupby('patient_id'):
        # For "timepoints" iterate over timepoints relevant for the current patient

        patient_data = data[data['patient_id'] == patient]
        unique_timepoints = patient_data['timepoint'].unique()


        # Filter the timepoints to include only those relevant for the patient id
        patient_timepoints = [tp for pid, tp in timepoints if pid == patient]
        missing_timepoints = list(set(patient_timepoints) - set(unique_timepoints))

        # Create a new DataFrame with missing timepoints
        missing_data = pd.DataFrame(missing_timepoints, columns=['timepoint'])
        missing_data['patient_id'] = patient
        missing_data['final_freq'] = 0

        # Merge the original DataFrame with the missing data
        data = pd.concat([data, missing_data], ignore_index=True)

    return data

def TRY_filter_out_suspected_problematic_by_mutation(df, report_df):

    # Extract all suspected problematic mutations from the report_df
    sus_mutations = report_df['mutation'].tolist()

    # Filter out the suspected problematic mutations from the main DataFrame
    filtered_out_df = df[df['mutation'].isin(sus_mutations)]
    filtered_df = df[~df['mutation'].isin(sus_mutations)]

    return filtered_df, filtered_out_df

##### Plotting Function #####

# def plot_divergence_rates(final_df, genome_length_without_non_coding, mutation_percentage_non_synonymous, mutation_percentage_synonymous):
#     """
#     Plot divergence rates for synonymous and non-synonymous mutations over time for each patient.
#     :param final_df: DataFrame containing the processed mutation data.
#     :param genome_length_without_non_coding: SARS-CoV-2 genome length without non-coding regions.
#     :param mutation_percentage_non_synonymous: The percentage of the genome where non-synonymous mutations can occur.
#     :param mutation_percentage_synonymous: The percentage of the genome where synonymous mutations can occur.
#     """
#
#     # plt.style.use('ggplot')
#
#     # Get a list of unique patient IDs
#     patients = final_df['patient_id'].unique()
#     num_patients = len(patients)
#
#     # Create subplots with 2 rows (one for synonymous and one for non-synonymous) and columns equal to number of patients
#     fig, axs = plt.subplots(2, num_patients, figsize=(15, 10), sharex=False, sharey=True)
#     fig.suptitle('Divergence Rates Over Time')
#
#     for i, patient in enumerate(patients):
#         patient_data = final_df[final_df['patient_id'] == patient]
#         all_time_points = list(set(patient_data['timepoint'].astype(int)))
#
#
#
#         # Calculate normalized divergence values
#         syn_data_normalized = calculate_divergence(patient_data, 'synonymous', mutation_percentage_synonymous, genome_length_without_non_coding)
#         non_syn_data_normalized = calculate_divergence(patient_data, 'non-synonymous', mutation_percentage_non_synonymous, genome_length_without_non_coding)
#
#         # Plot data for the current patient
#         axs[0, i].plot(syn_data_normalized.index, syn_data_normalized.values, marker='o', linestyle='-')
#         # axs[0, i].set_xticks(syn_data_normalized.index)
#         axs[0, i].set_title(f'Patient {patient}')
#         axs[0, i].set_ylim(0, None)  # Force start from 0
#
#         if i == 0:
#             axs[0, i].set_ylabel('Normalized Synonymous Frequencies')
#
#         axs[1, i].plot(non_syn_data_normalized.index, non_syn_data_normalized.values, marker='o', linestyle='-')
#         # axs[1, i].set_xticks(non_syn_data_normalized.index)  # Set x-axis limits based on patient data
#         axs[1, i].set_ylim(0, None)  # Force start from 0
#
#         if i == 0:
#             axs[1, i].set_ylabel('Normalized Non-Synonymous Frequencies')
#
#         plt.xticks(all_time_points, all_time_points)
#
#     # Add a single x-axis label for the whole figure
#     fig.text(0.5, 0.04, 'Time (days)', ha='center', va='center')
#
#     plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])
#     plt.show()


def linear_regression(x, y):
    model = LinearRegression(fit_intercept=False)
    model.fit(x.reshape(-1, 1), y)
    return model.coef_[0]


# def plot_normalized_divergence_rates_per_patient(final_df, genome_length_without_non_coding, mutation_percentage_non_synonymous, mutation_percentage_synonymous):
#     """
#     Plot divergence rates for synonymous and non-synonymous mutations over time for each patient.
#     :param final_df: DataFrame containing the processed mutation data.
#     :param genome_length_without_non_coding: SARS-CoV-2 genome length without non-coding regions.
#     :param mutation_percentage_non_synonymous: The percentage of the genome where non-synonymous mutations can occur.
#     :param mutation_percentage_synonymous: The percentage of the genome where synonymous mutations can occur.
#     """
#
#     # Get a list of unique patient IDs
#     patients = final_df['patient_id'].unique()
#     num_patients = len(patients)
#
#     # Set the font size for the plot
#     SMALL_SIZE = 10
#     MEDIUM_SIZE = 15
#     BIGGER_SIZE = 20
#     TITLE_SIZE = 25
#
#     plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
#     plt.rc('axes', titlesize=MEDIUM_SIZE)  # fontsize of the axes title
#     plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
#     plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
#     plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
#     plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
#     plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title
#
#     # # Filter final_df to contain only "ORF1ab" protein
#     # final_df = final_df[final_df['protein'] == 'ORF1ab']
#
#     # Create subplots with 2 rows (one for synonymous and one for non-synonymous) and columns equal to number of patients
#     fig, axs = plt.subplots(2, num_patients, figsize=(28, 13), sharex=False, sharey=True)
#     fig.suptitle('Divergence Rates Over Time')
#
#     for i, patient in enumerate(patients):
#         patient_data = final_df[final_df['patient_id'] == patient]
#         all_time_points = sorted(patient_data['timepoint'].astype(int).unique())
#
#         # Calculate normalized divergence values
#         syn_data_normalized = calculate_divergence(patient_data, 'synonymous', mutation_percentage_synonymous, genome_length_without_non_coding)
#         non_syn_data_normalized = calculate_divergence(patient_data, 'non-synonymous', mutation_percentage_non_synonymous, genome_length_without_non_coding)
#
#
#         # Plot data for the current patient
#         axs[0, i].scatter(syn_data_normalized.index, syn_data_normalized.values, marker='o', linestyle='-')
#         axs[0, i].set_title(f'Patient {patient}')
#         axs[0, i].set_ylim(0, None)  # Force start from 0
#
#         if i == 0:
#             axs[0, i].set_ylabel('Normalized Synonymous Frequencies')
#
#         axs[1, i].scatter(non_syn_data_normalized.index, non_syn_data_normalized.values, marker='o', linestyle='-')
#         axs[1, i].set_ylim(0, None)  # Force start from 0
#
#         # Calculate and plot linear regression line through the origin for non-synonymous
#         if len(non_syn_data_normalized) > 1:
#             try:
#                 slope_non_syn, _ = linear_regression(non_syn_data_normalized.index, non_syn_data_normalized.values)
#                 x_vals = np.array(axs[1, i].get_xlim())
#                 y_vals = slope_non_syn * x_vals
#                 axs[1, i].plot(x_vals, y_vals, color='r', linestyle='--')
#             except Exception as e:
#                 pass
#
#         if i == 0:
#             axs[1, i].set_ylabel('Normalized Non-Synonymous Frequencies')
#
#         # Set x-ticks and labels for both subplots
#         axs[0, i].set_xticks(all_time_points)
#         axs[1, i].set_xticks(all_time_points)
#
#     fig.text(0.5, 0.01, 'Time (days)', fontsize=MEDIUM_SIZE)
#     plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])
#     plt.savefig(
#         fr"C:\Users\natal\OneDrive - Open University of Israel\שולחן העבודה\N1_OLD_reps_iVar_results_(0.01_100_50)_COMBINED_w_TechnionResequencing_W_INDELS\try1.png",
#         bbox_inches="tight", dpi=600)
#     plt.show()
#
#     # Plot bar plot of regression slopes
#     fig, ax = plt.subplots(figsize=(10, 6))
#     bar_width = 0.35
#     index = np.arange(len(patients))
#
#     bar1 = ax.bar(index, syn_slopes, bar_width, label='Synonymous')
#     bar2 = ax.bar(index + bar_width, non_syn_slopes, bar_width, label='Non-Synonymous')
#
#     ax.set_xlabel('Patient ID')
#     ax.set_ylabel('Regression Slope')
#     ax.set_title('Regression Slopes by Patient')
#     ax.set_xticks(index + bar_width / 2)
#     ax.set_xticklabels(patients, rotation=45)
#     ax.legend()
#
#     plt.tight_layout()
#     plt.savefig(
#         fr"C:\Users\natal\OneDrive - Open University of Israel\שולחן העבודה\N1_OLD_reps_iVar_results_(0.01_100_50)_COMBINED_w_TechnionResequencing_W_INDELS\bar_plot_try1.png", bbox_inches="tight", dpi=600)
#     plt.show()

def plot_normalized_divergence_rates_per_patient_with_regression(final_df, genome_length_without_non_coding,
                                                 mutation_percentage_non_synonymous, mutation_percentage_synonymous):
    patients = final_df['patient_id'].unique()
    num_patients = len(patients)

    SMALL_SIZE = 10
    MEDIUM_SIZE = 15
    BIGGER_SIZE = 20
    TITLE_SIZE = 25

    plt.rc('font', size=SMALL_SIZE)
    plt.rc('axes', titlesize=MEDIUM_SIZE)
    plt.rc('axes', labelsize=MEDIUM_SIZE)
    plt.rc('xtick', labelsize=SMALL_SIZE)
    plt.rc('ytick', labelsize=SMALL_SIZE)
    plt.rc('legend', fontsize=SMALL_SIZE)
    plt.rc('figure', titlesize=MEDIUM_SIZE)



    fig, axs = plt.subplots(2, num_patients, figsize=(28, 13), sharex=False, sharey=True)
    fig.suptitle('Divergence Rates Over Time')

    syn_slopes = []
    non_syn_slopes = []
    summary_data = []

    for i, patient in enumerate(patients):
        patient_data = final_df[final_df['patient_id'] == patient]
        all_time_points = sorted(patient_data['timepoint'].astype(int).unique())

        syn_data_normalized = calculate_divergence(patient_data, 'synonymous', mutation_percentage_synonymous,
                                                   genome_length_without_non_coding)
        non_syn_data_normalized = calculate_divergence(patient_data, 'non-synonymous',
                                                       mutation_percentage_non_synonymous,
                                                       genome_length_without_non_coding)

        axs[0, i].scatter(syn_data_normalized.index, syn_data_normalized.values, marker='o', linestyle='-')
        axs[0, i].set_title(f'Patient {patient}')
        axs[0, i].set_ylim(0, None)
        if i == 0:
            axs[0, i].set_ylabel('Normalized Synonymous Frequencies')

        if len(syn_data_normalized) > 1:
            try:
                slope_syn = linear_regression(syn_data_normalized.index.values, syn_data_normalized.values)
                syn_slopes.append(slope_syn)
                x_vals = np.array(axs[0, i].get_xlim())
                y_vals = slope_syn * x_vals
                axs[0, i].plot(x_vals, y_vals, color='r', linestyle='--')
            except Exception as e:
                syn_slopes.append("NA")
        else:
            syn_slopes.append("NA")


        axs[1, i].scatter(non_syn_data_normalized.index, non_syn_data_normalized.values, marker='o', linestyle='-')
        axs[1, i].set_ylim(0, None)

        if len(non_syn_data_normalized) > 1:
            try:
                slope_non_syn = linear_regression(non_syn_data_normalized.index.values, non_syn_data_normalized.values)
                non_syn_slopes.append(slope_non_syn)
                x_vals = np.array(axs[1, i].get_xlim())
                y_vals = slope_non_syn * x_vals
                axs[1, i].plot(x_vals, y_vals, color='r', linestyle='--')
            except Exception as e:
                non_syn_slopes.append("NA")
        else:
            non_syn_slopes.append("NA")

        if i == 0:
            axs[1, i].set_ylabel('Normalized Non-Synonymous Frequencies')

        axs[0, i].set_xticks(all_time_points)
        axs[1, i].set_xticks(all_time_points)

        summary_data.append({
            'patient': patient,
            'syn_slope': syn_slopes[-1],
            'non_syn_slope': non_syn_slopes[-1]
        })

    fig.text(0.5, 0.01, 'Time (days)', fontsize=MEDIUM_SIZE)
    plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])
    plt.savefig("divergence_rates_V3.png", bbox_inches="tight", dpi=600)
    plt.show()

    fig, ax = plt.subplots(figsize=(10, 6))
    bar_width = 0.35
    index = np.arange(len(patients))

    # bar1 = ax.bar(index, syn_slopes, bar_width, label='Synonymous')
    # bar2 = ax.bar(index + bar_width, non_syn_slopes, bar_width, label='Non-Synonymous')

    ax.set_xlabel('Patient ID')
    ax.set_ylabel('Regression Slope')
    ax.set_title('Regression Slopes by Patient')
    ax.set_xticks(index + bar_width / 2)
    ax.set_xticklabels(patients, rotation=45)
    ax.legend()

    plt.tight_layout()
    # plt.savefig("regression_slopes.png", bbox_inches="tight", dpi=600)
    plt.show()

    summary_df = pd.DataFrame(summary_data)
    summary_df.to_excel("divergence_summary_V3.xlsx", index=False)

def plot_normalized_divergence_rates_per_patient_with_regression_P5_75_mod(final_df, genome_length_without_non_coding,
                                                    mutation_percentage_non_synonymous, mutation_percentage_synonymous):

        patients = final_df['patient_id'].unique()
        num_patients = len(patients)

        SMALL_SIZE = 10
        MEDIUM_SIZE = 15
        BIGGER_SIZE = 20
        TITLE_SIZE = 25

        plt.rc('font', size=SMALL_SIZE)
        plt.rc('axes', titlesize=MEDIUM_SIZE)
        plt.rc('axes', labelsize=MEDIUM_SIZE)
        plt.rc('xtick', labelsize=SMALL_SIZE)
        plt.rc('ytick', labelsize=SMALL_SIZE)
        plt.rc('legend', fontsize=SMALL_SIZE)
        plt.rc('figure', titlesize=MEDIUM_SIZE)

        fig, axs = plt.subplots(2, num_patients, figsize=(28, 13), sharex=False, sharey=True)
        fig.suptitle('Divergence Rates Over Time')

        syn_slopes = []
        non_syn_slopes = []
        summary_data = []

        for i, patient in enumerate(patients):
            patient_data = final_df[final_df['patient_id'] == patient]
            all_time_points = sorted(patient_data['timepoint'].astype(int).unique())

            # Modify the non_syn_data_normalized for P5 at timepoint 75
            if patient == "P5" and 75 in all_time_points:
                timepoint_68_value = patient_data.loc[patient_data['timepoint'] == 68, 'non-synonymous'].values
                if len(timepoint_68_value) > 0:
                    patient_data.loc[patient_data['timepoint'] == 75, 'non-synonymous'] = timepoint_68_value[0]

            syn_data_normalized = calculate_divergence(patient_data, 'synonymous', mutation_percentage_synonymous,
                                                       genome_length_without_non_coding)
            non_syn_data_normalized = calculate_divergence(patient_data, 'non-synonymous',
                                                           mutation_percentage_non_synonymous,
                                                           genome_length_without_non_coding)

            axs[0, i].scatter(syn_data_normalized.index, syn_data_normalized.values, marker='o', linestyle='-')
            axs[0, i].set_title(f'Patient {patient}')
            axs[0, i].set_ylim(0, None)
            if i == 0:
                axs[0, i].set_ylabel('Normalized Synonymous Frequencies')

            if len(syn_data_normalized) > 1:
                try:
                    slope_syn = linear_regression(syn_data_normalized.index.values, syn_data_normalized.values)
                    syn_slopes.append(slope_syn)
                    x_vals = np.array(axs[0, i].get_xlim())
                    y_vals = slope_syn * x_vals
                    axs[0, i].plot(x_vals, y_vals, color='r', linestyle='--')
                except Exception as e:
                    syn_slopes.append("NA")
            else:
                syn_slopes.append("NA")

            axs[1, i].scatter(non_syn_data_normalized.index, non_syn_data_normalized.values, marker='o', linestyle='-')
            axs[1, i].set_ylim(0, None)

            if len(non_syn_data_normalized) > 1:
                try:
                    slope_non_syn = linear_regression(non_syn_data_normalized.index.values,
                                                      non_syn_data_normalized.values)
                    non_syn_slopes.append(slope_non_syn)
                    x_vals = np.array(axs[1, i].get_xlim())
                    y_vals = slope_non_syn * x_vals
                    axs[1, i].plot(x_vals, y_vals, color='r', linestyle='--')
                except Exception as e:
                    non_syn_slopes.append("NA")
            else:
                non_syn_slopes.append("NA")

            if i == 0:
                axs[1, i].set_ylabel('Normalized Non-Synonymous Frequencies')

            axs[0, i].set_xticks(all_time_points)
            axs[1, i].set_xticks(all_time_points)

            summary_data.append({
                'patient': patient,
                'syn_slope': syn_slopes[-1],
                'non_syn_slope': non_syn_slopes[-1]
            })

        fig.text(0.5, 0.01, 'Time (days)', fontsize=MEDIUM_SIZE)
        plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])
        plt.savefig("divergence_rates.png", bbox_inches="tight", dpi=600)
        plt.show()

        fig, ax = plt.subplots(figsize=(10, 6))
        bar_width = 0.35
        index = np.arange(len(patients))

        # bar1 = ax.bar(index, syn_slopes, bar_width, label='Synonymous')
        # bar2 = ax.bar(index + bar_width, non_syn_slopes, bar_width, label='Non-Synonymous')

        ax.set_xlabel('Patient ID')
        ax.set_ylabel('Regression Slope')
        ax.set_title('Regression Slopes by Patient')
        ax.set_xticks(index + bar_width / 2)
        ax.set_xticklabels(patients, rotation=45)
        ax.legend()

        plt.tight_layout()
        # plt.savefig("regression_slopes.png", bbox_inches="tight", dpi=600)
        plt.show()

        summary_df = pd.DataFrame(summary_data)
        summary_df.to_excel("divergence_summary_V2.xlsx", index=False)


def plot_divergence_by_gene(final_df):
    """

    :param final_df:
    :return:
    """


    # Drop rows with protein = 'NA'
    final_df = final_df.dropna(subset=['protein'])


    # Set the font size for the plot
    SMALL_SIZE = 10
    MEDIUM_SIZE = 15
    BIGGER_SIZE = 20
    TITLE_SIZE = 25

    plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title
    # Get a list of unique gene types using numpy
    gene_types = np.unique(final_df['protein'])
    # gene_types = final_df['protein'].unique()
    num_genes = len(gene_types)

    # Create subplots with 2 rows (one for synonymous and one for non-synonymous) and columns equal to number of genes
    fig, axs = plt.subplots(2, num_genes, figsize=(15, 10), sharex=True, sharey=True)
    fig.suptitle('Divergence Rates Aggregated by Gene Type')

    syn_slopes = []
    non_syn_slopes = []


    for i, gene in enumerate(gene_types):
        gene_data = final_df[final_df['protein'] == gene]

        # Calculate normalized divergence values
        syn_data_aggregated = calculate_divergence_by_gene(gene_data, 'synonymous')
        non_syn_data_aggregated = calculate_divergence_by_gene(gene_data, 'non-synonymous')

        # Plot data for the current gene type
        axs[0, i].scatter(syn_data_aggregated['timepoint'], syn_data_aggregated["final_freq"])
        axs[0, i].set_title(f'{gene}')

        # Set x and y limits for the plots
        axs[0, i].set_xlim(0, 250)
        axs[0, i].set_ylim(0, 8)
        axs[0, i].set_xticks(np.arange(0, 251, 70))

        # Calculate and plot linear regression line through the origin for synonymous
        if len(syn_data_aggregated) > 1:
            try:
                slope_syn, _ = linear_regression(syn_data_aggregated['timepoint'], syn_data_aggregated["final_freq"])
                syn_slopes.append(slope_syn)
                x_vals = np.array(axs[0, i].get_xlim())
                y_vals = slope_syn * x_vals
                axs[0, i].plot(x_vals, y_vals, color='r', linestyle='--')

            except Exception as e:
                syn_slopes.append(np.nan)
        else:
            syn_slopes.append(np.nan)

        if i == 0:
            axs[0, i].set_ylabel('Sum of Synonymous Frequencies')

        axs[1, i].scatter(non_syn_data_aggregated["timepoint"], non_syn_data_aggregated["final_freq"])

        # Set x and y limits for the plots
        axs[1, i].set_xlim(0, 250)
        axs[1, i].set_ylim(0, 8)
        axs[1, i].set_xticks(np.arange(0, 251, 70))

        # Calculate and plot linear regression line through the origin for non-synonymous
        if len(non_syn_data_aggregated) > 1:
            try:
                slope_non_syn, _ = linear_regression(non_syn_data_aggregated["timepoint"],
                                                     non_syn_data_aggregated["final_freq"])
                non_syn_slopes.append(slope_non_syn)
                x_vals = np.array(axs[1, i].get_xlim())
                y_vals = slope_non_syn * x_vals
                axs[1, i].plot(x_vals, y_vals, color='r', linestyle='--')

            except Exception as e:
                non_syn_slopes.append(np.nan)

        else:
            non_syn_slopes.append(np.nan)
        if i == 0:
            axs[1, i].set_ylabel('Sum of Non-Synonymous Frequencies')




    # Set a single x-axis label for all subplots

    fig.text(0.5, 0.01, 'Gene Type', fontsize=MEDIUM_SIZE)
    plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])
    plt.savefig(
        fr"C:\Users\natal\OneDrive - Open University of Israel\שולחן העבודה\N1_OLD_reps_iVar_results_(0.01_100_50)_COMBINED_w_TechnionResequencing_W_INDELS\try3_linear_regression_DPI_600_LONG.png",
        bbox_inches="tight", dpi=600)
    plt.show()

 # Plot bar plot of regression slopes
    fig, ax = plt.subplots(figsize=(10, 6))
    bar_width = 0.35
    index = np.arange(len(gene_types))

    bar1 = ax.bar(index, syn_slopes, bar_width, label='Synonymous')
    bar2 = ax.bar(index + bar_width, non_syn_slopes, bar_width, label='Non-Synonymous')

    ax.set_xlabel('Gene Type')
    ax.set_ylabel('Regression Slope')
    ax.set_title('Regression Slopes by Gene Type')
    ax.set_xticks(index + bar_width / 2)
    ax.set_xticklabels(gene_types, rotation=45)
    ax.legend()

    plt.tight_layout()
    plt.savefig(fr"C:\Users\natal\OneDrive - Open University of Israel\שולחן העבודה\N1_OLD_reps_iVar_results_(0.01_100_50)_COMBINED_w_TechnionResequencing_W_INDELS\slopes_linear_regression_DPI_600_LONG.png", bbox_inches="tight", dpi=600)
    plt.show()

def plot_divergence_by_gene_per_patient(final_df):
    """

    :param final_df:
    :return:
    """


    # Drop rows with protein = 'NA'
    # final_df = final_df.dropna(subset=['protein'])


    # Keep only the 'ORF1ab' protein
    final_df = final_df[final_df['protein'] == 'ORF1ab']

    # Get a list of unique patient IDs
    patients = final_df['patient_id'].unique()
    num_patients = len(patients)

    # Set the font size for the plot
    SMALL_SIZE = 10
    MEDIUM_SIZE = 15
    BIGGER_SIZE = 20
    TITLE_SIZE = 25

    plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title
    # Get a list of unique gene types using numpy



    # Create subplots with 2 rows (one for synonymous and one for non-synonymous) and columns equal to number of genes
    fig, axs = plt.subplots(2, num_patients, figsize=(15, 10), sharex=False, sharey=False)
    fig.suptitle('ORF1ab')

    syn_slopes = []
    non_syn_slopes = []

    for i, patient in enumerate(patients):
        patient_data = final_df[final_df['patient_id'] == patient]
        all_time_points = sorted(patient_data['timepoint'].astype(int).unique())
        # Calculate normalized divergence values
        syn_data_aggregated = calculate_divergence_by_gene(patient_data, 'synonymous')
        non_syn_data_aggregated = calculate_divergence_by_gene(patient_data, 'non-synonymous')

        # Plot data for the current patient
        axs[0, i].scatter(syn_data_aggregated['timepoint'], syn_data_aggregated["final_freq"])
        axs[0, i].set_title(f'Patient {patient}')
        # Set x and y limits for the plots
        axs[1, i].set_xlim(0, 250)
        axs[1, i].set_ylim(0, 8)
        axs[1, i].set_xticks(np.arange(0, 251, 70))

        # Calculate and plot linear regression line through the origin for synonymous
        if len(syn_data_aggregated) > 1:
            try:
                slope_syn, _ = linear_regression(syn_data_aggregated['timepoint'], syn_data_aggregated["final_freq"])
                syn_slopes.append(slope_syn)
                x_vals = np.array(axs[0, i].get_xlim())
                y_vals = slope_syn * x_vals
                axs[0, i].plot(x_vals, y_vals, color='r', linestyle='--')
            except Exception as e:
                syn_slopes.append(np.nan)
        else:
            syn_slopes.append(np.nan)

        if i == 0:
            axs[0, i].set_ylabel('Sum of Synonymous Frequencies')

        axs[1, i].scatter(non_syn_data_aggregated["timepoint"], non_syn_data_aggregated["final_freq"])
        # Set x and y limits for the plots
        axs[1, i].set_xlim(0, 250)
        axs[1, i].set_ylim(0, 8)
        axs[1, i].set_xticks(np.arange(0, 251, 70))


        # Calculate and plot linear regression line through the origin for non-synonymous
        if len(non_syn_data_aggregated) > 1:
            try:
                slope_non_syn, _ = linear_regression(non_syn_data_aggregated["timepoint"],
                                                     non_syn_data_aggregated["final_freq"])
                non_syn_slopes.append(slope_non_syn)
                x_vals = np.array(axs[1, i].get_xlim())
                y_vals = slope_non_syn * x_vals
                axs[1, i].plot(x_vals, y_vals, color='r', linestyle='--')
            except Exception as e:
                non_syn_slopes.append(np.nan)
        else:
            non_syn_slopes.append(np.nan)

        if i == 0:
            axs[1, i].set_ylabel('Sum of Non-Synonymous Frequencies')



        # Set a single x-axis label for all subplots
    fig.text(0.5, 0.01, 'Time (days)', fontsize=MEDIUM_SIZE)
    plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])

    plt.show()

    # # Plot bar plot of regression slopes
    # fig, ax = plt.subplots(figsize=(10, 6))
    # bar_width = 0.35
    # index = np.arange(len(patient_ids))
    #
    # bar1 = ax.bar(index, syn_slopes, bar_width, label='Synonymous')
    # bar2 = ax.bar(index + bar_width, non_syn_slopes, bar_width, label='Non-Synonymous')
    #
    # ax.set_xlabel('Patient ID')
    # ax.set_ylabel('Regression Slope')
    # ax.set_title('Regression Slopes by Patient')
    # ax.set_xticks(index + bar_width / 2)
    # ax.set_xticklabels(patient_ids, rotation=45)
    # ax.legend()
    #
    # plt.tight_layout()
    #
    # plt.show()


# A scatter plot of all synonymous mutations aggregated per timepoint for each patient.
# One plot, where the x-axis represents the time and the y-axis represents the sum of synonymous frequencies in a normalized manner.
# Each dot is the normalized sum of frequencies of synonymous mutations at a specific timepoint for a specific patient.
def plot_non_synonymous_divergence_all_patients(final_df, genome_length_without_non_coding, mutation_percentage_non_synonymous):
    """
    Plot divergence rates for synonymous mutations over time for all patients in a single scatter plot.
    :param final_df:
    :param genome_length_without_non_coding:
    :param mutation_percentage_synonymous:
    :return:
    """

    # Get a list of unique patient IDs
    patients = final_df['patient_id'].unique()

    # Set the font size for the plot
    SMALL_SIZE = 30
    MEDIUM_SIZE = 40
    BIGGER_SIZE = 40
    TITLE_SIZE = 45


    plt.rc('font', size=30)  # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=20)  # legend fontsize
    plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title

    fig, ax = plt.subplots(figsize=(20, 15))


    # Filter df to contain only timepoint smaller than 40
    # final_df = final_df[final_df['timepoint'] <= 40]

    all_time_points = []
    all_syn_values = []

    for patient in patients:
        patient_data = final_df[final_df['patient_id'] == patient]

        # Calculate normalized divergence values
        syn_data_normalized = calculate_divergence(patient_data, 'non-synonymous', mutation_percentage_non_synonymous,
                                                   genome_length_without_non_coding)

        all_time_points.extend(syn_data_normalized.index)
        all_syn_values.extend(syn_data_normalized.values)
        # Plot data for the current patient
        ax.scatter(syn_data_normalized.index, syn_data_normalized.values, marker='o', linestyle='-', label=f'{patient}', s=150)

    # Perform linear regression through the origin using scikit-learn
    model = LinearRegression(fit_intercept=False)
    all_time_points = np.array(all_time_points).reshape(-1, 1)
    model.fit(all_time_points, all_syn_values)
    slope = model.coef_[0]
    intercept = model.intercept_


    regression_line = slope * np.array(all_time_points)

    # Calculate R^2
    r_squared = model.score(all_time_points, all_syn_values)



    # Perform linear regression using statsmodels to get the p-value
    all_time_points_sm = sm.add_constant(all_time_points, has_constant='add')
    sm_model = sm.OLS(all_syn_values, all_time_points_sm).fit()
    p_value = sm_model.pvalues[1]  # p-value for the slope

    # Plot the regression line
    ax.plot(all_time_points, regression_line, color='red', linestyle='-', linewidth=2,
            label=f'Regression Line $R^2$ = {r_squared:.2f}\np = {p_value:.2e}')

    # Debugging print statements
    print(f"Slope: {slope}, Intercept: {intercept}")
    print(f"R^2: {r_squared:.2f}")
    print(f"p-value: {p_value:.2e}")

    # Display the regression formula
    regression_formula = f'y = {slope:.2e}x'
    ax.text(0.95, 0.05, regression_formula, transform=ax.transAxes, fontsize=20, verticalalignment='bottom',
            horizontalalignment='right', bbox=dict(facecolor='white', alpha=0.5))

    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Non-Synonymous divergence')
    ax.set_ylim(0, None)  # Force start from 0
    ax.legend(title='Patients')
    plt.savefig(fr"C:\Users\natal\OneDrive - Open University of Israel\שולחן העבודה\N1_OLD_reps_iVar_results_(0.01_100_50)_COMBINED_w_TechnionResequencing_W_INDELS\non-syn_all_tp.png", bbox_inches="tight", dpi=600)
    plt.show()


def plot_synonymous_divergence_all_patients(final_df, genome_length_without_non_coding, mutation_percentage_synonymous):
    """
    Plot divergence rates for synonymous mutations over time for all patients in a single scatter plot.
    :param final_df:
    :param genome_length_without_non_coding:
    :param mutation_percentage_synonymous:
    :return:
    """

    # Get a list of unique patient IDs
    patients = final_df['patient_id'].unique()

    # Set the font size for the plot
    SMALL_SIZE = 30
    MEDIUM_SIZE = 40
    BIGGER_SIZE = 40
    TITLE_SIZE = 45


    plt.rc('font', size=30)  # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=20)  # legend fontsize
    plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title

    fig, ax = plt.subplots(figsize=(20, 15))


    # Filter df to contain only timepoint smaller than X
    final_df = final_df[final_df['timepoint'] <= 10]

    all_time_points = []
    all_syn_values = []

    for patient in patients:
        patient_data = final_df[final_df['patient_id'] == patient]

        # Calculate normalized divergence values
        syn_data_normalized = calculate_divergence(patient_data, 'synonymous', mutation_percentage_synonymous,
                                                   genome_length_without_non_coding)

        all_time_points.extend(syn_data_normalized.index)
        all_syn_values.extend(syn_data_normalized.values)
        # Plot data for the current patient
        ax.scatter(syn_data_normalized.index, syn_data_normalized.values, marker='o', linestyle='-', label=f'{patient}', s=150)

    # Perform linear regression through the origin using scikit-learn
    model = LinearRegression(fit_intercept=False)
    all_time_points = np.array(all_time_points).reshape(-1, 1)
    model.fit(all_time_points, all_syn_values)
    slope = model.coef_[0]
    intercept = model.intercept_


    regression_line = slope * np.array(all_time_points)

    # Calculate R^2
    r_squared = model.score(all_time_points, all_syn_values)



    # Perform linear regression using statsmodels to get the p-value
    all_time_points_sm = sm.add_constant(all_time_points, has_constant='add')
    sm_model = sm.OLS(all_syn_values, all_time_points_sm).fit()
    p_value = sm_model.pvalues[1]  # p-value for the slope

    # Plot the regression line
    ax.plot(all_time_points, regression_line, color='red', linestyle='-', linewidth=2,
            label=f'Regression Line $R^2$ = {r_squared:.2f}\np = {p_value:.2e}')

    # Debugging print statements
    print(f"Slope: {slope}, Intercept: {intercept}")
    print(f"R^2: {r_squared:.2f}")
    print(f"p-value: {p_value:.2e}")

    # Display the regression formula
    regression_formula = f'y = {slope:.2e}x'
    ax.text(0.95, 0.05, regression_formula, transform=ax.transAxes, fontsize=20, verticalalignment='bottom',
            horizontalalignment='right', bbox=dict(facecolor='white', alpha=0.5))

    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Non-Synonymous divergence')
    ax.set_ylim(0, None)  # Force start from 0
    ax.legend(title='Patients')
    plt.savefig(fr"C:\Users\natal\OneDrive - Open University of Israel\שולחן העבודה\N1_OLD_reps_iVar_results_(0.01_100_50)_COMBINED_w_TechnionResequencing_W_INDELS\syn_all_tp_10.png", bbox_inches="tight", dpi=600)
    plt.show()



if __name__ ==  "__main__":
    # Define the path to the base directory
    # base_directory = r"C:\Users\natal\OneDrive - Open University of Israel\שולחן העבודה\N1_OLD_reps_iVar_results_(0.01_100_50)_COMBINED_w_TechnionResequencing_W_INDELS"
    # # Define the non-coding regions of the SARS-CoV-2 genome
    non_coding_regions = [(1, 265 + 1), (29675, 29903 + 1)]

    # # # Define the number of possible non-synonymous and synonymous mutations
    genome_length_without_non_coding = 29264  # SARS-CoV-2 genome length without non-coding regions
    mutation_percentage_non_synonymous = 0.72
    mutation_percentage_synonymous = 0.22
    # report_df = pd.read_csv(
    #     r"Z:\home\volume1\natalie\projects\replicates_2023\suspected_problematic_by_voc\minor_alleles_disregarding_replicates_w_indels_0.01_f_0.2_mutations_summary_0.2_sample_frac.tsv",
    #     sep='\t')
    # Get the final DataFrame with divergence rates
    # final_df = get_final_df(base_directory, report_df, non_coding_regions, mutation_percentage_non_synonymous, mutation_percentage_synonymous)
    final_df = pd.read_csv(r"C:\Users\natal\OneDrive - Open University of Israel\שולחן העבודה\N1_OLD_reps_iVar_results_(0.01_100_50)_COMBINED_w_TechnionResequencing_W_INDELS\final_df_with_problematic_filter - Copy.csv")
    # Plot divergence rates for synonymous and non-synonymous mutations over time for each patient
    # plot_normalized_divergence_rates_per_patient(final_df, genome_length_without_non_coding, mutation_percentage_non_synonymous, mutation_percentage_synonymous)
    # plot_divergence_by_gene(final_df)
    # plot_divergence_by_gene_per_patient(final_df)
    # plot_non_synonymous_divergence_all_patients(final_df, genome_length_without_non_coding,
    #                                         mutation_percentage_non_synonymous)
    # plot_synonymous_divergence_all_patients(final_df, genome_length_without_non_coding, mutation_percentage_synonymous)
    plot_normalized_divergence_rates_per_patient_with_regression(final_df, genome_length_without_non_coding,mutation_percentage_non_synonymous, mutation_percentage_synonymous)
    # plot_normalized_divergence_rates_per_patient_with_regression_P5_75_mod(final_df, genome_length_without_non_coding,mutation_percentage_non_synonymous, mutation_percentage_synonymous)