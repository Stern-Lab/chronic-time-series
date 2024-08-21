############################################## 1. Import Libraries

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import colorcet as cc

############################################## 2. Utility Functions

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

def add_missing_timepoints(df, patient_mutations_dict):
    """
    Add missing time points with frequency 0 for mutations present in the patient mutations dictionary.
    :param df: DataFrame containing patient mutation data.
    :param patient_mutations_dict: Dictionary with patient names as keys and their mutations as values.
    :return: DataFrame with missing time points added.
    """
    mutated_df_list = []
    for patient, patient_data in df.groupby('patient_id'):
        if patient not in patient_mutations_dict:
            continue
        mutations = patient_mutations_dict[patient]
        timepoints = patient_data['timepoint'].unique()

        # Make a copy of timepoints and check if 0 is missing
        if '0' not in timepoints:
            timepoints = list(timepoints)  # Make a copy
            timepoints.append('0')  # Add 0 to the list

        print(timepoints)
        for mutation in mutations:
            for timepoint in timepoints:
                if mutation not in patient_data.loc[patient_data['timepoint'] == timepoint, 'mutation'].values:
                    # If mutation not present at this timepoint, add it with frequency 0
                    # Get mutation_type, protein, and POS_x corresponding to the mutation
                    mutation_row = df[(df['patient_id'] == patient) & (df['mutation'] == mutation)].iloc[0]
                    mutation_type = mutation_row['mutation_type']
                    protein = mutation_row['protein']
                    pos_x = mutation_row['POS_x']
                    mutated_df_list.append(pd.DataFrame(
                        {'patient_id': [patient], 'timepoint': [timepoint], 'mutation': [mutation],
                         'mutation_type': [mutation_type], 'protein': [protein], 'POS_x': [pos_x], 'final_freq': [0]}))
    if mutated_df_list:
        return pd.concat([df] + mutated_df_list, ignore_index=True)
    else:
        return df

def add_voc_to_df_from_metadata(df, path_to_metadata):

    metadata = pd.read_csv(path_to_metadata)
    # read only unique patients
    metadata = metadata.drop_duplicates(subset="patient")
    # add new column to the df
    df["voc"] = df["patient_id"].map(metadata.set_index("patient")["voc"])

    return df

def filter_out_suspected_problematic_by_voc(df, report_df):
    sus_mutations = {}
    filtered_out_df = pd.DataFrame()
    filtered_df = pd.DataFrame()
    for voc in df['voc'].unique():
        # put the suspected problematic mutations in a list
        sus_mutations[voc] = report_df[report_df[voc] > 0].mutation.tolist()

        # filer out the suspected problematic mutations
        voc_df = df[df['voc'] == voc]
        nonsus_voc = voc_df[~voc_df['mutation'].isin(sus_mutations[voc])]
        sus_voc = voc_df[voc_df['mutation'].isin(sus_mutations[voc])]
        filtered_out_df = pd.concat([filtered_out_df, sus_voc])
        if nonsus_voc.empty:
            continue
        filtered_df = pd.concat([filtered_df, nonsus_voc])

    return filtered_df, filtered_out_df

def filter_out_suspected_problematic_by_mutation(df, report_df):

    # Extract all suspected problematic mutations from the report_df
    sus_mutations = report_df['mutation'].tolist()

    # Filter out the suspected problematic mutations from the main DataFrame
    filtered_out_df = df[df['mutation'].isin(sus_mutations)]
    filtered_df = df[~df['mutation'].isin(sus_mutations)]

    return filtered_df, filtered_out_df

def filter_out_suspected_problematic_mutations_by_voc_and_indels(df, report_df):
    sus_mutations = {}
    filtered_out_df = pd.DataFrame()
    filtered_df = pd.DataFrame()

    for voc in df['voc'].unique():
        # put the suspected problematic mutations in a list
        sus_mutations[voc] = report_df[report_df[voc] > 0].mutation.tolist()

        # Filter out the suspected problematic mutations based on the voc
        voc_df = df[df['voc'] == voc]
        nonsus_voc = voc_df[~voc_df['mutation'].isin(sus_mutations[voc])]
        sus_voc = voc_df[voc_df['mutation'].isin(sus_mutations[voc])]
        filtered_out_df = pd.concat([filtered_out_df, sus_voc])
        if nonsus_voc.empty:
            continue
        filtered_df = pd.concat([filtered_df, nonsus_voc])

    # Additional filtering based on the presence of "+" or "-"
    sus_mutations_list = report_df['mutation'].tolist()
    mutation_with_sign = [mutation for mutation in sus_mutations_list if '+' in mutation or '-' in mutation]

    # Filter out mutations with "+" or "-" in them
    df_with_sign_filtered_out = df[df['mutation'].isin(mutation_with_sign)]
    df_with_sign_filtered_in = df[~df['mutation'].isin(mutation_with_sign)]

    # Combine the filtered DataFrames
    filtered_out_df = pd.concat([filtered_out_df, df_with_sign_filtered_out])
    filtered_df = df_with_sign_filtered_in

    return filtered_df, filtered_out_df

def update_final_freq(csv_file_path, df_to_update):

    # Read the Excel file into a DataFrame
    metadata_csv = pd.read_csv(csv_file_path)

    if df_to_update.empty:
        print("df_to_update is empty. Returning without any updates.")
        return df_to_update

    elif not df_to_update.empty:
        for index, row in metadata_csv.iterrows():
            patient = row["patient_id"]
            mutation = row["mutation"]
            timepoint = row["timepoint"]
            fixed_freq = row["fixed_freq"]


        # if patient, mutation, timepoint in df_to_update, update the final_freq value to the fixed_freq value
            mask = (df_to_update["patient_id"] == patient) & (df_to_update["mutation"] == mutation) & (df_to_update["timepoint"] == str(timepoint))
            matching_row = df_to_update[mask]
            if not matching_row.empty:
                df_to_update.loc[matching_row.index, "final_freq"] = fixed_freq

    return df_to_update

def per_patient_relevant_mutation(df):

    relevant_mutations_dict = {}
    patients = list(set(df["patient_id"]))
    for patient in patients:
        relevant_mutations_list = []
        patient_df = df[df["patient_id"] == patient].copy()
        all_patient_mutations = list(set(patient_df["mutation"]))
        for mutation in all_patient_mutations:
                df_patient_mut = patient_df[patient_df["mutation"] == mutation].copy()
                if any(df_patient_mut["final_freq"] > 0):
                    relevant_mutations_list.append(mutation)

        relevant_mutations_dict[patient] = relevant_mutations_list

    print(relevant_mutations_dict)

    return relevant_mutations_dict

# Function to exclude specified mutations and timepoints
def exclude_mutations(df, exclusion_list):
    """
    This function filters the DataFrame based on the exclusion list and exclusion positions.
    :param df: DataFrame containing patient mutation data.
    :param exclusion_list: List of dictionaries with patient_id and timepoint keys.
    :param exclusion_positions: List of mutation positions to exclude.
    :return: Filtered DataFrame.
    """
    # Convert exclusion list to a dictionary for easier access
    exclusion_dict = {entry['patient_id']: entry['timepoint'] for entry in exclusion_list}

    # Filter based on patient-specific timepoints
    for patient_id, timepoints in exclusion_dict.items():
        df = df[~((df['patient_id'] == patient_id) & (df['timepoint'].astype(int).isin(timepoints)))]

    # # Filter based on mutation positions
    # df = df[~df['POS_x'].astype(int).isin(exclusion_positions)]


    return df

############################################## 3. Data Processing Function (Main Function)

# Iterate over all patient_id directories
def get_final_df(path, non_coding_regions, csv_file_path, path_to_metadata, report_df, exclusion_list):

    # Initialize an empty DataFrame to store the final result
    # final_dataframe = pd.DataFrame()
    non_filtered_dfs_to_concat = []
    non_zeroes_dfs_to_concat = []
    filter_positions = []
    for patient_id_folder in os.listdir(base_directory):
        patient_id_path = os.path.join(base_directory, patient_id_folder)
        print(patient_id_path)
        # Check if it's a directory
        if os.path.isdir(patient_id_path):
            # Iterate over timepoint folders
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

                    # Extract patient_id and timepoint from the folder names
                    patient_id = patient_id_folder
                    timepoint = timepoint_folder

                    # Add patient_id and timepoint as new columns
                    df['patient_id'] = patient_id
                    df['timepoint'] = timepoint

                    # If mutation == T21990-TTA, change the mutation_type to 'Δ144', and the protein to surface glycoprotein
                    df.loc[df['mutation'] == 'T21990-TTA', 'mutation_type'] = 'Δ144'
                    df.loc[df['mutation'] == 'T21990-TTA', 'protein'] = 'surface glycoprotein'

                    # Label non-coding regions
                    df = label_non_coding_regions(df, non_coding_regions)

                    # Label synonymous, non-synonymous, and reversion mutations
                    df = label_syn_non_syn_reversion(df)

                    # Drop rows with mutation "T22204+GAGCCAGAA"
                    df = df[df['mutation'] != 'T22204+GAGCCAGAA']

                    # Replace -1 with NaN
                    df.loc[df['final_freq'] == -1, 'final_freq'] = np.nan

                    # Add the voc value to the DataFrame
                    df = add_voc_to_df_from_metadata(df, path_to_metadata)

                    # OPTION 1: Filter out suspected problematic mutations by voc
                    # df, filtered_out_df = filter_out_suspected_problematic_by_voc(df, report_df)
                    # filter_positions.append(filtered_out_df)

                    # OPTION 2: Filter out suspected problematic mutations by mutation
                    df, filtered_out_df = filter_out_suspected_problematic_by_mutation(df, report_df)

                    # OPTION 3: Filter out suspected problematic mutations by voc and specific mutations (indels)
                    # df, filtered_out_df = filter_out_suspected_problematic_mutations_by_voc_and_indels(df, report_df)

                    # Update the final_freq values from the manual fix Excel file
                    df_updated = update_final_freq(csv_file_path, df.copy())

                    # Check if the DataFrame is empty
                    if df_updated.empty:
                        print("df_updated is empty. Skipping this iteration.")
                        continue

                    df_with_excluded_positions = exclude_mutations(df_updated, exclusion_list)


                    # Filter rows with values > 0
                    df_get_non_zero = df_with_excluded_positions.loc[(df_with_excluded_positions['final_freq'] > 0)]

                    # Append DataFrame to the list
                    non_filtered_dfs_to_concat.append(df_with_excluded_positions)
                    non_zeroes_dfs_to_concat.append(df_get_non_zero)

                    # Concatenate all DataFrames in the list
    final_dataframe_non_filtered= pd.concat(non_filtered_dfs_to_concat, ignore_index=True)
    final_dataframe_non_zeros = pd.concat(non_zeroes_dfs_to_concat, ignore_index=True)
    # final_filtered_positions = pd.concat(filter_positions, ignore_index=True)
    # # Save filtered positions to a CSV file
    # final_filtered_positions.to_csv(
    #     fr"Z:\home\volume1\natalie\projects\replicates_2023\plots\Frequency_over_time\iVar\15-07-2024\filtered_positions_W_indels.csv", index=False)

    return final_dataframe_non_zeros, final_dataframe_non_filtered

############################################## 4. Plotting Functions

def find_duplicate_values(d):
    """
    This function takes a dictionary with patient names as keys and their mutations as values and returns a dictionary
    with mutations that are present in more than one patient and their corresponding color.
    :param d: Dictionary with patient names as keys and their mutations as values.
    :return: Dictionary with mutations that are present in more than one patient and their corresponding color.
    """

    value_counts = {}

    for patient, mutation_list in d.items():
        for mutation in mutation_list:
            if mutation in value_counts:
                value_counts[mutation].append(patient)
            else:
                value_counts[mutation] = [patient]

    duplicates_colors = cc.glasbey_cool
    duplicates_only = {}
    duplicate_values = {}
    duplicates_color_index = 0
    for mutation, patient_list in value_counts.items():
        if len(patient_list) > 1:
            duplicate_values[mutation] = duplicates_colors[duplicates_color_index]
            duplicates_only[mutation] = patient_list
            duplicates_color_index += 1
    return duplicate_values


def create_passage_graph(df, patient_mutations_dict, duplicate_mutations):
    """
    Groups DataFrame by patient and creates a plot of frequencies over time for each mutation
    present in their sequence. Because I had to decide which frequency of the two replicates I want to plot on a graph,
    I chose to plot the lower mutation frequency between the two.
    :param df: A DataFrame grouped by patient.
    :param patient_mutations_dict: Dictionary with patient names as keys and their mutations as values.
    """
    # df.loc[df['final_freq'] == -1, 'final_freq'] = 0
    df = df.dropna(subset=['final_freq'])

    SMALL_SIZE = 40
    MEDIUM_SIZE = 45
    BIGGER_SIZE = 40
    TITLE_SIZE = 45

    plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title

    grouped = df.groupby('patient_id')
    # df.sort_values(by=['timepoint'])
    # font = {'weight': 'normal', 'size': 30}
    # plt.rc('font', **font)
    plt.rcParams["figure.figsize"] = 26, 15
    plotted_values = []
    # Use 'glasbey_dark' colors from colorcet
    colors = cc.glasbey_light


    # Define a dictionary mapping full protein names to abbreviations
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

    # Create a line plot for each patient's mutation frequency over time
    for patient, patient_data in grouped:
        if patient not in patient_mutations_dict:
            continue  # Skip patients not in the dictionary

        print(patient)
        all_time_points = list(set(patient_data['timepoint'].astype(int)))
        per_mutation = patient_data.groupby('mutation')
        color_index = 0

        solid_legend_elements = []
        dashed_legend_elements = []


        for (mutation, mutation_data) in per_mutation:
            # Filter mutations based on the dictionary
            if mutation not in patient_mutations_dict[patient]:
                continue

            lw = 4
            mutation_data['timepoint'] = pd.to_numeric(mutation_data['timepoint'])
            sorted_mutation_data = mutation_data.sort_values('timepoint')  # Sort mutation data by timepoint
            abbreviated_proteins = sorted_mutation_data['protein'].map(protein_abbreviations).fillna('').infer_objects(copy=False)
            freq = np.ma.masked_invalid(sorted_mutation_data['final_freq'])
            time_points = sorted_mutation_data['timepoint'].astype(int).to_numpy()

            mutation_type = sorted_mutation_data['mutation_type'].iloc[0]  # Get mutation type for the first row
            linestyle = '-' if pd.notnull(mutation_type) else '--'  # Continuous line if mutation_type is not null, else dashed line

            # Check if the mutation was never at a frequency of 0.5 or above
            max_freq = sorted_mutation_data['final_freq'].max()
            should_have_label = True

            # Define specific styles for certain mutations
            if mutation in ['G23015C', 'G23012A', 'T23032G', 'T21990-TTA', 'A22893G']:
                linestyle = '-'
                color = 'black'
                marker = 'o'
                lw = 7

            elif mutation == 'T21990-TTA':
            # Change mutation label to 'S:'+'\u0394'+'144'
                linestyle = '-'
                color = 'black'
                marker = 'o'
                lw = 7


            elif mutation == 'T23010C':
                linestyle = '-'
                color = 'black'
                marker = 'o'
                lw = 7

            else:
                if mutation in duplicate_mutations:
                    color = duplicate_mutations[mutation]
                    marker = 'o'

            # if mutation == 'G23015C' or mutation == 'T23010C' or mutation == 'G23012A' or mutation == 'T23032G' or mutation == 'T21990-TTA':
            #     color = 'black'
            #     marker = 'o'
            #     # Make the line width thicker
            #     lw = 8



            # elif mutation == 'C15173A':  # Check if mutation is T23010C (S:V483A)
            #     color = 'black'  # Assign a specific color for T23010C mutation
            #     marker = 'o'
                else:
                    if max_freq < 0.5:
                        color = colors[color_index] + '20'  # Fade color if max frequency < 0.5
                        should_have_label = True  # Change to True if you want all the labels of the fade mutations to be present in the legend
                    else:
                        color = colors[color_index]
                marker = 'o'
                color_index += 1

            label = np.where(pd.notnull(sorted_mutation_data['mutation_type']),
                             abbreviated_proteins + ':' + sorted_mutation_data['mutation_type'],
                             mutation)[0]

            plotted_values.extend(
                list(zip(time_points, freq, [label] * len(time_points), [patient] * len(time_points), [mutation] * len(time_points))))
            plt.plot(time_points, freq, label=label if should_have_label else None, color=color, linestyle=linestyle, marker=marker, markersize=15,
                     lw=lw)

            if should_have_label:
                legend_element = plt.Line2D([0], [0], color=color, linestyle=linestyle, marker=marker, markersize=15,
                                            lw=lw, label=label)
                if linestyle == '-':
                    solid_legend_elements.append(legend_element)
                else:
                    dashed_legend_elements.append(legend_element)

        # Create custom legend handles for the linestyle types
        custom_solid_handle = plt.Line2D([0], [0], color='black', linestyle='-', lw=4, label='Non-Synonymous Mutations')
        custom_dashed_handle = plt.Line2D([0], [0], color='black', linestyle='--', lw=4, label='Synonymous Mutations')

        # Combine the custom handles with the actual mutation legend elements
        combined_legend_elements = solid_legend_elements + dashed_legend_elements + [custom_solid_handle] + [
            custom_dashed_handle]

        plt.xticks(all_time_points, all_time_points)
        plt.yticks((np.arange(0, 1.2, step=0.2)), [0, 0.2, 0.4, 0.6, 0.8, 1])
        plt.xlabel('')
        plt.ylabel('')
        plt.title(f'{patient}', fontsize=TITLE_SIZE, weight='bold')
        # plt.legend(handles=combined_legend_elements, bbox_to_anchor=(1.04, 1), loc="upper left", ncol=3))  # Change to ncol=3 (or more) if should_have_label = True
        plt.grid(True)
        plt.savefig(
            fr"Z:\home\volume1\natalie\projects\replicates_2023\plots\Frequency_over_time\iVar\13-08-2024\DPI_600_LONG_{patient}_frequency_over_time_TechnionResequencing_W_INDELS.png", bbox_inches="tight", dpi=600)
        plt.show()

        # Save the legend as a separate PNG file
        fig_legend, ax_legend = plt.subplots(figsize=(10, 5))
        ax_legend.legend(handles=combined_legend_elements, loc='center', ncol=3)
        ax_legend.axis('off')
        fig_legend.savefig(fr"Z:\home\volume1\natalie\projects\replicates_2023\plots\Frequency_over_time\iVar\13-08-2024\DPI_600_LONG_legend_{patient}_frequency_over_time_TechnionResequencing_W_INDELS.png", bbox_inches='tight', dpi=600)
        plt.close(fig_legend)

        # Convert the list to a DataFrame and save to CSV
        columns = ['Time', 'Frequency', 'Mutation', 'Patient', 'Mutation_Nt']
        plotted_values_df = pd.DataFrame(plotted_values, columns=columns)
        plotted_values_df.to_csv(fr"Z:\home\volume1\natalie\projects\replicates_2023\plots\Frequency_over_time\iVar\13-08-2024\SUMMARY_of_mutations_plotted_W.csv", index=False)


############################################## 5. Report Function
def generate_report(df):
    """
    Generates a summary report of the data processing steps.

    Parameters:
    df (pd.DataFrame): DataFrame containing the final processed data.

    Returns:
    None: Outputs a report.
    """

############################################## 6. Main Function
def main():
    # Specify the directory where your files are located
    base_directory = r"Z:\home\volume1\natalie\repositories\BN-SCRIPTS\Filter_Usecase\results\N1_OLD_reps_iVar_results_(0.01_100_50)_COMBINED_w_TechnionResequencing_W_INDELS"
    # Manual fix of mutation frequencies
    csv_file_path = r"Z:\home\volume1\natalie\projects\replicates_2023\Manual_correction_of_mutation_drops\FINAL_SUMMARY_hard_mutation_drops_w_tech_seq_p5.csv"
    path_to_metadata = r"Z:\nobackup\volume1\natalie\ichilov_chronic_omicron\libraries_analysis\replicates_2023\all_patients_global_content_initials_V5.csv"
    report_df = pd.read_csv(
        r"Z:\home\volume1\natalie\projects\replicates_2023\suspected_problematic_by_voc\minor_alleles_disregarding_replicates_w_indels_0.01_f_0.2_mutations_summary_0.2_sample_frac.tsv",
        sep='\t')

    # Define the exclusion list
    exclusion_list = [
        {'patient_id': 'P5', 'timepoint': [27, 44]},
        {'patient_id': 'P4', 'timepoint': [13, 28]},
    ]
    # Define the non-coding regions of the SARS-CoV-2 genome
    non_coding_regions = [(1, 265 + 1), (29675, 29903 + 1)]

    non_zero, non_filtered = get_final_df(base_directory, non_coding_regions, csv_file_path, path_to_metadata, report_df, exclusion_list)

    relevant_mut_dict = per_patient_relevant_mutation(non_zero)

    # Adding missing time points with frequency 0 for relevant mutations
    df = add_missing_timepoints(non_filtered, relevant_mut_dict)
    # Example usage:
    duplicates = find_duplicate_values(relevant_mut_dict)



    # Define the range of mutation positions to exclude
    # exclusion_positions = list(range(15168, 15174))
    # df = exclude_mutations(df, exclusion_list)


    # Save the filtered DataFrame to a CSV file - for Shir
    # output_file_path = fr"Z:\home\volume1\natalie\projects\replicates_2023\N2+P3\filtered_dataframe.csv"
    # df.to_csv(output_file_path, index=False)
    # print(f"Filtered DataFrame saved as {output_file_path}")
    create_passage_graph(df, relevant_mut_dict, duplicates)

if __name__ == "__main__":
    main()
