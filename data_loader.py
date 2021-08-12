
import pandas as pd
import numpy as np
import os
import pyteomics.mzml
import spectrum_utils.spectrum as sus



#how we're going to format the peptide modifictations:
# +15.995
# +57.021
def format_oxidation(row, column, to_replace):
    peptide = row[column]
#     print(to_replace)
    replace_with = "+15.995"
    if pd.isna(peptide):
        new_pep = peptide
    else:
        if to_replace in peptide:
            new_pep = peptide.replace(to_replace, replace_with)
        else:
            new_pep = peptide
    return new_pep
def format_carbamidomethyl(row, column, to_replace):
    peptide = row[column]
#     print(to_replace)
    replace_with = ""
    if pd.isna(peptide):
        new_pep = peptide
    else:
        if to_replace in peptide:
            new_pep = peptide.replace(to_replace, replace_with)
        else:
            new_pep = peptide
    return new_pep


def download_parsed_psm():

    list_of_file_names = [ 'bulk_rep1',"bulk_rep2","bulk_rep3",
                            '2ng_rep1','2ng_rep2','2ng_rep3','2ng_rep4','2ng_rep5','2ng_rep6',
                            '0.2ng_rep1','0.2ng_rep2','0.2ng_rep3','0.2ng_rep4','0.2ng_rep5','0.2ng_rep6']


    mm_files = {}
    #bulk
    mm_files["bulk_rep1"] = "data/Project_PXD011163/Metamorpheus_output/OR11_20160122_PG_HeLa_CVB3_CT_A-calib_PSMs.psmtsv"
    mm_files["bulk_rep2"] = "data/Project_PXD011163/Metamorpheus_output/OR11_20160122_PG_HeLa_CVB3_CT_B-calib_PSMs.psmtsv"
    mm_files["bulk_rep3"] = "data/Project_PXD011163/Metamorpheus_output/OR11_20160122_PG_HeLa_CVB3_CT_C-calib_PSMs.psmtsv"

    #2ng
    mm_files["2ng_rep1"] = "data/2ng_Files_March_19_2021/Metamorpheus_output/Ex_Auto_J3_30umTB_2ngQC_60m_1-calib_PSMs.psmtsv"
    mm_files["2ng_rep2"] = "data/2ng_Files_March_19_2021/Metamorpheus_output/Ex_Auto_J3_30umTB_2ngQC_60m_2-calib_PSMs.psmtsv"
    mm_files["2ng_rep3"] = "data/2ng_Files_March_19_2021/Metamorpheus_output/Ex_Auto_K13_30umTA_2ngQC_60m_1-calib_PSMs.psmtsv"
    mm_files["2ng_rep4"] = "data/2ng_Files_March_19_2021/Metamorpheus_output/Ex_Auto_K13_30umTA_2ngQC_60m_2-calib_PSMs.psmtsv"
    mm_files["2ng_rep5"] = "data/2ng_Files_March_19_2021/Metamorpheus_output/Ex_Auto_W17_30umTB_2ngQC_60m_1-calib_PSMs.psmtsv"
    mm_files["2ng_rep6"] = "data/2ng_Files_March_19_2021/Metamorpheus_output/Ex_Auto_W17_30umTB_2ngQC_60m_2-calib_PSMs.psmtsv"

    #.2ng
    mm_files["0.2ng_rep1"] = "data/02ng_Files_March_19_2021/Meatmorpheus_output/Ex_Auto_J3_30umTB_02ngQC_60m_1-calib_PSMs.psmtsv"
    mm_files["0.2ng_rep2"] = "data/02ng_Files_March_19_2021/Meatmorpheus_output/Ex_Auto_J3_30umTB_02ngQC_60m_2-calib_PSMs.psmtsv"
    mm_files["0.2ng_rep3"] = "data/02ng_Files_March_19_2021/Meatmorpheus_output/Ex_Auto_K13_30umTA_02ngQC_60m_1-calib_PSMs.psmtsv"
    mm_files["0.2ng_rep4"] = "data/02ng_Files_March_19_2021/Meatmorpheus_output/Ex_Auto_K13_30umTA_02ngQC_60m_2-calib_PSMs.psmtsv"
    mm_files["0.2ng_rep5"] = "data/02ng_Files_March_19_2021/Meatmorpheus_output/Ex_Auto_W17_30umTA_02ngQC_60m_3-calib_PSMs.psmtsv"
    mm_files["0.2ng_rep6"] = "data/02ng_Files_March_19_2021/Meatmorpheus_output/Ex_Auto_W17_30umTA_02ngQC_60m_4-calib_PSMs.psmtsv"
    #
    #
    #
    for file_name in list_of_file_names:
        #get path to file
        path_to_data_loader = os.path.abspath(os.path.dirname(__file__)) # This gets the absolute path to the location of the data_loader.py file
        complete_path_to_data = os.path.join(path_to_data_loader, mm_files.get(file_name)) # We then append the relative path to the data files

        data = pd.read_csv(complete_path_to_data, sep = '\t', dtype={'Base Sequence': str, 'Missed Cleavages': str, 'Peptide Monoisotopic Mass':str, 'Mass Diff (Da)':str,'Mass Diff (ppm)':str,'	Protein Accession':str,'Peptide Description':str, 'Notch':str, 'Num Variable Mods':str, 'Decoy': str})

        #make a new col that includes modifide peptides
        data['temp_peptide'] = data.apply(lambda row: format_oxidation(row, "Full Sequence", "[Common Variable:Oxidation on M]"), axis=1)
        data["temp2"] = data.apply(lambda row: format_carbamidomethyl(row, "temp_peptide", "[Common Fixed:Carbamidomethyl on C]"), axis=1)


        #uniform naming
        data_new = data.rename({"Decoy": "decoy", "Scan Number": "scan", "temp2": "peptide", 'QValue': 'probability'}, axis=1)

        #remove duplicate scans
        data_new = data_new.sort_values("probability")
        data_new = data_new.drop_duplicates(subset=["scan"], keep="first")

        write_file_path = "data/parsed_psm/" + file_name + ".csv"
        data_new.to_csv(write_file_path, sep="\t")
def load_psm(file_name):
    #return the joined table (we need the start_scan time)
    mm_files = {}
    mm_files["bulk_rep1"] = "data/parsed_psm/bulk_rep1.csv"
    mm_files["bulk_rep2"] = "data/parsed_psm/bulk_rep2.csv"
    mm_files["bulk_rep3"] = "data/parsed_psm/bulk_rep3.csv"

    mm_files["2ng_rep1"] = "data/parsed_psm/2ng_rep1.csv"
    mm_files["2ng_rep2"] = "data/parsed_psm/2ng_rep2.csv"
    mm_files["2ng_rep3"] = "data/parsed_psm/2ng_rep3.csv"
    mm_files["2ng_rep4"] = "data/parsed_psm/2ng_rep4.csv"
    mm_files["2ng_rep5"] = "data/parsed_psm/2ng_rep5.csv"
    mm_files["2ng_rep6"] = "data/parsed_psm/2ng_rep6.csv"

    mm_files["0.2ng_rep1"] = "data/parsed_psm/0.2ng_rep1.csv"
    mm_files["0.2ng_rep2"] = "data/parsed_psm/0.2ng_rep2.csv"
    mm_files["0.2ng_rep3"] = "data/parsed_psm/0.2ng_rep3.csv"
    mm_files["0.2ng_rep4"] = "data/parsed_psm/0.2ng_rep4.csv"
    mm_files["0.2ng_rep5"] = "data/parsed_psm/0.2ng_rep5.csv"
    mm_files["0.2ng_rep6"] = "data/parsed_psm/0.2ng_rep6.csv"

    mm_files["Bcell_rep1"] = "data/parsed_psm/Bcell_rep1.csv"
    mm_files["Bcell_rep2"] = "data/parsed_psm/Bcell_rep2.csv"
    mm_files["Tcell_rep1"] = "data/parsed_psm/Tcell_rep1.csv"
    mm_files["Tcell_rep2"] = "data/parsed_psm/Tcell_rep2.csv"

    path_to_data_loader = os.path.abspath(os.path.dirname(__file__)) # This gets the absolute path to the location of the data_loader.py file
    complete_path_to_data = os.path.join(path_to_data_loader, mm_files.get(file_name))


    df = pd.read_csv(complete_path_to_data, sep='\t', index_col=0, dtype={'Base Sequence': str, 'Missed Cleavages': str, 'Peptide Monoisotopic Mass':str, 'Mass Diff (Da)':str,'Mass Diff (ppm)':str,'	Protein Accession':str,'Peptide Description':str, 'Notch':str, 'Num Variable Mods':str, 'Peptide Description':str, 'Decoy': str})


    return(df)

def download_mzml():
    list_of_file_names = [ "bulk_rep1", 'bulk_rep2', 'bulk_rep3',
                            '2ng_rep1','2ng_rep2','2ng_rep3','2ng_rep4','2ng_rep5','2ng_rep6',
                            '0.2ng_rep1','0.2ng_rep2','0.2ng_rep3','0.2ng_rep4','0.2ng_rep5','0.2ng_rep6']


    mzml = {}
    mzml["bulk_rep1"] = "data/Project_PXD011163/mzMLs/OR11_20160122_PG_HeLa_CVB3_CT_A.mzML"
    mzml["bulk_rep2"] = "data/Project_PXD011163/mzMLs/OR11_20160122_PG_HeLa_CVB3_CT_B.mzML"
    mzml["bulk_rep3"] = "data/Project_PXD011163/mzMLs/OR11_20160122_PG_HeLa_CVB3_CT_C.mzML"

    mzml["2ng_rep1"] = "data/2ng_Files_March_19_2021/mzMLs/Ex_Auto_J3_30umTB_2ngQC_60m_1.mzML"
    mzml["2ng_rep2"] = "data/2ng_Files_March_19_2021/mzMLs/Ex_Auto_J3_30umTB_2ngQC_60m_2.mzML"
    mzml["2ng_rep3"] = "data/2ng_Files_March_19_2021/mzMLs/Ex_Auto_K13_30umTA_2ngQC_60m_1.mzML"
    mzml["2ng_rep4"] = "data/2ng_Files_March_19_2021/mzMLs/Ex_Auto_K13_30umTA_2ngQC_60m_2.mzML"
    mzml["2ng_rep5"] = "data/2ng_Files_March_19_2021/mzMLs/Ex_Auto_W17_30umTB_2ngQC_60m_1.mzML"
    mzml["2ng_rep6"] = "data/2ng_Files_March_19_2021/mzMLs/Ex_Auto_W17_30umTB_2ngQC_60m_2.mzML"
    #
    mzml["0.2ng_rep1"] = "data/02ng_Files_March_19_2021/mzMLs/Ex_Auto_J3_30umTB_02ngQC_60m_1.mzML"
    mzml["0.2ng_rep2"] = "data/02ng_Files_March_19_2021/mzMLs/Ex_Auto_J3_30umTB_02ngQC_60m_2.mzML"
    mzml["0.2ng_rep3"] = "data/02ng_Files_March_19_2021/mzMLs/Ex_Auto_K13_30umTA_02ngQC_60m_1.mzML"
    mzml["0.2ng_rep4"] = "data/02ng_Files_March_19_2021/mzMLs/Ex_Auto_K13_30umTA_02ngQC_60m_2.mzML"
    mzml["0.2ng_rep5"] = "data/02ng_Files_March_19_2021/mzMLs/Ex_Auto_W17_30umTA_02ngQC_60m_3.mzML"
    mzml["0.2ng_rep6"] = "data/02ng_Files_March_19_2021/mzMLs/Ex_Auto_W17_30umTA_02ngQC_60m_4.mzML"


    for file_name in list_of_file_names:
        #get path to file
        path_to_data_loader = os.path.abspath(os.path.dirname(__file__)) # This gets the absolute path to the location of the data_loader.py file
        complete_path_to_data = os.path.join(path_to_data_loader, mzml.get(file_name))

        columns = ["scan",'ms_level',"total_ion_curr","ion_time",'mz_array', 'scan_start_time', 'intensity_array', 'retention_time','precursor_charge']
        lst = []
        with pyteomics.mzml.MzML(complete_path_to_data) as f_in:
            for spectrum_dict in f_in:
                scan_num = spectrum_dict["id"][spectrum_dict["id"].find('scan=') + 5:]
                ms_level = spectrum_dict["ms level"]
                total_ion_curr = spectrum_dict["total ion current"]
                ion_time = spectrum_dict["scanList"]['scan'][0]["ion injection time"]
                mz_array = list(spectrum_dict['m/z array'])
                scan_start_time = spectrum_dict['scanList']['scan'][0]['scan start time']

                #added for graphing purposes
                if ms_level == 2:
                    intensity_array = list(spectrum_dict['intensity array'])
                    retention_time = (spectrum_dict['scanList']['scan'][0].get('scan start time', -1))
                    precursor = spectrum_dict['precursorList']['precursor'][0]
                    precursor_ion = precursor['selectedIonList']['selectedIon'][0]
                    precursor_mz = precursor_ion['selected ion m/z']
                    if 'charge state' in precursor_ion:
                        precursor_charge = int(precursor_ion['charge state'])
                    elif 'possible charge state' in precursor_ion:
                        precursor_charge = int(precursor_ion['possible charge state'])
                    else:
                        precursor_charge = 'NAN'
                        # raise ValueError('Unknown precursor charge')


                    lst.append([scan_num, ms_level, total_ion_curr, ion_time, mz_array, scan_start_time, intensity_array, retention_time, precursor_charge])

                else:
                    lst.append([scan_num, ms_level, total_ion_curr, ion_time, mz_array, scan_start_time, 'NAN', 'NAN','NAN'])


        mzml_df = pd.DataFrame(lst, columns=columns)

        mzml_df[["scan"]] = mzml_df[["scan"]].apply(pd.to_numeric)
        mzml_df["scan_start_time"] = mzml_df["scan_start_time"].astype(str)

        #split the time column so you get one that's just minutes
        mzml_df["temp_minute"] = mzml_df["scan_start_time"].str.split("\.")
        mzml_df.loc[:, 'minute'] = mzml_df['temp_minute'].map(lambda x: x[0])

        mzml_df[["minute"]] = mzml_df[["minute"]].apply(pd.to_numeric)
        mzml_df= mzml_df.sort_values("minute")

        write_file_path = "data/parsed_mzml/" + file_name + ".csv"
        mzml_df.to_csv(write_file_path, sep="\t")
def load_mzml(file_name):
    mzml = {}
    mzml["bulk_rep1"] = "data/parsed_mzml/bulk_rep1.csv"
    mzml["bulk_rep2"] = "data/parsed_mzml/bulk_rep2.csv"
    mzml["bulk_rep3"] = "data/parsed_mzml/bulk_rep3.csv"

    mzml["2ng_rep1"] = "data/parsed_mzml/2ng_rep1.csv"
    mzml["2ng_rep2"] = "data/parsed_mzml/2ng_rep2.csv"
    mzml["2ng_rep3"] = "data/parsed_mzml/2ng_rep3.csv"
    mzml["2ng_rep4"] = "data/parsed_mzml/2ng_rep4.csv"
    mzml["2ng_rep5"] = "data/parsed_mzml/2ng_rep5.csv"
    mzml["2ng_rep6"] = "data/parsed_mzml/2ng_rep6.csv"

    mzml["0.2ng_rep1"] = "data/parsed_mzml/0.2ng_rep1.csv"
    mzml["0.2ng_rep2"] = "data/parsed_mzml/0.2ng_rep2.csv"
    mzml["0.2ng_rep3"] = "data/parsed_mzml/0.2ng_rep3.csv"
    mzml["0.2ng_rep4"] = "data/parsed_mzml/0.2ng_rep4.csv"
    mzml["0.2ng_rep5"] = "data/parsed_mzml/0.2ng_rep5.csv"
    mzml["0.2ng_rep6"] = "data/parsed_mzml/0.2ng_rep6.csv"

    path_to_data_loader = os.path.abspath(os.path.dirname(__file__)) # This gets the absolute path to the location of the data_loader.py file
    complete_path_to_data = os.path.join(path_to_data_loader, mzml.get(file_name))


    df = pd.read_csv(complete_path_to_data, sep='\t')

    return(df)

def download_joined_psm_mzml():
    list_of_file_names = [ "bulk_rep1", 'bulk_rep2', 'bulk_rep3',
                            '2ng_rep1','2ng_rep2','2ng_rep3','2ng_rep4','2ng_rep5','2ng_rep6',
                            '0.2ng_rep1','0.2ng_rep2','0.2ng_rep3','0.2ng_rep4','0.2ng_rep5','0.2ng_rep6']

    for file_name in list_of_file_names:
        df_mm = load_psm(file_name)
        df_mz = load_mzml(file_name)

        #join the two together
        joined = pd.merge(df_mm, df_mz, on="scan").sort_values("scan")

        write_file_path = "data/joined_psm_and_mzml/" + file_name + ".csv"
        joined.to_csv(write_file_path, sep="\t")
def load_joined_psm_mzml(file_name):
    joined_files = {}
    joined_files["bulk_rep1"] = "data/joined_psm_and_mzml/bulk_rep1.csv"
    joined_files["bulk_rep2"] = "data/joined_psm_and_mzml/bulk_rep2.csv"
    joined_files["bulk_rep3"] = "data/joined_psm_and_mzml/bulk_rep3.csv"

    joined_files["2ng_rep1"] = "data/joined_psm_and_mzml/2ng_rep1.csv"
    joined_files["2ng_rep2"] = "data/joined_psm_and_mzml/2ng_rep2.csv"
    joined_files["2ng_rep3"] = "data/joined_psm_and_mzml/2ng_rep3.csv"
    joined_files["2ng_rep4"] = "data/joined_psm_and_mzml/2ng_rep4.csv"
    joined_files["2ng_rep5"] = "data/joined_psm_and_mzml/2ng_rep5.csv"
    joined_files["2ng_rep6"] = "data/joined_psm_and_mzml/2ng_rep6.csv"

    joined_files["0.2ng_rep1"] = "data/joined_psm_and_mzml/0.2ng_rep1.csv"
    joined_files["0.2ng_rep2"] = "data/joined_psm_and_mzml/0.2ng_rep2.csv"
    joined_files["0.2ng_rep3"] = "data/joined_psm_and_mzml/0.2ng_rep3.csv"
    joined_files["0.2ng_rep4"] = "data/joined_psm_and_mzml/0.2ng_rep4.csv"
    joined_files["0.2ng_rep5"] = "data/joined_psm_and_mzml/0.2ng_rep5.csv"
    joined_files["0.2ng_rep6"] = "data/joined_psm_and_mzml/0.2ng_rep6.csv"

    path_to_data_loader = os.path.abspath(os.path.dirname(__file__)) # This gets the absolute path to the location of the data_loader.py file
    complete_path_to_data = os.path.join(path_to_data_loader, joined_files.get(file_name))

    df = pd.read_csv(complete_path_to_data, sep='\t', index_col=0, dtype={'Base Sequence': str, 'Missed Cleavages': str, 'Peptide Monoisotopic Mass':str, 'Mass Diff (Da)':str,'Mass Diff (ppm)':str,'	Protein Accession':str,'Peptide Description':str, 'Notch':str, 'Num Variable Mods':str, 'Peptide Description':str, 'Decoy': str})

    return(df)
