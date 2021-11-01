
import pandas as pd
import numpy as np
import os
import pyteomics.mzml
import spectrum_utils.spectrum as sus
from pathlib import Path



def download_parsed_psm():

    list_of_file_names = [ 'bulk_rep1',"bulk_rep2","bulk_rep3",
                            '2ng_rep1','2ng_rep2','2ng_rep3','2ng_rep4','2ng_rep5','2ng_rep6',
                            '0.2ng_rep1','0.2ng_rep2','0.2ng_rep3','0.2ng_rep4','0.2ng_rep5','0.2ng_rep6',
                            'sc_rep1','sc_rep2','sc_rep3','sc_rep4','sc_rep5']

    # list_of_file_names = ['2ng_rep1']
    mm_files = {}
    #bulk
    mm_files["bulk_rep1"] = "data/MetaMorpheus_output/Project_PXD011163/OR11_20160122_PG_HeLa_CVB3_CT_A-calib_PSMs.psmtsv.gz"
    mm_files["bulk_rep2"] = "data/MetaMorpheus_output/Project_PXD011163/OR11_20160122_PG_HeLa_CVB3_CT_B-calib_PSMs.psmtsv.gz"
    mm_files["bulk_rep3"] = "data/MetaMorpheus_output/Project_PXD011163/OR11_20160122_PG_HeLa_CVB3_CT_C-calib_PSMs.psmtsv.gz"

    #2ng
    mm_files["2ng_rep1"] = "data/MetaMorpheus_output/2ng/Ex_Auto_J3_30umTB_2ngQC_60m_1-calib_PSMs.psmtsv.gz"
    mm_files["2ng_rep2"] = "data/MetaMorpheus_output/2ng/Ex_Auto_J3_30umTB_2ngQC_60m_2-calib_PSMs.psmtsv.gz"
    mm_files["2ng_rep3"] = "data/MetaMorpheus_output/2ng/Ex_Auto_K13_30umTA_2ngQC_60m_1-calib_PSMs.psmtsv.gz"
    mm_files["2ng_rep4"] = "data/MetaMorpheus_output/2ng/Ex_Auto_K13_30umTA_2ngQC_60m_2-calib_PSMs.psmtsv.gz"
    mm_files["2ng_rep5"] = "data/MetaMorpheus_output/2ng/Ex_Auto_W17_30umTB_2ngQC_60m_1-calib_PSMs.psmtsv.gz"
    mm_files["2ng_rep6"] = "data/MetaMorpheus_output/2ng/Ex_Auto_W17_30umTB_2ngQC_60m_2-calib_PSMs.psmtsv.gz"

    #.2ng
    mm_files["0.2ng_rep1"] = "data/MetaMorpheus_output/02ng/Ex_Auto_J3_30umTB_02ngQC_60m_1-calib_PSMs.psmtsv.gz"
    mm_files["0.2ng_rep2"] = "data/MetaMorpheus_output/02ng/Ex_Auto_J3_30umTB_02ngQC_60m_2-calib_PSMs.psmtsv.gz"
    mm_files["0.2ng_rep3"] = "data/MetaMorpheus_output/02ng/Ex_Auto_K13_30umTA_02ngQC_60m_1-calib_PSMs.psmtsv.gz"
    mm_files["0.2ng_rep4"] = "data/MetaMorpheus_output/02ng/Ex_Auto_K13_30umTA_02ngQC_60m_2-calib_PSMs.psmtsv.gz"
    mm_files["0.2ng_rep5"] = "data/MetaMorpheus_output/02ng/Ex_Auto_W17_30umTA_02ngQC_60m_3-calib_PSMs.psmtsv.gz"
    mm_files["0.2ng_rep6"] = "data/MetaMorpheus_output/02ng/Ex_Auto_W17_30umTA_02ngQC_60m_4-calib_PSMs.psmtsv.gz"
    #
    #True SC
    mm_files["sc_rep1"] = "data/MetaMorpheus_output/true_SC/D19_15um30cm_SC1-calib_PSMs.psmtsv.gz"
    mm_files["sc_rep2"] = "data/MetaMorpheus_output/true_SC/D19_15um30cm_SC2-calib_PSMs.psmtsv.gz"
    mm_files["sc_rep3"] = "data/MetaMorpheus_output/true_SC/D19_15um30cm_SC3-calib_PSMs.psmtsv.gz"
    mm_files["sc_rep4"] = "data/MetaMorpheus_output/true_SC/D19_15um30cm_SC4-calib_PSMs.psmtsv.gz"
    mm_files["sc_rep5"] = "data/MetaMorpheus_output/true_SC/D19_15um30cm_SC5-calib_PSMs.psmtsv.gz"

    mzml_files = {}
    mzml_files["bulk_rep1"] = "data/mzMLs/OR11_20160122_PG_HeLa_CVB3_CT_A.mzML"
    mzml_files["bulk_rep2"] = "data/mzMLs/OR11_20160122_PG_HeLa_CVB3_CT_B.mzML"
    mzml_files["bulk_rep3"] = "data/mzMLs/OR11_20160122_PG_HeLa_CVB3_CT_C.mzML"

    mzml_files["2ng_rep1"] = "data/mzMLs/Ex_Auto_J3_30umTB_2ngQC_60m_1.mzML"
    mzml_files["2ng_rep2"] = "data/mzMLs/Ex_Auto_J3_30umTB_2ngQC_60m_2.mzML"
    mzml_files["2ng_rep3"] = "data/mzMLs/Ex_Auto_K13_30umTA_2ngQC_60m_1.mzML"
    mzml_files["2ng_rep4"] = "data/mzMLs/Ex_Auto_K13_30umTA_2ngQC_60m_2.mzML"
    mzml_files["2ng_rep5"] = "data/mzMLs/Ex_Auto_W17_30umTB_2ngQC_60m_1.mzML"
    mzml_files["2ng_rep6"] = "data/mzMLs/Ex_Auto_W17_30umTB_2ngQC_60m_2.mzML"
    #
    mzml_files["0.2ng_rep1"] = "data/mzMLs/Ex_Auto_J3_30umTB_02ngQC_60m_1.mzML"
    mzml_files["0.2ng_rep2"] = "data/mzMLs/Ex_Auto_J3_30umTB_02ngQC_60m_2.mzML"
    mzml_files["0.2ng_rep3"] = "data/mzMLs/Ex_Auto_K13_30umTA_02ngQC_60m_1.mzML"
    mzml_files["0.2ng_rep4"] = "data/mzMLs/Ex_Auto_K13_30umTA_02ngQC_60m_2.mzML"
    mzml_files["0.2ng_rep5"] = "data/mzMLs/Ex_Auto_W17_30umTA_02ngQC_60m_3.mzML"
    mzml_files["0.2ng_rep6"] = "data/mzMLs/Ex_Auto_W17_30umTA_02ngQC_60m_4.mzML"

    mzml_files["sc_rep1"] = "data/mzMLs/D19_15um30cm_SC1.mzML"
    mzml_files["sc_rep2"] = "data/mzMLs/D19_15um30cm_SC2.mzML"
    mzml_files["sc_rep3"] = "data/mzMLs/D19_15um30cm_SC3.mzML"
    mzml_files["sc_rep4"] = "data/mzMLs/D19_15um30cm_SC4.mzML"
    mzml_files["sc_rep5"] = "data/mzMLs/D19_15um30cm_SC5.mzML"




    Path('data/parsed_psm/').mkdir(parents=True, exist_ok=True)#make the folder that we'll store all parsed psms in

    for file_name in list_of_file_names:
        path_to_data_loader = os.path.abspath(os.path.dirname(__file__))
        complete_path_to_psm = os.path.join(path_to_data_loader, mm_files.get(file_name)) # We then append the relative path to the data file
        complete_path_to_mzml = os.path.join(path_to_data_loader, mzml_files.get(file_name))

        psm = pd.read_csv(complete_path_to_psm, sep = '\t', dtype={'Base Sequence': str, 'Missed Cleavages': str, 'Peptide Monoisotopic Mass':str, 'Mass Diff (Da)':str,'Mass Diff (ppm)':str,'	Protein Accession':str,'Peptide Description':str, 'Notch':str, 'Num Variable Mods':str, 'Decoy': str})

        #drop columns that we don't use
        unused = ['Score', 'Delta Score', 'Notch','Base Sequence','Essential Sequence',
        'Mods','Mods Chemical Formulas', 'Mods Combined Chemical Formula',
         'Num Variable Mods', 'Missed Cleavages','Organism Name','Identified Sequence Variations',
         'Splice Sites', 'Contaminant','Peptide Description', 'Start and End Residues In Protein',
          'Previous Amino Acid', 'Next Amino Acid', 'Theoreticals Searched', 'Decoy/Contaminant/Target',
          'Normalized Spectral Angle', 'Localized Scores',
       'Improvement Possible', 'Cumulative Target', 'Cumulative Decoy',
       'Cumulative Target Notch', 'Cumulative Decoy Notch',
       'QValue Notch', 'PEP', 'PEP_QValue']
        psm = psm.drop(columns = unused)

        psm = psm.rename({"Scan Number": "scan", "Full Sequence": "peptide"}, axis=1)

        #remove duplicate scans
        psm = psm.sort_values("QValue")
        psm = psm.drop_duplicates(subset=["scan"], keep="first")

        psm[["scan"]] = psm[["scan"]].apply(pd.to_numeric)
        psm["Scan Retention Time"] = psm["Scan Retention Time"].astype(str)

        #split the time column so you get one that's just minutes
        psm["temp_minute"] = psm["Scan Retention Time"].str.split("\.")
        psm.loc[:, 'minute'] = psm['temp_minute'].map(lambda x: x[0])

        psm[["minute"]] = psm[["minute"]].apply(pd.to_numeric)
        psm = psm.drop(columns='temp_minute')

        #append info from mzML
        # import pdb;pdb.set_trace()
        mzml = pyteomics.mzml.MzML(complete_path_to_mzml)
        psm['mzml_info'] = psm.apply(lambda row: _mzml_helper(row, mzml),axis=1)
        psm[['mz_array', 'intensity_array', 'precursor_intenisty']] = pd.DataFrame(psm['mzml_info'].tolist(), index= psm.index)
        psm=psm.drop(columns='mzml_info')
        write_file_path = "data/parsed_psm/" + file_name + ".csv"
        psm.to_csv(write_file_path, sep="\t")

def load_joined_psm_mzml(file_name):
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

    mm_files["sc_rep1"] = "data/parsed_psm/sc_rep1.csv"
    mm_files["sc_rep2"] = "data/parsed_psm/sc_rep2.csv"
    mm_files["sc_rep3"] = "data/parsed_psm/sc_rep3.csv"
    mm_files["sc_rep4"] = "data/parsed_psm/sc_rep4.csv"
    mm_files["sc_rep5"] = "data/parsed_psm/sc_rep5.csv"

    path_to_data_loader = os.path.abspath(os.path.dirname(__file__)) # This gets the absolute path to the location of the data_loader.py file
    complete_path_to_data = os.path.join(path_to_data_loader, mm_files.get(file_name))


    df = pd.read_csv(complete_path_to_data, sep='\t', index_col=0, dtype={'Base Sequence': str, 'Missed Cleavages': str, 'Peptide Monoisotopic Mass':str, 'Mass Diff (Da)':str,'Mass Diff (ppm)':str,'	Protein Accession':str,'Peptide Description':str, 'Notch':str, 'Num Variable Mods':str, 'Peptide Description':str, 'Decoy': str})


    return(df)

def _mzml_helper(row, mzml):
    scan = str(row['scan'])
    my_id = 'controllerType=0 controllerNumber=1 scan='+ str(scan)
    spectrum_dict = mzml.get_by_id(my_id)

    spectrum_id = spectrum_dict['id']
    retention_time = (spectrum_dict['scanList']['scan'][0].get('scan start time', -1))
    mz_array = list(spectrum_dict['m/z array'])
    intensity_array = list(spectrum_dict['intensity array'])

    #precursor information
    precursor = spectrum_dict['precursorList']['precursor'][0]
    precursor_ion = precursor['selectedIonList']['selectedIon'][0]
    precursor_mz = precursor_ion['selected ion m/z']
    if 'peak intensity' in precursor_ion:
        precursor_intenisty =  precursor_ion['peak intensity']
    else:
        precursor_intenisty = None
    if 'charge state' in precursor_ion:
        precursor_charge = int(precursor_ion['charge state'])
    elif 'possible charge state' in precursor_ion:
        precursor_charge = int(precursor_ion['possible charge state'])
    else:
        precursor_charge = 'NAN'

    all_info = [mz_array,intensity_array,precursor_intenisty]

    return all_info






# def append_mzml_info():
#     list_of_file_names = [ "bulk_rep1", 'bulk_rep2', 'bulk_rep3',
#                             '2ng_rep1','2ng_rep2','2ng_rep3','2ng_rep4','2ng_rep5','2ng_rep6',
#                             '0.2ng_rep1','0.2ng_rep2','0.2ng_rep3','0.2ng_rep4','0.2ng_rep5','0.2ng_rep6',
#                             'sc_rep1','sc_rep2','sc_rep3','sc_rep4','sc_rep5']
#
#
#     mzml = {}
#     mzml["bulk_rep1"] = "data/mzMLs/OR11_20160122_PG_HeLa_CVB3_CT_A.mzML"
#     mzml["bulk_rep2"] = "data/mzMLs/OR11_20160122_PG_HeLa_CVB3_CT_B.mzML"
#     mzml["bulk_rep3"] = "data/mzMLs/OR11_20160122_PG_HeLa_CVB3_CT_C.mzML"
#
#     mzml["2ng_rep1"] = "data/mzMLs/Ex_Auto_J3_30umTB_2ngQC_60m_1.mzML"
#     mzml["2ng_rep2"] = "data/mzMLs/Ex_Auto_J3_30umTB_2ngQC_60m_2.mzML"
#     mzml["2ng_rep3"] = "data/mzMLs/Ex_Auto_K13_30umTA_2ngQC_60m_1.mzML"
#     mzml["2ng_rep4"] = "data/mzMLs/Ex_Auto_K13_30umTA_2ngQC_60m_2.mzML"
#     mzml["2ng_rep5"] = "data/mzMLs/Ex_Auto_W17_30umTB_2ngQC_60m_1.mzML"
#     mzml["2ng_rep6"] = "data/mzMLs/Ex_Auto_W17_30umTB_2ngQC_60m_2.mzML"
#     #
#     mzml["0.2ng_rep1"] = "data/mzMLs/Ex_Auto_J3_30umTB_02ngQC_60m_1.mzML"
#     mzml["0.2ng_rep2"] = "data/mzMLs/Ex_Auto_J3_30umTB_02ngQC_60m_2.mzML"
#     mzml["0.2ng_rep3"] = "data/mzMLs/Ex_Auto_K13_30umTA_02ngQC_60m_1.mzML"
#     mzml["0.2ng_rep4"] = "data/mzMLs/Ex_Auto_K13_30umTA_02ngQC_60m_2.mzML"
#     mzml["0.2ng_rep5"] = "data/mzMLs/Ex_Auto_W17_30umTA_02ngQC_60m_3.mzML"
#     mzml["0.2ng_rep6"] = "data/mzMLs/Ex_Auto_W17_30umTA_02ngQC_60m_4.mzML"
#
#     mzml["sc_rep1"] = "data/mzMLs/D19_15um30cm_SC1.mzML"
#     mzml["sc_rep2"] = "data/mzMLs/D19_15um30cm_SC2.mzML"
#     mzml["sc_rep3"] = "data/mzMLs/D19_15um30cm_SC3.mzML"
#     mzml["sc_rep4"] = "data/mzMLs/D19_15um30cm_SC4.mzML"
#     mzml["sc_rep5"] = "data/mzMLs/D19_15um30cm_SC5.mzML"
#
#
#     mm_files = {}
#     mm_files["bulk_rep1"] = "data/parsed_psm/bulk_rep1.csv"
#     mm_files["bulk_rep2"] = "data/parsed_psm/bulk_rep2.csv"
#     mm_files["bulk_rep3"] = "data/parsed_psm/bulk_rep3.csv"
#
#     mm_files["2ng_rep1"] = "data/parsed_psm/2ng_rep1.csv"
#     mm_files["2ng_rep2"] = "data/parsed_psm/2ng_rep2.csv"
#     mm_files["2ng_rep3"] = "data/parsed_psm/2ng_rep3.csv"
#     mm_files["2ng_rep4"] = "data/parsed_psm/2ng_rep4.csv"
#     mm_files["2ng_rep5"] = "data/parsed_psm/2ng_rep5.csv"
#     mm_files["2ng_rep6"] = "data/parsed_psm/2ng_rep6.csv"
#
#     mm_files["0.2ng_rep1"] = "data/parsed_psm/0.2ng_rep1.csv"
#     mm_files["0.2ng_rep2"] = "data/parsed_psm/0.2ng_rep2.csv"
#     mm_files["0.2ng_rep3"] = "data/parsed_psm/0.2ng_rep3.csv"
#     mm_files["0.2ng_rep4"] = "data/parsed_psm/0.2ng_rep4.csv"
#     mm_files["0.2ng_rep5"] = "data/parsed_psm/0.2ng_rep5.csv"
#     mm_files["0.2ng_rep6"] = "data/parsed_psm/0.2ng_rep6.csv"
#
#     mm_files["sc_rep1"] = "data/parsed_psm/sc_rep1.csv"
#     mm_files["sc_rep2"] = "data/parsed_psm/sc_rep2.csv"
#     mm_files["sc_rep3"] = "data/parsed_psm/sc_rep3.csv"
#     mm_files["sc_rep4"] = "data/parsed_psm/sc_rep4.csv"
#     mm_files["sc_rep5"] = "data/parsed_psm/sc_rep5.csv"
#
#     for file_name in list_of_file_names:
#         #get path to file
#         path_to_data_loader = os.path.abspath(os.path.dirname(__file__)) # This gets the absolute path to the location of the data_loader.py file
#         complete_path_to_mzml = os.path.join(path_to_data_loader, mzml.get(file_name))
#         complete_path_to_psm = os.path.join(path_to_data_loader, mm_files.get(file_name))
#
#         mzml = pyteomics.mzml.MzML('complete_path_to_mzml')
#
#         psm = pd.read_csv(complete_path_to_psm, sep = '\t', dtype={'Base Sequence': str, 'Missed Cleavages': str, 'Peptide Monoisotopic Mass':str, 'Mass Diff (Da)':str,'Mass Diff (ppm)':str,'	Protein Accession':str,'Peptide Description':str, 'Notch':str, 'Num Variable Mods':str, 'Decoy': str})
#
# def download_mzml():
#     list_of_file_names = [ "bulk_rep1", 'bulk_rep2', 'bulk_rep3',
#                             '2ng_rep1','2ng_rep2','2ng_rep3','2ng_rep4','2ng_rep5','2ng_rep6',
#                             '0.2ng_rep1','0.2ng_rep2','0.2ng_rep3','0.2ng_rep4','0.2ng_rep5','0.2ng_rep6']
#
#
#     mzml = {}
#     mzml["bulk_rep1"] = "data/mzMLs/OR11_20160122_PG_HeLa_CVB3_CT_A.mzML"
#     mzml["bulk_rep2"] = "data/mzMLs/OR11_20160122_PG_HeLa_CVB3_CT_B.mzML"
#     mzml["bulk_rep3"] = "data/mzMLs/OR11_20160122_PG_HeLa_CVB3_CT_C.mzML"
#
#     mzml["2ng_rep1"] = "data/mzMLs/Ex_Auto_J3_30umTB_2ngQC_60m_1.mzML"
#     mzml["2ng_rep2"] = "data/mzMLs/Ex_Auto_J3_30umTB_2ngQC_60m_2.mzML"
#     mzml["2ng_rep3"] = "data/mzMLs/Ex_Auto_K13_30umTA_2ngQC_60m_1.mzML"
#     mzml["2ng_rep4"] = "data/mzMLs/Ex_Auto_K13_30umTA_2ngQC_60m_2.mzML"
#     mzml["2ng_rep5"] = "data/mzMLs/Ex_Auto_W17_30umTB_2ngQC_60m_1.mzML"
#     mzml["2ng_rep6"] = "data/mzMLs/Ex_Auto_W17_30umTB_2ngQC_60m_2.mzML"
#     #
#     mzml["0.2ng_rep1"] = "data/mzMLs/Ex_Auto_J3_30umTB_02ngQC_60m_1.mzML"
#     mzml["0.2ng_rep2"] = "data/mzMLs/Ex_Auto_J3_30umTB_02ngQC_60m_2.mzML"
#     mzml["0.2ng_rep3"] = "data/mzMLs/Ex_Auto_K13_30umTA_02ngQC_60m_1.mzML"
#     mzml["0.2ng_rep4"] = "data/mzMLs/Ex_Auto_K13_30umTA_02ngQC_60m_2.mzML"
#     mzml["0.2ng_rep5"] = "data/mzMLs/Ex_Auto_W17_30umTA_02ngQC_60m_3.mzML"
#     mzml["0.2ng_rep6"] = "data/mzMLs/Ex_Auto_W17_30umTA_02ngQC_60m_4.mzML"
#
#
#     for file_name in list_of_file_names:
#         #get path to file
#         path_to_data_loader = os.path.abspath(os.path.dirname(__file__)) # This gets the absolute path to the location of the data_loader.py file
#         complete_path_to_data = os.path.join(path_to_data_loader, mzml.get(file_name))
#
#         columns = ["scan",'ms_level',"total_ion_curr","ion_time",'mz_array', 'scan_start_time', 'intensity_array', 'retention_time','precursor_charge']
#         lst = []
#         with pyteomics.mzml.MzML(complete_path_to_data) as f_in:
#             for spectrum_dict in f_in:
#                 scan_num = spectrum_dict["id"][spectrum_dict["id"].find('scan=') + 5:]
#                 ms_level = spectrum_dict["ms level"]
#                 total_ion_curr = spectrum_dict["total ion current"]
#                 ion_time = spectrum_dict["scanList"]['scan'][0]["ion injection time"]
#                 mz_array = list(spectrum_dict['m/z array'])
#                 scan_start_time = spectrum_dict['scanList']['scan'][0]['scan start time']
#
#                 #added for graphing purposes
#                 if ms_level == 2:
#                     intensity_array = list(spectrum_dict['intensity array'])
#                     retention_time = (spectrum_dict['scanList']['scan'][0].get('scan start time', -1))
#                     precursor = spectrum_dict['precursorList']['precursor'][0]
#                     precursor_ion = precursor['selectedIonList']['selectedIon'][0]
#                     precursor_mz = precursor_ion['selected ion m/z']
#                     if 'charge state' in precursor_ion:
#                         precursor_charge = int(precursor_ion['charge state'])
#                     elif 'possible charge state' in precursor_ion:
#                         precursor_charge = int(precursor_ion['possible charge state'])
#                     else:
#                         precursor_charge = 'NAN'
#                         # raise ValueError('Unknown precursor charge')
#
#
#                     lst.append([scan_num, ms_level, total_ion_curr, ion_time, mz_array, scan_start_time, intensity_array, retention_time, precursor_charge])
#
#                 else:
#                     lst.append([scan_num, ms_level, total_ion_curr, ion_time, mz_array, scan_start_time, 'NAN', 'NAN','NAN'])
#
#
#         mzml_df = pd.DataFrame(lst, columns=columns)
#
#         mzml_df[["scan"]] = mzml_df[["scan"]].apply(pd.to_numeric)
#         mzml_df["scan_start_time"] = mzml_df["scan_start_time"].astype(str)
#
#         #split the time column so you get one that's just minutes
#         mzml_df["temp_minute"] = mzml_df["scan_start_time"].str.split("\.")
#         mzml_df.loc[:, 'minute'] = mzml_df['temp_minute'].map(lambda x: x[0])
#
#         mzml_df[["minute"]] = mzml_df[["minute"]].apply(pd.to_numeric)
#         mzml_df= mzml_df.sort_values("minute")
#
#         write_file_path = "data/parsed_mzml/" + file_name + ".csv"
#         mzml_df.to_csv(write_file_path, sep="\t")
# def load_mzml(file_name):
#     mzml = {}
#     mzml["bulk_rep1"] = "data/parsed_mzml/bulk_rep1.csv"
#     mzml["bulk_rep2"] = "data/parsed_mzml/bulk_rep2.csv"
#     mzml["bulk_rep3"] = "data/parsed_mzml/bulk_rep3.csv"
#
#     mzml["2ng_rep1"] = "data/parsed_mzml/2ng_rep1.csv"
#     mzml["2ng_rep2"] = "data/parsed_mzml/2ng_rep2.csv"
#     mzml["2ng_rep3"] = "data/parsed_mzml/2ng_rep3.csv"
#     mzml["2ng_rep4"] = "data/parsed_mzml/2ng_rep4.csv"
#     mzml["2ng_rep5"] = "data/parsed_mzml/2ng_rep5.csv"
#     mzml["2ng_rep6"] = "data/parsed_mzml/2ng_rep6.csv"
#
#     mzml["0.2ng_rep1"] = "data/parsed_mzml/0.2ng_rep1.csv"
#     mzml["0.2ng_rep2"] = "data/parsed_mzml/0.2ng_rep2.csv"
#     mzml["0.2ng_rep3"] = "data/parsed_mzml/0.2ng_rep3.csv"
#     mzml["0.2ng_rep4"] = "data/parsed_mzml/0.2ng_rep4.csv"
#     mzml["0.2ng_rep5"] = "data/parsed_mzml/0.2ng_rep5.csv"
#     mzml["0.2ng_rep6"] = "data/parsed_mzml/0.2ng_rep6.csv"
#
#     path_to_data_loader = os.path.abspath(os.path.dirname(__file__)) # This gets the absolute path to the location of the data_loader.py file
#     complete_path_to_data = os.path.join(path_to_data_loader, mzml.get(file_name))
#
#
#     df = pd.read_csv(complete_path_to_data, sep='\t')
#
#     return(df)
#
# def download_joined_psm_mzml():
#     list_of_file_names = [ "bulk_rep1", 'bulk_rep2', 'bulk_rep3',
#                             '2ng_rep1','2ng_rep2','2ng_rep3','2ng_rep4','2ng_rep5','2ng_rep6',
#                             '0.2ng_rep1','0.2ng_rep2','0.2ng_rep3','0.2ng_rep4','0.2ng_rep5','0.2ng_rep6']
#
#     for file_name in list_of_file_names:
#         df_mm = load_psm(file_name)
#         df_mz = load_mzml(file_name)
#
#         #join the two together
#         joined = pd.merge(df_mm, df_mz, on="scan").sort_values("scan")
#
#         write_file_path = "data/joined_psm_and_mzml/" + file_name + ".csv"
#         joined.to_csv(write_file_path, sep="\t")
# def load_joined_psm_mzml(file_name):
#     joined_files = {}
#     joined_files["bulk_rep1"] = "data/joined_psm_and_mzml/bulk_rep1.csv"
#     joined_files["bulk_rep2"] = "data/joined_psm_and_mzml/bulk_rep2.csv"
#     joined_files["bulk_rep3"] = "data/joined_psm_and_mzml/bulk_rep3.csv"
#
#     joined_files["2ng_rep1"] = "data/joined_psm_and_mzml/2ng_rep1.csv"
#     joined_files["2ng_rep2"] = "data/joined_psm_and_mzml/2ng_rep2.csv"
#     joined_files["2ng_rep3"] = "data/joined_psm_and_mzml/2ng_rep3.csv"
#     joined_files["2ng_rep4"] = "data/joined_psm_and_mzml/2ng_rep4.csv"
#     joined_files["2ng_rep5"] = "data/joined_psm_and_mzml/2ng_rep5.csv"
#     joined_files["2ng_rep6"] = "data/joined_psm_and_mzml/2ng_rep6.csv"
#
#     joined_files["0.2ng_rep1"] = "data/joined_psm_and_mzml/0.2ng_rep1.csv"
#     joined_files["0.2ng_rep2"] = "data/joined_psm_and_mzml/0.2ng_rep2.csv"
#     joined_files["0.2ng_rep3"] = "data/joined_psm_and_mzml/0.2ng_rep3.csv"
#     joined_files["0.2ng_rep4"] = "data/joined_psm_and_mzml/0.2ng_rep4.csv"
#     joined_files["0.2ng_rep5"] = "data/joined_psm_and_mzml/0.2ng_rep5.csv"
#     joined_files["0.2ng_rep6"] = "data/joined_psm_and_mzml/0.2ng_rep6.csv"
#
#     path_to_data_loader = os.path.abspath(os.path.dirname(__file__)) # This gets the absolute path to the location of the data_loader.py file
#     complete_path_to_data = os.path.join(path_to_data_loader, joined_files.get(file_name))
#
#     df = pd.read_csv(complete_path_to_data, sep='\t', index_col=0, dtype={'Base Sequence': str, 'Missed Cleavages': str, 'Peptide Monoisotopic Mass':str, 'Mass Diff (Da)':str,'Mass Diff (ppm)':str,'	Protein Accession':str,'Peptide Description':str, 'Notch':str, 'Num Variable Mods':str, 'Peptide Description':str, 'Decoy': str})
#
#     return(df)
