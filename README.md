# SingleCellSpectra

The data files are too large to be uploaded to github, but are publicly available on MassIVE. The single cell data can be found on MassIVE (MSV000087524) and the bulk data can be found on MassIVE (MSV000087689). 

Before any notebook can be run, the data must be downloaded and parsed. The following steps decribe how.
1) download mzMLs and put into folder /data/mzMLs.  
2) Run the script join_psm_and_mzml_files.ipynb, which will generate the data tables used in all notebooks. 
3) All data can now be accessed through data_loader.py. 
