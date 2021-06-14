# OBaNK
OBaNK takes as input a single matrix (CSV) of all normalized multi-omics data (feature x sample) generated from a two-group (case and control) experimental design. Each experimental condition (no more than two conditions) should be sequentially grouped in the input matrix. OBaNK returns three files: graphical_rendering.pdf, Strength.csv, and summary_generated_BN.txt. 

INPUT MATRIX:

Column 1:
•	Header: ‘Symbol’
•	Values: An identifying symbol or name for the feature. The symbol can be, but is not limited to, a common name (ex: plasminogen_activator), a protein/gene symbol (ex: PLAT), or a user-defined name (ex: secreted_protease_1). 
•	The symbol must be unique for each feature and can be a mixture of numbers and letters. Avoid names that include special characters other than a dash, dot, or underscore. An example would be to change the lipid PE(18:0/18:0) to PE_18.0_18.0 or something similar.

Column 2:
•	Header: ‘ID’
•	Values: The KEGG (ex: hsa5327, C00047), SwissLipid (ex: SLM000094711), Lipid Maps (ex: LMFA01010002, LMST01010057) or HMDB (ex: HMDB00237) identifier where all special characters and spaces are omitted. These values are used to query the External Data Bank.
•	If an ID is unavailable, use the value from column 1 in its place.  

Columns 3-n:
•	Headers: Group 1 names/ Group 2 names
•	Each sample should have a different header name (ex: control_1, control_2, control_3, case1a, case2a, case3a). Names can be a mixture of numbers, letters, and underscores. 
•	Values: The experimental data for the features 	

Row 2:
•	Position 1: ‘Group’
•	Position 2: ‘GR2’. The trailing number represents the number of experimental groups. 
•	Values: Place a ‘1’ for the first phenotype (ex: case samples) and ‘2’ for the second phenotype (ex: control samples).  


OUTPUT FILES:	
graphical_rendering.pdf
•	OBaNK generates two networks. The links are the same, but the node naming convention is different. One includes the names provided in column 1 (symbols/names) and the other with IDs from column 2. 
Strength.csv
•	The strength file provides the node from/to information and the associated link’s strength value. This file can be converted to a SIF for use with graphical visualization tools such as Cytoscape. 
summary_generated_BN.txt
•	This file includes a summary of the learned model and associated significance threshold (strength) used to define valid links.
