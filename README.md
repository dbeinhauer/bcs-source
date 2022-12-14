# About
This repository contains the source code of my bachelor thesis 
with the name "Optimization of the placement of the charging 
stations for the electric vehicles". The thesis was written
during my studies in the Faculty of Mathematics and Physics in the 
Charles University in Prague. The text of the thesis could be found
in the repository [bcs-thesis]() (note that the text is written in 
Czech).

# Structure of The Repository
The repository contains multiple files and directories with the 
specific content. The repository contains:

* `analysis_results/` - directory which contains information about the 
experiments and its results
* `experiment_analysis/` - directory which contains programs for 
preparation of the scripts to run the experiments and display the
output data
* `map_preprocessing/` - directory which contains program for
preparation of the input data for the traffic simulator
* `prepared_grapg/` - directory which contains the prepared
input data which represents the road network of the Czech Republic
* `source_code/` - directory which contains the source code of the
traffic simulator which is the main program of the thesis
* `Traffic_Simulator_Manual.pdf` - file with the development 
documentation generated from the comments in the source code
(Note that more detailed documentation can be found in the text of the 
thesis [here]() (only in Czech))


# Abstract (EN)
As the number of electric vehicles grows, so does the need to create a suitable
network of charging stations. A solution of this problem can be significantly
improved by the usage of suitable optimization techniques.
We implement a simplified traffic simulator serving as a suitable tool for 
their analysis.
We also analyze optimization techniques using the so-called greedy algorithm,
genetic algorithm and k-means algorithm. Based on the experiments, the optimizations 
using the genetic algorithm and the greedy algorithm showed noticeably better results.
The k-means method did not show signs of results better than a
random approach.


# Abstract (CZ)
S rostouc??m po??tem elektrick??ch vozidel roste i pot??eba vytvo??it vhodnou 
infrastrukturu pro jejich nab??jen??. K ??e??en?? tohoto probl??mu m????e v??razn?? 
napomoci pou??it?? vhodn??ch optimaliza??n??ch metod. V pr??ci jsme implementovali 
zjednodu??en?? simul??tor dopravy slou????c?? jako vhodn?? n??stroj pro jejich anal??zu.
Analyzovali jsme tak?? optimaliza??n?? metody tzv. hladov??m algoritmem, 
genetick??m algoritmem a algoritmem k-means. Na z??klad?? experiment?? 
vykazovala prokazateln?? lep???? v??sledky optimalizace za vyu??it?? 
genetick??ho algoritmu a hladov?? optimalizace. K-means optimalizace 
nevykazovala zn??mky lep????ch v??sledk?? oproti n??hodn??mu p????stupu.