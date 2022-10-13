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
S rostoucím počtem elektrických vozidel roste i potřeba vytvořit vhodnou 
infrastrukturu pro jejich nabíjení. K řešení tohoto problému může výrazně 
napomoci použití vhodných optimalizačních metod. V práci jsme implementovali 
zjednodušený simulátor dopravy sloužící jako vhodný nástroj pro jejich analýzu.
Analyzovali jsme také optimalizační metody tzv. hladovým algoritmem, 
genetickým algoritmem a algoritmem k-means. Na základě experimentů 
vykazovala prokazatelně lepší výsledky optimalizace za využití 
genetického algoritmu a hladová optimalizace. K-means optimalizace 
nevykazovala známky lepších výsledků oproti náhodnému přístupu.