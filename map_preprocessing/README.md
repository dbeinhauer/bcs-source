# Předpříprava mapy

Tento adresář obsahuje program pro přípravu mapy do vstupního formátu
programu pro simulování dopravy a optimalizaci rozmístění nabíjecích stanic,
jehož zdrojový kód se nachází v adresáři `source_code`.

## Požadavky
Program vyžaduje pro správné fungování knihovnu Pyrosm a Numpy.

## Příklad použití

```
./map_reader.py \
--input_map=$(road_network_input_file_OSM_PBF_format) \
--input_cities=$(cities_input_file_OSM_PBF_format) \
--output_nodes=$(nodes_output_in_simulator_format) \
--output_edges=$(edges_output_in_simulator_format) \
--output_cities=$(cities_output_in_simulator_format) \
--output_combined=$(nodes_cities_output_in_simulator_format) \
```

## Redukce mapy pomocí nástroje Osmium
V případě rozsáhlé dopravní sítě může být program pro předpřípravu mapy
velmi paměťově náročný. Pro velké mapy proto doporučujeme zredukovat
původní mapu pouze na mapu dopravní sítě, či mapu měst pomocí nástroje
Osmium. Níže uvádíme možné příklady použití.

Redukce mapy na mapu dopravní sítě se silnicemi typu `motorway`,
`trunk` a `primary`:

```
osmium tags-filter czech-republic-latest.osm.pbf \
w/highway=motorway,trunk,primary,motorway_link,trunk_link, primary_link
-o top-roads-czech-republic.osm.pbf
```

Redukce mapy na mapu osídlených oblastí typů `city` a `town`:

```
osmium tags-filter maps/czech-republic-latest.osm.pbf \
place=city,town \
-o maps/czech_cities.osm.pbf
```