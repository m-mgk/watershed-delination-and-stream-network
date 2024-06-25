# watershed delineation and stream network 

This project involves analyzing a Digital Elevation Model (DEM) to delineate catchment areas, extract river networks, and compute various hydrological metrics. The analysis is performed using Python and several geospatial libraries.

## Table of Contents

- [Introduction](#introduction)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Outputs](#outputs)
- [License](#license)

## Introduction

The purpose of this project is to:
1. Load and process a DEM.
2. Delineate catchment areas.
3. Extract and analyze river networks.
4. Compute hydrological metrics such as the principal channel length, average slope, and watershed area.

## Requirements

The following Python packages are required:

- `pysheds`
- `matplotlib`
- `pyproj`
- `shapely`
- `geopandas`
- `numpy`
- `networkx`
- `geopy`
- `folium`

## Installation

To install the necessary packages, run the following command:

```bash
pip install pysheds matplotlib pyproj shapely geopandas numpy networkx geopy folium
```

## Usage

1. Place the DEM file (`example_dem_name.tif`) in the same directory as the script.
2. Edit the script to set the pour point location (latitude and longitude) and the stream flow accumulation threshold.

### Running the Script

Execute the script using Python:

```bash
python dem_analysis.py
```

The script will:

1. Load the DEM.
2. Plot and save the DEM.
3. Process the DEM to fill pits and depressions and resolve flats.
4. Calculate flow direction and flow accumulation.
5. Snap the pour point to the nearest high accumulation cell.
6. Delineate the catchment area.
7. Extract and save the river network and catchment boundary.
8. Compute the principal channel length and average slope.
9. Generate and save plots for various stages of the analysis.
10. Create an interactive map with the catchment boundary and river network.

## Outputs

The results are saved in a directory named after the DEM file and the current timestamp. The outputs include:

- `Digital_elevation_map.png`: Plot of the DEM.
- `flow_direction_grid.png`: Plot of the flow direction grid.
- `flow_accumulation_map.png`: Plot of the flow accumulation.
- `rivers.shp`: Shapefile of the extracted river network.
- `catchment_boundary.geojson`: GeoJSON file of the catchment boundary.
- `stream_network_plot.png`: Plot of the stream network with the principal channel and pour points.
- `principal_channel_elevation_profile.png`: Elevation profile of the principal channel.
- `catchment_area_plot.png`: Plot of the delineated catchment area with the stream network.
- `catchment_map.html`: Interactive map of the catchment area.
- `analysis_output.txt`: Text file summarizing the principal channel length, elevations, average slope, and watershed area.

### Example Output

```text
DEM filename: example_dem_name.tif
Location: latitude xxxxxxx, longitude xxxxxxx
Principal Channel Length: X.XX km
Stream Node 1 Elevation: XXX.XXX meters
Stream Node 2 Elevation: XXX.XXX meters
Weighted Average Slope: X.XXXXX or X.XXX%
Watershed Area: XX.XXX square kilometers
```

## License

This project is licensed under the MIT License. See the LICENSE file for details.

