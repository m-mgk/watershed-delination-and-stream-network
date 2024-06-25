from pysheds.grid import Grid
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from pyproj import Proj, transform
from shapely.geometry import LineString, Point, Polygon
import geopandas as gpd
import numpy as np
import networkx as nx
from geopy.distance import geodesic
import warnings
warnings.filterwarnings('ignore')
import os
from datetime import datetime
import folium

# DEM MUST BE IN WGS 1984 COORDINATES
##########################################
# INPUTS
# Enter the name of the file that must be in the same directory
dem_file_name = 'example_dem_filename.tif'
# pour point location in decimal format
lat = 32.455136  
long = 71.921049
# stream flow accumulation lower threshold
accumulation_value = 50
# stream flow accumulation lower threshold for the pour point to snap to closest stream
accumulation_value_for_snap_pour_point=1000 #generally 1000 values works ok.
##########################################

#creating floder to save the results
current_time = datetime.now().strftime('%Y%m%d_%H%M%S')
folder_name = f"{dem_file_name}_{current_time}"
# Create the folder
os.makedirs(folder_name, exist_ok=True)
##########################################

# loading the grid and dem
grid = Grid.from_raster(dem_file_name)
dem = grid.read_raster(dem_file_name)

# plotting DEM
fig, ax = plt.subplots(figsize=(8, 6))
fig.patch.set_alpha(0)
plt.imshow(dem, extent=grid.extent, cmap='terrain', zorder=1)
plt.colorbar(label='Elevation (m)')
plt.grid(zorder=0)
plt.title('Digital elevation map', size=14)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.tight_layout()
# Define the file path for saving the plot
plot_file_path = os.path.join(folder_name, 'Digital_elevation_map')
# Save the plot
plt.savefig(plot_file_path, bbox_inches='tight')

# Fill pits in DEM
pit_filled_dem = grid.fill_pits(dem)

# Fill depressions in DEM
flooded_dem = grid.fill_depressions(pit_filled_dem)

# Resolve flats in DEM
inflated_dem = grid.resolve_flats(flooded_dem)

# Specify directional mapping
dirmap = (64, 128, 1, 2, 4, 8, 16, 32)
fdir = grid.flowdir(inflated_dem, dirmap=dirmap)

# plotting flow direction

fig = plt.figure(figsize=(8, 6))
fig.patch.set_alpha(0)
plt.imshow(fdir, extent=grid.extent, cmap='viridis', zorder=2)
boundaries = ([0] + sorted(list(dirmap)))
plt.colorbar(boundaries=boundaries,
             values=sorted(dirmap))
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Flow direction grid', size=14)
plt.grid(zorder=-1)
plt.tight_layout()
# Define the file path for saving the plot
plot_file_path = os.path.join(folder_name, 'flow_direction_grid')
# Save the plot
plt.savefig(plot_file_path, bbox_inches='tight')


# flow accumulation code:
acc = grid.accumulation(fdir, dirmap=dirmap)
# plotting flow accumulation
fig, ax = plt.subplots(figsize=(8, 6))
fig.patch.set_alpha(0)
plt.grid('on', zorder=0)
im = ax.imshow(acc, extent=grid.extent, zorder=2,
               cmap='cubehelix',
               norm=colors.LogNorm(1, acc.max()),
               interpolation='bilinear')
plt.colorbar(im, ax=ax, label='Upstream Cells')
plt.title('Flow Accumulation', size=14)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.tight_layout()
# Define the file path for saving the plot
plot_file_path = os.path.join(folder_name, 'flow_accumulation_map')
# Save the plot
plt.savefig(plot_file_path, bbox_inches='tight')


# Export to tiff
# grid.to_raster(acc, 'accumulation_image.tif')

# Snap pour point to high accumulation cell
x, y = long, lat
x_snap, y_snap = grid.snap_to_mask(acc > accumulation_value_for_snap_pour_point, (x, y))

# Delineate the catchment
catch = grid.catchment(x=x_snap, y=y_snap, fdir=fdir, dirmap=dirmap, xytype='coordinate')

# Clip the bounding box to the catchment
grid.clip_to(catch)
clipped_catch = grid.view(catch)

# Extract river network
branches = grid.extract_river_network(fdir, acc > accumulation_value, dirmap=dirmap)

# Plotting the streams network with directions
line_strings = []

for branch in branches['features']:
    line_coords = branch['geometry']['coordinates']
    line = LineString(line_coords)
    line_strings.append(line)

gdf = gpd.GeoDataFrame(geometry=line_strings)


# Define the file path
file_path = os.path.join(folder_name, 'rivers.shp')

# Save the shapefile
gdf = gpd.GeoDataFrame(geometry=line_strings)
gdf.to_file(file_path)


#calculating the catchment area boundary and saving as geojason

# Ensure the catchment array is of type uint8
catch_uint8 = catch.astype(np.uint8)

# Extract the catchment boundary
shapes = list(grid.polygonize(catch_uint8))
catch_boundary = [Polygon(shape['coordinates'][0]) for shape, value in shapes if value]

# Create a GeoDataFrame for the catchment boundary
catchment_gdf = gpd.GeoDataFrame(geometry=catch_boundary)

# Create the folder (if not already created)
os.makedirs(folder_name, exist_ok=True)

# Define the file path for GeoJSON
geojson_file_path = os.path.join(folder_name, 'catchment_boundary.geojson')

# Save the catchment boundary as a GeoJSON file
catchment_gdf.to_file(geojson_file_path, driver='GeoJSON')



# Reload the DEM to ensure fresh data
grid = Grid.from_raster(dem_file_name)
dem = grid.read_raster(dem_file_name)

# Calculate the principal channel length and slope considering the pour point

def create_network_from_branches(branches):
    G = nx.Graph()
    for branch in branches['features']:
        line = np.asarray(branch['geometry']['coordinates'])
        for start, end in zip(line[:-1], line[1:]):
            dist = geodesic(start[::-1], end[::-1]).meters
            G.add_edge(tuple(start), tuple(end), weight=dist)
    return G

def snap_to_nearest_node(G, point):
    nearest_node = min(G.nodes, key=lambda node: geodesic((point[1], point[0]), (node[1], node[0])).meters)
    return nearest_node

def find_principal_channel(G, pour_point):
    all_pairs_lengths = nx.single_source_dijkstra_path_length(G, pour_point, weight='weight')
    max_length = 0
    principal_channel = None
    end_point = None

    for target, length in all_pairs_lengths.items():
        if length > max_length:
            max_length = length
            end_point = target
    
    if end_point:
        principal_channel = nx.shortest_path(G, pour_point, end_point, weight='weight')
    
    return principal_channel, max_length

def get_elevation(dem, point, transform):
    col, row = ~transform * point
    col, row = int(col), int(row)
    return dem[row, col]

def calculate_slope_per_section(principal_channel, dem, transform):
    slopes = []
    elevations = [get_elevation(dem, (x, y), transform) for x, y in principal_channel]
    distances = [geodesic((principal_channel[i][1], principal_channel[i][0]), 
                          (principal_channel[i+1][1], principal_channel[i+1][0])).meters 
                 for i in range(len(principal_channel) - 1)]
    for i in range(len(distances)):
        slope = (-elevations[i] + elevations[i+1]) / distances[i]
        slopes.append(slope)
    return slopes, elevations, distances

def calculate_weighted_average_slope(slopes, distances):
    weighted_slope_sum = sum(slope * distance for slope, distance in zip(slopes, distances))
    total_distance = sum(distances)
    return weighted_slope_sum / total_distance

G = create_network_from_branches(branches)

# Snap pour point to the nearest stream point in the river network
pour_point = (x_snap, y_snap)
pour_point = snap_to_nearest_node(G, pour_point)

principal_channel, max_length = find_principal_channel(G, pour_point)
principal_channel_length_km = max_length / 1000
slopes, elevations, distances = calculate_slope_per_section(principal_channel, dem, grid.affine)

# Calculate the weighted average slope considering varying slopes
average_slope = calculate_weighted_average_slope(slopes, distances)

# Plot the principal channel and the stream network
fig, ax = plt.subplots(figsize=(8, 6))
fig.patch.set_alpha(0)

# Plot the stream network
gdf.plot(ax=ax, color='blue', label='Stream Network')

# Plot the principal channel
principal_channel_line = LineString(principal_channel)
principal_channel_gdf = gpd.GeoDataFrame(geometry=[principal_channel_line])
principal_channel_gdf.plot(ax=ax, color='red', linewidth=2, label='Principal Channel')

# Plot the pour points
provided_pour_point = gpd.GeoDataFrame(geometry=[Point(long, lat)])
snapped_pour_point = gpd.GeoDataFrame(geometry=[Point(x_snap, y_snap)])
provided_pour_point.plot(ax=ax, color='green', markersize=50, label='Provided Pour Point')
snapped_pour_point.plot(ax=ax, color='yellow', markersize=50, label='Snapped Pour Point')

# Customize title and labels
plt.title('Stream Network with Principal Channel and Pour Points', size=14)
plt.xlabel('Longitude')
plt.ylabel('Latitude')

# Create a custom legend and place it at the bottom of the plot
handles, labels = ax.get_legend_handles_labels()
legend = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, -0.08), ncol=4, frameon=False)

# Customize grid lines to be dotted
plt.grid(True, linestyle='dotted', linewidth=0.5, zorder=0)

# Adjust the layout
plt.tight_layout()

# Define the file path for saving the plot
plot_file_path = os.path.join(folder_name, 'stream_network_plot.png')

# Save the plot
plt.savefig(plot_file_path, bbox_inches='tight')

# Show the plot
plt.show()


# Plot the elevation profile of the principal channel

# Calculate distances along the principal channel
distances = [0] + np.cumsum([geodesic((principal_channel[i][1], principal_channel[i][0]), 
                                      (principal_channel[i+1][1], principal_channel[i+1][0])).meters 
                             for i in range(len(principal_channel) - 1)]).tolist()

# Assuming `elevations` is already defined from previous calculations
# For demonstration purposes, let's define it with sample data
elevations = [get_elevation(dem, (x, y), grid.affine) for x, y in principal_channel]

# Plot the elevation profile of the principal channel
fig, ax = plt.subplots(figsize=(12, 8))
fig.patch.set_alpha(0)

# Plot the elevation profile
ax.plot(distances, elevations, linestyle='-', color='brown', label='Elevation')
ax.fill_between(distances, elevations, color='brown', alpha=0.3)

# Set title and labels with customized font sizes
ax.set_title('Elevation Profile of the Principal Channel', fontsize=16)
ax.set_xlabel('Distance along the channel (meters)', fontsize=14)
ax.set_ylabel('Elevation (meters)', fontsize=14)

# Customize ticks
ax.tick_params(axis='both', which='major', labelsize=12)

# Add grid lines with customization
ax.grid(True, linestyle='dashed', linewidth=0.9, alpha=1)


# Remove the top and right spines for a cleaner look
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Adjust the layout to make room for labels and title
plt.tight_layout()

# Define the file path for saving the plot
plot_file_path = os.path.join(folder_name, 'principal_channel_elevation_profile')
# Save the plot
plt.savefig(plot_file_path, bbox_inches='tight')

# Show the plot
plt.show()


# Catchment area plot and calculations

def approximate_cell_area(lat, cell_size_deg):
    R = 6378137
    cell_size_rad = np.radians(cell_size_deg)
    return (R ** 2) * cell_size_rad ** 2 * np.cos(np.radians(lat))


#calculating the number of cells and catchment area
if catch is not None and np.any(catch):
    cell_size_deg = abs(grid.affine[0])

    lats, lons = np.where(catch)
    latitudes = grid.affine[5] + lats * grid.affine[4]
    avg_lat = np.mean(latitudes)

    cell_area_m2 = approximate_cell_area(avg_lat, cell_size_deg)
    num_cells = np.sum(catch > 0)
    # catchment area in sq.km
    watershed_area = (num_cells * cell_area_m2) / 10**6

else:
    print("Catchment mask is empty or invalid")


# Plot the catchment

# Convert DEM to float32
dem = dem.astype(np.float32)


# Create a mask for the catchment area
catchment_mask = catch.astype(bool)

# Mask the DEM outside the catchment area
masked_dem = np.where(catchment_mask, dem, np.nan)

# Determine the bounding box of the catchment area
minx, miny, maxx, maxy = catchment_gdf.total_bounds

# Calculate a buffer as a fraction of the catchment size
buffer_x = (maxx - minx) * 0.1
buffer_y = (maxy - miny) * 0.1

# Plot the catchment
fig, ax = plt.subplots(figsize=(10, 8))
fig.patch.set_alpha(0)

# Plot masked DEM with transparency
plt.imshow(masked_dem, extent=grid.extent, cmap='terrain', zorder=0, alpha=0.5)

# Adjust plot limits to the catchment area with a buffer
plt.xlim(minx - buffer_x, maxx + buffer_x)
plt.ylim(miny - buffer_y, maxy + buffer_y)

# Plot catchment boundary
catchment_gdf.boundary.plot(ax=ax, color='red', linewidth=2, label='Catchment Boundary')

# Plot stream network
gdf.plot(ax=ax, color='blue', linewidth=1, zorder=1, label='Stream Network')

# Plot principal stream
principal_channel_gdf.plot(ax=ax, color='purple', linewidth=2, zorder=2, label='Principal Stream')

# Plot pour points
provided_pour_point.plot(ax=ax, color='green', markersize=50, label='Provided Pour Point', zorder=2)
snapped_pour_point.plot(ax=ax, color='yellow', markersize=50, label='Snapped Pour Point', zorder=2)

# Create a custom legend
handles, labels = ax.get_legend_handles_labels()
legend = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, -0.08), ncol=4, frameon=False)

# Add color bar
cbar = plt.colorbar(ax.images[0], ax=ax, orientation='vertical', fraction=0.046, pad=0.04)
cbar.set_label('Elevation (m)')

# Turn off the grid
ax.grid(False)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Delineated Catchment with Stream Network', size=14)
plt.tight_layout()
# Define the file path for saving the plot
plot_file_path = os.path.join(folder_name, 'catchment_area_plot')
# Save the plot
plt.savefig(plot_file_path, bbox_inches='tight')

plt.show()

#interactive map

# Set CRS for the catchment boundary GeoDataFrame
catchment_gdf.crs = {'init': 'epsg:4326'}

# Set CRS for the river network GeoDataFrame
gdf.crs = {'init': 'epsg:4326'}

# Set CRS for the principal channel GeoDataFrame
principal_channel_gdf.crs = {'init': 'epsg:4326'}

# Create a folium map centered at the pour point
folium_map = folium.Map(location=[lat, long], zoom_start=12, control_scale=True)

# Add the catchment boundary to the map
catchment_gdf_folium = folium.GeoJson(
    catchment_gdf,
    name="Watershed Boundary",
    style_function=lambda x: {
        'color': 'red',
        'weight': 2,
        'fill': False
    }
).add_to(folium_map)

# Add the river network to the map
rivers_gdf = gdf.copy()
rivers_gdf['linewidth'] = 2  # Add a column with a constant linewidth

rivers_gdf_folium = folium.GeoJson(
    rivers_gdf,
    name="River Centerlines",
    style_function=lambda x: {
        'color': 'blue',
        'weight': x['properties']['linewidth']
    }
).add_to(folium_map)

# Add the principal channel to the map
principal_channel_gdf_folium = folium.GeoJson(
    principal_channel_gdf,
    name="Principal Channel",
    style_function=lambda x: {
        'color': 'purple',  # Change the color to your preference
        'weight': 3
    }
).add_to(folium_map)

# Add satellite imagery layer to the map
satellite_tiles = folium.TileLayer(
    tiles='https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
    attr='Esri',
    name='Esri Satellite',
    overlay=False,
    control=True
).add_to(folium_map)

# Add layer control to the map
folium.LayerControl().add_to(folium_map)

# Save the map to an HTML file
map_file_path = os.path.join(folder_name, 'catchment_map.html')
folium_map.save(map_file_path)

# Optionally display the map
# folium_map  # Uncomment this line if you are running in a Jupyter Notebook and want to display the map inline





############################################################################
# printing outputs
# Define the file path for the text file
output_text_file_path = os.path.join(folder_name, 'analysis_output.txt')

# Define the outputs
outputs = [
    f"DEM filename: {dem_file_name}",
    f"Location: latitude {lat}, longitude {long}",
    f"Principal Channel Length: {principal_channel_length_km:.2f} km",
    f"Stream Node 1 Elevation: {elevations[0]:.3f} meters",
    f"Stream Node 2 Elevation: {elevations[-1]:.3f} meters",
    f"Weighted Average Slope: {average_slope:.5f} or {100*average_slope:.3f}%",
    f"Watershed Area: {watershed_area:.3f} square kilometers"
]

# Write the outputs to the text file
with open(output_text_file_path, 'w') as file:
    file.write("\n".join(outputs))

# Printing outputs to the console
print("\n".join(outputs))
