import pandas as pd
import numpy as np
from modules.mesh import RectangularMesh
from modules.point import Point
from modules.planar import PlanarSurface

from matplotlib import pyplot
from mpl_toolkits.basemap import Basemap
import time, math
from threading import Thread
import os  # Import the os module for file existence check
# import customtkinter
import geopandas as gpd


# Your processing function
def process_eqrm_data_and_find_closest(row, sites_df):
    # Create PlanarSurface from the row data
    surf = PlanarSurface(
        strike=row['StrikeAzimuth'], dip=row['DipAngle'],
        top_left=Point(row['top_left_lon'], row['top_left_lat'], row['top_left_depth']),
        top_right=Point(row['top_right_lon'], row['top_right_lat'], row['top_right_depth']),
        bottom_left=Point(row['bottom_left_lon'], row['bottom_left_lat'], row['bottom_left_depth']),
        bottom_right=Point(row['bottom_right_lon'], row['bottom_right_lat'], row['bottom_right_depth'])
    )
    
    # Define buffer and delta for the grid
    buf = 1.8
    delta = 0.001

    # Get bounding box
    min_lon, max_lon, max_lat, min_lat = surf.get_bounding_box()
    min_lon -= buf
    max_lon += buf
    min_lat -= buf
    max_lat += buf

    # Create grid of points
    lons = np.arange(min_lon, max_lon + delta, delta)
    lats = np.arange(min_lat, max_lat + delta, delta)
    lons, lats = np.meshgrid(lons, lats)
    mesh = RectangularMesh(lons=lons, lats=lats, depths=None)

    
    r_rup = surf.get_min_distance(mesh)
    r_jb = surf.get_joyner_boore_distance(mesh)
    r_x = surf.get_rx_distance(mesh)
    r_y0 = surf.get_ry0_distance(mesh)



    # Flatten the arrays
    lats_flattened = np.ravel(mesh.lats)
    lons_flattened = np.ravel(mesh.lons)
    r_rup_flattened = np.ravel(r_rup)
    r_jb_flattened = np.ravel(r_jb)
    r_x_flattened = np.ravel(r_x)
    r_y0_flattened = np.ravel(r_y0)

    # Convert site DataFrame to numpy arrays for efficient distance computation
    site_lons = sites_df['SMLongitude'].values
    site_lats = sites_df['SMLatitude'].values

    # Initialize list for closest points
    closest_points = []

    # Compute distances between all grid points and all sites
    for _, site in sites_df.iterrows():
        # Calculate Euclidean distances
        distances = np.sqrt((lats_flattened - site['SMLatitude'])**2 + (lons_flattened - site['SMLongitude'])**2)
        
        # Find the index of the closest point
        closest_idx = np.argmin(distances)
        
        # Append closest point info to list
        closest_points.append({
            'EventName': row['EventName'],  
            'EpicenterLongitude': row['EpicenterLongitude'],  
            'EpicenterLatitude': row['EpicenterLatitude'],  
            'SMStation': site['SMStation'],
            'SMLongitude': site['SMLongitude'],
            'SMLatitude': site['SMLatitude'],
            'Vs30': site['SMVs30'],
            'PGA': None,
            'Magnitude': row['Magnitude'],
            'ModelDepth': row['ModelDepth'],
            'StrikeAzimuth': row['StrikeAzimuth'],
            'DipAngle': row['DipAngle'],
            'RakeAngle': row['RakeAngle'],
            'FaultMechanism': row['FaultMechanism'],
            'rrup': r_rup_flattened[closest_idx],
            'rjb': r_jb_flattened[closest_idx],
            'rx': r_x_flattened[closest_idx],
            'ry0': r_y0_flattened[closest_idx]
        })

        # # Version 1 cols
        # closest_points.append({
        #     'eq_event_id': row['EventName'],  
        #     'EpicenterLongitude': row['EpicenterLongitude'],  
        #     'EpicenterLatitude': row['EpicenterLatitude'],  
        #     'station': site['SMStation'],
        #     'SMLongitude': site['SMLongitude'],
        #     'SMLatitude': site['SMLatitude'],
        #     'vs30': site['SMVs30'],
        #     'pga': None,
        #     'magnitude': row['Magnitude'],
        #     'depth': row['ModelDepth'],
        #     'strike': row['StrikeAzimuth'],
        #     'dip': row['DipAngle'],
        #     'rake': row['RakeAngle'],
        #     'FaultMechanism': row['FaultMechanism'],
        #     'rrup': r_rup_flattened[closest_idx],
        #     'rjb': r_jb_flattened[closest_idx],
        #     'rx': r_x_flattened[closest_idx],
        #     'ry0': r_y0_flattened[closest_idx]
        # })

    return closest_points


def _compute_width(row) -> float:
    fault_type = row['FaultMechanism'].lower()  # Convert to lowercase for case-insensitivity
    magnitude = row['Magnitude']
    
    if 'strike-slip' in fault_type:
        return 10**(-0.76 + (0.27 * magnitude))
    elif 'reverse' in fault_type:
        return 10**(-1.61 + (0.41 * magnitude))
    elif 'normal' in fault_type:
        return 10**(-1.14 + (0.35 * magnitude))
    else:
        return 1
    
def _compute_Ztor(row) -> float:
    Zhyp = float(row['ModelDepth']) # Actual Value
    Zhyp = float(row['FaultWidth'])
    dip = float(row['DipAngle'])

    return (Zhyp - (0.6 * math.sin(dip)))

def vs30_z1pt0_AbrahamsonSilva2008(row) -> float:
    """Extracts a depth of 1.0 km/s velocity layer using the relationship proposed in Abrahamson & Silva 2008"""
    vs30 = row['Vs30']
    if vs30 < 180:
        return np.exp(6.745)                                # vs30 < 180
    elif vs30 > 500:
        return np.exp(5.3945 - 4.48 * np.log(vs30/500))     # vs30 > 500
    else:
        return np.exp(5.3945 - 1.35 * np.log(vs30/180))     # 180 <= vs30 <= 500

def vs30_z1pt0_ChiouYoungs2008(row) -> float:
    vs30 = row['Vs30']
    return np.exp(28.5 - (3.82/8) * np.log(vs30**8 + 378.7**8))

def z1pt0_z2pt5_CampbellBozorgnia2007(row) -> float:
    """ z2pt5 or z2.5 (Minimum depth (km) at which vs30 ≥ 2.5 km/s) - Vertical distance from earth surface to 
    layer where seismic waves start to propagate with a speed above 2.5 km/sec, in km.

    To estimate z2.5, Campbell and Bozorgnia (2007) offer guidelines for extrapolating z2.5 from z1.0 or z1.5"""
    
    z1pt0 = row['Z1pt0']
    return 0.519 + 3.595 * (z1pt0/1000) # in kilometers
    #return 519 + 3.595 * z1pt0          # in meters
    
def vs30_z1pt0_ChiouYoungs2014(row) -> float:
    vs30 = row['Vs30']
    
    Japan = False
    if Japan:
        c1 = 412 ** 2
        c2 = 1360 ** 2
        return np.exp((-5.23/2) * np.log((vs30**2 + c1) / (c2 + c1)))
    else:
        c1 = 571 ** 4
        c2 = 1360 ** 4
        return np.exp((-7.14/4) * np.log((vs30**4 + c1) / (c2 + c1)))
    
def vs30_z2pt5_CampbellBozorgnia2014(row):
    vs30 = row['Vs30']
    
    Japan = False
    if Japan:
        return np.exp(5.359 - 1.102 * np.log(vs30))
    else:
        return np.exp(7.089 - 1.144 * np.log(vs30))
    

def _get_distance_Repi(row):
    """
    This function 
    """
    from shapely.geometry import Point
    # try: 
    point1 = Point(float(row['EpicenterLongitude']), float(row['EpicenterLatitude']))
    point2 = Point(float(row['SMLongitude']), float(row['SMLatitude']))
    geoseries = gpd.GeoSeries([point1, point2])

    # Distance of the 2 points in degrees
    distance = geoseries[0].distance(geoseries[1])

    # return the distance which is converted to kilometers (Multiply by 111)
    return distance*111
    
    # except (TypeError, ValueError) as e:
    #     print(f"Error calculating distance: {e}")
    #     return None  # Return a default value or handle the error as needed
    
    

def complete_rupture_smtk(df: pd.DataFrame) -> pd.DataFrame:
    """
    This computes for additional parameters based on the user input. 
    This includes Width, Ztor, Z1pt0, and Z2pt5.
    """
    # df['Repi']  = df.apply(_get_distance_Repi, axis=1)
    # df['FaultWidth']  = df.apply(_compute_width, axis=1)
    # df['ztor']   = df.apply(_compute_Ztor, axis=1)
    # df['z1pt0']  = df.apply(_z1pt0_estimationAS, axis=1)
    # df['z2pt5']  = df.apply(_z2pt5_estimation, axis=1)
    
    # old Version
    df['repi']  = df.apply(_get_distance_Repi, axis=1)
    df['FaultWidth']  = df.apply(_compute_width, axis=1)
    df['Ztor']   = df.apply(_compute_Ztor, axis=1)
    df['Z1pt0']  = df.apply(vs30_z1pt0_ChiouYoungs2014, axis=1)
    df['Z2pt5']  = df.apply(vs30_z2pt5_CampbellBozorgnia2014, axis=1)

    return df


# vs30_z1pt0_AbrahamsonSilva2008