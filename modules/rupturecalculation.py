# Define a function to extract rupture data
def extract_rupture_data(row):
    hypocenter = {"lon": row['EpicenterLongitude'], "lat": row['EpicenterLatitude'], "depth": row['ModelDepth']}
    rupture_data = get_rupture_surface(row['Magnitude'], hypocenter, row['StrikeAzimuth'], row['DipAngle'], row['RakeAngle'])
    vertices = rupture_data
    print(rupture_data)
    return {
        'top_left_lon': vertices['topLeft']['lon'],
        'top_left_lat': vertices['topLeft']['lat'],
        'top_left_depth': vertices['topLeft']['depth'],
        'top_right_lon': vertices['topRight']['lon'],
        'top_right_lat': vertices['topRight']['lat'],
        'top_right_depth': vertices['topRight']['depth'],
        'bottom_left_lon': vertices['bottomLeft']['lon'],
        'bottom_left_lat': vertices['bottomLeft']['lat'],
        'bottom_left_depth': vertices['bottomLeft']['depth'],
        'bottom_right_lon': vertices['bottomRight']['lon'],
        'bottom_right_lat': vertices['bottomRight']['lat'],
        'bottom_right_depth': vertices['bottomRight']['depth'],
    }
    
#------------------------------------------------------------------------------------------------    
from datetime import datetime
import zipfile, random, logging
import math, os, numpy
import xml.dom.minidom
import xml.etree.ElementTree as ET


#: Earth radius in km.
EARTH_RADIUS = 6371.0


def _point_at(origin, horizontal_distance, vertical_increment, azimuth):
    lon, lat = numpy.radians(origin["lon"]), numpy.radians(origin["lat"])
    tc = numpy.radians(360 - azimuth)
    sin_dists = numpy.sin(horizontal_distance / EARTH_RADIUS)
    cos_dists = numpy.cos(horizontal_distance / EARTH_RADIUS)
    sin_lat = numpy.sin(lat)
    cos_lat = numpy.cos(lat)

    sin_lats = sin_lat * cos_dists + cos_lat * sin_dists * numpy.cos(tc)
    sin_lats = sin_lats.clip(-1., 1.)
    lats = numpy.degrees(numpy.arcsin(sin_lats))

    dlon = numpy.arctan2(numpy.sin(tc) * sin_dists * cos_lat,
                        cos_dists - sin_lat * sin_lats)
    lons = numpy.mod(lon - dlon + numpy.pi, 2 * numpy.pi) - numpy.pi
    lons = numpy.degrees(lons)

    deps = origin["depth"] + vertical_increment

    target = {"lon": lons, "lat": lats, "depth": deps}

    return target

def _get_rupture_length_subsurface(mag, rake):
    assert rake is None or -180 <= rake <= 180
    if rake is None:
        # their "All" case
        return 10.0 ** (-2.44 + 0.59 * mag)
    elif (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
        # strike slip
        return 10.0 ** (-2.57 + 0.62 * mag)
    elif rake > 0:
        # thrust/reverse
        return 10.0 ** (-2.42 + 0.58 * mag)
    else:
        # normal
        return 10.0 ** (-1.88 + 0.50 * mag)

def _get_rupture_width(mag, rake):
    assert rake is None or -180 <= rake <= 180
    if rake is None:
        # their "All" case
        return 10.0 ** (-1.01 + 0.32 * mag)
    elif (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
        # strike slip
        return 10.0 ** (-0.76 + 0.27 * mag)
    elif rake > 0:
        # thrust/reverse
        return 10.0 ** (-1.61 + 0.41 * mag)
    else:
        # normal
        return 10.0 ** (-1.14 + 0.35 * mag)
    
def get_rupture_surface(mag, hypocenter, strike, dip, rake):
    
    rdip = math.radians(dip)

    azimuth_right = strike
    azimuth_down = (azimuth_right + 90) % 360
    azimuth_left = (azimuth_down + 90) % 360
    azimuth_up = (azimuth_left + 90) % 360

    rup_length = _get_rupture_length_subsurface(mag, rake)
    rup_width = _get_rupture_width(mag, rake)
    # calculate the height of the rupture being projected
    # on the vertical plane:
    rup_proj_height = rup_width * math.sin(rdip)
    # and its width being projected on the horizontal one:
    rup_proj_width = rup_width * math.cos(rdip)

    # half height of the vertical component of rupture width
    # is the vertical distance between the rupture geometrical
    # center and it's upper and lower borders:
    hheight = rup_proj_height / 2.
    # calculate how much shallower the upper border of the rupture
    # is than the upper seismogenic depth:
    vshift = hheight - hypocenter["depth"]
    # if it is shallower (vshift > 0) than we need to move the rupture
    # by that value vertically.

    rupture_center = hypocenter

    if vshift > 0:
        # we need to move the rupture center to make the rupture plane
        # lie below the surface
        hshift = abs(vshift / math.tan(rdip))
        rupture_center = _point_at(
            hypocenter,
            horizontal_distance=hshift, vertical_increment=vshift,
            azimuth=azimuth_down)

    theta = math.degrees(
        math.atan((rup_proj_width / 2.) / (rup_length / 2.))
    )
    hor_dist = math.sqrt(
        (rup_length / 2.) ** 2 + (rup_proj_width / 2.) ** 2
    )

    left_top = _point_at(
        rupture_center,
        horizontal_distance=hor_dist,
        vertical_increment=-rup_proj_height / 2.,
        azimuth=(strike + 180 + theta) % 360
    )
    right_top = _point_at(
        rupture_center,
        horizontal_distance=hor_dist,
        vertical_increment=-rup_proj_height / 2.,
        azimuth=(strike - theta) % 360
    )
    left_bottom = _point_at(
        rupture_center,
        horizontal_distance=hor_dist,
        vertical_increment=rup_proj_height / 2.,
        azimuth=(strike + 180 - theta) % 360
    )
    right_bottom = _point_at(
        rupture_center,
        horizontal_distance=hor_dist,
        vertical_increment=rup_proj_height / 2.,
        azimuth=(strike + theta) % 360
    )
    rupture_plane = {"topLeft": left_top,
                    "topRight": right_top,
                    "bottomLeft": left_bottom,
                    "bottomRight": right_bottom}
    return rupture_plane
    
def get_rupture_surface_round(mag, hypocenter, strike, dip, rake):
    rupture_plane = get_rupture_surface(mag, hypocenter, strike, dip, rake)
    for corner in ["topLeft", "bottomLeft", "topRight", "bottomRight"]:
        for comp in ["lat", "lon", "depth"]:
            rupture_plane[corner][comp] = "%.5f" % rupture_plane[corner][comp]

    return rupture_plane