# Plot a country outline world map with Bokeh.
# The Bokeh map appears has linear axes.
#   Will need to use bokeh.models.map_plots to get true Mercator y-axis.
#   But bokeh.models.map_plots uses Google Maps API and these maps are very slow to load and the Bokeh
#   zoom controls hardly work at all.
# https://groups.google.com/a/continuum.io/forum/#!topic/bokeh/cJ7tzsm1aLc
# Apparently with bokeh.plotting.figure, the default toolbar has Pan, BoxZoom, WheelZoom, Save and Reset.
# BoxZoomTool works with XY figure and WheelZoomTool is smooth - not so with map_plot objects.

# Here is some more documentation on Bokeh tools: http://bokeh.pydata.org/en/latest/docs/user_guide/tools.html#userguide-tools-clicktap

# Nice - side-by-side plots: zoom inside one plot, the other one remains full scale, but shows box of zoom region
# http://bokeh.pydata.org/en/latest/docs/user_guide/interaction/callbacks.html#userguide-interaction-callbacks

# Instructions for outputting both JS script tag and DIV tag for the plot:
# http://bokeh.pydata.org/en/latest/docs/user_guide/embed.html#userguide-embed

# GDN - this program modifies cntr_bokeh_2.py to output script and div tags instead of complete HTML.

from __future__ import print_function 

import requests 
import numpy as np 

import bokeh 
print(bokeh.__version__) 

from bokeh.io import output_file, show 
from bokeh.models import (GeoJSONDataSource, CustomJS,Circle, Line, Rect, DataRange1d, HoverTool, PanTool, WheelZoomTool, BoxSelectTool, BoxZoomTool, SaveTool, CrosshairTool, TapTool)
from bokeh.plotting import figure, ColumnDataSource
from bokeh.embed import components
import pyproj

# Define map center (not origin) and extent
lon0 = 0.0
lat0 = 0.0
lon_extent = 360
lat_extent = 160    # -80 to +80
    
def geo_conv(xy, coord_from=pyproj.Proj(init='epsg:4326'), coord_to=pyproj.Proj(init='epsg:3857')):
    x0,y0 = xy
    x1,y1 = pyproj.transform(coord_from, coord_to, x0, y0)
    return x1,y1
    
# geo_inv gives correct longitude values, because this guy says so:
# http://gis.stackexchange.com/questions/169171/why-do-pyproj-and-matplotlib-basemap-produce-different-results-when-they-convert
# base_ll_proj is obtained from m=Basemap(projection='mill',lon_0=0)
# base_ll_proj = m.proj4string   
base_ll_proj =  '+lon_0=0 +y_0=14675034.4036 +R=6370997.0 +proj=mill +x_0=-0.0 +units=m'
map_trans = pyproj.Proj(base_ll_proj)   # transform from LL to XY
def geo_inv(xy):
    x0,y0 = xy
    x1,y1 = map_trans(x0, y0, inverse=True)
    return x1,y1

# GDN: coordinates are longitude, latitude
url = 'https://raw.githubusercontent.com/johan/world.geo.json/master/countries.geo.json' 

r = requests.get(url) 
geo_json_data = r.json() 

# Read contours created by matplotlib and plot onto Bokeh 'fig'
def plot_contours(fig):
    new_cntr = False
    xmin = xmax = ymin = ymax = 10000000.0
    lonmin = lonmax = latmin = latmax = 0.0
    cntr = []
    # colors are defined by SVG names: http://www.december.com/html/spec/colorsvg.html
    colors = ['red', 'green', 'slateblue']
    line_color = 'black'
    with open('cntr_paths.out','r') as fd:
        for line in fd:
            line = line.strip()
            if line.startswith('['):
                if len(cntr) > 0:
                    # plot it
                    x,y = zip(*cntr)
                    fig.line(x,y,line_color=line_color)
                    #print 'plotted a line'
                ii = line.find('cs=')+3
                jj = line.find(',')
                cntr_idx = int(line[ii:jj])
                line_color = colors[cntr_idx % len(colors)]
                new_cntr = True
                cntr = []
            else:
                #xy = [float(num)-20000000 for num in line.split(',')]
                xy = [float(num) for num in line.split(',')]
                xmin = min(xmin,xy[0])
                xmax = max(xmax,xy[0])
                ymin = min(ymin,xy[1])
                ymax = max(ymax,xy[1])
                xy[0] -= 20000000.0
                #xy[1] -= (29400000.0/2.0)   # 294 is the y-range from Basemap
                # yrange of contours: 4852325.85 to 24748119.58
                #lonlat = geo_conv(xy,coord_from=pyproj.Proj(init='epsg:3857'), coord_to=pyproj.Proj(init='epsg:4326'))
                lonlat = geo_inv(xy)
                lonmin = min(lonmin,lonlat[0])
                lonmax = max(lonmax,lonlat[0])
                latmin = min(latmin,lonlat[1])
                latmax = max(latmax,lonlat[1])
                cntr.append(lonlat)
    print('x:   %.2f, %.2f, y:   %.2f, %.2f' %(xmin,xmax,ymin,ymax))
    print('lon: %.2f, %.2f, lat: %.2f, %.2f' %(lonmin,lonmax,latmin,latmax))

# Convert world map features into long-lat coordinate arrays
def get_coordinates(features):
    # GDN: why +1 for depth?
    # GDN: just what does this function do?
    depth = lambda L: isinstance(L, list) and max(map(depth, L))+1 
    # GDN: xs and ys are list of lists
    xs = [] 
    ys = [] 
    for feature in features: 
        coords = feature['geometry']['coordinates'] 
        nbdims = depth(coords)    # GDN: num border dims
        # one border 
        if nbdims == 3: 
            # GDN: feature has a single polygon border
            pts = np.array(coords[0], 'f') 
            xs.append(pts[:, 0]) 
            ys.append(pts[:, 1]) 
        # several borders 
        else: 
            # GDN: feature is multiple polygons
            for shape in coords: 
                pts = np.array(shape[0], 'f') 
                xs.append(pts[:, 0]) 
                ys.append(pts[:, 1]) 
    return xs, ys 
   
# callback to process muose movement so we can get long-lat
s = ColumnDataSource(data = dict(x=[0,0.5,1],y=[0,0.5,1])) # points of the line
hover_callback = CustomJS(code="""
        var geometry = cb_data['geometry'];
        var y_data = geometry.y; // current mouse y position in plot coordinates
        var x_data = geometry.x; // current mouse x position in plot coordinates
        console.log("(x,y)=" + x_data+","+y_data); //monitors values in Javascript console
    """)
#hover_tool = HoverTool(callback=hover_callback)
# The callback is active everywhere on the canvas. The tooltips only appear
# when the HoverTool is above a glyph, such as country or ensemble field contour lines.
'''
hover_tool = HoverTool(callback=callback,
            tooltips=[
                ("index", "$index"),
                ("data (using $) (x,y)", "($x, $y)"),
                ("data (using @) (x,y)", "(@x, @y)"),
                ("canvas (x,y)", "($sx, $sy)")
                ])
'''
base_callback = CustomJS(code="""
        console.log("tap= " + JSON.stringify(cb_data)); //monitors values in Javascript console
    """)
tap_callback = CustomJS(code="""
        var geometry = cb_data['geometries'][0];
        var x = geometry['x'];
        var y = geometry['y'];
        console.log("tap: (x,y)=" + x+","+y); //monitors values in Javascript console
    """)
tap_tool = TapTool(behavior='inspect',callback=tap_callback)

# This is an XY figure with linear scales.
#p = figure(plot_width=1000, plot_height=600, x_range=(-180,180), y_range=(-80,80),tools= [hover_tool,                 "crosshair,box_zoom,wheel_zoom,pan,reset"]) 
p = figure(plot_width=1000, plot_height=600, x_range=(-lon_extent/2,lon_extent/2), y_range=(-lat_extent/2,lat_extent/2)) 
#p.add_tools(hover_tool, CrosshairTool(), tap_tool)
p.add_tools(tap_tool)

# customize the BoxZoomTool appearance
zoom_overlay = p.select_one(BoxZoomTool).overlay
zoom_overlay.line_color = "olive"
zoom_overlay.line_width = 8
zoom_overlay.line_dash = "solid"
zoom_overlay.fill_color = None

# GDN: cannot use next line because geojson does not have subpackage 'utils'
#country_xs, country_ys = zip(*geojson.utils.coords(geo_json_data['features']))
# GDN: so instead someone created 'get_coordinates'
xs, ys = get_coordinates(geo_json_data['features']) 
p.patches(xs, ys, fill_color="#F1EEF6", fill_alpha=0.7, line_width=2)
plot_contours(p)

def draw_boxes(plt):
    # define two center points
    boxes = ColumnDataSource(
        data=dict(
            x=[-125.0, 0.0],
            y=[0.0, 0.0],
        )
    )
    # Create rectangle object of specified width, height, color
    rect = Rect(x="x", y="y", width=60, height=20, fill_color="blue", fill_alpha=0.5, line_color=None)
    # Draw rectangles at the centers specified in 'boxes'
    plt.add_glyph(boxes, rect)

    # Now create another rectangle, this time bypass the Rect class
    plt.rect([0.0],[60.0],width=40,height=20,color="green",alpha=0.5)

#draw_boxes(p)

# Cover entire map with a box, but make it transparent.
# We need this box so that TapTool will register a click everywhere on the map.
p.rect([lon0,],[lat0,],width=lon_extent,height=lat_extent,color="LightGray",alpha=0.2)

if False:
    # this code outputs complete HTML page
    output_file("cntr_bokeh_2.html") 
    show(p) 
else:
    from string import Template
    script,div = components(p)
    page = '''
<!DOCTYPE html>
<html lang="en">
    <head>
        <meta charset="utf-8">
        <title>Bokeh Scatter Plots</title>

        <link rel="stylesheet" href="http://cdn.pydata.org/bokeh/release/bokeh-0.12.3.min.css" type="text/css" />
        <script type="text/javascript" src="http://cdn.pydata.org/bokeh/release/bokeh-0.12.3.min.js"></script>
        <link rel="stylesheet" href="http://cdn.pydata.org/bokeh/release/bokeh-widgets-0.12.3.min.css" type="text/css">
        <script type="text/javascript" src="http://cdn.pydata.org/bokeh/release/bokeh-widgets-0.12.3.min.js"></script>

        <!-- COPY/PASTE SCRIPT HERE -->
        $script

    </head>
    <body>
        <h1>Here's Hoping!</h1>
        <!-- INSERT DIVS HERE -->
        $divs
    </body>
</html>
    '''
    pageTemplate = Template(page)
    outStr = pageTemplate.safe_substitute(script=script,divs=div)
    #print(outStr)
    fl = open('cntr_21.html','w')
    fl.write(outStr)
    fl.close()
