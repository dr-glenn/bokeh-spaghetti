#!/usr/bin/env python

'''read ISIS grid files, display with color contour lines'''

# spaghetti plot. Choose 3 values for wind or whatever field. Each value gets a unique color.
# Plot contours for all ensemble members.

# Example filename from ISIS:
# ENSEMBLE:2015123000:global_360x181:wnd_spd:isbr_lvl:08500000:00000000:fcst_et001:0240

'''
NOTES:
- How to handle user settable contour values.
Contour values are displayed in GUI.
Checkbox enables user to override default values.
  If checkbox is off, then values should be fetched from server and will be default values.
  If checkbox is on, then values are whatever is in the number text boxes.
Contour image files must have contour values in the filename.
Javascript must know the currently selected contour values in order to create hidden img filenames.
We can use session cookie to store user overrides, but ultimately must use identity.

'''

import os
import sys
from datetime import datetime
import numpy as np
import re       # regular expression parser
import traceback    # for debugging disasters

#import efsm_config as config
#import efsv_isis_fnmoc as isis  # read ISIS files, convert to object ISIS_GRID
import pyproj
    
def geo_conv(xy, coord_from=pyproj.Proj(init='epsg:3857'), coord_to=pyproj.Proj(init='epsg:4326')):
    x0,y0 = xy
    x1,y1 = pyproj.transform(coord_from, coord_to, x0, y0)
    return x1,y1

cgi_app = False     # default assumes run from command line (which is actually unlikely)

isis_params = None

ifiles = ['2016041400/ENSEMBLE:2016041400:global_360x181:wnd_spd:ht_sfc:00100000:00000000:fcst_et001:0000','2016041400/ENSEMBLE:2016041400:global_360x181:wnd_spd:ht_sfc:00100000:00000000:fcst_et002:0000']
ifiles = ['2016041400/ENSEMBLE:2016041400:global_360x181:wnd_spd:isbr_lvl:02500000:00000000:fcst_et001:0000','2016041400/ENSEMBLE:2016041400:global_360x181:wnd_spd:isbr_lvl:02500000:00000000:fcst_et002:0000',
'2016041400/ENSEMBLE:2016041400:global_360x181:wnd_spd:isbr_lvl:02500000:00000000:fcst_et003:0000',
'2016041400/ENSEMBLE:2016041400:global_360x181:wnd_spd:isbr_lvl:02500000:00000000:fcst_et004:0000',
]
contours = []

def httpSendCache(cacheFile):
    '''
    Map image files are cached by app. If the cacheFile exists, output to HTTP stream and return True.
    If cacheFile does not exist, simply return False.
    :param cacheFile: test for existence of cacheFile and output to HTTP
    :return: True if cache found, False otherwise
    '''
    if os.path.exists(cacheFile):
        logger.debug('output http: cache: %s' %(cacheFile))
        print 'Content-Type: image/png'
        # Output of "Content-Length" header is optional, but if not used, then Apache will
        # not send output until Python CGI completely exits. When "Content-Length" is present,
        # Apache will send output as soon as its buffer is filled.
        fstat = os.stat(cacheFile)
        print 'Content-Length: '+str(fstat.st_size)
        print ''    # extra newline required to terminate HTTP header block
        # TODO: maybe should also send Expires header?
        infile = open( cacheFile, 'rb' )
        while True:
            data = infile.read(4096)    # buffered read
            if not data:
                break   # no more input
            # TODO - when changed to Python 3.x, I believe this statement needs to be
            #    changed to sys.stdout.buffer.write(data) and maybe even wrap with
            #    UTF-8 encoder.
            sys.stdout.write(data)
        infile.close()
        # flush stdout - shouldn't be needed, but screen draw seems slow.
        sys.stdout.flush()
        return True
    else:
        return False

def main(argv):
    '''
    main does all the plotting
    '''
    print('main 0.0')
    from mpl_toolkits.basemap import Basemap
    print('main 0.1')
    import matplotlib.pyplot as plt
    #from matplotlib import rcParams
    #rcParams['figure.figsize'] = (800,800)
    fig = plt.figure(figsize=(8,8), dpi=100)    # this method to set figsize might not always work
    #fig = plt.gcf() # get current figure - the entire plotting area
    #fig = plt.figure(1) # get current figure - the entire plotting area
    print('main 0.2')
    # TODO: not sure why these next 3 lines are here - I think I only wanted them for command line debugging
    fig.set_frameon(True)
    fig.patch.set_visible(True)
    fig.patch.set_facecolor('green')  # just want to know what part of display is "figure"

    #fig.set_dpi(100)
    global isis_params
    global cgi_app
    global contours
    upper_limit = 1.0e8
    bSetmax = False
    lower_limit = 1.0e8
    bSetmin = False
    hasHeader = True
    # TODO: gridx, gridy should not be fixed here
    gridx = 360
    gridy = 181

    ######### Parse command line args ###########
    # there are better ways to parse args, but this will do for this simple program
    if len(argv) > 1:
        iarg = 1
        while iarg < len(argv):
            arg = argv[iarg]
            if arg[0] == '-':
                if arg == '-max':
                    bSetmax = True
                    iarg += 1
                    upper_limit = float(argv[iarg])
                elif arg == '-min':
                    bSetmin = True
                    iarg += 1
                    lower_limit = float(argv[iarg])
                elif arg == '-gridx':
                    iarg += 1
                    gridx = int(argv[iarg])
            else:
                infiles = argv[iarg:]
                break
            iarg += 1
    else:
        pass
    gridy = gridx/2 + 1     # program only works for global grid
    grid_incr = 1.0/float(gridx/360)
    
    # TODO: do something here to possibly display multiple panels if more than one infiles
    #files = allFiles(infiles[0])    # from single filename, generate list of all enesmeble member files
    files = ifiles
    #files = []
    #files.append(ifiles[0])

    # if all graphs should be same range
    # 10m winds: 5,10,20
    # 2m air temp: -10,0,10
    # 250mb winds: 40,60,80
    zbreaks = (40,60,80)
    zcolors = ('red','green','blue')     # too garish
    zcolors = ('maroon','seagreen','steelblue')   # steelblue is too dark
    zcolors = ('maroon','seagreen','lightskyblue')

    ######### Open grid file and plot it ###########
    iplot = 0    # current plot number, used to flag initial setup
    
    mpl99 = False       # special code may be needed for maplotlib version 0.99
    n_cols = n_rows = 1

    # 'ENSEMBLE-2016041400-global_360x181-wnd_spd-ht_sfc-00100000-00000000-fcst_et001-0000'
    for infile in files:
        if not os.path.exists(infile):
            print('file not found: %s' %(infile))
            continue
        # Get field name to put into map title
        file_fields = infile.split(':')
        fld_title = file_fields[3]
        print 'File = %s'%(infile)
        ####################################################################################
        # Read the map grid and sanitize it by setting max and NaN values to something sane
        try:
            #logger.debug('call ISIS_GRID.iread: %d' %(iplot+1))
            fd = open(infile, "rb")
            # read the 512 byte ISIS header
            header = np.fromfile(fd, dtype='<f', count=128)
            grid = np.fromfile(file=fd, dtype='<f')   # little-Endian 4 byte float

            map_data = np.reshape(grid,(181,360),order='C')

            valid_data = True
            dtg = file_fields[1]
            tau = int(file_fields[8])
            level = int(file_fields[6][:4])
            fld_title = file_fields[3]
        except:
            valid_data = False
            e = sys.exc_info()
            print('ex0=%s' %(e[0]))
            print('ex1=%s' %(e[1]))
            for line in traceback.extract_tb(e[2]):
                print line
        finally:
            pass

        if not valid_data:
            continue        # bad file, skip to next

        # Unique PNG filename for caching
        str_contours = '%d' %(int(zbreaks[0]))
        for z in zbreaks[1:]:
            str_contours += '-%d' %(int(z))

        # If we get here, then the image was not found in cache
        # TODO: remaing code in main() should be a separate function and should be called by httpSendCache
        if iplot == 0:
            # Setting up plot figure first time through
            # my_axes is a numpy array (n_rows x n_cols) - neat!
            plt.close('all')
            try:
                # squeeze=False guarantees that returned object 'my_axes' is always a 2-D numpy array.
                # When squeeze=True (default), it will return a 2-D array only if nrows>1 and ncols>1.
                # TODO: spaghetti not expected to have subplots
                # NOTE: subplots not valid for mpl 0.99
                fig,my_axes = plt.subplots(nrows=n_rows, ncols=n_cols, squeeze=False)
                #plt.subplots_adjust(bottom=0.10,left=0.10)
            except:
                pass
            if True:    # instead of using subplots functionality, just pretend this is mpl 0.99 and do it the old-fashioned way
                # old matplotlib 0.99 does not have method "subplots"
                mpl99 = True
                fig = None
                #my_axes = np.array((n_rows,n_cols))
                #my_axes = np.array(dtype=object,ndmin=2)
                # TODO: should make my_axes into a numpy 2-D array, but how?
                my_axes = []
                for r in range(n_rows):
                    my_axes.append([])
                    for c in range(n_cols):
                        my_axes[r].append(plt.subplot(1,1,(c+1)+(r*n_rows)))
                    
        nmax = 0
        nmin = 0
        num_nan = 0
    
        try:
            # TODO: convert degree K to F
            #map_data = (9.0/5.0)*(map_data - 273.15) + 32.0
            # set upper_limit on the grid
            maxval = np.nanmax(map_data)
            minval = np.nanmin(map_data)
            print 'min=%f, max=%f' %(minval,maxval)
            # Count the NaN values and establish min and max for grid
            #logger.debug('shape(map_data)=%s' %str(map_data.shape))
            # TODO: There are numpy methods that could eliminate the loop here
            for i in range(map_data.shape[0]):
                for j in range(map_data.shape[1]):
                    if np.isnan(map_data[i,j]):
                        num_nan += 1
                    if bSetmax:
                        # if user specified upper_limit, then set values greater than limit
                        # to a value greater than upper_limit.
                        # And if value is NaN, set it to an even larger value so we can see
                        # it easily on the map.
                        if map_data[i,j] > upper_limit:
                            map_data[i,j] = upper_limit * 1.2
                            nmax += 1
                        elif np.isnan(map_data[i,j]):
                            map_data[i,j] = upper_limit * 1.5
                    if bSetmin:
                        # if user specified lower_limit, then set values less than limit
                        # to a value less than lower_limit.
                        if map_data[i,j] < lower_limit:
                            # TODO: need a better way to alter map_data at lower_limit
                            map_data[i,j] = lower_limit - 1.0
                            nmin += 1
            #print len(map_data)
            if not bSetmax:
                upper_limit = maxval
            if not bSetmin:
                lower_limit = minval
            #logger.debug('main: finished create map_data: valid=%s' %str(map_data.valid))
            #logger.debug('-- min=%f, max=%f, num_nan=%d' %(minval,maxval,num_nan))
        except:
            e = sys.exc_info()
            #logger.error('ex0=%s' %(e[0]))
            #logger.error('ex1=%s' %(e[1]))
            for line in traceback.extract_tb(e[2]):
                print(line)

        ####################################################################################
        # Plot the map_data
    
        try:
            # Which subplot?
            irow = icol = 0
            if mpl99 == False:
                sub_axis = my_axes[irow,icol]
            else:
                sub_axis = my_axes[irow][icol]
                plt.subplot(n_rows,n_cols,(icol+1)+(irow*n_cols)) # make sure to select current subplot
            # TODO: next section should be only done once, but what if multiplot display?
            lonmax = 180
            if iplot == 0:
                # GDN: I tried projection='merc' and found that (a) the map was skinny and vertically stretched,
                #   and (b) the projection coordinate Y values went up to above 66,000,000, even though
                #   total earth should not be > 40,000,000.
                if lonmax == 360:
                    # Must have center longitude = 180, because ISIS grids go from 0 to 360
                    map_ax = Basemap(ax=sub_axis, projection='mill', resolution='c', lon_0=180)
                    lons0 = np.arange(0.0,359.9,grid_incr)
                else:
                    # Bokeh and other maps default to range of -180 to 180.
                    # I have not figured out how to change Bokeh map range to 0 to 360.
                    map_ax = Basemap(ax=sub_axis, projection='mill', resolution='c', lon_0=0)
                    # Setting lons0 to this range does not work, because we are falsely stating that
                    # this is the range for ISIS grid. If I could rearrange map_data, it should work.
                    lons0 = np.arange(-180.0,179.9,grid_incr)
                print 'basemap proj: '+map_ax.proj4string
                # plotting grid data on top of basemap
                lats0 = np.arange(-90.0,90.1,grid_incr)
                lons,lats = np.meshgrid(lons0, lats0)   # create 2D grid of all (lon,lat) points
                #print '%d,%d' %(len(lons),len(lats))
                x,y = map_ax(lons, lats)    # convert (lon,lat) to x,y - projection coordinates
                print('contour shapes: map_data=%s, x=%s, y=%s' %(str(map_data.shape),str(x.shape),str(y.shape)))
                print 'x range: %f, %f' %(np.amin(x), np.amax(x))
                print 'y range: %f, %f' %(np.amin(y), np.amax(y))
            if lonmax == 360:
                pass
            else:
                # rearrange map_data: 180 to 360 and then 0 to 180
                def _swap180(xy):
                    xy1 = np.zeros((181,360))
                    xy1[:,range(0,180)] = xy[:,range(180,360)]
                    xy1[:,range(180,360)] = xy[:,range(0,180)]
                    return xy1
                    
                map_data = _swap180(map_data)
            if zbreaks:
                css = sub_axis.contour(x,y,map_data,levels=zbreaks,colors=zcolors)
                #cb = plt.colorbar(cs, shrink=0.8, extend='both')
            else:
                css = sub_axis.contourf(x,y,map_data)

            
            # Get the actual contours so that we can use them somewhere else
            if iplot == 0:
                file_stat = 'w'
            else:
                file_stat = 'a'
            path_file = open('cntr_paths.out',file_stat)
            ll_file = open('cntr_ll.out',file_stat)
            cs_out = css.collections
            print 'cs_out: %s' %(str(cs_out))
            cs_idx = 0
            xmin = xmax = ymin = ymax = 20000000.0    # 0 longitude, 0 latitude
            offx = 20000000.0
            offy = 29400000.0 / 2.0   # 294 is the y-range from Basemap

            for cs in cs_out:
                print '-- cs: %s' %(str(cs))
                #print '---- cs.get_paths: %s' %(str(cs.get_paths()))
                path_idx = 0
                for path in cs.get_paths():
                    xcoords = path.vertices.transpose()[0]
                    ycoords = path.vertices.transpose()[1]
                    path_file.write('[cs=%d, path=%d]\n' %(cs_idx,path_idx))
                    ll_file.write('[cs=%d, path=%d]\n' %(cs_idx,path_idx))
                    for xx,yy in zip(xcoords,ycoords):
                        xmin = min(xmin,xx)
                        xmax = max(xmax,xx)
                        ymin = min(ymin,yy)
                        ymax = max(ymax,yy)
                        path_file.write('  %f, %f\n' %(xx,yy))
                        #lon,lat = geo_conv((xx,yy-offy))   # this gives wrong loongitude
                        lon,lat = map_ax(xx,yy,inverse=True) # this is correct
                        ll_file.write('  %f, %f\n' %(lon,lat))
                    path_idx += 1
                cs_idx += 1
            path_file.close()
            ll_file.close()
            
            offx = 20000000.0
            offy = 29400000.0 / 2.0
            print 'x: %.2f, %.2f, y: %.2f, %.2f' %(xmin,xmax,ymin,ymax)
            print 'x: %.2f, %.2f, y: %.2f, %.2f' %(xmin-offx,xmax-offx,ymin-offy,ymax-offy)
            '''
            # get_paths returns a list of 'Path'. Each 'Path' contains a tuple of ('array',None).
            # Each 'array' contains a single 'list'; the 'list' contains lists of [X,Y].
            # GDN: I suppose each 'array' could contain more than one list, and maybe a 'Path'
            # could contain more than one 'array'.
            paths = output.get_paths()[i] # for the various contours

            # the x,y coordinates can then be accessed as
            xcoords = paths.vertices.transpose()[0]
            ycoords = paths.vertices.transpose()[1]
            '''

            if iplot == 0:
                # plot coastlines, draw label meridians and parallels.
                map_ax.drawcoastlines()
                map_ax.drawparallels(np.arange(-90,90,10),labels=[1,0,0,0])
                #map_ax.drawmeridians(np.arange(map_ax.lonmin,map_ax.lonmax+30,60),labels=[0,0,0,1])
                map_ax.drawmeridians(np.arange(lons0[0],lons0[-1]+30,60),labels=[0,0,0,1])
                # fill continents 'coral' (with zorder=0), color wet areas 'aqua'
                map_ax.drawmapboundary(fill_color='grey')
                #map_ax.fillcontinents(color='coral',lake_color='aqua')
                #plt.colorbar(cs,ax=sub_axis,shrink=0.8)
                #plt.clabel(cs, inline=1, fontsize=10)
                labels = ['%.1f' %(z) for z in zbreaks]
                for i in range(len(labels)):
                    css.collections[i].set_label(labels[i])
                plt.legend()

                sub_axis.set_title('%s, level=%.1f, tau=%03d, dtg=%s\nmin/max = %g/%g' \
                                   %(fld_title,level,tau,dtg,minval,maxval))

        except:
            e = sys.exc_info()
            for line in traceback.extract_tb(e[2]):
                print(line)
        iplot += 1
        #logger.debug('iplot=%d: memory use=%d'%(iplot,resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    
    if iplot > 0:
        # Show map in command line viewer
        #logger.debug('output show')
        #plt.tight_layout(pad=0.4)
        print('about to call show')
        plt.show()

if __name__ == '__main__':
    '''
    If run from command line, we are just testing and will get a GTK window for display.
    argv is expected to be a filename. We need to break it apart, so that it seems like it
    was called from web.
    When run from cmd line, be sure to give full path to data file. Also you cannot run from
    web server dir because it attempts to write to log files that are owned by apache.
    So for cmd line testing you should run from dirs that you own.
    '''
    #import Tkinter
    import matplotlib
    #matplotlib.use('Qt4Agg')
    matplotlib.use('TKAgg')
    #matplotlib.use('GTKAgg')
    #matplotlib.use('Agg')
    try:
        print('start')
        main(sys.argv)
    except:
        e = sys.exc_info()
        print('ex0=%s' %(e[0]))
        print('ex1=%s' %(e[1]))
        for line in traceback.extract_tb(e[2]):
            print line
