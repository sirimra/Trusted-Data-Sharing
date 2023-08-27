import numpy as np
from shapely.ops import nearest_points
import shapefile
import pandas as pd
import math
import geopy.distance
import seaborn as sns
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from tqdm import tqdm
import matplotlib.pyplot as plt
from minerals import pointClass


# those magic values comes from this particular paper
#Porwal, A., Gonzalez-Alvarez, I., Markwitz, V., McCuaig, T., Mamuse, A., 2010. Weights-of-evidence and logistic regression modeling of magmatic nickel sulfide prospectivity in the yilgarn craton, western australia. Ore Geology Reviews 38, 184–196.
weight_model = {'komatiite':2.672682,'ultramafic':2.371440,' mafic':0.826142,'crust':2.288389, 'Cr':0.144607, 'Cu':0.082966,'Ni':0.437637}

def class_weight(z):
    '''
    # those magic values comes also from this particular paper
    #Porwal, A., Gonzalez-Alvarez, I., Markwitz, V., McCuaig, T., Mamuse, A., 2010. Weights-of-evidence and logistic regression modeling of magmatic nickel sulfide prospectivity in the yilgarn craton, western australia. Ore Geology Reviews 38, 184–196.

    '''
    ll =[]
    for distance in z :
        if distance  == 0:
            ll.append( 1)
        elif distance <= 1:
            ll.append( 0.9)
        elif distance <=2:
            ll.append( 0.4)
        elif distance <=3:
            ll.append( 0.3)
        elif distance <=4:
            ll.append( 0.2)
        elif distance <=5:
            ll.append( 0.1)
        elif distance <=10:
            ll.append( 0.01)
        else:
            ll.append( 0)
    return ll


def calculate_closest_geology(df, df_minerals_grid):
    '''
       The geology attributes are selected according to the following paper
       This is just a exmaple and not reflect the real life mineral exploration.
       #Porwal, A., Gonzalez-Alvarez, I., Markwitz, V., McCuaig, T., Mamuse, A., 2010. Weights-of-evidence and logistic regression modeling of magmatic nickel sulfide prospectivity in the yilgarn craton, western australia. Ore Geology Reviews 38, 184–196.

    '''
    for col in ['komatiite', 'ultramafic', ' mafic']:
        df_minerals_grid[col + '_distance'] = ''

    for id, row in tqdm(df_minerals_grid.iterrows()):
        dico = {}


        for geo_str in eval(row.geology_points):
            geo = eval(geo_str)
            dico[geo[0]] = geo[1]
            dico_sorted = dict(sorted(dico.items(), key=lambda item: item[1]))

        df_minerals_grid['crust_distance'] = list(dico_sorted.values())[0]
        for col in [' mafic', 'komatiite', 'ultramafic']:
            for geo in dico_sorted.keys():
                val = df['geology'][df['geology']['gid'] == geo]['descriptn'].values[0].lower()
                # print (val)
                # 'komatiite',
                if col in val:
                    df_minerals_grid.at[id, col + '_distance'] = dico_sorted[geo]

                    break
    return df_minerals_grid

def calculate_closest_chemistry(df, df_minerals_grid):
    '''
           The geology attributes are selected according to the following paper
           This is just a exmaple and not reflect the real life mineral exploration.
           #Porwal, A., Gonzalez-Alvarez, I., Markwitz, V., McCuaig, T., Mamuse, A., 2010. Weights-of-evidence and logistic regression modeling of magmatic nickel sulfide prospectivity in the yilgarn craton, western australia. Ore Geology Reviews 38, 184–196.

        '''
    for col in ['Cr', 'Cu', 'Ni']:
        df_minerals_grid[col + '_distance'] = ''


    for id, row in tqdm(df_minerals_grid.iterrows()):
        dico, dico_sorted = {}, {}
        for geo_str in eval(row.chemistry_points):
            geo = eval(geo_str)
            dico[geo[0]] = geo[1]
            dico_sorted = dict(sorted(dico.items(), key=lambda item: item[1]))

        for col in ['Cr', 'Cu', 'Ni']:

            for geo in dico_sorted.keys():
                val = df['sampling'][df['sampling'].index == geo][col].values[0]
                if val.dtype == float and np.isnan(val) == False and val > 0:
                    df_minerals_grid.at[id, col + '_distance'] = str((geo, val, dico_sorted[geo]))

                    break
    return df_minerals_grid

def calculate_closerchem( df_sampling, df_min, list_sampl_into_tenement):
    dico_sorted = {}
    for id, row in tqdm(df_min.iterrows()):
        dico = {}
        dico[9999] = 9999
        #print ('eval(row.chemistry_points', eval(row.chemistry_points), type(eval(row.chemistry_points)))
        for geo_str in list(eval(row.chemistry_points)):
            geo = eval(geo_str)
            if not geo[0] in list_sampl_into_tenement:
                dico[geo[0]] = geo[1]

        dico_sorted = dict(sorted(dico.items(), key=lambda item: item[1]))

        for col in ['Cr', 'Cu', 'Ni']:

            for geo in dico_sorted.keys():
                if geo == 9999:
                    df_min.at[id, col + '_distance'] = (geo, 9999, dico_sorted[geo])
                else:
                    val = df_sampling[df_sampling.index == geo][col].values[0]
                    if val.dtype == float and np.isnan(val) == False and val > 0:
                        df_min.at[id, col + '_distance'] = (geo, val, dico_sorted[geo])

                        break

    return df_min




def display_graph(axi, df, tenement_id, tit,innercells=[],  lls=[5], conf_loss_count=1, conf_loss_tnmt_count=0):
    tenid = tenement_id.replace(' ', '')
    # display (df.head())
    meanerr={}
    for l in lls:
        klargest_errors = df[~df['tenements'].str.contains(tenement_id)]['abs_error_' + tenid].nlargest(l).values

        klargest_errors_index = df[(~df['tenements'].str.contains(tenement_id)) & (df['abs_error_' + tenid] > 0.00000001)][
        'abs_error_' + tenid].nlargest(l).index.tolist()

        meanerr[l] = str(round(np.mean(klargest_errors), 3))


    df['proba' + tenid] = df['proba_prospectivity_' + tenid] * 2

    a = []
    text_a = []
    for y in sorted(np.unique(df.coords_y.tolist()), reverse=True):
        ll = []
        text_ll = []
        for x in sorted(np.unique(df.coords_x.tolist())):
            ll.append(df[(df.coords_x == x) & (df.coords_y == y)]['proba' + tenid].values[0])
            idd = df[(df.coords_x == x) & (df.coords_y == y)].index.values[0]
            ten = (df[(df.coords_x == x) & (df.coords_y == y)]['tenements'].values[0])
            # print (ten, tenid)
            if len(innercells)>0:
                if idd in innercells:
                    text_ll.append('X')
                else:
                    # text_ll.append(idd)
                    #if idd in klargest_errors_index:
                    #    text_ll.append(idd)
                    #else:
                        #
                    text_ll.append(' ')
            else:
                if tenement_id in ten:
                    text_ll.append('X')
                else:
                    # text_ll.append(idd)
                    #if idd in klargest_errors_index:
                    #    text_ll.append(idd)
                    #else:
                        #
                    text_ll.append(' ')
        a.append(ll)
        text_a.append(text_ll)

    # print (text_a)
    if len(df)>1000:
        step =5
    else:
        step=4

    q = sns.heatmap(a, annot=text_a, fmt="", linewidth=0.5, center=1, ax=axi, cmap='bwr')
    xlab = sorted([round(x, 2) for x in np.unique(df.coords_x.tolist())])
    lx = np.arange(len(xlab))

    axi.tick_params(axis='x', labelsize=14, rotation=30)

    axi.set(xticks=lx[0::step], xticklabels=xlab[0::step])#)
    axi.set_xlabel('Longitude',fontsize=14)
    axi.set_ylabel('Latitude',fontsize=14)

    ylab = sorted([round(x, 2) for x in np.unique(df.coords_y.tolist())])
    ly = np.arange(len(ylab))

    axi.tick_params(axis='y', labelsize=14, rotation=0)

    axi.set(yticks=ly[0::step], yticklabels=ylab[0::step])  # , rotation=30)





    df_counts = pd.DataFrame(
        df[~df['tenements'].str.contains(tenement_id)]['proba_prospectivity_' + tenid].round(decimals=1).value_counts())
    conf_loss_weight = (df_counts.index * df_counts['proba_prospectivity_' + tenid]).sum() / len(df)

    
    if len(innercells) == 0 :
        conf_loss_tnmt_weight = 0
         
    else:

        allcels = (df[df['tenements'].str.contains(tenement_id)].index.tolist())
        if len(allcels) >0:
            df_counts_tnmt_border = pd.DataFrame(
                df.iloc[list(set(allcels) - set(innercells))]['proba_prospectivity_' + tenid].round(
                    decimals=1).value_counts())
            df_counts_tnmt_border
            conf_loss_tot_border = (
                        df_counts_tnmt_border.index * df_counts_tnmt_border['proba_prospectivity_' + tenid]).sum()
            conf_loss_tnmt_weight = conf_loss_tot_border / len(allcels)
            
        else:
            conf_loss_tnmt_weight=0


    dico_metrics = {}
    if tenement_id == 'orig':
        dico_metrics['conf_loss_count'] = conf_loss_count
        dico_metrics['conf_loss_weight']=conf_loss_weight
    else:
        dico_metrics['conf_loss_count'] = conf_loss_count
        dico_metrics['conf_loss_weight'] = conf_loss_weight
        dico_metrics['conf_loss_tnmt_count'] = conf_loss_tnmt_count
        dico_metrics['conf_loss_tnmt_weight'] = conf_loss_tnmt_weight
        dico_metrics['meanerr']=meanerr
    if tenement_id=='orig':
        axi.set_title( tit +' \nConf_loss_count =' + \
        str(round(conf_loss_count, 3)) \
        + '   Conf_loss_weight=' + str(round(conf_loss_weight, 3)),\
         fontsize=22)
    else:
        tit1 = ''
        for l in lls:
            tit1 += "loss_MAE-"+str(l)+"=" + meanerr[l]+ " "
        axi.set_title(\
            tit + "\n" +tit1 +\
            ' \nConf_loss_count =' + str(round(conf_loss_count, 3)) \
            + '   Conf_loss_weight=' + str(round(conf_loss_weight, 3))
            + ' \n Conf_loss_tnmt_count =' + str(round(conf_loss_tnmt_count, 3)) \
            + '   Conf_loss_tnmt_weight=' + str(round(conf_loss_tnmt_weight, 3)),  fontsize=20)


    q.collections[0].colorbar.set_ticks([0, 0.5, 1, 1.5, 2])
    q.collections[0].colorbar.set_ticklabels(['0', '0.25', '0.5', '0.75', '1'])
    q.collections[0].colorbar.ax.tick_params(labelsize=20)
    axi.collections[0].colorbar.set_label("Prospectivity probability", size=20)
    return dico_metrics

def getboder_cells(df_orig, tenement_id):
    from shapely.ops import cascaded_union
    list_cells_of = df_orig[df_orig.tenements.str.contains(tenement_id)].index.tolist()
    polygon_cell = {}
    for id, cell in (df_orig.iloc[df_orig[df_orig.tenements.str.contains(tenement_id)].index.tolist()].iterrows()):
        try:
            cell_coordinates = eval(cell.coords)
        except:
            cell_coordinates = (cell.coords)
        polygon_cell[id] = Polygon([cell_coordinates[0], cell_coordinates[1], cell_coordinates[2], cell_coordinates[3]])

    polygons = [polygon_cell[x] for x in polygon_cell]

    u = cascaded_union(polygons)

    try:
        
        cells_in_border = []
        if not u.type == 'MultiPolygon':
            borders = u.exterior.xy


            for id, cell in (df_orig.iterrows()):
                legend = False

                cell_coordinates = eval(cell.coords)
                polygon_cell[id] = Polygon([cell_coordinates[0], cell_coordinates[1], cell_coordinates[2], cell_coordinates[3]])

                for i in range(len(borders[0])):

                    if polygon_cell[id].intersects(Point(borders[0][i], borders[1][i])) == True:
                        cells_in_border.append(id)

            cellsborders = set(intersection(cells_in_border, list_cells_of))
            innercells = set(list_cells_of) - cellsborders
        else:
            cellsborders =set()
            innercells=set()
            for u2 in u.geoms:
                borders = u2.exterior.xy

                for id, cell in (df_orig.iterrows()):
                    legend = False

                    cell_coordinates = eval(cell.coords)
                    polygon_cell[id] = Polygon([cell_coordinates[0], cell_coordinates[1], cell_coordinates[2], cell_coordinates[3]])

                    for i in range(len(borders[0])):

                        if polygon_cell[id].intersects(Point(borders[0][i], borders[1][i])) == True:
                            cells_in_border.append(id)

                cellsborders0 = set(intersection(cells_in_border, list_cells_of))
                innercells0 = set(list_cells_of) - cellsborders0
            cellsborders.update(cellsborders0)
            innercells.update(innercells0)
        return list_cells_of, list(cellsborders), list(innercells)
    except:
        print(tenement_id, 'no object has no attribute exterior')
        return list_cells_of, [], []



def connect_tenement_tosample(df):
    df_sampling_tenement = pd.DataFrame()

    # Iterate through each row in the 'tenement_large' column of the input DataFrame
    for id_ten, ten in tqdm(df['tenement_large'].iterrows()):
        # Create a polygon from the coordinates of the tenement
        polygon_ten = Polygon(ten.coords)

        # Iterate through each row in the 'sampling' column of the input DataFrame
        for id, samp in (df['sampling'].iterrows()):

            # Create a point from the coordinates of the sampling point
            point = Point(samp.coords)
            # Check if the sampling point is contained within the tenement polygon
            if polygon_ten.contains(point):
                # Create a new row for the DataFrame with the tenement and sample IDs
                df_new_row = pd.DataFrame({'TENID': [ten.TENID], 'sample': [id]})
                # Concatenate the new row to the results DataFrame
                df_sampling_tenement = pd.concat([df_sampling_tenement, df_new_row])


    df_sampling_tenement['sample'] = df_sampling_tenement['sample'].astype(int)

    df_sampling_tenement.reset_index(drop=True, inplace=True)
    # Return the results DataFrame
    return df_sampling_tenement



def create_grid_dataframe(df, side_len):
    # Add new columns 'coords_x' and 'coords_y' to the 'sampling' DataFrame
    df['sampling']['coords_x'] = [x[0] for x in df['sampling'].coords]
    df['sampling']['coords_y'] = [x[1] for x in df['sampling'].coords]

    # Calculate minimum and maximum x and y coordinates
    coords_x_min = min(df['sampling']['coords_x'])
    coords_x_min_round = round(min(df['sampling']['coords_x']), 2) - .01
    coords_x_max = max(df['sampling']['coords_x'])
    coords_x_max_round = round(max(df['sampling']['coords_x']), 2)
    coords_y_min = min(df['sampling']['coords_y'])
    coords_y_min_round = round(min(df['sampling']['coords_y']), 2) - .01
    coords_y_max = max(df['sampling']['coords_y'])
    coords_y_max_round = round(max(df['sampling']['coords_y']), 2)

    df_minerals_grid = pd.DataFrame(columns=['coords', 'coords_x', 'coords_y', 'chemistry_points', 'geology_points'])
    # Loop through x and y coordinates to create grid cells
    for coord_x in np.arange(coords_x_min_round, coords_x_max_round, (coords_x_max_round - coords_x_min_round) / side_len):


        for coord_y in np.arange(coords_y_min_round, coords_y_max_round,
                                 (coords_y_max_round - coords_y_min_round) / side_len):
            # Calculate the coordinates of the grid cell
            coords_x1 = coord_x
            coords_x2 = coord_x + (coords_x_max_round - coords_x_min_round) / side_len
            coords_y1 = coord_y
            coords_y2 = coord_y + (coords_y_max_round - coords_y_min_round) / side_len

            # Create a new row for the grid DataFrame
            df_new_row = pd.DataFrame({'coords': [[(coords_x1, coords_y1),(coords_x2, coords_y1),
                                                                   (coords_x2, coords_y2),(coords_x1, coords_y2)]],
                                                        'coords_x': [coords_x1 + (coords_x2 - coords_x1) / 2],
                                                        'coords_y': [coords_y1 + (coords_y2 - coords_y1) / 2]})
            df_minerals_grid = pd.concat([df_minerals_grid, df_new_row])

    # Sort the grid DataFrame
    df_minerals_grid.sort_values(by=['coords_x','coords_y'] )
    df_minerals_grid.reset_index(inplace=True, drop=True)
    # Perform multiple calculations and connections to tenements and togeology_and_chemistry_points
    df_minerals_grid = connect_grid_totenement(df,df_minerals_grid)
    df_minerals_grid.reset_index(inplace=True, drop=True)
    df_minerals_grid = connect_grid_togeology_and_chemistry_points(df, df_minerals_grid)
    df_minerals_grid.reset_index(inplace=True, drop=True)
    df_minerals_grid = calculate_closest_geology(df, df_minerals_grid)
    df_minerals_grid.reset_index(inplace=True, drop=True)
    df_minerals_grid = calculate_closest_chemistry(df, df_minerals_grid)
    df_minerals_grid.reset_index(inplace=True, drop=True)
    df_minerals_grid = calculate_weight(df_minerals_grid)

    return df_minerals_grid

def connect_grid_totenement( df,df_minerals_grid):
    df_minerals_grid['tenements'] = ''

    for id, cell in tqdm(df_minerals_grid.iterrows()):

        polygon_cell = Polygon((cell.coords))
        list_tens = []

        for id_ten, ten in df['tenement_large'].iterrows():

            # Create a polygon from the tenement's coordinates
            polygon_ten = Polygon(ten.coords)

            # Check if the cell and tenement polygons intersect
            if polygon_ten.intersects(polygon_cell) == True or polygon_cell.intersects(polygon_ten):
                list_tens.append(ten.TENID)

            # Update the 'tenements' column in the grid DataFrame
            df_minerals_grid.at[id, 'tenements'] = str(list_tens)
    return df_minerals_grid


def connect_tenement_togrid(df_sampling_tenement, df, df_minerals_grid):
    df_sampling_tenement['cell'] = None
    for id, row in tqdm(df_sampling_tenement.iterrows()):

        sample_p = Point(df['sampling'].iloc[row['sample']].coords)
        for id_cell, cell in df_minerals_grid.iterrows():
            try:
                coords = (cell['coords'])
            except:
                coords = eval(cell['coords'])

            cell_coord = Polygon((coords[0], coords[1], coords[2], coords[3]))

            if cell_coord.contains(sample_p):
                df_sampling_tenement.at[id, 'cell'] = id_cell
                df_sampling_tenement.at[id, 'cell_coords_x'] = cell.coords_x
                df_sampling_tenement.at[id, 'cell_coords_y'] = cell.coords_y
                break


    return df_sampling_tenement


def connect_grid_togeology_and_chemistry_points(df, df_minerals_grid):
    # Iterate through each cell in the grid DataFrame
    for i in tqdm(range(len(df_minerals_grid))):
        coords_cell_i = df_minerals_grid.iloc[i].coords
        coords_cell_i_point = (df_minerals_grid.iloc[i].coords_x, df_minerals_grid.iloc[i].coords_y)
        polygon_cell_i = Polygon(coords_cell_i)

        # Extract chemistry points within the cell
        chemistry_points = []
        for id, row in (df['sampling'].iterrows()):
            point = Point(row.coords)
            p1, p2 = nearest_points(point, polygon_cell_i)
            distance = geopy.distance.geodesic((p1.y, p1.x), (p2.y, p2.x))
            if (distance < 10):

                chemistry_points.append(pointClass.PointClass(id, distance).to_str())

        # Update the 'chemistry_points' column in the grid DataFrame
        df_minerals_grid['chemistry_points'][i] = str(chemistry_points)

        # Extract geology zones close to the cell
        geology_zones = []
        for id, row in (df['geology'].iterrows()):

            polygon_geology = Polygon(row.coords)

            # Find nearest points between cell and geology area polygons
            p1, p2 = nearest_points(polygon_geology, polygon_cell_i)  # Point(coords_cell_i_point))
            # calculate distance between these 2 points
            distance = geopy.distance.geodesic((p1.y, p1.x), (p2.y, p2.x))

            # if distance less than 10km, store the information
            if (distance < 10):
                geology_zones.append(pointClass.PointClass(row.gid, distance).to_str())

        # Update the 'geology_points' column in the grid DataFrame
        df_minerals_grid['geology_points'][i] = str(geology_zones)
    return df_minerals_grid

def calculate_weight(df_min2, alpha = -40):
    '''
        The employed machine learning model is based on the following paper
        Porwal, A., Gonzalez-Alvarez, I., Markwitz, V., McCuaig, T., Mamuse, A., 2010. Weights-of-evidence and logistic regression modeling of magmatic nickel sulfide prospectivity in the yilgarn craton, western australia. Ore Geology Reviews 38, 184–196.
    '''
    try:
        df_min2.reset_index(inplace=True, drop=True)
    except:
        print('Error ')

    # Prepare x and y data for class weight model
    xdata = np.arange(10)
    ydata = class_weight(xdata)
    # Create a polynomial model for class weight
    mymodel = np.poly1d(np.polyfit(xdata, ydata, 3))

    for id, row in tqdm(df_min2.iterrows()):
        # Calculate weights based on distance and weight model for specific columns
        for col in ['komatiite', 'ultramafic', ' mafic', 'crust']:  # ,'Cr','Cu', 'Ni']:
            df_min2.at[id, col + '_weight'] = mymodel(row[col + '_distance']) * weight_model[col]

        for col in ['Cr', 'Cu', 'Ni']:
            try:
                wwi = eval(row[col + '_distance'])
            except:
                wwi = (row[col + '_distance'])
            if wwi == '': wwi = (0, 0, 0)
            qt = wwi[1]
            dist = mymodel([wwi[2]])[0]


            # Calculate weights for specific columns based on quantity and distance
            df_min2.at[id, col + '_weight'] = (qt * dist * weight_model[col])

            if df_min2.at[id, col + '_weight'] < 0:
                df_min2.at[id, col + '_weight'] = 0

    # Calculate total weights and adjust if needed
    column_list = []
    for col in ['komatiite', 'ultramafic', ' mafic', 'crust', 'Cr', 'Cu', 'Ni']:
        column_list.append(col + '_weight')

    # Calculate prospectivity probability using specified alpha value

    df_min2['weight_sum'] = df_min2[column_list].sum(axis=1)
    df_min2.loc[df_min2.weight_sum > 700, 'weight_sum'] = 700
    df_min2['proba_prospectivity' + str(alpha)] = 0.0
    df_min2['proba_prospectivity' + str(alpha)] = [math.exp(alpha + x) / (1 + math.exp(alpha + x)) for x in
                                                   df_min2['weight_sum']]
    return df_min2


def plot_2cells(ax, df_min, tit,  colormap='bwr'):
    """
    Plot heatmap for two cells of prospectivity probability.

    Args:
        ax (matplotlib.axes.Axes): The axis to plot the heatmap on.
        df_min (pd.DataFrame): The dataframe containing prospectivity probability values.
        tit (str): The title for the plot.
        colormap (str, optional): The colormap to use for the heatmap. Defaults to 'bwr'.
    """

    alpha = -40

    df_min['proba' + str(alpha)] = df_min['proba_prospectivity'+  str(alpha)] * 2
    a = []

    for y in sorted(np.unique(df_min.coords_y.tolist()), reverse=True):
        ll = []

        for x in sorted(np.unique(df_min.coords_x.tolist())):
            ll.append(df_min[(df_min.coords_x == x) & (df_min.coords_y == y)][
                          'proba' + str(alpha)].values[0])
            ten = eval(df_min[(df_min.coords_x == x) & (df_min.coords_y == y)]['tenements'].values[0])

        a.append(ll)


    # Create heatmap
    q = sns.heatmap(a, fmt="", linewidth=0.5, center=1, ax=ax, cmap=colormap)
    xlab = sorted([round(x, 3) for x in np.unique(df_min.coords_x.tolist())])
    lx = np.arange(len(xlab))

    if len(df_min)>1000:
        step =4
    else:
        step=3
    # Customize x-axis and y-axis labels
    ax.tick_params(axis='x', labelsize=14, rotation=30)

    ax.set(xticks=lx[0::step], xticklabels=xlab[0::step])#)
    ylab = sorted([round(x, 3) for x in np.unique(df_min.coords_y.tolist())])
    ly = np.arange(len(ylab))

    ax.tick_params(axis='y', labelsize=14, rotation=0)

    ax.set(yticks=ly[0::step], yticklabels=ylab[0::step])  # , rotation=30)
    ax.set_title(tit+'\n', fontsize=22)

    # Customize colorbar
    q.collections[0].colorbar.set_ticks([0, 0.5, 1, 1.5, 2])
    q.collections[0].colorbar.set_ticklabels(['0', '0.25', '0.5', '0.75', '1'])
    q.collections[0].colorbar.ax.tick_params(labelsize=20)
    ax.collections[0].colorbar.set_label("Prospectivity probability", size=20)



def read_shapefile(sf_shape):
    """
    Read a shapefile into a Pandas dataframe with a 'coords'
    column holding the geometry information. This uses the pyshp
    package
    """
    fields = [x[0] for x in sf_shape.fields][1:]
    records = [y[:] for y in sf_shape.records()]
    # records = sf_shape.records()
    shps = [s.points for s in sf_shape.shapes()]
    df = pd.DataFrame(columns=fields, data=records)
    df = df.assign(coords=shps)
    return df


def allpoints_into_zone(df, polygon, datatype, key):
    list_points = []

    for id, row in df[datatype].iterrows():
        coord_point = row.coords[0]
        point = Point(coord_point[0], coord_point[1])
        if polygon.contains(point):
            list_points.append(row[key])
    return list_points


def remove_and_recalculate_prospectivy(df, df_minerals_grid2, tenid, path, list_sampl_into_tenement):
    # select all cells linked to that list1, t2
    df_tmp = df_minerals_grid2[['Cr_distance', 'Cu_distance', 'Ni_distance']]

    list_converned_cells = []
    for symb in ['Cr', 'Cu', 'Ni']:
        try:
            df_tmp[symb + '_sample'] = [eval(x)[0] for x in df_tmp[symb + '_distance']]
        except:
            df_tmp[symb + '_sample'] = [x[0] for x in df_tmp[symb + '_distance']]
        list_converned_cells += df_tmp[df_tmp[symb + '_sample'].isin(list_sampl_into_tenement)].index.tolist()
    set_converned_cells = set(list_converned_cells)

    # Calculate proba prospertivity before deleteting tenement
    df_minerals_grid3 = df_minerals_grid2.iloc[list(set_converned_cells)]
    df_minerals_grid6 = df_minerals_grid2[~df_minerals_grid2.index.isin(df_minerals_grid3.index)]
    #print ('subpart', len(df_minerals_grid3))

    # Calculate proba prospertivity  after deleteting tenement
    df_minerals_grid4 = calculate_closerchem(df['sampling'], df_minerals_grid3, list_sampl_into_tenement)
    #print ('closest cgen',df_minerals_grid4)
    df_minerals_grid5 = calculate_weight(df_minerals_grid4)
    #l = set(df_minerals_grid2.index.tolist()) - set(df_minerals_grid4.index.tolist())
    #df_minerals_grid6 = df_minerals_grid2.iloc[list(l)]
    alpha = -40
    df_minerals_grid6['proba_prospectivity' + str(alpha)]= df_minerals_grid6['proba_prospectivity_orig'+ str(alpha)]
    df_minerals_grid6 = df_minerals_grid6.sort_values(by=['coords_x', 'coords_y'])

    pd.concat([df_minerals_grid5, df_minerals_grid6]).to_csv(
        path + '/tenement_no_' + tenid.replace(' ', '_') + '.csv', index=False)

    return pd.concat([df_minerals_grid5, df_minerals_grid6])


def get_sample_list_into_cell(df_sampling_tenement, cellid):
    return df_sampling_tenement[df_sampling_tenement.cell.isin([cellid])]['sample'].tolist()

def remove_and_recalculate_tenement(df, df_minerals_grid2, tenid, path, list_sampl_into_tenement):


    df_res= remove_and_recalculate_prospectivy(df, df_minerals_grid2, tenid, path, list_sampl_into_tenement)

    # Plot two maps before and after deleteing tenement
    plot_2cells([df_minerals_grid2, df_res],
                ['Cells before deleting tenement '+tenid, 'Cells after deleting tenement '+tenid], tenid)

    # plot CDF absolute Error
    ddf0 = df_minerals_grid2[['coords_x', 'coords_y', 'proba_prospectivity-40']]
    ddf1 = df_res[['coords_x', 'coords_y', 'proba_prospectivity-40']]
    ddf2 = pd.merge(ddf0, ddf1, on=['coords_x', 'coords_y'])
    ddf2['error'] = ddf2['proba_prospectivity-40_x'] - ddf2['proba_prospectivity-40_y']
    ddf2['absolute_error'] = abs(ddf2['proba_prospectivity-40_x'] - ddf2['proba_prospectivity-40_y'])
    df_absolute_error = ddf2.sort_values(by='absolute_error').tail( int(len(ddf2)*.2))

    plt.plot(sorted(ddf2['absolute_error'].tolist()))
    plt.xlabel('cells')
    plt.ylabel('absolute_error')
    plt.title('CDF probability prospectiviry difference before and after deleting tenement'+tenid)
    plt.show()
    plt.close()

    # Plot highest errors into barplot
    df_absolute_error = df_absolute_error.sort_values(by='absolute_error', ascending=False)
    fig, ax = plt.subplots(figsize=(20, 10))
    all_nm = [y for y in df_absolute_error.error.tolist()]
    color = ['r' if y < 0 else 'g' for y in all_nm]
    ax.bar(range(len(all_nm)), all_nm, color=color)
    plt.xlabel('cells')
    plt.ylabel('Error of prospectivity probability')
    ax.set_xticks(range(len(all_nm)))
    ax.set_xticklabels(['Cell #' + str(x) for x in df_absolute_error.index.tolist()], rotation='vertical')
    plt.title('probability prospectiviry difference before and after deleting tenement '+tenid)
    plt.show()




def read_files(path):
    '''
    This function takes a 'path' parameter and reads multiple shapefiles and an Excel file based on the specified paths. It returns a dictionary of dataframes, each corresponding to a specific shapefile or sheet within the Excel file. The code also handles different encodings for reading the shapefiles and performs error handling in case the specified encoding fails.
    The following files are used :
    *   Western Australia State Interpreted Bedrock Geology:
        Source: This dataset can be accessed at https://catalogue.data.wa.gov.au/dataset/1-500-000-state-interpreted-bedrock-geology-dmirs-016. Purpose: It provides a comprehensive interpretation of the state's bedrock geology at a 1:500,000 scale. Utilization: The function employs this data for analysis within the context of the larger project.

    *  Tenements List with Size and Coordinates:
        Content: This file contains a compilation of tenements, including their respective sizes and coordinates. Importance: It plays a crucial role in understanding land tenure and resource distribution within the study area. Integration: The function incorporates this data to enhance the overall dataset's depth and accuracy.

    *   Study Area Surface Auger Sampling:
    Scope: This dataset pertains to Geochemical samples from the near-surface auger measurements data: Ni anomalies, Cu anomalies, and Cr298 anomalies Significance: It offers insights into soil composition and mineral content, contributing to the project's geological assessment. Inclusion: The function integrates this dataset to encompass a comprehensive overview of the study area's geological characteristics.

    !! Due to concerns regarding confidentiality, we are unable to provide access to the final two files: the tenement list and Auger Sampling.

    '''
    # List of files to be read
    list_of_files = ['geology', 'tenement', 'tenement_large']

    # Dictionary to store paths for different shapefiles
    sf_path = {
        'geology': path + '/data/State_interpreted_bedrock_geology_1_500000_DMIRS_016_WA_GDA2020_Public_Shapefile/State_interpreted_bedrock_geology_1_500000_DMIRS_016.shp',
        'tenement': path + '/data/XXXX.shp',
        'tenement_large': path + '/data/XXX.shp',

    }
    # Dictionary to store shapefile objects and dataframes
    sf = {}
    df = {}
    for file in list_of_files:
        try:
            # Try reading shapefiles with specified encoding
            sf[file] = shapefile.Reader(sf_path[file], encoding='cp1252', errors='ignore')
            df[file] = read_shapefile(sf[file])


        except:
            # If the specified encoding fails, try another encoding
            # shift_jisx0213
            sf[file] = shapefile.Reader(sf_path[file], encoding='Shift-JIS', errors='ignore')
            df[file] = read_shapefile(sf[file])

    # Read Excel file for 'sampling'
    sf_path['sampling'] = path+'/data/XXX.XLSX'
    df['sampling'] = pd.read_excel(sf_path['sampling'], sheet_name=None)['Sheet 1']
    pp = path+'/data/XXX2.shp'
    ss = shapefile.Reader(pp, encoding='cp1252', errors='ignore')
    sampling = read_shapefile(ss)

    # Add 'coords' column to the 'sampling' dataframe
    df['sampling']['coords'] = [x[0] for x in sampling.coords]

    return df