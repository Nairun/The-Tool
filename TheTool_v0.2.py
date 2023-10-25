# -*- coding: utf-8 -*-
"""
Created on Fri Jul 01 09:21:45 2022

@author: lukas Broich
"""
# ------------------------------- IMPORTS

import PySimpleGUI as sg
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import trackpy as tp
import pims
from scipy import ndimage
import os

# ------------------------------- FUNCTIONS
# ---------- This is to include a matplotlib figure in a Tkinter canvas

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk


def draw_figure_w_toolbar(canvas, fig, canvas_toolbar):
    if canvas.children:
        for child in canvas.winfo_children():
            child.destroy()
    if canvas_toolbar.children:
        for child in canvas_toolbar.winfo_children():
            child.destroy()
    figure_canvas_agg = FigureCanvasTkAgg(fig, master=canvas)
    figure_canvas_agg.draw()
    toolbar = Toolbar(figure_canvas_agg, canvas_toolbar)
    toolbar.update()
    figure_canvas_agg.get_tk_widget().pack(side='right', fill='both', expand=1)


class Toolbar(NavigationToolbar2Tk):
    def __init__(self, *args, **kwargs):
        super(Toolbar, self).__init__(*args, **kwargs)

#%% Load and filter the dataframe, calculate trajectories and apply driftcorrection

def Calculate_Traj(localizations, trackdist, trackmemory, stubs, maxlength, driftbool, driftsmooth):
    sg.popup_non_blocking('Linking, this might take a few minutes!',title = 'Linking', keep_on_top=True)
## Load data
    try:
        if decodelocs == False:
            data_raw = pd.read_csv(localizations)
            # print('Linking localizations...')
            all_traj = tp.link(data_raw, pos_columns = ['x', 'y'], t_column = 'frame', search_range = trackdist, memory = trackmemory, adaptive_stop = 0.1, adaptive_step = 0.9)
        else:
            try:
                data_raw = pd.read_csv(localizations, skiprows = 3)
                data_1 = data_raw.rename(columns={'frame_ix': 'frame'}) ## Rename 'frame' column so Trackpy recognizes it
                data_2 = data_1[['x','y','phot','frame','prob','bg','phot_sig','x_sig','y_sig']] ## Only keep the columns with relevant information to reduce dataframe
            except:
                try:
                    data_raw = pd.read_csv(localizations, skiprows = 0)
                    data_1 = data_raw.rename(columns={'frame_ix': 'frame'}) ## Rename 'frame' column so Trackpy recognizes it
                    data_2 = data_1[['x','y','phot','frame','prob','bg','phot_sig','x_sig','y_sig']] ## Only keep the columns with relevant information to reduce dataframe
                except:
                    sg.Popup('Something went wrong!',title = 'Error', keep_on_top=True)
                    t = pd.DataFrame()
                    return t
        
## Link localizations into trajectories
            # print('Linking localizations...')
            all_traj = tp.link(data_2, search_range = trackdist, memory = trackmemory, adaptive_stop = 0.1, adaptive_step = 0.9)

        ## Save the unfiltered trajectories in a .csv
        # print('Saving unfiltered trajectories...')
        all_traj.to_csv(folder + '/Unfiltered trajectories.csv', index = False)
        
## Filter stubs
        # print('Filtering stubs...')
        non_stub_traj = tp.filter_stubs(all_traj, stubs) ## Filter out tracks shorter than the defined stub value
        filt_traj = non_stub_traj.groupby('particle').filter(lambda x : len(x)<=maxlength) ## Filter tracks longer than the defined maxlength value

        ## Print the before and after values
        print('Trajectories before filtering:', all_traj['particle'].nunique()) ## Prints trajectory number before stubfiltering
        print('Trajectories after filtering:', filt_traj['particle'].nunique()) ## Prints trajectory number after stubfiltering

        ## Save the filtered trajectories
        # print('Saving stubfiltered trajectories...')
        filt_traj.to_csv(folder + '/Stubfiltered trajectories.csv', index = False)
        
## Drift correction
        ## Apply Driftcorrection
        if driftbool == True:
            print('Applying Drift correction...')
            drift = tp.compute_drift(filt_traj, smoothing = driftsmooth) ## Calculated the drift
            filtrajsmoothed = tp.subtract_drift(filt_traj.copy(), drift) ## Substracts the drift
            print('Saving drift corrected trajectories...')
            filt_traj.to_csv(folder + '/Drift corrected trajectories.csv', index = False)
            
            ## Plot the graph
            fig, ax = plt.subplots() ## Create fig
            ax.plot(drift[drift.columns[0]], label='y') ## Plot y
            ax.plot(drift[drift.columns[1]], label='x') ## Plot x
            ax.set(ylabel='Drift [px]', xlabel='Frame') ## Set axis labels
            plt.title('Drift correction') ## Set Plot title
            plt.legend(); ## Actually plot
            fig = plt.gcf() ## Hook to canvas
            draw_figure_w_toolbar(mainwindow['fig_cv'].TKCanvas, fig, mainwindow['controls_cv'].TKCanvas) ## Show on canvas
            
            ## Save the plot
            # print('Saving drift correction plot...')
            plt.savefig(folder + '/Drift correction.png')
            plt.close(); ## Close the plot
            sg.Popup('Done! \nData has been saved at ' + folder + '/', title = 'Done', keep_on_top=True)
            return filtrajsmoothed
        
        ## Don't apply drift correction
        else:
            sg.Popup('Done! \nData has been saved at ' + folder + '/', title = 'Done', keep_on_top=True)
            return filt_traj
    except:
        sg.Popup('Something went wrong!',title = 'Error', keep_on_top=True)
        t = pd.DataFrame()
        return t
    
#%% Plotting
def Plot(filtraj, decodelocs):
    if not filtraj.empty:
        
        ##################### If there is an image loaded, the plot will created from this code.
        try:
            frame = pims.open(image) ## Open Image
            sg.popup_non_blocking('Plotting Trajectories, this might take a while!', title = 'Plotting', keep_on_top=True)
            if decodelocs == True:
                      frame = ndimage.rotate(frame[0], 90) ## Rotate image if you use decode localizations
                      frame = np.flipud(frame) ## Flip image if you use decode localizations
            else:
                    frame = frame[0] ## Just load the image when you don't use decode localizations
            #################################
            
            fig, ax = plt.subplots() ## Create the Plot

            ## Making the plot
            for traj_id, traj_data in filtraj.groupby('particle'): ## Iterate through each trajectory in your DataFrame
                ax.plot(traj_data['x'], traj_data['y'], lw=1) ## Plot the trajectories

            ## Set plot parameters
            ax.imshow(frame, cmap='gray') ## Superimpose on the image, adjust colormap as needed
            ax.set_xlabel('X-coordinate') ## X-Axis title
            ax.set_ylabel('Y-coordinate') ## Y-Axis title
            
            ## Save and display
            fig.savefig(folder +'/Trajectory_plot.png', dpi=600) ## Save plot with 600 DPI
            DPI = fig.get_dpi() ## Get DPI to adjust offset
        # ------------------------------- you have to play with this size to reduce the movement error when the mouse hovers over the figure, it's close to canvas size
            fig.set_size_inches(404 * 2 / float(DPI), 404 / float(DPI))
        # -------------------------------
            draw_figure_w_toolbar(mainwindow['fig_cv'].TKCanvas, fig, mainwindow['controls_cv'].TKCanvas) ## Send to main window
            
            ## Finished!
            sg.Popup('Finished Plotting', title = 'Finished', keep_on_top=True)
            return
        
        ##################### If there is no image loaded, the plot will created from this code.
        except:
            fig, ax = plt.subplots() ## Create the Plot

            ## Making the plot
            for traj_id, traj_data in filtraj.groupby('particle'): ## Iterate through each trajectory in your DataFrame
                ax.plot(traj_data['x'], traj_data['y'], lw=1) ## Plot the trajectories

            ## Set plot parameters
            ax.set_aspect('equal') ## Ensures equal axis
            ax.set_xlabel('X-coordinate') ## X-Axis title
            ax.set_ylabel('Y-coordinate') ## Y-Axis title
            
            ## Save and display
            fig.savefig(folder +'/Trajectory_plot.png', dpi=600) ## Save plot with 600 DPI
            DPI = fig.get_dpi() ## Get DPI to adjust offset
        # ------------------------------- you have to play with this size to reduce the movement error when the mouse hovers over the figure, it's close to canvas size
            fig.set_size_inches(404 * 2 / float(DPI), 404 / float(DPI))
        # -------------------------------
            draw_figure_w_toolbar(mainwindow['fig_cv'].TKCanvas, fig, mainwindow['controls_cv'].TKCanvas)  ## Send to main window
            
            ## Finished!
            sg.Popup('Finished Plotting', title = 'Finished', keep_on_top=True)
            return

    else:
        sg.Popup('You need to calculate Trajectories first!', title = 'Error', keep_on_top=True)
        return
        
        
#%% Calculate Diffusions       
def calculate_diffusion(filt_traj, campix, framerate, lagtime, conf_tresh_D):
    try:
        ## Popup
        sg.popup_non_blocking('Calculating, this might take a few minutes! Do not interrupt!', title = 'Calculating', keep_on_top=True)
        
        ## Calculate individual MSDs
        campix = campix/1000
        iMSD = tp.imsd(filt_traj, mpp = campix, fps = framerate)
         
        ## Save the IMSDs in a .csv
        iMSD.to_csv(folder + '/iMSDs.csv')
        
        ## Plot the iMSD graph
        fig, ax = plt.subplots() ## Create fig
        ax.plot(iMSD.index, iMSD, 'k-', alpha = 0.1)  ## Black lines, semitransparent
        ax.set(ylabel=r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]', xlabel='lag time $t$ [s]') ## Label axis
        ax.set_xscale('log') ## Log x-axis
        ax.set_yscale('log') ## Log y-axis  
        plt.title('Individual mean squared displacements') ## Set Plot title
        fig = plt.gcf() ## Hook to canvas
        draw_figure_w_toolbar(mainwindow['fig_cv'].TKCanvas, fig, mainwindow['controls_cv'].TKCanvas) ## Show on canvas
        
        ## Save the iMSD graph
        plt.savefig(folder + '/iMSD plot.png')
        plt.close(); ## Close the plot
        
        ## Calculating the diffusion coefficients
        taustart = 1/framerate
        tauend = (1/framerate)*lagtime
        
        ## Calculate all diffusion and power-law exponents
        filt_traj['Init D'] = 0
        filt_traj['Confined D'] = 0
        count_all = 0
        count_conf = 0
    
        filt_traj['Power-law exponent'] = 0
        
        ## Calculate all diffusion coefficients and power law exponents.
        for i in iMSD:
            iCalc = iMSD.loc[taustart:tauend,i] ## Access MSD of  the first x taus
            iFIT = tp.utils.fit_powerlaw(iCalc, plot = False) ## Do a fit
            iDiffu = iFIT.iloc[0]['A'] / 4 ## Recieve Diffusion coefficient
            iAlpha = iFIT.iloc[0]['n'] ## Recieve Powerlawexponent
            filt_traj.loc[filt_traj['particle'] == i,'Init D'] = iDiffu ## Add Diffusion coefficent to the table
            filt_traj.loc[filt_traj['particle'] == i,'Power-law exponent'] = iAlpha ## Add Powerlawexponent to the table
            count_all = count_all + 1 ## Count the number of analyzed particles
            
            ## Check if a particle is considered to be confined 
            if iDiffu <= conf_tresh_D: 
                filt_traj.loc[filt_traj['particle'] == i,'Confined D'] = 1 ## Add a 1 to the confined box in the table
                count_conf = count_conf + 1 ## Count overall number of immobile particles for % analysis
              
            ## Calculate percent of immobile particles
            p_d = round((count_conf/count_all)*100 , 2)    
         
        ## Save the table with the diffusion coefficients
        filt_traj_no_dupes = filt_traj.drop_duplicates(subset=['particle']) ## Remove all duplicate entries from particle column
        filt_traj_no_dupes_2 = filt_traj_no_dupes[['x','y','frame','particle','Init D','Confined D','Power-law exponent']] ## Only keep the columns with relevant information to reduce dataframe
        filt_traj_no_dupes_2 = filt_traj_no_dupes_2[['particle','x','y','frame','Init D', 'Power-law exponent', 'Confined D']] ## Reorder columns
        filt_traj_no_dupes_2 = filt_traj_no_dupes_2.rename(columns={'frame': 'First occurence frame','particle': 'Track ID', 'x': 'x-origin', 'y': 'y-origin', 'Init D': 'Diffusion coefficient', 'Confined D': 'Confined Y/N'}) ## Rename columns
        filt_traj_no_dupes_2.to_csv(folder + '/Trajectories with Diffusion and Powerlaw.csv', index=False) ## Save
    
        ## Write and save a textfile with the percentage of immobile particles
        out_str_d = str(f"{p_d}" + '% of your particles are considered immobile.\nTreshold for being immobile was: ' + f"{conf_tresh_D}" + ' µm²/s\n The total number of particles is: ' + f"{count_all}")
        conf_perc_name = folder + '/Fraction of immobile particles by diffusion.txt' ## Adds directory and base name of localizations
        text_file = open(conf_perc_name, "w") ## Opens a new textfile
        text_file.write(out_str_d) ## Write the textfile
        text_file.close() ## Close the textfile
        
        ## Popup that the calculation finished
        outputstring = str('Calculation finished,\n' + f"{count_all}" + ' trajectories were analyzed!\nData can be found at:\n' + folder) ## Create string
        sg.Popup(outputstring, title = 'Done!', keep_on_top=True) ## Create pop up
        
        return filt_traj_no_dupes_2 ## Return the data
    except:
        sg.Popup('Something went wrong!',title = 'Error', keep_on_top=True) ## Create error popup
        return
        
#%% Calculate Trajectories from coordinates
def Calculate_Coords(filt_traj_analyzed, xcoord, ycoord,framerate, campix, squaresize, lagtime):
    try:
        df = filt_traj_analyzed ## Rename df for easier programming
    
        ## Make subdirectory
        p1 = 'x' + str(xcoord) ## Str for x coordinate
        p2 = 'y' + str(ycoord) ## Str for y coordinate
    
        outfolder_path = folder.replace("\\", "/") + '/Squares/' + p1 + p2 ## Assign folder
        if not os.path.exists(outfolder_path):
            os.makedirs(outfolder_path) ## Create folder if it does not exist
    
        xlimup = xcoord + squaresize/2 ## Upper X
        xlimdn = xcoord - squaresize/2 ## Lower X
        ylimup = ycoord + squaresize/2 ## Upper Y
        ylimdn = ycoord - squaresize/2 ## Lower Y
    
        ## Subset the dataframe and extract values
        spot_traj = df[(df['x-origin'] > xlimdn) & (df['x-origin'] < xlimup) & (df['y-origin'] > ylimdn) & (df['y-origin'] < ylimup)] ## Subset the dataframe
        particles_to_select = spot_traj['Track ID'].unique().tolist() ## Create list with Track-IDs
        spot_values = df[df['Track ID'].isin(particles_to_select)] ## Grab Rows with matching Track-IDs
    
        ## Save the values
        spot_values.to_csv(outfolder_path + '/' + p1 + p2 + '_Trajectories.csv', sep = ',', index = False)
    
        ## Pop up
        length = len(particles_to_select) ## Get number of Tracks
        outputstring = ('Done! \nWithin your square are ' + f"{length}" + ' tracks.\nData can be found at:\n' + outfolder_path)
        sg.Popup(outputstring, title = 'Done!', keep_on_top=True)
        
        return

    except:
        sg.Popup('Something went wrong!',title = 'Error', keep_on_top=True) ## Error message popup
        return
        
# ----------------------Menu Window
def menuwindow(framerate, decodelocs, stubs, trackdist, trackmemory, maxlength, driftbool, driftsmooth, campix, squaresize, lagtime, conf_tresh_D):
    
    ## Layout of the Options menu
    menulayout = [
        [sg.Text("General", font=('Any 15'))],
        [sg.Checkbox('Use DECODE localizations', default=decodelocs, change_submits=False, key='-DECODELOCS-')],
        [sg.Text("Camera pixel size [nm]:", size=(32, 1)), sg.In(size=(5, 1), default_text=campix, change_submits=False, key='-CAMPIX-')],
        [sg.Text("Framerate [FPS]:", size=(32, 1)), sg.In(size=(5, 1), default_text=framerate, change_submits=False, key='-FRAMERATE-')],
        [sg.Text("Linking parameters", font=('Any 15'))],
        [sg.Checkbox('Apply drift correction', size=(14, 1), default=driftbool, change_submits=False, key='-DRIFTBOOL-'),
         sg.Text("Smoothing factor: "), sg.In(size=(5, 1), default_text=driftsmooth, change_submits=False, key='-DRIFTSMOOTH-')],
        [sg.Text("Max. linking distance [px]:", size=(32, 1)), sg.In(size=(5, 1), default_text=trackdist, change_submits=False, key='-TRACKINGDISTANCE-')],
        [sg.Text("Linking memory [Frames]:", size=(32, 1)), sg.In(size=(5, 1), default_text=trackmemory, change_submits=False, key='-TRACKINGMEMORY-')],
        [sg.Text("Min. Track length [Frames]:", size=(32, 1)), sg.In(size=(5, 1), default_text=stubs, change_submits=False, key='-STUBS-')],
        [sg.Text("Max. Track length [Frames]:", size=(32, 1)), sg.In(size=(5, 1), default_text=maxlength, change_submits=False, key='-MAXLENGTH-')],     
        [sg.Text("Calculation parameters", font=('Any 15'))],
        [sg.Text("Square size [px]:", size=(32, 1)), sg.In(size=(5, 1), default_text=squaresize, change_submits=False, key='-SQUARESIZE-')],
        [sg.Text("Max. lag time [Frames]:", size=(32, 1)), sg.In(size=(5, 1), default_text=lagtime, change_submits=False, key='-LAGTIME-')],
        [sg.Text("Immobile particle threshold [µm²/s]:", size=(32, 1), key='-CONF_D_TEXT-'), sg.In(size=(5, 1), default_text=conf_tresh_D, change_submits=False, key='-CONF_D-')],
        [sg.Button("Submit"), sg.Button('Exit')],
    ]
    
    window = sg.Window("Set Parameter", menulayout, modal=True) ## Create Window

    ## Functions in the options menu
    ##Keep Window open
    while True:
        event, values = window.read()
        if event in (sg.WINDOW_CLOSED, "Cancel"):
            break
        ## Submit
        elif event == 'Submit':
            
            ## Framerate
            try:
                framerate = float(values['-FRAMERATE-'])
                if framerate <= 0:
                    sg.Popup('Framerate has to be a positive value \nValue is set to default', title = 'Error', keep_on_top=True)
                    framerate = float(31.5)
                    continue
            except ValueError:
                sg.Popup('Framerate must be float!', keep_on_top=True)
                continue
            
            ## Trackingdistance
            try:
                trackdist = float(values['-TRACKINGDISTANCE-'])
                if trackdist <= 0:
                    sg.Popup('Trackingdistance has to be a postive value \nValue is set to default', title = 'Error', keep_on_top=True)
                    trackdist = float(0.8)
                    continue
            except ValueError:
                sg.Popup('Trackingdistance must be float!', keep_on_top=True)
                continue
            
            ## Trackingmemory
            try:
                trackmemory = int(values['-TRACKINGMEMORY-'])
                if trackmemory < 0:
                    sg.Popup('Trackingmemory may not be negative, \nValue is set to default', title = 'Error', keep_on_top=True)
                    trackmemory = int(1)
                    continue
            except ValueError:
                sg.Popup('Trackingmemory must be int!', keep_on_top=True)
                continue
            
            ## Maxlength
            try:
                maxlength = int(values['-MAXLENGTH-'])
                if maxlength <= 0:
                    sg.Popup('Maximum tracklength has to be a postive value \nValue is set to default', title = 'Error', keep_on_top=True)
                    maxlength = int(100)
                    continue
            except ValueError:
                sg.Popup('Maximum tracklength must be int!', keep_on_top=True)
                continue
            
            ## Stubs
            try:
                stubs = int(values['-STUBS-'])
                if stubs <= 0:
                    sg.Popup('Minumum tracklength has to be a postive value \nValue is set to default', title = 'Error', keep_on_top=True)
                    stubs = int(10)
                    continue
            except ValueError:
                sg.Popup('Stubs must be int!', keep_on_top=True)
                continue
            
            ## Driftsmoothing
            try:
                driftsmooth = int(values['-DRIFTSMOOTH-'])
                if driftsmooth < 0:
                    sg.Popup('Smoothing factor may not be negative \nValue is set to default', title = 'Error', keep_on_top=True)
                    driftsmooth = int(0)
                    continue
            except ValueError:
                sg.Popup('Smoothing factor must be int!', keep_on_top=True)
                continue
            
            ## Camerapixelsize
            try:
                campix = float(values['-CAMPIX-'])
                if campix <= 0:
                    sg.Popup('Camera pixel size has to be a positive value \nValue is set to default', title = 'Error', keep_on_top=True)
                    campix = float(160)
                    continue
            except ValueError:
                sg.Popup('Camera pixel size must be float!', keep_on_top=True)
                continue
            
            ## Squaresize
            try:
                squaresize = float(values['-SQUARESIZE-'])
                if squaresize <= 0:
                    sg.Popup('Squaresize has to be a positive value \nValue is set to default', title = 'Error', keep_on_top=True)
                    squaresize = float(5)
                    continue
            except ValueError:
                sg.Popup('Squaresize must be float!', keep_on_top=True)
                continue
            
            ## Lagtime
            try:
                lagtime = int(values['-LAGTIME-'])
                if lagtime <= 0:
                    sg.Popup('Lagtime has to be a postive value \nValue is set to default', title = 'Error', keep_on_top=True)
                    lagtime = int(10)
                    continue
            except ValueError:
                sg.Popup('Lagtime must be int!', keep_on_top=True)
                continue
            
            ## Treshold for confinement
            try:
                conf_tresh_D = float(values['-CONF_D-'])
                if conf_tresh_D <= 0:
                    sg.Popup('Lagtime has to be a postive value \nValue is set to default', title = 'Error', keep_on_top=True)
                    conf_tresh_D = float(0.01)
                    continue
            except ValueError:
                sg.Popup('Threshold must be float!', keep_on_top=True)
                
            ## Decodelocs Bool
            decodelocs = values['-DECODELOCS-']
            
            ### Driftbool Bool
            driftbool = values['-DRIFTBOOL-']
            
            
        break
    window.close()
    return framerate, decodelocs, stubs, trackdist, trackmemory, maxlength, driftbool, driftsmooth, campix, squaresize, lagtime, conf_tresh_D

###### End of Menu Window


# ------------------------------- Main Window
## Set all variables default values
colorscheme = '#5c7ea1'
decodelocs = bool(True)
framerate = float(31.5)
stubs = int(10)
trackdist = float(0.8)
trackmemory = int(0)
maxlength = int(100)
filtraj = pd.DataFrame()
localizations = str("")
driftbool = bool(True)
driftsmooth = int(0)
campix = float(160)
folder = os.getcwd()
squaresize = float(5)
image = str()
lagtime = int(3)
conf_tresh_D = float(0.01)
filt_traj_analyzed = pd.DataFrame()

# ---------- Menubar Layout

menu_def = [['&Options...', ['&Set Parameter', 'E&xit']],
            ['&Help...',['&Parameter info', '&About']], ]


# ---------- Main Window Layout
menuline = [sg.Menu(menu_def, pad=(0,0), k='-CUST MENUBAR-')]

firstline = [sg.Text('Outputfolder', size=(10,1), background_color=colorscheme), sg.In(size=(15,1), 
             change_submits=True, enable_events=True ,key='-FOLDER-'), sg.FolderBrowse(), 
             sg.Text('Imagefile', background_color=colorscheme), sg.In(size=(15,1), change_submits=True, enable_events=True ,key='-IMAGE-'), 
             sg.FileBrowse(file_types=(("*.png *.jpg *.tif", "*.png *.jpg *.tif"),))]

secondline = [sg.Text('Localizations', size=(10,1), background_color=colorscheme), sg.In(size=(15,1), change_submits=True, 
              enable_events=True ,key='-LOCALIZATIONS-'), sg.FileBrowse(file_types=(("*.csv", "*.csv"),))]

thirdline = [sg.B('Link trajectories', key ='-LINKTRAJECTORIES-'), sg.B('Calculate diffusions', key='-CALCDIFFU-'), 
             sg.B('Plot Trajectories', key='-PLOTTRAJECTORIES-')]
             

fourthline = [sg.T('Controls:', background_color=colorscheme), sg.Canvas(key='controls_cv'),
              sg.Canvas(key='controls_cv')]

plotline =    [sg.Column(
        layout=[
            [sg.Canvas(key='fig_cv',
                       size=(400 * 2, 400) # it's important that you set this size
                       )]
        ],
        background_color='#DAE0E6',
        pad=(0, 0)
    )]

lastline = [sg.Text('Coordinates:', size=(10, 1), background_color=colorscheme),
            sg.Text('X', size=(1, 1), background_color=colorscheme), sg.InputText('', size= (5,1), key="-X-"), 
            sg.Text('Y', size=(1, 1), background_color=colorscheme), sg.InputText('', size= (5,1), key="-Y-"), 
            sg.B('Calculate spot', key = '-CALCULATECOORDS-')]
######

layout = menuline, firstline, secondline, thirdline, fourthline, plotline, [sg.B('Exit'), lastline]
mainwindow = sg.Window('The Tool', layout, background_color=colorscheme)

## Texts
Parameter_info_text = '''Use decode localizations:
Check this box if you use Decode (2020, Artur Speiser, Lucas-Raphael Mueller et al.) to fit your data. CSV formatted localization files created in DECODE can be directly loaded into the program without any prior processing. However, it is recommended to filter the localizations before tracking. 
If you DON’T use Decode localizations your CSV file needs to have at least an ‘x’, an ‘y’ and a ‘frame’ column.
NOTE: Checking this box also rotates and mirrors the underlying image, this is necessary because of different conventions on the coordinate system in Decode.

Camera pixel size:
The pixel size of the camera in nanometers that is used to capture the video data. It can be found in the manufacturer's specifications. For correct results it is mandatory that this parameter is set correctly.

Framerate: 
The framerate used during acquisition of the video data in frames per second. For correct results it is mandatory that this parameter is set correctly.

Apply drift correction:
Calculate the ensemble drift x,y(t) and substract it from the trajectories. 
The smoothing factor determines over how many frames the drift should be smoothed using a forward-looking rolling mean.

Max. linking distance:
Sets the maximum distance in pixels, over that two localizations in consecutive frames will be linked.

Linking memory:
The number of frames before a localization is ‘forgotten’. Useful if your localizations disappear for a few frames, because of a bad signal or fitting.

Min. Track length:
A filter function to remove tracks that are shorter than the set amount of frames.

Max. Track length:
A filter function to remove tracks that are longer than the set amount of frames.

Square size:
Sets the side length of the square in pixel used when analyzing tracks at specific coordinates. 

Max. lag time:
Sets the number of frames used to calculate the diffusion coefficient for each particle. For example, if set to ‘3’ the diffusion coefficients will only be calculated from the first three data points of each track.

Immobile particle threshold:
Particles with a diffusion coefficient slower than this value will be considered immobile. '''

About_text = ''' 
This software was created by Lukas Broich in 2023 at the nanoinfectionlab from the HZI Brunswick (AG Sieben).
For further information visit https://nanoinfection.org/

Feel free to modify and change the software for any non-commercial needs.

This software makes use of the 
TrackPy library v0.5.0 (10.5281/zenodo.4682814)
For more information check out https://github.com/soft-matter/trackpy/tree/v0.5.0
'''

# ---------- Main Window Loop
while True:
    event, values = mainwindow.read()
    print(event, values)
    
    if event in (sg.WIN_CLOSED, 'Exit'):  ## Always,  always give a way out!
        break
    
    elif event == '-PLOTTRAJECTORIES-':
            Plot(filtraj, decodelocs)

    elif event == '-FOLDER-':
            folder = values['-FOLDER-']  
            
    elif event == '-IMAGE-':
            image = values['-IMAGE-'] 
            
    elif event == '-LOCALIZATIONS-':
            localizations = values['-LOCALIZATIONS-']
            
    elif event == '-LINKTRAJECTORIES-':
        if localizations:
           filtraj = Calculate_Traj(localizations, trackdist, trackmemory, stubs, maxlength, driftbool, driftsmooth)
        else:
            sg.Popup('You need to select a file first!', keep_on_top=True)
            
    elif event == 'Set Parameter':
            framerate, decodelocs, stubs, trackdist, trackmemory, maxlength, driftbool, driftsmooth, campix, squaresize, lagtime, conf_tresh_D = menuwindow(framerate, decodelocs, stubs, trackdist, trackmemory, maxlength, driftbool, driftsmooth, campix, squaresize, lagtime, conf_tresh_D)
    
    elif event == '-CALCDIFFU-':
            filt_traj_analyzed = calculate_diffusion(filtraj, campix, framerate, lagtime, conf_tresh_D)
    
    elif event == '-CALCULATECOORDS-':
        try:
            xcoord = float(values['-X-'])
            ycoord = float(values['-Y-'])
            if not filt_traj_analyzed.empty:
                Calculate_Coords(filt_traj_analyzed, xcoord, ycoord, framerate, campix, squaresize, lagtime)
            else:
                sg.Popup('You need to calculate diffusion coefficients first!', keep_on_top=True)
                
        except ValueError:
            sg.Popup('Your coordinates need to be floats', keep_on_top=True)
    
    elif event == 'About':
            sg.Popup(About_text, title = 'About', custom_text = 'Close', keep_on_top=True)
            
    elif event == 'Parameter info':
            sg.Popup(Parameter_info_text, title = 'About', custom_text = 'Close', keep_on_top=True)
mainwindow.close()