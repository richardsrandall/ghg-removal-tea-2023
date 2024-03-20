
# This module contains a bunch of helper functions for quickly creating and modifying plots of model outputs.

# Imports
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.ticker as mtick
import numpy
from math import floor
import scipy
#import ghg_conversions
from labellines import labelLines
import scienceplots
from matplotlib import rcParams

# Some convenient fixed formatting values and conversion functions
percent_formatter = mtick.FuncFormatter(lambda y, _: ('{:g}').format(y*100)+'%')
dollar_formatter = mtick.FuncFormatter(lambda y, _: ('$'+'{:g}').format(y))
default_formatter = mtick.FuncFormatter(lambda y, _: ('{:g}').format(y))

# Class representing a custom plot with some labor-saving elements. r_plot just refers to my name to disambiguate from plt.
class r_plot:

    # Constructor
    def __init__(self,style="Large & Bright",fig=None,ax=None,figsize=None,**kwargs):
        if fig is None:
            if figsize is None:
                self.fig,self.ax = plt.subplots()
            else:
                self.fig,self.ax = plt.subplots(figsize=figsize)
        else:
            self.fig = fig
            self.ax = ax
        self.style=style
        if 'force_cycler' in kwargs.keys():
            self.force_cycler = kwargs['force_cycler']
        else:
            self.force_cycler = 0

    # Flexible function to generate many of the plots we want with $ on Y-axis and ppm, diameter, or QY on X-axis
    # which_is_x and which_is_lines should be fields of params, e.g. params.wind_speed or params.DQY
    # function is the actual cost solver function for a particular system
    def generate_plot(self, params, function, which_is_x, x_values, which_is_lines=(None), line_values=(None), line_labels=("Series 1"),zorder=0,**kwargs):
        fig = self.fig
        ax = self.ax
        z = zorder
        for l, label in zip(line_values,line_labels): # Step through each line you want to draw (e.g. wind speed)
            if which_is_lines is not None:
                setattr(params, which_is_lines, l)
            y_values = []
            for x in x_values: # Step through the x values to create ordered pairs to plot
                setattr(params, which_is_x, x)
                y_values.append(function(params)) # Modify this line once solvers return dict's instead of floats
            z=z+(-10 if which_is_lines=='air_velocity' else 10)#Kludge to reverse Z order
            ax.plot(x_values, y_values, label=label,zorder=z)
        # At this point, the data are plotted on lin-lin axes but are not formatted in any way.

    # Format left and/or right axes with some converting function(s). Set the figure size as needed.
    def format_axes_log(self,x_formatter='{:g}', y_formatter='{:g}', x_lims=None, y_lims=None, minor_tick_labels=(), rotate_x=0, which_grids=(), sec_axis_functions=None, **kwargs):
        fig = self.fig
        ax = self.ax
        # Start formatting the x axis...
        ax.set_xscale('log')
        ax.xaxis.set_major_formatter(x_formatter)
        # Start formatting the primary y axis...
        ax.set_yscale('log')
        ax.yaxis.set_major_formatter(y_formatter)
        # Start formatting the secondary y axis....
        if sec_axis_functions is not None:
            secax = ax.secondary_yaxis('right',functions=sec_axis_functions)
            secax.yaxis.set_major_formatter(y_formatter)
            self.secax = secax
        else:
            secax = None
            self.secax = None
        # Format limits of plot
        if x_lims is not None:
            ax.set_xlim(x_lims[0],x_lims[1])
        if y_lims is not None:
            ax.set_ylim(y_lims[0],y_lims[1])      
        # Do grids all at once:
        if 'major x' in which_grids:
            ax.grid(visible=True,axis='x',which='major')
        if 'both x' in which_grids:
            ax.grid(visible=True,axis='x',which='major')
            ax.grid(visible=True,axis='x',which='minor')
        if 'major y' in which_grids:
            ax.grid(visible=True,axis='y',which='major')
        if 'both y' in which_grids:
            ax.grid(visible=True,axis='y',which='minor')
            ax.grid(visible=True,axis='y',which='major')   
        # Do minor tickmarks:
        if ('x' in minor_tick_labels) or ('both' in minor_tick_labels):
            ax.xaxis.set_minor_formatter(x_formatter)
            for label in ax.get_xticklabels(minor=True)[1::2]:
                label.set_visible(False)
            for label in ax.get_xticklabels(minor=True)[3::4]:#Deletes 8's because they look bad on log plots   
                label.set_visible(False)
            label.set_visible(False)
        if ('y' in minor_tick_labels) or ('both' in minor_tick_labels):
            ax.yaxis.set_minor_formatter(y_formatter)
            for label in ax.get_yticklabels(minor=True)[1::2]:
                label.set_visible(False)
            for label in ax.get_xticklabels(minor=True)[3::4]:#Deletes 8's because they look bad on log plots
                label.set_visible(False)
            if secax is not None:
                secax.yaxis.set_minor_formatter(y_formatter)
                for label in secax.get_yticklabels(minor=True)[1::2]:
                    label.set_visible(False)
                for label in ax.get_xticklabels(minor=True)[3::4]:#Deletes 8's because they look bad on log plots
                    label.set_visible(False)
        # Rotate ticks if desired
        if rotate_x!=0:
            for label in ax.get_xticklabels(): # Make the labels line up with the ticks
                label.set_rotation(rotate_x)
                if rotate_x == 90:
                    pass # Turns out default settings give the desired alignment
                else:
                    label.set_ha('right')
                    label.set_rotation_mode('anchor')
            for label in ax.get_xticklabels(minor=True):
                label.set_rotation(rotate_x)
                if rotate_x == 90:
                    pass # Turns out default settings give the desired alignment
                else:
                    label.set_ha('right') # Make the labels line up with the ticks
                    label.set_rotation_mode('anchor')

    # Format left and/or right axes with some converting function(s). Set the figure size as needed.
    def format_axes_linear(self,x_formatter='{:g}', y_formatter='{:g}', x_lims=None, y_lims=None, rotate_x=0, which_grids=(),**kwargs):
        fig = self.fig
        ax = self.ax
        self.secax = None
        # Start formatting the x axis...
        ax.xaxis.set_major_formatter(x_formatter)
        # Start formatting the primary y axis...
        ax.yaxis.set_major_formatter(y_formatter)
        # Format limits of plot
        if x_lims is not None:
            ax.set_xlim(x_lims[0],x_lims[1])
        if y_lims is not None:
            ax.set_ylim(y_lims[0],y_lims[1])      
        # Do grids all at once:
        if ('x' in which_grids) or ('both' in which_grids):
            ax.grid(visible=True,axis='x',which='major')
        if ('y' in which_grids) or ('both' in which_grids):
            ax.grid(visible=True,axis='y',which='major')
        # Rotate ticks if desired
        if rotate_x!=0:
            for label in ax.get_xticklabels():
                label.set_rotation(rotate_x)
                label.set_ha('right')
            for label in ax.get_xticklabels(minor=True):
                label.set_rotation(rotate_x)
                label.set_ha('right')
        # Return values for future functions to use; could have used gca() too.
        return (fig, ax, None)

    # Add labels and titles to the plot
    def add_labels(self, x_label=None, y_label=None, sec_y_label=None, title=None,**kwargs):
        fig = self.fig
        ax = self.ax
        secax = self.secax
        if x_label is not None:
            ax.set(xlabel=x_label)
        if y_label is not None:
            ax.set(ylabel=y_label)
        if title is not None:
            ax.set(title=title)
        if (secax is not None) and (sec_y_label is not None):
            secax.set(ylabel=sec_y_label)


    def do_cosmetics_setup(self):
            # Colors available to me can be found here: https://github.com/garrettj403/SciencePlots/wiki/Gallery#color-cycles
            # Styles available can be found here: https://github.com/garrettj403/SciencePlots
        if self.style=="1 Column":
            rcParams.update({'figure.autolayout': True})
            plt.style.use(['science','nature','no-latex','vibrant'])
        if self.style=="Large & Bright":
            rcParams.update({'figure.autolayout': True})
            plt.style.use(['high-vis'])
        # Delete this
        from cycler import cycler
        #print('formatting')
        #plt.style.use(['bright'])
        plt.style.use(['high-vis'])
        plt.rcParams['lines.linewidth'] = 1.3  
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        colors.append('#0504aa')
        colors = colors[self.force_cycler:]+colors[:self.force_cycler]
        linestyles = ['-', '--', '-.','-','--','-.','-']
        linestyles = linestyles[:len(colors)]
        self.ax.set_prop_cycle((cycler('color',colors) + cycler('linestyle', linestyles)))

    # Make any cosmetic changes desired
    def do_cosmetics_post(self,no_labels=False):
        fig = self.fig
        ax = self.ax
        secax = self.secax
        if self.style=="1 Column":
            if not no_labels:
                ax.legend()
        if self.style=="Large & Bright":
            xl = ax.get_xlim()
            lines = ax.get_lines()
            xv = numpy.geomspace(xl[0],xl[len(xl)-1],2+len(lines))
            xv = xv[1:(len(xv)-1)]
            if not no_labels:
                labelLines(lines,xvals=xv)

    # Finish the plot
    def show_plot(self):
        plt.show()

    # Save the plot
    def save_plot(self,savename,filetype):
        plt.savefig((savename+filetype))














