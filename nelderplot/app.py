"""
GUI interface module for Nelder Plot Designer

@author: Degi Harja Asmara
"""

from tkinter import Tk, Label, Radiobutton, StringVar, IntVar, Checkbutton, Scrollbar, filedialog,\
    PhotoImage, simpledialog, messagebox
from tkinter.ttk import Notebook, Treeview
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Wedge, Circle
import pkg_resources
import tkinter as tk

from nelderplot.core import NelderPlot, Point 
    
class PlotCanvas(tk.Frame):
    """ the canvas class for drawing the Nelder plot design """
    
    borderpoint_color = "black"
    datapoint_color = ("red", "blue", "green")
    wheel_color = "black"
    spoke_color = "black"
    defaul_dpi = 100
        
    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        self.figure = Figure(figsize=(5,5), dpi=self.defaul_dpi, tight_layout = True)
        self.plot = self.figure.add_subplot(111)
        self.plot.axis('scaled')
        self.plot.set_axis_off()
        self.canvas = FigureCanvasTkAgg(self.figure, self)
        self.canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.nelder_plot = None
            
    def updatePlotSize(self):
        """ redraw the plot when the size changes """
        if self.nelder_plot.width > 0 and self.nelder_plot.height > 0:
            self.plot.axis([0,self.nelder_plot.width,0,self.nelder_plot.height])
            self.plot.set_axis_on()
        else:
            self.plot.set_axis_off()
                 
    def setNelderPlot(self, nelder_plot: NelderPlot):
        """ define the Nelder plot object to be drawn """
        self.plot.clear()
        self.nelder_plot = nelder_plot
        self.updatePlotSize()
        if not nelder_plot.isValid(): 
            self.canvas.draw()
            return
        for w in nelder_plot.wheels.values():
            wheel = Wedge((w.center_point.x, w.center_point.y), w.radius, 0, 360, linewidth=0.5, width = 0, 
                           edgecolor = self.wheel_color, alpha = 0.5, linestyle = (0, (4, 8)))
            self.plot.add_patch(wheel)
        for l in nelder_plot.spokes:
            line = Line2D([l.intersect_point.x, l.end_point.x], [l.intersect_point.y, l.end_point.y], linewidth=0.5,
                              color = self.spoke_color, alpha = 0.5, linestyle = ":")
            self.plot.add_line(line) 
        p_size = min(self.nelder_plot.width, self.nelder_plot.height)/120        
        for p in nelder_plot.datapoints:
            point = Circle((p.x, p.y), p_size, facecolor = self.datapoint_color[p.group_id])
            if p.is_border: 
                point.set_edgecolor(self.borderpoint_color)
                point.set_radius(p_size * 0.9)
            self.plot.add_patch(point)
        self.canvas.draw()
    
    def saveImage(self, fileName, dpi = defaul_dpi):
        """ save the image to the file and size (DPI) defined """
        self.figure.set_dpi(dpi)
        self.figure.savefig(fileName)
        self.figure.set_dpi(self.defaul_dpi)

class NelderGUI:
    """ The main GUI for designing the plot """
    
    title = "Nelder Plot Designer"
    param_type = ()
    field_center = ("X", "Y")
    center_options = ("Center of plot", "Bottom left corner", "Half width", "Half height", "Customize")
    label_alt_spoke = "Alternate spokes"
    output_page = ("Plot Design", "Variables", "Data Points")
    font_header = "default 9 bold"
    txt_inp_save_image = "Resolution size (DPI):"
    txt_inp_save_image_title = "Input image size"
    output_text_options = ("Nelder parameters and summary", "Plot variables", "Data points")
    toolbar_text = ("Save Data", "Export Plot Image", "Export KML file", "About")
    toolbar_images = ("save.png", "save_image.png", "kml.png", "about.png")

    def __init__(self, master):
        self.master = master
        master.title(self.title)
        form_main = tk.Frame(master)

        # create the interface for the parameters input
        self.fieldEnt = {}
        self.nelder_plot = NelderPlot()
        self.field_summary = self.nelder_plot.summary_info
        self.param_type = self.nelder_plot.param_type_label
        self.opt_param = StringVar()
        self.opt_param.set(self.param_type[0])
        r1 = Radiobutton(form_main, text=self.param_type[0], variable=self.opt_param, 
                         value=self.param_type[0], command=self.selectParams, font=self.font_header)
        r2 = Radiobutton(form_main, text=self.param_type[1], variable=self.opt_param, 
                         value=self.param_type[1], command=self.selectParams, font=self.font_header)
        form_density = self.makeForm(form_main, self.nelder_plot.param_practical_label, "", r1)
        form_param = self.makeForm(form_main, self.nelder_plot.param_basic_label, "", r2)
        form_density.pack(padx=5, pady=2, fill=tk.X)
        form_param.pack(padx=5, pady=2, fill=tk.X)
        form_main.pack(side=tk.LEFT, padx=5, pady=5)
        r1.invoke()
        
        #options for center of rotation
        rotPanel = tk.LabelFrame(form_main, bd = 2, text = self.param_type[2], font=self.font_header)
        self.opt_rot_center = StringVar()
        self.opt_rot_center.set(self.nelder_plot.center_config[0])
        n = 0
        for t, v in zip(self.center_options, self.nelder_plot.center_config):
            r, c = divmod(n, 2) 
            n += 1
            Radiobutton(rotPanel, text=t, variable=self.opt_rot_center, value=v, 
                             command=self.selectCenter).grid(row = r, column = c, sticky=tk.W)
        centerxyF = self.makeFormRow(rotPanel, self.field_center)        
        centerxyF.grid(row = 2, column = 1, sticky=tk.W)
        rotPanel.pack(padx=5, pady=2, fill=tk.X)
        self.enableEntry(self.field_center, False)

        # input for number of species
        self.opt_species = IntVar()
        self.opt_species.set(1)
        speciesPanel = tk.LabelFrame(form_main, bd = 2, text = self.param_type[3], font=self.font_header)
        for i in range(3):
            Radiobutton(speciesPanel, text=str(i+1), variable=self.opt_species, value=i+1, 
                             command=self.selectSpeciesOptions).pack(side=tk.LEFT)
        self.opt_alt = IntVar()
        self.cb_alt = Checkbutton(speciesPanel, text=self.label_alt_spoke, variable=self.opt_alt,
            command=self.selectSpeciesOptions, state = 'disabled')
        self.cb_alt.pack(side=tk.LEFT)
        speciesPanel.pack(padx=5, pady=2, fill=tk.X)
                             
        # create summary info panel
        form_summary = self.makeLabel(form_main, self.field_summary, self.param_type[4])
        form_summary.pack(padx=5, pady=2, fill=tk.X)
        form_canvas = tk.Frame(master)
        self.makeToolbar(form_canvas)
        nb = Notebook(form_canvas)
        page1 = tk.Frame(nb)
        page2 = tk.Frame(nb)
        page3 = tk.Frame(nb)
        nb.add(page1, text= self.output_page[0])
        nb.add(page2, text= self.output_page[1])
        nb.add(page3, text= self.output_page[2])
        nb.pack(fill=tk.BOTH, expand=True)
        form_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        # create tabbed panel for the canvas dan data table
        self.plot_canvas = PlotCanvas(page1)
        self.plot_canvas.pack(fill=tk.BOTH, expand=True)
        self.spacingTable = self.makeTable(page2, self.nelder_plot.variable_table_header)
        self.datapointsTable = self.makeTable(page3, self.nelder_plot.datapoint_table_header)
        # get the supported image formats 
        ft = self.plot_canvas.figure.canvas.get_supported_filetypes()
        self.image_file_types = [(v, "." + k) for k, v in ft.items()]
        png = ('Portable Network Graphics', '.png')
        self.image_file_types.remove(png)
        self.image_file_types.insert(0, png)
       
        
    def makeToolbar(self, master):
        """ create toolbar menu """
        toolbar = tk.Frame(master)
        self.img_save = self.makeIcon(toolbar, self.toolbar_text[0], self.toolbar_images[0], self.saveText)
        self.img_image = self.makeIcon(toolbar, self.toolbar_text[1], self.toolbar_images[1], self.saveImage)
        self.img_kml = self.makeIcon(toolbar, self.toolbar_text[2], self.toolbar_images[2], self.saveKML)
        self.img_about = self.makeIcon(toolbar, self.toolbar_text[3], self.toolbar_images[3], self.about, tk.RIGHT)
        toolbar.pack(side=tk.TOP, fill=tk.X)
    
    def makeIcon(self, master, label, image, cmd, side = tk.LEFT):
        """ create icon image """
        resource_package = __name__  # Could be any module/package name
        resource_path = '/'.join(('images', image))  # Do not use os.path.join()
        imagefile = pkg_resources.resource_filename(resource_package, resource_path)
        img=PhotoImage(file=imagefile)
        b = tk.Button(master, image = img, command = cmd, text = label,
                      relief="flat", overrelief="raised", compound="left", padx=5)
        b.pack(side=side, padx=2, pady=2)
        return img
        
    def makeTable(self, master, headers):
        """ create table for data point and variables"""    
        tree = Treeview(columns=headers, show="headings")
        vsb = Scrollbar(orient="vertical", command=tree.yview)
        hsb = Scrollbar(orient="horizontal", command=tree.xview)
        tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        tree.grid(column=0, row=0, sticky='nsew', in_=master)
        vsb.grid(column=1, row=0, sticky='ns', in_=master)
        hsb.grid(column=0, row=1, sticky='ew', in_=master)
        master.grid_columnconfigure(0, weight=1)
        master.grid_rowconfigure(0, weight=1)
        for col in headers:
            tree.heading(col, text=col, command=lambda c=col: self.sortby(tree, c, 0))
            tree.column(col, width= 30, anchor=tk.E)
        tree.column(headers[0], width= 10, anchor=tk.W)
        return tree
    
    def sortby(self, tree, col, descending):
        """ Sort tree contents when a column is clicked on."""
        # grab values to sort
        data = [(tree.set(child, col), child) for child in tree.get_children('')]
        # reorder data
        data.sort(reverse=descending)
        for indx, item in enumerate(data):
            tree.move(item[1], '', indx)
        # switch the heading so that it will sort in the opposite direction
        tree.heading(col, command=lambda col=col: self.sortby(tree, col, int(not descending)))
    
    def makeForm(self, master, fields, label, opt_label = None):
        """ create form for data input """
        root = tk.LabelFrame(master, bd = 2, text = label, font="default 9 bold")
        if not opt_label == None:
            root['labelwidget'] = opt_label 
        for field in fields:
            row = tk.Frame(root)
            lab = Label(row, text=field, width=24, anchor='w')
            ent = tk.Entry(row, width=8, justify='right', disabledforeground='black')
            ent.bind('<FocusOut>', (lambda _: self.entryUpdated(fields)))
            row.pack(side=tk.TOP, fill=tk.X, padx=5, pady=3)
            lab.pack(side=tk.LEFT)
            ent.pack(side=tk.RIGHT, expand=tk.YES, fill=tk.X)
            self.fieldEnt[field] = ent
        return root
    
    def makeFormRow(self, master, fields):
        """ create form for data input in row format """
        row = tk.Frame(master)
        for field in fields:
            lab = Label(row, text=field, width=2, anchor='e')
            ent = tk.Entry(row, width=6, justify='right', disabledforeground='black')
            ent.bind('<FocusOut>', (lambda _: self.entryUpdated(fields)))
            lab.pack(side=tk.LEFT, padx=2 )
            ent.pack(side=tk.LEFT, expand=tk.YES, fill=tk.X)
            self.fieldEnt[field] = ent
        return row
    
    def makeLabel(self, master, fields, label):
        """ create label for information data """
        root = tk.LabelFrame(master, bd = 2, text = label, font="default 9 bold")
        for field in fields:
            row = tk.Frame(root)
            lab = Label(row, text=field + " :", width=24, anchor='w')
            ent = Label(row, width=8, anchor='e')
            row.pack(side=tk.TOP, fill=tk.X, padx=5, pady=2)
            lab.pack(side=tk.LEFT)
            ent.pack(side=tk.RIGHT, expand=tk.YES, fill=tk.X)
            self.fieldEnt[field] = ent
        return root
    
    def selectParams(self):
        """ enable/disable the selected set of parameterization """
        d = (self.opt_param.get() == self.param_type[0])
        self.enableEntry(self.nelder_plot.param_practical_label, d)
        self.enableEntry(self.nelder_plot.param_basic_label, not d)
    
    def selectCenter(self):
        """ enable/disable the center of rotation coordinates input """
        self.enableEntry(self.field_center, 
                         self.opt_rot_center.get() == self.nelder_plot.center_config[4])
        self.entryUpdated(self.field_center)
        
    def selectSpeciesOptions(self):
        """ update for number of species selection """
        self.nelder_plot.number_of_species = self.opt_species.get()
        if self.nelder_plot.number_of_species == 1:
            self.cb_alt['state'] = 'disabled'
        else:
            self.cb_alt['state'] = 'normal'
        self.nelder_plot.is_alternate_spokes = self.opt_alt.get() == 1
        self.entryUpdated([])
    
    def enableEntry(self, fields, is_enable):
        """ enable/disable the input field """
        if is_enable: s = 'normal'  
        else: s ='disabled'
        for field in fields:
            self.fieldEnt.get(field)['state'] = s
            
    def validateValues(self, fields, isAllowNegative = False):
        """ validate the input value """
        for field in fields:
            try:
                f = float(self.fieldEnt.get(field).get())
                if f < 0 and not isAllowNegative:
                    e = self.fieldEnt.get(field)
                    e.delete(0, tk.END)
                    e.insert(0, abs(f))
            except ValueError:
                e = self.fieldEnt.get(field)
                e.delete(0, tk.END)
                e.insert(0, "")
            
    def getValueList(self, fields):
        """ return the value from defined list of input fields """
        vals = []
        for field in fields:
            vals.append(self.getValue(field))
        return tuple(vals)
    
    def getValue(self, fieldLabel):
        """ return the value from defined input field """
        try:
            return float(self.fieldEnt.get(fieldLabel).get())
        except ValueError:
            return 0
    
    def setValue(self, fieldLabel, value):
        """ set the value to the defined input fields """
        e = self.fieldEnt.get(fieldLabel)
        s = e['state']
        if not s == 'normal': e['state'] = 'normal'
        e.delete(0, tk.END)
        e.insert(0, '{:.5g}'.format(value))
        e['state'] = s
    
    def setValueList(self, fields, value_list):
        """ set the list of values to the list of fields """
        for i in range(0, len(fields)):
            self.setValue(fields[i], value_list[i])
    
    def entryUpdated(self, fields):
        """ do something when some parameter entry is updated """
        if fields == self.field_center:
            self.validateValues(fields, True)
        else:
            self.validateValues(fields)
        self.nelder_plot.clearData()
        cp = Point(*self.getValueList(self.field_center))
        cOpt = self.opt_rot_center.get()
        if self.opt_param.get() == self.param_type[0]:
            self.nelder_plot.setPracticalParameters(*self.getValueList(self.nelder_plot.param_practical_label),  cOpt, cp)
            self.setValueList(self.nelder_plot.param_basic_label, self.nelder_plot.getBasicParameters())
            if self.fieldEnt.get(self.nelder_plot.param_practical_label[3]).get() == '': 
                self.setValue(self.nelder_plot.param_practical_label[3], self.nelder_plot.rectangularity)
        elif self.opt_param.get() == self.param_type[1]:
            self.nelder_plot.setBasicParameters(*self.getValueList(self.nelder_plot.param_basic_label), cOpt, cp)
            self.setValueList(self.nelder_plot.param_practical_label, self.nelder_plot.getPracticalParameters())
        if not self.nelder_plot.center_point == None:
            self.setValue("X", self.nelder_plot.center_point.x)
            self.setValue("Y", self.nelder_plot.center_point.y)
        self.summary = self.nelder_plot.getSummaryInfo()    
        for l in self.field_summary:
            self.fieldEnt.get(l)["text"] = '{:.5g}'.format(self.summary[l])
           
        self.plot_canvas.setNelderPlot(self.nelder_plot)
        self.updateDataTable()

    def updateDataTable(self):
        """ update the data list interface """
        self.spacingTable.delete(*self.spacingTable.get_children())
        sp = self.nelder_plot.getVariableTable()
        for item in sp: self.spacingTable.insert('', 'end', values=item)
        
        self.datapointsTable.delete(*self.datapointsTable.get_children())
        dp = self.nelder_plot.getDatapointsTable()
        for item in dp: self.datapointsTable.insert('', 'end', values=item)
        self.sortby(self.datapointsTable, self.nelder_plot.datapoint_table_header[0], 0)
    
    def saveImage(self):
        """ save the plot image """
        answer = simpledialog.askinteger(self.txt_inp_save_image_title, self.txt_inp_save_image,
                                         initialvalue = self.plot_canvas.defaul_dpi,
                                         parent=self.plot_canvas, minvalue=50, maxvalue=1000)
        if answer is None: return
        save_to = filedialog.asksaveasfilename(**dict(
            filetypes=self.image_file_types,
            defaultextension='.png',
            title='Select file'))
        if not save_to: return
        self.plot_canvas.saveImage(save_to, answer)

    class SaveTextOptionDialog(simpledialog.Dialog):
        """ a dialog for asking user input on text saving options """
        output_text_options = ("Nelder parameters and summary", "Plot variables", "Data points")
            
        def body(self, master):
            self.optOuts = []
            self.cbs = []
            for optText in self.output_text_options:
                v = IntVar()
                cb = Checkbutton(master, text=optText, variable = v)
                cb.select()
                cb.pack(anchor=tk.W)
                self.cbs.append(cb)
                self.optOuts.append(v) 
            return self.cbs[0]
            
        def apply(self):
            self.result = [cb.get() for cb in self.optOuts]
        
    def saveText(self):
        """ Saving the data and parameters on text file """
        options = self.SaveTextOptionDialog(self.plot_canvas,
                                            "Select the data to be saved").result
        if options is None: return
        save_to = filedialog.asksaveasfilename(**dict(
            defaultextension='.txt',
            title='Select file'))
        if not save_to: return
        file = open(save_to, "w") 
        if options[0] == 1:
            file.write("# " + self.param_type[0] + "\n\n")
            self.writeFields(file, self.nelder_plot.param_practical_label)
            file.write("\n# " + self.param_type[1] + "\n\n")
            self.writeFields(file, self.nelder_plot.param_basic_label)
            file.write("\n# " + self.param_type[2] + "\n\n")
            self.writeFields(file, self.field_center)
            file.write("\n# " + self.param_type[4] + "\n\n")
            file.write("{} = {:.5g}\n".format(self.param_type[3], self.opt_species.get()))
            for l in self.field_summary:
                file.write("{} = {:.5g}\n".format(l, self.summary[l]))
            file.write("\n")
        if options[1] == 1:
            file.write("# Plot variables\n\n")
            self.writeTab(file, self.nelder_plot.variable_table_header)
            vt = self.nelder_plot.getVariableTable()
            for r in vt: self.writeTab(file, r)
            file.write("\n")
        if options[2] == 1:
            file.write("# Data points\n\n")
            self.writeTab(file, self.nelder_plot.datapoint_table_header)
            dt = self.nelder_plot.getDatapointsTable()
            for r in dt: self.writeTab(file, r)
        file.close()
    
    def writeFields(self, file, fields):
        """ write the data to file """
        for f in fields:
            v = self.getValue(f)
            file.write("{} = {:.5g}\n".format(f, v))
    
    def writeTab(self, file, fields):
        """ write the tab format character """
        for f in fields: file.write("{}\t".format(f))
        file.write("\n")  
    
    class CoordInputDialog(simpledialog.Dialog):
        """ a dialog for latitude and longitude coordinates input"""
            
        def body(self, master):
            Label(master, text="Latitude:").grid(row=0, sticky=tk.W)
            Label(master, text="Longitude:").grid(row=1, sticky=tk.W)
            Label(master, text="Rotation:").grid(row=2, sticky=tk.W)
            self.lat = tk.Entry(master)
            self.long = tk.Entry(master)
            self.rot = tk.Entry(master)
            self.rot.insert(tk.END, "0")
            self.lat.grid(row=0, column=1)
            self.long.grid(row=1, column=1)
            self.rot.grid(row=2, column=1)
            return self.lat # initial focus
        
        def validate(self):
            try:
                latitude = float(self.lat.get())
                longitude = float(self.long.get())
                rotation = float(self.rot.get())
                self.result =  (latitude, longitude, rotation) 
                return 1
            except ValueError:
                messagebox.showwarning(
                    "Bad input",
                    "Illegal values, please try again"
                )
                return 0

    def saveKML(self):
        """ save output to KML format file """
        coord = self.CoordInputDialog(self.plot_canvas, 
                                      title = "Input plot coordinates").result
        if coord == None: return
        save_to = filedialog.asksaveasfilename(**dict(
            defaultextension='.kml',
            title='Select file'))
        if not save_to: return
        file = open(save_to, "w") 
        file.write(self.nelder_plot.getKML(*coord))
        file.close()

    def about(self):
        info = "Nelder Plot Designer\n\n@author: Degi Harja Asmara"
        messagebox.showinfo("Information", info)
        
                
root = Tk()
my_gui = NelderGUI(root)
resource_package = __name__  
resource_path = '/'.join(('images', 'nelder_app.ico')) 
imagefile = pkg_resources.resource_filename(resource_package, resource_path)
root.iconbitmap(imagefile)
root.mainloop()