# managing athena files
from larch.io import create_athena, read_athena, extract_athenagroup

# calculate pre-edge and post edge for normalisation
from larch.xafs import pre_edge
# perform background removal
from larch.xafs import autobk
# calculate fourier transform
from larch.xafs import xftf

# plotting library
import matplotlib.pyplot as plt



 #######################################################
# |              initialise log file                  | #
# V                                                   V #
 #######################################################
def set_logger(log_file):
    logger = logging.getLogger()
    fhandler = logging.FileHandler(filename=log_file, mode='a')
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fhandler.setFormatter(formatter)
    logger.addHandler(fhandler)
    # prevent matplotlib font manager from writing to log
    logging.getLogger('matplotlib.font_manager').disabled = True
    logger.setLevel(logging.DEBUG)


 #######################################################
# |    Build a list of all files matching the given   | #
# V                     pattern                       V #
 #######################################################
def get_files_list(source_dir, f_pattern):
    i_counter = 0
    files_list = []
    for filepath in sorted(source_dir.glob(f_pattern)):
        i_counter += 1
        files_list.append(filepath)
    return files_list


 #######################################################
# |    Read mu data from a text file, use labels      | #
# V         parameter to name columns                 V #
 #######################################################
def read_text(a_file, labels = "energy mu"):
    a_group = read_ascii(a_file, labels=use_labels)
    return a_group


 #######################################################
# |          Create mu on energy plot                 | #
# V                                                   V #
 #######################################################
def plot_mu(xafs_group):
        plt.plot(xafs_group.energy, xafs_group.mu, label=xafs_group.filename) # plot mu in blue
        plt.grid(color='r', linestyle=':', linewidth=1) #show and format grid
        plt.xlabel('Energy (eV)') # label y graph
        plt.ylabel('x$\mu$(E)') # label y axis
        plt.title("$\mu$" + xafs_group.filename)
        plt.legend() # show legend
        return plt
    
 #######################################################
# |    Create pre-edge and post-edge fitting plot     | #
# V                                                   V #
 #######################################################      
def plot_edge_fit(xafs_group):
        plt.plot(xafs_group.energy, xafs_group.pre_edge, 'g', label='pre-edge') # plot pre-edge in green
        plt.plot(xafs_group.energy, xafs_group.post_edge, 'r', label='post-edge')# plot post-edge in red
        plt.plot(xafs_group.energy, xafs_group.mu, 'b', label=xafs_group.filename) # plot mu in blue
        plt.grid(color='r', linestyle=':', linewidth=1) #show and format grid
        plt.xlabel('Energy (eV)') # label y graph
        plt.ylabel('x$\mu$(E)') # label y axis
        plt.title("pre-edge and post_edge fitting to $\mu$")
        plt.legend() # show legend
        return plt
 #######################################################
# |    Create pre-edge and post-edge fitting plot     | #
# V                                                   V #
 #######################################################      
    
# show plot of normalised data
def plot_normalised(xafs_group):
        #plt.plot(xafs_group.energy, xafs_group.pre_edge, 'g', label='pre-edge') # plot pre-edge in green
        #plt.plot(xafs_group.energy, xafs_group.post_edge, 'r', label='post-edge')# plot post-edge in red
        #plt.plot(xafs_group.energy, xafs_group.mu, 'b', label=xafs_group.filename) # plot mu in blue
        plt.plot(xafs_group.energy, xafs_group.flat, label=xafs_group.filename)
        plt.grid(color='r', linestyle=':', linewidth=1) #show and format grid
        plt.xlabel('Energy (eV)') # label y graph
        plt.ylabel('x$\mu$(E)') # label y axis
        plt.title("pre-edge and post_edge fitting to $\mu$")
        plt.legend() # show legend
        return plt

 #######################################################
# |         Read data from Athena project file        | #
# V              returns a project object             V #
 #######################################################
def read_project(project_file):
    return read_athena(project_file)

 #######################################################
# |         Extract groups from Athena project        | #
# V              returns a list of groups             V #
 #######################################################
def get_groups(athena_project):
    athena_groups = []
    group_keys=list(athena_project._athena_groups.keys())
    for group_key in group_keys:
        gr_0 = extract_athenagroup(athena_project._athena_groups[group_key])
        athena_groups.append(gr_0)
    return athena_groups

 #######################################################
# |         Athena recalculates everything so we      | #
# |      need to create a function that calculates    | #
# V               all for each new group              V #
 #######################################################
def calc_with_defaults(xafs_group):
    # calculate mu and normalise with background extraction
    # should let the user specify the colums for i0, it, mu, iR. 
    if not hasattr(xafs_group, 'mu'):
        xafs_group = get_mu(xafs_group)
    # calculate pre-edge and post edge and add them to group
    # need to read parameters for pre-edge before background calculation with  
    # defaul values undo the work of previous step (setting pre-edge limits).
    pre_edge(xafs_group, pre1=xafs_group.athena_params.bkg.pre1, pre2=xafs_group.athena_params.bkg.pre2)
    #pre_edge(xafs_group)
    # perform background removal
    autobk(xafs_group) # using defaults so no additional parameters are passed
    # calculate fourier transform
    xftf(xafs_group)#, kweight=0.5, kmin=3.0, kmax=12.871, dk=1, kwindow='Hanning')
    return xafs_group

 #######################################################
# |       The code for plotting Nmu vs E repeats      | #
# |   so it is useful to have a plotting function     | #
# V            to reduce duplicated code              V #
 #######################################################
def plot_normalised(xafs_group):
        plt.plot(xafs_group.energy, xafs_group.mu, label=xafs_group.filename) # plot mu in blue
        plt.grid(color='r', linestyle=':', linewidth=1) #show and format grid
        plt.xlabel('Energy (eV)') # label y graph
        plt.ylabel('x$\mu$(E)') # label y axis
        plt.title("pre-edge and post_edge fitting to $\mu$")
        plt.legend() # show legend
        return plt
    
 #######################################################
# |       The code for plotting Nmu vs E repeats      | #
# |   so it is useful to have a plotting function     | #
# V            to reduce duplicated code              V #
 #######################################################
# plot flat normalised mu vs Energy for selected groups
def plot_Nxmu_E(athena_project, group_keys, group_names,
                            title = "Normalised $\mu$ vs E", xlimits = None,
                            ylimits = None):    
    # plot mu vs flat normalised mu for selected groups
    for group_key in group_keys:
        gr_0 = extract_athenagroup(athena_project._athena_groups[group_key])
        # recalculate normalisation
        calc_with_defaults(gr_0)
        plt.plot(gr_0.energy, gr_0.flat, label=group_names[group_key])

    # set plot format
    plt.xlabel("Energy")
    plt.ylabel("Normalised $\mu$" )
    plt.title(title)
    plt.grid(linestyle=':', linewidth=1) #show and format grid
    if xlimits != None:
        plt.xlim(xlimits[0],xlimits[1])
    if ylimits != None:
        plt.ylim(ylimits[0],ylimits[1])
    plt.legend()
    return plt, gr_0

    
    