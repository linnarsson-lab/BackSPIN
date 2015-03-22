#!/usr/bin/env python
from __future__ import division
#from matplotlib.pylab import *
from numpy import *
import getopt
import sys
import os
import time
import pandas as pd
from Cef_tools import CEF_obj


def _calc_weights_matrix(mat_size, wid):
    '''Calculate Weight Matrix
    Parameters
    ----------
    mat_size: int
        dimension of the distance matrix
    wid: int
        parameter that controls the width of the neighbourood
    Returns
    -------
    weights_mat: 2-D array
        the weights matrix to multiply with the distance matrix

    '''
    #calculate square distance from the diagonal
    sqd = (arange(1,mat_size+1)[newaxis,:] - arange(1,mat_size+1)[:,newaxis])**2
    #make the distance relative to the mat_size
    norm_sqd = sqd/wid
    #evaluate a normal pdf
    weights_mat = exp(-norm_sqd/mat_size)
    #avoid useless precision that would slow down the matrix multiplication
    weights_mat -= 1e-6
    weights_mat[weights_mat<0] = 0
    #normalize row and column sum
    weights_mat /= sum(weights_mat,0)[newaxis,:]
    weights_mat /= sum(weights_mat,1)[:, newaxis]
    #fix asimmetries
    weights_mat = (weights_mat + weights_mat.T) / 2.
    return weights_mat


def _sort_neighbourhood( dist_matrix, wid ):
    '''Perform a single iteration of SPIN
    Parameters
    ----------
    dist_matrix: 2-D array
        distance matrix
    wid: int
        parameter that controls the width of the neighbourood
    Returns
    -------
    sorted_ind: 1-D array
        indexes that order the matrix

    '''
    assert wid > 0, 'Parameter wid < 0 is not allowed'
    mat_size = dist_matrix.shape[0]
    #assert mat_size>2, 'Matrix is too small to be sorted'
    weights_mat = _calc_weights_matrix(mat_size, wid)
    #Calculate the dot product (can be very slow for big mat_size)
    mismatch_score = dot(dist_matrix, weights_mat)
    energy, target_permutation = mismatch_score.min(1), mismatch_score.argmin(1)
    max_energy = max(energy)
    #Avoid points that have the same target_permutation value
    sort_score = target_permutation - 0.1 * sign( (mat_size/2 - target_permutation) ) * energy/max_energy
    #sort_score = target_permutation - 0.1 * sign( 1-2*(int(1000*energy/max_energy) % 2) ) * energy/max_energy # Alternative
    # Sorting the matrix
    sorted_ind = sort_score.argsort(0)[::-1]
    return sorted_ind


def sort_mat_by_neighborhood(dist_matrix, wid, times):
    '''Perform several iterations of SPIN using a fixed wid parameter
    Parameters
    ----------
    dist_matrix: 2-D array
        distance matrix
    wid: int
        parameter that controls the width of the neighbourood
    times: int
        number of repetitions
    verbose: bool
        print the progress
    Returns
    -------
    indexes: 1-D array
        indexes that order the matrix

    '''
    # original indexes
    indexes = arange(dist_matrix.shape[0])
    for i in range(times):
        #sort the sitance matrix according the previous iteration
        tmpmat = dist_matrix[indexes,:] 
        tmpmat = tmpmat[:,indexes]
        sorted_ind = _sort_neighbourhood(tmpmat, wid);
        #resort the original indexes
        indexes = indexes[sorted_ind]
    return indexes


def _generate_widlist(data, axis=1, step=0.6):
    '''Generate a list of wid parameters to execute sort_mat_by_neighborhood
    Parameters
    ----------
    data: 2-D array
        the data matrix
    axis: int
        the axis to take in consideration
    step: float
        the increment between two successive wid parameters
    Returns
    -------
    wid_list: list of int
        list of wid parameters to run SPIN

    '''
    max_wid = data.shape[axis]*0.6
    new_wid = 1
    wid_list = []
    while new_wid < (1+step)*max_wid:
        wid_list.append( new_wid )
        new_wid = int(ceil( new_wid + new_wid*(step) +1))
    return wid_list[::-1]



def SPIN(dt, widlist=[10,1], iters=30, axis='both', verbose=False):
    """Run the original SPIN algorithm
    Parameters
    ----------
    dt: 2-D array
        the data matrix
    widlist: float or list of int
        If float is passed, it is used as step parameted of _generate_widlist, 
        and widlist is generated to run SPIN.
        If list is passed it is used directly to run SPIN.
    iters: int
        number of repetitions for every wid in widlist
    axis: int
        the axis to take in consideration (must be 0, 1 or 'both')
    step: float
        the increment between two successive wid parameters
    Returns
    -------
    indexes: 1-D array (if axis in [0,1]) or tuple of 1-D array (if axis = 'both')
        indexes that sort the data matrix
    Notes
    -----
    Typical usage
    sorted_dt0 = SPIN(dt, iters=30, axis=0)
    sorted_dt1 = SPIN(dt, iters=30, axis=1)
    dt = dt[sorted_dt0,:]
    dt = dt[:,sorted_dt1]
    """
    IXc = arange(dt.shape[1])
    IXr = arange(dt.shape[0])
    assert axis in ['both', 0,1], 'axis must be 0, 1 or \'both\' '
    #Sort both axis
    if axis == 'both':
        CCc = 1 - corrcoef(dt.T) 
        CCr = 1 - corrcoef(dt) 
        if type(widlist) != list:
            widlist_r = _generate_widlist(dt, axis=0, step=widlist)
            widlist_c = _generate_widlist(dt, axis=1, step=widlist)
        if verbose:
                print '\nSorting genes.'
                print 'Neighbourood=',
        for wid in widlist_r:
            if verbose:
                print ('%i, ' % wid),
                sys.stdout.flush()
            INDr = sort_mat_by_neighborhood(CCr, wid, iters)
            CCr = CCr[INDr,:][:,INDr]
            IXr = IXr[INDr]
        if verbose:
                print '\nSorting cells.'
                print 'Neighbourood=',
        for wid in widlist_c:
            if verbose:
                print ('%i, ' % wid),
                sys.stdout.flush()
            INDc = sort_mat_by_neighborhood(CCc, wid, iters)
            CCc = CCc[:,INDc][INDc,:]
            IXc= IXc[INDc]
        return IXr, IXc
    #Sort rows
    elif axis == 0:
        CCr = 1 - corrcoef(dt)
        if type(widlist) != list:
            widlist = _generate_widlist(dt, axis=0, step=widlist)
        if verbose:
                print '\nSorting genes.\nNeighbourood=',
        for wid in widlist:
            if verbose:
                print '%i, ' % wid,
                sys.stdout.flush()
            INDr = sort_mat_by_neighborhood(CCr, wid, iters)
            CCr = CCr[INDr,:][:,INDr]
            IXr = IXr[INDr]
        return IXr
    #Sort columns
    elif axis == 1:
        CCc = 1 - corrcoef(dt.T)
        if type(widlist) != list:
            widlist = _generate_widlist(dt, axis=1, step=widlist)
        if verbose:
            print '\nSorting cells.\nNeighbourood=',
        for wid in widlist:
            if verbose:
                print '%i, ' % wid,
                sys.stdout.flush()
            INDc = sort_mat_by_neighborhood(CCc, wid, iters)
            CCc = CCc[:,INDc][INDc,:]
            IXc = IXc[INDc]
        return IXc


def backSPIN(data, numLevels=2, first_run_iters=10, first_run_step=0.05, runs_iters=8 ,runs_step=0.25,\
    split_limit_g=2, split_limit_c=2, stop_const = 1.15, low_thrs=0.2, verbose=False):
    '''Run the backSPIN algorithm
    Parameters
    ----------
    data: 2-D array
        the data matrix, rows should be genes and columns single cells/samples
    numLevels: int
        the number of splits that will be tried
    first_run_iters: float
        the iterations of the preparatory SPIN
    first_run_step: float
        the step parameter passed to _generate_widlist for the preparatory SPIN
    runs_iters: int
        the iterations parameter passed to the _divide_to_2and_resort.
        influences all the SPIN iterations except the first
    runs_step: float
        the step parameter passed to the _divide_to_2and_resort.
        influences all the SPIN iterations except the first
    wid: float
        the wid of every iteration of the splitting and resorting
    split_limit_g: int
        If the number of specific genes in a subgroup is smaller than this number
         splitting of that subgrup is not allowed
    split_limit_c: int
        If the number cells in a subgroup is smaller than this number splitting of
        that subgrup is not allowed
    stop_const: float
        minimum score that a breaking point has to reach to be suitable for splitting
    low_thrs: float
        genes with average lower than this threshold are assigned to either of the 
        splitting group reling on genes that are higly correlated with them

    Returns
    -------
    results: Result object
        The results object contain the following attributes
        dataout_sorted:

        genes_order:

        cells_order:

        genes_gr_level:

        cells_gr_level:

        cells_gr_level_sc:

        genes_bor_level:

        cells_bor_level:


    Notes
    -----
    Typical usage
    
    '''
    assert numLevels>0, '0 is not an available depth for backSPIN, use SPIN instead'
    #initialize some varaibles
    genes_bor_level = [[] for i in range(numLevels)] 
    cells_bor_level = [[] for i in range(numLevels)] 
    N,M = data.shape
    genes_order = arange(N)
    cells_order = arange(M)
    genes_gr_level = zeros((N,numLevels+1))
    cells_gr_level = zeros((M,numLevels+1))
    cells_gr_level_sc = zeros((M,numLevels+1))

    # Do a Preparatory SPIN on cells
    if verbose:
        print '\nPreparatory SPIN'
    ix1 = SPIN(data, widlist=_generate_widlist(data, axis=1, step=first_run_step), iters=first_run_iters, axis=1, verbose=verbose)
    cells_order = cells_order[ix1]

    #For every level of depth DO:
    for i in range(numLevels): 
        k=0 # initialize group id counter
        # For every group generated at the parent level DO:
        for j in range( len( set(cells_gr_level[:,i]) ) ): 
            # Extract the a data matrix of the genes at that level
            g_settmp = nonzero(genes_gr_level[:,i]==j)[0] #indexes of genes in the level j
            c_settmp = nonzero(cells_gr_level[:,i]==j)[0] #indexes of cells in the level j
            datatmp = data[ ix_(genes_order[g_settmp], cells_order[c_settmp]) ]
            # If we are not below the splitting limit for both genes and cells DO:
            if (len(g_settmp)>split_limit_g) & (len(c_settmp)>split_limit_c): 
                # Split and SPINsort the two halves
                divided = _divide_to_2and_resort(datatmp, wid=runs_step, iters_spin=runs_iters,\
                    stop_const=stop_const, low_thrs=low_thrs, verbose=verbose)
                # _divide_to_2and_resort retruns an empty array in gr2 if the splitting condition was not satisfied
                if divided:
                    sorted_data_resort1, genes_resort1, cells_resort1,\
                    gr1, gr2, genesgr1, genesgr2, score1, score2 = divided
                    # Resort from the previous level
                    genes_order[g_settmp] = genes_order[g_settmp[genes_resort1]]
                    cells_order[c_settmp] = cells_order[c_settmp[cells_resort1]]
                    # Assign a numerical identifier to the groups
                    genes_gr_level[g_settmp[genesgr1],i+1] = k
                    genes_gr_level[g_settmp[genesgr2],i+1] = k+1
                    cells_gr_level[c_settmp[gr1],i+1] = k
                    cells_gr_level[c_settmp[gr2],i+1] = k+1
                    # Not really clear what sc is
                    cells_gr_level_sc[c_settmp[gr1],i+1] = score1
                    cells_gr_level_sc[c_settmp[gr2],i+1] = score2
                    # Augment the counter of 2 becouse two groups were generated from one
                    k = k+2
                else:
                    # The split is not convenient, keep everithing the same
                    genes_gr_level[g_settmp,i+1] = k
                    cells_gr_level[c_settmp,i+1] = k
                    cells_gr_level_sc[c_settmp,i+1] = cells_gr_level_sc[c_settmp,i]
                    # Augment of 1 becouse no new group was generated
                    k = k+1
            else:
                # The split is not convenient, keep everithing the same
                genes_gr_level[g_settmp,i+1] = k
                cells_gr_level[c_settmp,i+1] = k
                cells_gr_level_sc[c_settmp,i+1] = cells_gr_level_sc[c_settmp,i]
                # Augment of 1 becouse no new group was generated
                k = k+1
        
        # Find boundaries
        genes_bor_level[i] = r_[0, nonzero(diff(genes_gr_level[:,i+1])>0)[0]+1, data.shape[0] ]
        cells_bor_level[i] = r_[0, nonzero(diff(cells_gr_level[:,i+1])>0)[0]+1, data.shape[1] ]

    #dataout_sorted = data[ ix_(genes_order,cells_order) ]

    class Results:
        pass

    results = Results()
    results.genes_order = genes_order
    results.cells_order = cells_order
    results.genes_gr_level = genes_gr_level
    results.cells_gr_level = cells_gr_level
    results.cells_gr_level_sc = cells_gr_level_sc
    results.genes_bor_level = genes_bor_level
    results.cells_bor_level = cells_bor_level

    return results
    
    

def _divide_to_2and_resort(sorted_data, wid, iters_spin=8, stop_const = 1.15, low_thrs=0.2 ,verbose=False):
    '''Core function of backSPIN: split the datamatrix in two and resort the two halves

    Parameters
    ----------
    sorted_data: 2-D array
        the data matrix, rows should be genes and columns single cells/samples
    wid: float
        wid parameter to give to widlist parameter of th SPIN fucntion
    stop_const: float
        minimum score that a breaking point has to reach to be suitable for splitting
    low_thrs: float
        if the difference between the average expression of two groups is lower than threshold the algorythm 
        uses higly correlated gens to assign the gene to one of the two groups
    verbose: bool
        information about the split is printed

    Returns
    -------
    '''
    
    # Calculate correlation matrix for cells and genes
    Rcells = corrcoef(sorted_data.T)
    Rgenes = corrcoef(sorted_data)
    # Look for the optimal breaking point
    N = Rcells.shape[0]
    score = zeros(N)
    for i in range(2,N-2):
        tmp1 = Rcells[:i,:i]
        tmp2 = Rcells[i:,i:]
        score[i] = (sum(tmp1)+sum(tmp2)) / float(i**2 + (N-i)**2)
    
    breakp1 = argmax(score)
    score1 = Rcells[:breakp1,:breakp1]
    score1 = triu(score1)
    score1 = mean( score1[score1 != 0] ) # WHY ???????????????????????????????
    score2 = Rcells[breakp1:,breakp1:]
    score2 = triu(score2)
    score2 = mean( score2[score2 != 0] ) # WHY ???????????????????????????????
    avg_tot = triu(Rcells)
    avg_tot = mean( avg_tot[avg_tot != 0] ) # WHY ???????????????????????????????

    # If it is convenient to break
    if (max([score1,score2])/avg_tot) > stop_const:
        # Divide in two groups
        gr1 = arange(N)[:breakp1]
        gr2 = arange(N)[breakp1:]
        # and assign the genes into the two groups on the basis of the mean
        mean_gr1 = mean( sorted_data[:,gr1],1 )
        mean_gr2 = mean( sorted_data[:,gr2],1 )
        d = abs( mean_gr1 - mean_gr2 )
        # Deal with low variance genes using correlation with other genes to assign them to one of the groups
        # This is  considered reliable if the original group contained more than 20 genes 
        if len(d) > 20:
            # For every difference lower than a threshold 
            for i in range(len(d)): 
                if d[i] < low_thrs:
                    IN = Rgenes[i,:] > percentile(Rgenes[i,:], 100 - 100*(20/len(d)))
                    mean_gr1[i] = sorted_data[ix_(IN,gr1)].sum(0).mean() #the mean of the sum of the columns
                    mean_gr2[i] = sorted_data[ix_(IN,gr2)].sum(0).mean()
                    
        bigger_gr1 = (mean_gr1 - mean_gr2) > 0 # boolean vector
        
        # Avoid group of cells with no genes to be formed by adding the highest 
        # expressed gene to the gene-empty group 
        genesgr1 = nonzero(bigger_gr1)[0]
        genesgr2 = nonzero(~bigger_gr1)[0]
        if size(genesgr1) == 0:
            IN = argmax(mean_gr1)
            genesgr1 = array([IN])
            genesgr2 = setdiff1d(genesgr2, IN)
        elif size(genesgr2) == 0:
            IN = argmax(mean_gr2)
            genesgr2 = array([IN])
            genesgr1 = setdiff1d(genesgr1, IN)
        
        if verbose:
            print '\nSplitting (%i, %i) ' %  sorted_data.shape
            print 'in (%i,%i) ' % (genesgr1.shape[0],gr1.shape[0])
            print 'and (%i,%i)' % (genesgr2.shape[0],gr2.shape[0]),
            sys.stdout.flush()

        # Data of group1
        datagr1 = sorted_data[ix_(genesgr1,gr1)]
        # zero center
        datagr1 = datagr1 - datagr1.mean(1)[:,newaxis]
        # Resort group1
        if min( datagr1.shape ) > 1:
            genesorder1,cellorder1 = SPIN(datagr1, widlist=wid, iters=iters_spin, axis='both', verbose=verbose)
        elif len(genesgr1) == 1:
            genesorder1 = 0
            cellorder1 = argsort( datagr1[0,:] )
        elif len(gr1) == 1:
            cellorder1 = 0
            genesorder1 = argsort( datagr1[:,0] )

        # Data of group2
        datagr2 = sorted_data[ix_(genesgr2,gr2)]
        # zero center
        datagr2 = datagr2 - datagr2.mean(1)[:,newaxis]
        # Resort group2
        if min( datagr2.shape )>1:
            genesorder2, cellorder2 = SPIN(datagr2, widlist=wid, iters=iters_spin, axis='both',verbose=verbose)
        elif len(genesgr2) == 1:
            genesorder2 = 0
            cellorder2 = argsort(datagr2[0,:])
        elif len(gr2) == 1:
            cellorder2 = 0
            genesorder2 = argsort(datagr2[:,0])
        
        # contcatenate cells and genes indexes
        genes_resort1 = r_[genesgr1[genesorder1], genesgr2[genesorder2] ]
        cells_resort1 = r_[gr1[cellorder1], gr2[cellorder2] ]
        genesgr1 = arange(len(genesgr1))
        genesgr2 = arange(len(genesgr1), len(sorted_data[:,0]))
        # resort
        sorted_data_resort1 = sorted_data[ix_(genes_resort1,cells_resort1)]

        return sorted_data_resort1, genes_resort1, cells_resort1, gr1, gr2, genesgr1, genesgr2, score1, score2

    else:
        if verbose:
            print 'No splitting score was : %.4f' % (max([score1,score2])/avg_tot)
        return False


def usage_quick():

    message ='''usage: backSPIN [-hbv] [-i inputfile] [-o outputfolder] [-d int] [-t int] [-s float] [-T int] [-S float] [-g int] [-c int] [-k float] [-r float]
    manual: backSPIN -h
    '''
    print message

def usage():

    message='''
       backSPIN commandline tool
       -------------------------

       The options are as follows:

       -i [inputfile]
       --input=[inputfile]
              Path of the tab delimited file.
              Rows should be genes and columns single cells/samples

       -o [outputfolder]
       --output=[outputfolder]
              The name of the folder where the output will be written (output will be a 
              file named results_ddMMyyhhmmss.cef)

       -d [int]
              Depth/Number of levels: The number of nested splits that will be tried by the algorythm
       -t [int]
              Number of the iterations used in the preparatory SPIN.
              Defaults to 10
       -s [float]
              Controls the decrease rate of the wid parameter used in the preparatory SPIN.
              Smaller values will increase the number of SPIN iterations and result in higher 
              precision in the first step but longer execution time.
              Defaults to 0.05
       -T [int]
              Number of the iterations used for every wid parameter.
              Does not apply on the first run (use -t instead)
              Defaults to 8
       -S [float]
              Controls the decrease rate of the wid parameter.
              Smaller values will increase the number of SPIN iterations and result in higher 
              precision but longer execution time.
              Does not apply on the first run (use -s isntead)
              Defaults to 0.25
       -g [int]
              Minimal number of genes that a group must contain for splitting to be allowed.
              Defaults to 2
       -c [int]
              Minimal number of cells that a group must contain for splitting to be allowed.
              Defaults to 2
       -k [float]
              Minimum score that a breaking point has to reach to be suitable for splitting.
              Defaults to 1.15
       -r [float]
              If the difference between the average expression of two groups is lower than threshold the algorythm 
              uses higly correlated gens to assign the gene to one of the two groups
              Defaults to 0.2
       -b [[axisvalue]]
              Run normal SPIN instead of backSPIN.
              Normal spin accepts the parameters -T -S
              optionally one can pass an axis value 0 to only sort genes (rows), 1 to only sort cells (columns)
       -v  
              Verbose. Print extra details of what is happening to the stdoutput 

    '''

    print message



if __name__ == '__main__':
    print 
    #defaults arguments
    input_path = None
    outputfolder_path = None
    numLevels=2 # -d
    first_run_iters=10 # -t
    first_run_step=0.1 # -s
    runs_iters=8 # -T
    runs_step=0.3 # -S
    split_limit_g=2 # -g
    split_limit_c=2 # -c
    stop_const = 1.15 # -k
    low_thrs=0.2 # -r
    normal_spin = False #-b
    normal_spin_axis = 'both'
    verbose=False # -v

    optlist, args = getopt.getopt(sys.argv[1:], "hvi:o:d:t:s:T:S:g:c:k:r:b:", ["help", "input=","output="])

    if optlist== [] and args == []:
        usage_quick()
        sys.exit()
    for opt, a in optlist:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ('-i', '--input'):
            input_path = a
        elif opt in ("-o", "--output"):
            if os.path.exists(a):
                outputfolder_path = a
            else:
                'Output folder %s was not found' % os.path.abspath(a)
        elif opt == '-d':
            numLevels = int(a)
        elif opt == '-t':
            first_run_iters = int(a)
        elif opt == '-s':
            first_run_step = float(a)
        elif opt == '-T':
            runs_iters = int(a)
        elif opt == '-S':
            runs_step = float(a)
        elif opt == '-g':
            split_limit_g = int(a)
        elif opt == '-c':
            split_limit_c = int(a)
        elif opt == '-k':
            stop_const = float(a)
        elif opt == '-r':
            low_thrs = float(a)
        elif opt == '-v':
            verbose = True
        elif opt == '-b':
            normal_spin = True
            normal_spin_axis = a
        else:
            assert False, "%s option is not supported" % opt

    if input_path == None:
        print 'No input file was provided.\nYou need to specify an input file\n(e.g. backSPIN -i path/to/your/file/foo.txt)\n'
        sys.exit()
    if outputfolder_path == None:
        outputfolder_path = os.getcwd()
    
    outfiles_path = os.path.join( outputfolder_path, 'results_%s.cef' % time.strftime('%d%m%y%H%M%S') )

    try:
        if verbose:
            print 'Loading file.'
        input_cef = CEF_obj()
        input_cef.readCEF(input_path)
        data = array(input_cef.matrix)
        data = log2(data+1)
        data = data - data.mean(1)[:,newaxis]
        if data.shape[0] <= 3 and data.shape[1] <= 3:
            print 'Inputfile is not correctly formatted.\n'
            sys.exit()
    except:
        print 'Error occurred in parsing the input file.'
        print 'Plase check that your input file is a correctly formatted tab separated file.\n'
        sys.exit()

    if normal_spin == False:

        print 'backSPIN started\n----------------\n'
        print 'Input file seleceted:\n%s\n' % input_path
        print 'Output files will be saved in:\n%s\n' % outfiles_path
        print 'numLevels: %i\nfirst_run_iters: %i\nfirst_run_step: %.3f\nruns_iters: %i\nruns_step: %.3f\nsplit_limit_g: %i\nsplit_limit_c: %i\nstop_const: %.3f\nlow_thrs: %.3f\n' % (numLevels, first_run_iters, first_run_step, runs_iters,\
            runs_step, split_limit_g, split_limit_c, stop_const, low_thrs)

        results = backSPIN(data, numLevels, first_run_iters, first_run_step, runs_iters, runs_step,\
            split_limit_g, split_limit_c, stop_const, low_thrs, verbose)

        sys.stdout.flush()
        print '\nWriting output files.\n'

        output_cef = CEF_obj()

        for h_name, h_val in zip( input_cef.header_names, input_cef.header_values):
            output_cef.add_header(h_name, h_val )
        for c_name, c_val in zip( input_cef.col_attr_names, input_cef.col_attr_values):
            output_cef.add_col_attr(c_name, array(c_val)[results.cells_order])
        for r_name, r_val in zip( input_cef.row_attr_names, input_cef.row_attr_values):
            output_cef.add_row_attr(r_name, array(r_val)[results.genes_order])

        for level, groups in enumerate( results.cells_gr_level.T ):
            output_cef.add_row_attr('Groups_lvl%i' % level, [int(el) for el in groups])
        for level, groups in enumerate( results.cells_gr_level.T ):
            output_cef.add_col_attr('Groups_lvl%i' % level, [int(el) for el in groups])

        output_cef.set_matrix(array(input_cef.matrix)[results.genes_order,:][:,results.cells_order])

        output_cef.writeCEF( outfiles_path )

    else:

        print 'normal SPIN started\n----------------\n'
        print 'Input file seleceted:\n%s\n' % input_path
        print 'Output files will be saved in:\n%s\n' % outfiles_path

        results = SPIN(dt, widlist=runs_step, iters=runs_iters, axis=normal_spin_axis, verbose=verbose)

        print '\nWriting output files.\n'

        output_cef = CEF_obj()

        for h_name, h_val in zip( input_cef.header_names, input_cef.header_values):
            output_cef.add_header(h_name, h_val )

        if normal_spin_axis == 'both':
            for c_name, c_val in zip( input_cef.col_attr_names, input_cef.col_attr_values):
                output_cef.add_col_attr(c_name, array(c_val)[results[1]])
            for r_name, r_val in zip( input_cef.row_attr_names, input_cef.row_attr_values):
                output_cef.add_row_attr(r_name, array(r_val)[results[0]])

        if normal_spin_axis == 0:
            for r_name, r_val in zip( input_cef.row_attr_names, input_cef.row_attr_values):
                output_cef.add_row_attr(r_name, array(r_val)[results])

        if normal_spin_axis == 1:
            for c_name, c_val in zip( input_cef.col_attr_names, input_cef.col_attr_values):
                output_cef.add_col_attr(c_name, array(c_val)[results])

        output_cef.set_matrix(data[results[0],:][:,results[1]])

        output_cef.writeCEF( outfiles_path )







