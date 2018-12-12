import numpy as np
from numpy import genfromtxt, log, array, exp, save
from networkx import DiGraph, shortest_path_length, adjacency_matrix, shortest_path, all_simple_paths, strongly_connected_component_subgraphs
from scipy.sparse.linalg import inv
from tqdm import tqdm
from scipy.sparse import diags, eye

###############################################################################

class EffectiveDistances:
    def __init__(self, fname="", verbose=True, **kwargs):
        assert isinstance(fname,str)
        self.graph = None
        self.nodes = 0
        self.dominant_path_distance = None
        self.multiple_path_distance = None
        self.random_walk_distance = None
        
        if len(fname) > 0:
            self.graph = self.load(fname, verbose=verbose, **kwargs)
    
    def __str__():
        pass
    
    def load(self,fname, verbose=True, **kwargs):
        """
        Load a data file. The expected data format is three columns 
        (comma seperated by default) with source, target, flux.
        No header should be included and the node IDs have to run contuously 
        from 0 to Number_of_nodes-1.

        Parameters
        ----------
            fname : str
                Path to the file
            
            verbose : bool
                Print information about the data. True by Default
                
            kwargs : dict
                Default parameters can be changed here. Supported key words are
                    dtype     : float (default)
                    delimiter : ","   (default)
        
            return_graph : bool
                If True, the graph is returned (False by default).
                
        Returns:
        --------
            The graph is saved internally in self.graph.
                
        """
        delimiter = kwargs["delimiter"]      if "delimiter"      in kwargs.keys() else " "
        
        data = np.genfromtxt(fname, delimiter=delimiter, dtype=int, unpack=False)
        source, target = data[:,0], data[:,1]
        if data.shape[1] > 2:
            flux = data[:,2]
        else:
            flux = np.ones_like(source)
        nodes  = set(source) | set(target)
        self.nodes = len(nodes)
        lines  = len(flux)
        if set(range(self.nodes)) != nodes:
            new_node_ID = {old:new for new,old in enumerate(nodes)}
            map_new_node_ID = np.vectorize(new_node_ID.__getitem__)
            source = map_new_node_ID(source)
            target = map_new_node_ID(target)
            if verbose:
                print "\nThe node IDs have to run continuously from 0 to Number_of_nodes-1."
                print "Node IDs have been changed according to the requirement.\n-----------------------------------\n"
                
        
            print 'Lines: ',lines , ', Nodes: ', self.nodes
            print '-----------------------------------\nData Structure:\n\nsource,    target,    weight \n'
            for ii in range(7):            
                print "%i,       %i,       %1.2e" %(source[ii], target[ii], flux[ii])
            print '-----------------------------------\n'
        
        
        G = DiGraph()         # Empty, directed Graph
        G.add_nodes_from(range(self.nodes))
        for ii in xrange(lines):
            u, v, w = int(source[ii]), int(target[ii]), float(flux[ii])
            if u != v: # ignore self loops
                assert not G.has_edge(u,v), "Edge appeared twice - not supported"                    
                G.add_edge(u,v,weight=w)
            else:
                if verbose:
                    print "ignore self loop at node", u
        
        symmetric = True
        for s,t,w in G.edges(data=True):
            w1 = G[s][t]["weight"]
            try:
                w2 = G[t][s]["weight"]
            except KeyError:
                symmetric = False
                G.add_edge(t,s,weight=w1)
                w2 = w1
            if w1 != w2:
                symmetric = False
                G[s][t]["weight"] += G[t][s]["weight"]
                G[s][t]["weight"] /= 2
                G[t][s]["weight"]  = G[s][t]["weight"]
        if verbose:
            if not symmetric:
                print "The network has been symmetricised."
        
        
        ccs = strongly_connected_component_subgraphs(G)
        ccs = sorted(ccs, key=len, reverse=True)
        
        G_GSCC = ccs[0]
        if G_GSCC.number_of_nodes() != G.number_of_nodes():
            G = G_GSCC
            if verbose:
                print "\n--------------------------------------------------------------------------"
                print "The network has been restricted to the giant strongly connected component."
        self.nodes = G.number_of_nodes()
        
        
        
        
        for u, v, data in G.edges(data=True):
            weight = G.out_degree(u,weight='weight')
            data['transition_rate'] = 1.*data['weight']/weight
        
        
        for u, v, data in G.edges(data=True):
            data['effective_distance'] = 1. - log(data['transition_rate'])
        
        if verbose:
            print "\n--------------------------------------------------------------------------"
            print "\nnode ID, out-weight,   normalized out-weight,  sum of effective distances \n "
            for ii in range(7):
                out_edges = G.out_edges(ii, data=True)
                out_weight, effective_distance, transition_rate = 0, 0, 0
                for u, v, data in out_edges:
                    out_weight          += data["weight"]
                    effective_distance  += data["effective_distance"]
                    transition_rate     += data["transition_rate"]
                print "  %i       %1.2e           %2.3f                 %1.2e " %(ii,out_weight, transition_rate, effective_distance)
            print "\n ... graph is saved in self.graph"
        return G
        
    ###############################################################################
    
    def get_dominant_path_distance(self, source=None, target=None, parameter=1, saveto=""):
        """
        Compute the (dominant path) effective distance:
        Gautreau, A, Barrat, A. and Barthelemy, M., Global disease spread: 
        Statistics and estimation of arrival times, Journal of Theoretical Biology 251, 3, 509 (2008)

        Parameters
        ----------
             source : int or None
                If source is None, the distances from all nodes to the target is calculated
                Otherwise the integer has to correspond to a node index
            
            target : int or None
                If target is None, the distances from the source to all other nodes is calculated
                Otherwise the integer has to correspond to a node index
                
            parameter : float
                compound parameter which includes the infection and recovery rate alpha and beta, respectively, 
                the mobility rate kappa and the Euler-Mascheroni constant lambda:
                    log[ (alpha-beta)/kappa - lambda ]
        
            saveto : string
                If empty, the result is saved internally in self.dominant_path_distance           
                
        Returns:
        --------
            dominant_path_distance : ndarray or float
                If source and target are specified, a float value is returned that specifies the distance.
                
                If either source or target is None a numpy array is returned.
                The position corresponds to the node ID. 
                shape = (Nnodes,)
                
                If both are None a numpy array is returned.
                Each row corresponds to the node ID. 
                shape = (Nnodes,Nnodes)
        """    
        
        assert (isinstance(parameter,float) or isinstance(parameter,int)) and parameter > 0
        assert isinstance(saveto,str)
        assert self.graph != None, "Load graph first."        
        
        if parameter != 1:
            for n,m,data in self.graph.edges(data=True):
                data['effective_distance'] = parameter - log(data['transition_rate'])
        
        DPED_dic = shortest_path_length(self.graph, source=source, target=target, weight="effective_distance")
        
        if source is None and target is None:
            DPED = np.array([DPED_dic[s].values() for s in xrange(self.nodes)]).transpose()
        elif (source is None) != (target is None):
            DPED = array(DPED_dic.values())
        else:
            DPED = DPED_dic
        if saveto is not "":
            save( saveto, DPED )
        else:
            return DPED
        
        
    ###############################################################################
    
    def get_random_walk_distance(self, source=None, target=None, parameter=1, saveto=""):
        """
        Compute the random walk effective distance:
        F. Iannelli, A. Koher, P. Hoevel, I.M. Sokolov (in preparation)
        
        Parameters
        ----------
             source : int or None
                If source is None, the distances from all nodes to the target is calculated
                Otherwise the integer has to correspond to a node index
            
            target : int or None
                If target is None, the distances from the source to all other nodes is calculated
                Otherwise the integer has to correspond to a node index
                
            parameter : float
                compound parameter which includes the infection and recovery rate alpha and beta, respectively, 
                the mobility rate kappa and the Euler-Mascheroni constant lambda:
                    log[ (alpha-beta)/kappa - lambda ]
        
            saveto : string
                If empty, the result is saved internally in self.dominant_path_distance           
                
        Returns:
        --------
            random_walk_distance : ndarray or float
                If source and target are specified, a float value is returned that specifies the distance.
                
                If either source or target is None a numpy array is returned.
                The position corresponds to the node ID. 
                shape = (Nnodes,)
                
                If both are None a numpy array is returned.
                Each row corresponds to the node ID. 
                shape = (Nnodes,Nnodes)
        """           
        
        assert (isinstance(parameter,float) or isinstance(parameter,int)) and parameter > 0
        assert isinstance(saveto,str)
        assert self.graph != None, "Load graph first."        
        
        P = adjacency_matrix(self.graph, weight="transition_rate").tocsc()
        assert np.all(np.isclose(P.sum(axis=1), 1, rtol=1e-15)), "The transition matrix has to be row normalized"
        
        one = eye( self.nodes, format="csc")
        Z = inv( one - P * np.exp(-parameter))
        D = diags(1./Z.diagonal(), format="csc")
        RWED = -np.log( Z.dot(D).toarray() )

        if source is not None:
            if target is not None:
                RWED = RWED[source, target]
            else:
                RWED = RWED[source,:]
        elif target is not None:
            RWED = RWED[:,target]

        if saveto is not "":
            save( saveto, RWED )
        
        return RWED
    
    ###############################################################################
    
    def get_multiple_path_distance(self, source=None, target=None, parameter=1, cutoff=0, saveto="", verbose=False):
        """
        Compute the multiple path effective distance:
        Gautreau, A, Barrat, A. and Barthelemy, M., Global disease spread:
        Statistics and estimation of arrival times, Journal of Theoretical Biology 251, 3, 509 (2008)
        
        Parameters
        ----------
             source : int or None
                If source is None, the distances from all nodes to the target is calculated
                Otherwise the integer has to correspond to a node index
            
            target : int or None
                If target is None, the distances from the source to all other nodes is calculated
                Otherwise the integer has to correspond to a node index
                
            parameter : float
                compound parameter which includes the infection and recovery rate alpha and beta, respectively, 
                the mobility rate kappa and the Euler-Mascheroni constant lambda:
                    log[ (alpha-beta)/kappa - lambda ]
        
            save : string
                If empty, the result is saved internally in self.dominant_path_distance      
            
            verbose : bool
                Print progress if set as True (False by default)
                
        Returns:
        --------
            dominant_path_distance : ndarray or float
                If source and target are specified, a float value is returned that specifies the distance.
                
                If either source or target is None a numpy array is returned.
                The position corresponds to the node ID. 
                shape = (Nnodes,)
                
                If both are None a numpy array is returned.
                Each row corresponds to the node ID. 
                shape = (Nnodes,Nnodes)
        """
        
        assert (isinstance(parameter,float) or isinstance(parameter,int)) and parameter > 0
        assert isinstance(saveto,str)
        assert self.graph != None, "Load graph first."  
        
        if source is None:
            sources = range(self.nodes)
        else:
            sources = [source]
        if target is None:
            targets = range(self.nodes)
        else:
            targets = [target]
        
        MPED_dic = {}
        for s in sources:
            if verbose:
                print s, "out of", self.nodes
            MPED_dic[s] = np.zeros((self.nodes,))
            for t in targets:
                if s != t:
                    shortest = len(shortest_path(self.graph, source=s, target=t, weight=None))-1
                    paths = all_simple_paths(self.graph, source=s, target=t, cutoff=shortest+cutoff)
                    psum = 0
                    for path in paths:
                        n = len(path)-1
                        prob = np.prod([self.graph[path[ii]][path[ii+1]]["transition_rate"] for ii in range(n)])
                        psum += prob*exp(-parameter*n)
                    MPED_dic[s][t] = -np.log(psum)
        
        if source is not None and target is not None:
            MPED = MPED_dic[source][target]
        elif source is None and target is None:
            MPED = np.array([MPED_dic[s] for s in xrange(self.nodes)]).transpose()
        else:
            if target is None:
                MPED = MPED_dic[source]
            else:
                MPED = np.array([MPED_dic[s][target] for s in xrange(self.nodes)])
        
        if saveto is not "":
            save( saveto, MPED )
        else:
            return MPED
    
    ###############################################################################
    
    def get_shortest_path_distance(self, source=None, target=None, saveto=""):
        """
        Compute the (topological) shortest path distance:

        Parameters
        ----------
             source : int or None
                If source is None, the distances from all nodes to the target is calculated
                Otherwise the integer has to correspond to a node index
            
            target : int or None
                If target is None, the distances from the source to all other nodes is calculated
                Otherwise the integer has to correspond to a node index
                
            saveto : string
                If empty, the result is saved internally in self.dominant_path_distance           
                
        Returns:
        --------
            dominant_path_distance : ndarray or float
                If source and target are specified, a float value is returned that specifies the distance.
                
                If either source or target is None a numpy array is returned.
                The position corresponds to the node ID. 
                shape = (Nnodes,)
                
                If both are None a numpy array is returned.
                Each row corresponds to the node ID. 
                shape = (Nnodes,Nnodes)
        """    
        
        assert isinstance(saveto,str)
        assert self.graph != None, "Load graph first."        
        
        SPD_dic = shortest_path_length(self.graph, source=source, target=target)
        
        if source is None and target is None:
            SPD = np.array([SPD_dic[s].values() for s in xrange(self.nodes)]).transpose()
        elif (source is None) != (target is None):
            SPD = array(SPD_dic.values())
        else:
            SPD = SPD_dic
        if saveto is not "":
            save( saveto, SPD )
        else:
            return SPD
            
###############################################################################
    
if __name__ == "__main__":
    
    # CHECK VERSIONS 
    vers_python0 = '2.7.11'
    vers_numpy0  = '1.11.0'
    vers_scipy0  = '0.17.0'
    vers_netx0   = '1.9.1'
    
    from sys import version_info
    from scipy import __version__ as vers_scipy    
    from networkx import __version__ as vers_netx
    vers_python = '%s.%s.%s' % version_info[:3]
    vers_numpy  = np.__version__
    
    
    print '\n---------'
    print '---------'
    print '---------'
    print '---------'
    print '------------------- Effective Distances ---------------------------'
    print '---------'
    print 'Required modules:'
    print 'Python:   tested for: %s.  Yours: %s'    % (vers_python0, vers_python)
    print 'numpy:    tested for: %s.  Yours: %s'    % (vers_numpy0, vers_numpy)
    print 'scipy:    tested for: %s.  Yours: %s'    % (vers_scipy0, vers_scipy)
    print 'networkx: tested for: %s.   Yours: %s'    % (vers_netx0, vers_netx)
    print '--------'
    print '--------'
    print '--------'
    print '--------\n'
    
    G = EffectiveDistances("/home/andreasko/effective_distance/line_graph.txt", verbose=False, dtype=int)
    #distances = G.get_dominant_path_distance(source=1, target=100, parameter=1, saveto="")
    #print "Dominant Path Effective Distance :", distances
    #distances = G.get_dominant_path_distance(source=0, target=None, parameter=3., saveto="/home/andreasko/effective_distance/line_graph_DPED")
    distances = G.get_multiple_path_distance(source=0, target=None, parameter=3., saveto="/home/andreasko/effective_distance/line_graph_MPED")
    print "Multiple Path Effective Distance :", distances
    #distances = G.get_random_walk_distance(source=1, target=100, parameter=1, saveto="")
    #print "Random Walk Effective Distance   :", distances
    
    
