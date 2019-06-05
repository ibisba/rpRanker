import csv
import os
import itertools
import logging
import collections
import pickle
import logging
import gzip
import sys

## @package InputReader
#
# Documentation for the input files reader of rpFBA

## \brief Class to read all the input files
#
# Contains all the functions that read the cache files and input files to reconstruct the heterologous pathways
# To include in the input directory the following files are required:
# - chemicals.csv (MetaNetX)
# - compartments.csv (MetaNetX)
# - compounds.txt (RP2paths output)
# - out_paths.csv (RP2paths output)
# - scope.csv (RetroPath2 output)
# - xref.csv (User, if needed)
# - sink.csv (User, if needed)
class rpReader:
    """ WARNING: if you define inputPath, then all the files must have specific names to
        make sure that it can find the appropriate files
    """
    ## InputReader constructor
    # 
    #  @param self The object pointer
    #  @param inputPath The path to the folder that contains all the input/output files required
    #  @param Database The database name of the user's xref
    def __init__(self):
        #cache files
        self.deprecatedMNXM_mnxm = None
        #input files
        self.rp_paths = None
        self.rp_smiles = None
        #TODO: not the best strategy to have this set in this fashion -- change it
        #self.model_chemicals = None
        #self.model_compartments = None
        #TODO: not the best strategy to have this set in this fashion -- change it
        self.rp_transformation = None
        self.rp_smiles_inchi = None
        self.smiles_inchi = None
        self.in_xref = None
        self.in_inchi = None
        if not self._loadCache():
            raise ValueError


    #######################################################################
    ############################# PRIVATE FUNCTIONS ####################### 
    #######################################################################


    ## Private function to load the required cache parameters
    #
    #
    def _loadCache(self):
        dirname = os.path.dirname(os.path.abspath( __file__ ))
        try:
            self.deprecatedMNXM_mnxm = pickle.load(open(dirname+'/cache/deprecatedMNXM_mnxm.pickle', 'rb'))
        except FileNotFoundError as e:
            logging.error(e)
            return False
        try:
            self.smiles_inchi = pickle.load(gzip.open(dirname+'/cache/smiles_inchi.pickle.gz', 'rb'))
        except FileNotFoundError as e:
            logging.error(e)
            return False
        return True


    ###############################################################
    ############################# PUBLIC FUNCTIONS ################
    ############################################################### 


    ## Function to parse a given cross references file by the user
    #
    #  Extract the identifiers for each compound
    #
    #  @param self Object pointer
    #  @param path The input file path
    #  @return xref The dictionnary of cross references
    def user_xref(self, path):
        """The file has to be of the form : ID   DatabaseID  Database
        """
        self.in_xref = {}
        try:
            with open(path) as f:
                reader = csv.reader(f, delimiter='\t')
                next(reader)
                for row in reader:
                    for i in range(1,len(row)-1):
                        self.in_xref[row[0]] = {'{}' .format(str(row[i+1])): row[i]}
                        i += 2
        except(FileNotFoundError, TypeError):
            logging.warning('You did not provide an xref file, or the path is not correct')


    ## Function to extract the inchi for each compound in he model of the user
    #
    #  @param self The Object pointer
    #  @return inchi_sink Dictionnary of inchis for the user's ids  
    def user_sink(self, path):
        try:
            self.in_inchi = {}
            with open(path) as f:
                    reader = csv.reader(f, delimiter='\t')
                    next(reader)
                    for row in reader:
                        self.in_inchi[row[0]]=row[1]
        except(FileNotFoundError, TypeError):
            logging.warning('You did not provide a sink file, or the path is not correct')


    ## Function to parse the compounds.txt file
    #
    #  Extract the smile and the structure of each compounds of RP2Path output
    #
    #  @param self Object pointer
    #  @param path The compounds.txt file path
    #  @return rp_compounds Dictionnary of smile and structure for each compound
    def compounds(self, path):
        """ Method to parse all the RP output compounds.
        """
        self.rp_smiles = {}
        try:
            with open(path) as f:
                reader = csv.reader(f, delimiter='\t')
                next(reader)
                for row in reader:
                    self.rp_smiles[row[0]] = {'smiles': row[1]}  #, 'structure':row[1].replace('[','').replace(']','')
        except (TypeError, FileNotFoundError) as e:
            logging.error('Could not read the compounds file ('+str(path)+')')


    ## Function to parse the scope.csv file
    #
    #  Extract the reaction rules from the retroPath2.0 output using the scope.csv file
    #
    #  @param self Object pointer
    #  @param path The scope.csv file path
    #  @return rp_transformation Dictionnary discribing each transformation
    #  @return smiles_inchi dictionnary describing the inchi for each smile
    def transformation(self, path):
        self.rp_transformation = {}
        self.rp_smiles_inchi = {}
        try:
            with open(path) as f:
                reader = csv.reader(f, delimiter=',')
                next(reader)
                for row in reader:
                    self.rp_transformation[row[1]] = row[2]
                    if not row[3] in self.smiles_inchi:
                        self.smiles_inchi[row[3]] = row[4]
                    if not row[5] in self.smiles_inchi:
                        self.smiles_inchi[row[5]] = row[6]
        except (TypeError, FileNotFoundError) as e:
            logging.error('Could not read the compounds file ('+str(path)+')')


    ''' DEPRECATED
    ## Function to parse the compartments.csv file
    #
    #  Parse the compartments.csv file to extract the full name and short name of the different compartments
    #
    #  @param self Object pointer 
    #  @param path The compartments.csv file path
    #  @return model_compartments Dictionnary of compartments names
    def compartments(self, path):
        """ Open the compartments.tsv file from Metanetx that gives the common name to the MNX ID's
        TODO: make this optional
        """
        self.model_compartments = {}
        try:
            with open(path) as f:
                reader = csv.reader(f, delimiter='\t')
                next(reader)
                for row in reader:
                    self.model_compartments[row[0]] = {'full_name': row[1], 'short_name': row[2]}
        except (TypeError, FileNotFoundError) as e:
            logging.error('Could not read the compartments file ('+str(path)+')')
    '''


    ''' DEPRECATED
    ## Function to parse chemicals.csv file
    #
    #  Extract different information about components
    #
    #  @param self Object pointer
    #  @param The chemicals.csv file path
    #  @return model_chemicals Dictionnary of component information
    def chemicals(self, path):
        """ Open the chemicals.tsv file from MetaNetX that describes the sink from a model with InChI
            TODO: Replace this with a method that scans an SBML document and extracts all the chemical
            species
        """
        self.model_chemicals = {}
        ################# open the chemicals file #############
        try:
            with open(path) as f:
                reader = csv.reader(f, delimiter='\t')
                for row in reader:
                    ######### chemical formula #############
                    chem_formula = None
                    if not row[3]=='':
                        chem_formula = row[3]
                    ########## mass #####################
                    mass = None
                    if not row[4]=='':
                        try:
                            mass = float(row[4])
                        except ValueError:
                            logging.error('Could not convert the mass to float ('+str(row[4])+')')
                    ########## charge ##################
                    charge = None
                    if not row[5]=='':
                        try:
                            charge = int(row[5])
                        except ValueError:
                            logging.error('Could not convert charge to int ('+str(row[5])+')')
                    ######### xref #####################
                    xref = {} #construct xref dict
                    for i in list(set([i.split(':')[0] for i in row[6].split(';')])): #unique xref db names
                        xref[i] = []
                    for i in [i.split(':') for i in row[6].split(';')]:
                        if len(i)==2:
                            xref[i[0]].append(i[1])
                    self.model_chemicals[row[0]] = {'name': row[1],
                            'names': row[2].split(';'),
                            'chem_formula': chem_formula,
                            'mass': mass,
                            'charge': charge,
                            'xref': xref}
        except (TypeError, FileNotFoundError) as e:
            logging.error('Could not read the chemicals file ('+str(path)+')')
    '''


    ## Function to parse the out_paths.csv file
    #
    #  Reading the RP2path output and extract all the information for each pathway
    #
    #  @param self Object pointer
    #  @param path The out_path.csv file path
    #  @return toRet_rp_paths Pathway object
    def outPaths(self, path):
        """RP2path Metabolic pathways from out_paths.csv
        create all the different values for heterologous paths from the RP2path out_paths.csv file
        Note that path_step are in reverse order here
        """
        ########## open either the global path or the local defined path ############
        #### (with priority with the local path)
        try:
            rp_paths = {}
            #reactions = self.rr_reactions
            with open(path) as f:
                reader = csv.reader(f)
                next(reader)
                current_path_id = 0
                path_step = 0
                for row in reader:
                    if not int(row[0])==current_path_id:
                        path_step = 0
                    else:
                        path_step += 1
                    current_path_id = int(row[0])
                    ################################################################
                    # WARNING: we are using ONLY the first rule and not the others #
                    ################################################################
                    #for singleRule in row[2].split(','):
                    tmpReac = {'rule_id': row[2].split(',')[0],
                            'right': {},
                            'left': {},
                            'step': path_step,
                            'path_id': int(row[0]),
                            'transformation_id': row[1][:-2]}
                    for l in row[3].split(':'):
                        tmp_l = l.split('.')
                        try:
                            #tmpReac['left'].append({'stoichio': int(tmp_l[0]), 'name': tmp_l[1]})
                            mnxm = ''
                            if tmp_l[1] in self.deprecatedMNXM_mnxm:
                                mnxm = self.deprecatedMNXM_mnxm[tmp_l[1]]
                            else:
                                mnxm = tmp_l[1]
                            tmpReac['left'][mnxm] = int(tmp_l[0])
                        except ValueError:
                            logging.error('Cannot convert tmp_l[0] to int ('+str(tmp_l[0])+')')
                            return {}
                    for r in row[4].split(':'):
                        tmp_r = r.split('.')
                        try:
                            #tmpReac['right'].append({'stoichio': int(tmp_r[0]), 'name': tmp_r[1]})
                            mnxm = ''
                            if tmp_r[1] in self.deprecatedMNXM_mnxm:
                                mnxm = self.deprecatedMNXM_mnxm[tmp_r[1]]  #+':'+self.rr_reactions[tmpReac['rule_id']]['left']
                            else:
                                mnxm = tmp_r[1]  #+':'+self.rr_reactions[tmpReac['rule_id']]['left']
                            tmpReac['right'][mnxm] = int(tmp_r[0])
                        except ValueError:
                            logging.error('Cannot convert tmp_r[0] to int ('+str(tmp_r[0])+')')
                            return {}
                    '''##associate an MNX id for compounds annotated as CMPD ##Don't work when more than 1 compound in the left side
                    to_remove = []
                    for i in self.rr_reactions[tmpReac['rule_id']]['right'].split('.'):   
                        if not i in tmpReac['left']:
                            for n in list(tmpReac['left']):
                                if n[:4] == 'CMPD':
                                    tmpReac['left'][n+':'+i] = tmpReac['left'][n]
                                    to_remove.append(n)
                    for i in to_remove:
                        del tmpReac['left'][i]'''
                    try:
                        if not int(row[0]) in rp_paths:
                            rp_paths[int(row[0])] = []
                        rp_paths[int(row[0])].insert(0, tmpReac)
                    except ValueError:
                        logging.error('Cannot convert path_id to int ('+str(row[0])+')')
                        return {}
            ####### now check where are the duplicates path_steps in each path and duplicate if yes ###
            self.rp_paths = [] # we make this into a list instead of a dict 
            #to find index positions in an array: usage find([1,3,4,5],[2,3])
            find = lambda searchList, elem: [[i for i, x in enumerate(searchList) if x == e] for e in elem]
            #loop through all path metabolic steps in order
            for path in rp_paths:
                dupli_items = [item for item, count in collections.Counter([i['step'] for i in rp_paths[path]]).items() if count>1]
                dupli_index = find([i['step'] for i in rp_paths[path]], dupli_items)
                flat_dupli_index = [item for sublist in dupli_index for item in sublist]
                if not dupli_items:
                    self.rp_paths.append(rp_paths[path])
                else:
                    keep_always_index = [i for i in [y for y in range(len(rp_paths[path]))] if i not in flat_dupli_index]
                    for dupli_include in list(itertools.product(*dupli_index)):
                        toAdd_index = list(keep_always_index+list(dupli_include))
                        new_path = []
                        for ta_i in toAdd_index:
                            new_path.append(rp_paths[path][ta_i])
                        new_path = sorted(new_path, key=lambda k: k['step'], reverse=True)
                        self.rp_paths.append(new_path)
        except (TypeError, FileNotFoundError) as e:
            logging.error('Could not read the out_paths file ('+str(path)+')')


    ## Generate the sink from a given model and the 
    #
    # NOTE: this only works for MNX models, since we are parsing the id
    # TODO: change this to read the annotations and extract the MNX id's
    #
    def genSink(self, rpsbml, file_out, compartment_id='MNXC3'):
        ### open the cache ###
        cytoplasm_species = []
        for i in rpsbml.model.getListOfSpecies():
            if i.getCompartment()==compartment_id:
                cytoplasm_species.append(i)
        with open(file_out, mode='w') as f:
            writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
            writer.writerow(['Name','InChI'])
            for i in cytoplasm_species:
                res = rpsbml.readAnnotation(i.getAnnotation())
                #extract the MNX id's
                try:
                    mnx = res['metanetx.chemical'][0]
                except KeyError:
                    continue
                #mnx = i.getId().split('__')[0]
                try:
                    inchi = self.smiles_inchi[mnx]['inchi']
                except KeyError:
                    inchi = None
                if mnx and inchi:
                    writer.writerow([mnx,inchi])

