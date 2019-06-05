import libsbml
from hashlib import md5
import logging
import os
import pickle
import gzip

## @package RetroPath SBML writer
# Documentation for SBML representation of the different model
#
# To exchange between the different workflow nodes, the SBML (XML) format is used. This
# implies using the libSBML library to create the standard definitions of species, reactions, etc...
# Here we also define our own annotations that are used internally in that we call IBISBA nodes.
# The object holds an SBML object and a series of methods to write and access IBISBA related annotations


##################################################################
############################### rpSBML ###########################
##################################################################


## libSBML reader for RetroPath
# Converts an SBML object (or file) into the internal format
#
class rpSBML:
    ## Constructor
    #
    # @param model libSBML model object
    # @param docModel libSBML Document object
    # @param nameSpaceModel libSBML name space (not required)
    def __init__(self, modelName, document=None, path=None):
        self.modelName = modelName
        self.document = document
        if self.document==None:
            self.model = None
        else:
            self.model = self.document.getModel()
        self.path = path
        #need to scan for these if we are passing mode and documents
        self.hetero_group = None
        self.compartmentName = None
        self.compartmentId = None
        #if self.model==None:
        #    logging.warning('rpSBML object was initiated as empty. Please call createModel() to initalise the model')

    #######################################################################
    ############################# PRIVATE FUNCTIONS ####################### 
    #######################################################################


    ## Check the libSBML calls
    #
    # Check that the libSBML python calls do not return error INT and if so, display the error. Taken from: http://sbml.org/Software/libSBML/docs/python-api/create_simple_model_8py-example.html
    #
    # @param value The SBML call
    # @param message The string that describes the call
    def _checklibSBML(self, value, message):
        if value is None:
            logging.error('LibSBML returned a null value trying to ' + message + '.')
            raise SystemExit('LibSBML returned a null value trying to ' + message + '.')
        elif type(value) is int:
            if value==libsbml.LIBSBML_OPERATION_SUCCESS:
                return
            else:
                err_msg = 'Error encountered trying to ' + message + '.' \
                        + 'LibSBML returned error code ' + str(value) + ': "' \
                        + libsbml.OperationReturnValue_toString(value).strip() + '"'
                logging.error(err_msg)
                raise SystemExit(err_msg)
        else:
            logging.info(message)
            return None


    ## String to SBML ID
    #
    # Convert any String to one that is compatible with the SBML metaID formatting requirements
    #
    # @param name The input string
    def _nameToSbmlId(self, name):
        IdStream = []
        count = 0
        end = len(name)
        if '0' <= name[count] and name[count] <= '9':
            IdStream.append('_')
        for count in range(0, end):
            if (('0' <= name[count] and name[count] <= '9') or
                    ('a' <= name[count] and name[count] <= 'z') or
                    ('A' <= name[count] and name[count] <= 'Z')):
                IdStream.append(name[count])
            else:
                IdStream.append('_')
        Id = ''.join(IdStream)
        if Id[len(Id) - 1] != '_':
            return Id
        return Id[:-1]


    ## String to hashed ID
    #
    # Hash an input string and then pass it to _nameToSbmlId()
    #
    # @param input string
    def _genMetaID(self, name):
        return self._nameToSbmlId(md5(str(name).encode('utf-8')).hexdigest())        


    #####################################################################
    ########################## READ/WRITE ###############################
    #####################################################################

    
    ## Open an SBML using libSBML 
    #
    # Situation where an SBML is passed to add the heterologous pathway
    #
    # @param inFile String Path to the input SBML file
    def readSBML(self, inFile):
        if not os.path.isfile(inFile):
            logging.error('Invalid input file')
            raise FileNotFoundError
        document = libsbml.readSBML(inFile)
        self._checklibSBML(document, 'readinf input file')
        errors = document.getNumErrors()
        #display the errors in the log accordning to the severity
        for err in [document.getError(i) for i in range(document.getNumErrors())]:
            if err.isFatal:
                logging.error('libSBML reading error: '+str(err.getShortMessage()))
                raise FileNotFoundError
            else:
                logging.warning('libSBML reading warning: '+str(err.getShortMessage()))
        model = document.getModel()
        if not model:
            loging.error('Either the file was not read correctly or the SBML is empty')
            raise FileNotFoundError
        self.document = document
        self.model = model


    ## Export a libSBML model to file
    #
    # Export the libSBML model to an SBML file
    # 
    # @param model libSBML model to be saved to file
    # @param model_id model id, note that the name of the file will be that
    # @param path Non required parameter that will define the path where the model will be saved
    def writeSBML(self, path):
        ####### check the path #########
        #need to determine where are the path id's coming from
        p = None
        if path:
            if path[-1:]=='/':
                path = path[:-1]
            if not os.path.isdir(path):
                if self.path:
                    p = self.path
                else:
                    logging.error('The output path is not a directory: '+str(path))
                    return False
            else:
                p = path
        else:
            p = self.path
        ########## check and create folder #####
        if not os.path.exists(p):
            os.makedirs(p)
        libsbml.writeSBMLToFile(self.document, p+'/'+str(self.modelName)+'.sbml')
        return True


    ## Return the reaction ID's and the pathway annotation
    #
    #
    def readRPpathway(self, pathId='rp_pathway'):
        groups = self.model.getPlugin('groups')
        rp_pathway = groups.getGroup(pathId)
        toRet = {}
        toRet['annotation'] = rp_pathway.getAnnotation()
        toRet['members'] = []
        for member in rp_pathway.getListOfMembers():
            toRet['members'].append(member.getIdRef())
        return toRet


    ##
    #
    #
    def readRPrules(self, path_id='rp_pathway'):
        toRet = {}
        for reacId in self.readRPpathway(path_id)['members']:
            reac = rpsbml.model.getReaction(reacId)
            ibibsa_annot = rpsbml.readIBISBAAnnotation(reac.getAnnotation())
            toRet[ibibsa_annot['rule_id']] = ibibsa_annot['smiles'].replace('&gt;', '>')
        return toRet

    ## Return the the species annitations 
    #
    #
    def readRPspeciesAnnotation(self, path_id='rp_pathway'):
        #reacMembers = {'reactants': {}, 'products': {}}
        for reacId in self.readRPpathway(path_id):
            reacMembers[reacId] = {}
            reac = self.model.getReaction(reacId)
            for pro in reac.getListOfProducts():
                reacMembers['products'][pro.getSpecies()] = pro.getStoichiometry()
            for rea in reac.getListOfReactants():
                reacMembers['reactants'][rea.getSpecies()] = rea.getStoichiometry() 
        return reacMembers
    

    ## Return the MIRIAM annotations of species
    #
    #
    def readMIRIAMSpeciesAnnotation(self, annot, cid):
        id_annotId = {'bigg': 'bigg.metabolite', 
                'mnx': 'metanetx.chemical', 
                'chebi': 'chebi', 
                'hmdb': 'hmdb', 
                'kegg': 'kegg.compound',
                'seed': 'seed.compound'}
        toRet = []
        if not cid in id_annotId:
            logging.warning('Cannnot find '+str(cid)+' in id_annotId')
            return toRet
        bag = annot.getChild('RDF').getChild('Description').getChild('is').getChild('Bag')
        for i in range(bag.getNumChildren()):
            str_annot = bag.getChild(i).getAttrValue(0)
            if str_annot=='':
                logging.error('This contains no attributes: '+str(bag.getChild(i).toXMLString()))
                continue
            if str_annot.split('/')[-2]==id_annotId[cid]:
                toRet.append(str_annot.split('/')[-1])
        return toRet


    ## Takes a libSBML Reactions or Species object and returns a dictionnary for all its elements
    #
    def readAnnotation(self, annot):
        toRet = {}
        bag = annot.getChild('RDF').getChild('Description').getChild('is').getChild('Bag')
        for i in range(bag.getNumChildren()):
            str_annot = bag.getChild(i).getAttrValue(0)
            if str_annot=='':
                logging.error('This contains no attributes: '+str(bag.getChild(i).toXMLString()))
                continue
            if not str_annot.split('/')[-2] in toRet:
                toRet[str_annot.split('/')[-2]] = []
            #TODO: check if the MNXM and if so check against depreceated
            toRet[str_annot.split('/')[-2]].append(str_annot.split('/')[-1])
        return toRet


    ## Takes for input a libSBML annotatio object and returns a dictionnary of the annotations
    #
    #
    def readIBISBAAnnotation(self, annot):
        toRet = {}
        bag = annot.getChild('RDF').getChild('Ibisba').getChild('ibisba')
        for i in range(bag.getNumChildren()):
            ann = bag.getChild(i)
            if ann=='':
                logging.error('This contains no attributes: '+str(ann.toXMLString()))
                continue
            if not ann.getName() in toRet:
                if ann.getName()=='dG_prime_m' or ann.getName()=='dG_uncert' or ann.getName()=='dG_prime_o':
                    toRet[ann.getName()] = {
                            'units': ann.getAttrValue('units'), 
                            'value': ann.getAttrValue('value')}
                else:
                    toRet[ann.getName()] = ann.getChild(0).toXMLString()
        return toRet


    ## Find out if two libSBML Species or Reactions come from the same species
    #
    # Compare two dictionnaries and if any of the values of any of the same keys are the same then the 
    # function return True, and if none are found then return False
    #
    # @param libSBML Annotation object for one of the 
    # @return Boolean to determine if they are the same
    def compareAnnotations(self, source_annot, target_annot):
        source_dict = self.readAnnotation(source_annot)
        target_dict = self.readAnnotation(target_annot)
        #list the common keys between the two
        for com_key in set(list(source_dict.keys()))-(set(list(source_dict.keys()))-set(list(target_dict.keys()))):
            #compare the keys and if same is non-empty means that there 
            #are at least one instance of the key that is the same
            if bool(set(source_dict[com_key]) & set(target_dict[com_key])):
                return True
        return False
 
    
    #########################################################################
    ############################# MODEL APPEND ##############################
    #########################################################################


    ## Merge two models species and reactions using the annotations to recognise the same species and reactions
    #
    # The source mode has to have both the GROUPS and FBC packages enabled in its SBML. The course must have a groups
    #called rp_pathway.
    #
    def mergeModels(self, target_model):
        #target_model = target_document.getModel()
        #Find the ID's of the similar target_model species
        ################ UNITDEFINITIONS ######
        #return the list of unit definitions id's for the target to avoid overwritting
        #WARNING: this means that the original unit definitions will be prefered over the new one
        target_unitDefID = [i.getId() for i in target_model.getListOfUnitDefinitions()]
        for source_unitDef in self.model.getListOfUnitDefinitions():
            if not source_unitDef.getId() in target_unitDefID: #have to compare by ID since no annotation
                #create a new unitDef in the target
                target_unitDef = target_model.createUnitDefinition()
                self._checklibSBML(target_unitDef, 'fetching target unit definition')
                #copy unitDef info to the target
                self._checklibSBML(target_unitDef.setId(source_unitDef.getId()), 
                    'setting target unit definition ID')
                #self._checklibSBML(target_unitDef.setMetaId(source_unitDef.getMetaId()), 
                #    'setting target unit definition MetaId')
                self._checklibSBML(target_unitDef.setAnnotation(source_unitDef.getAnnotation()), 
                    'setting target unit definition Annotation')
                for source_unit in source_unitDef.getListOfUnits():
                    #copy unit info to the target unitDef
                    target_unit = target_unitDef.createUnit()
                    self._checklibSBML(target_unit, 'creating target unit')
                    self._checklibSBML(target_unit.setKind(source_unit.getKind()), 
                        'setting target unit kind')
                    self._checklibSBML(target_unit.setExponent(source_unit.getExponent()), 
                        'setting target unit exponent')
                    self._checklibSBML(target_unit.setScale(source_unit.getScale()), 
                        'setting target unit scale')
                    self._checklibSBML(target_unit.setMultiplier(source_unit.getMultiplier()), 
                        'setting target unit multiplier')
                target_unitDefID.append(source_unitDef.getId()) #add to the list to make sure its not added twice
        ################ COMPARTMENTS ###############
        sourceCompartmentID_targetCompartmentID = {}
        toAddNum = []
        ####### compare the annotations to find the same ones #######
        for i in range(self.model.getNumCompartments()):
            found = False
            source_compartment = self.model.getCompartment(i)
            self._checklibSBML(source_compartment, 'Getting target compartment')
            source_annotation = source_compartment.getAnnotation()
            #self._checklibSBML(source_annotation, 'Getting compartment target annotation')
            if not source_annotation:
                logging.warning('No annotation for the source of compartment '+str(source_compartment.getId()))
                continue
            for y in range(target_model.getNumCompartments()):
                target_compartment = target_model.getCompartment(y)
                self._checklibSBML(target_compartment, 'Getting target compartment')
                target_annotation = target_compartment.getAnnotation()
                #self._checklibSBML(target_annotation, 'Getting target annotation')
                if not target_annotation:    
                    logging.warning('No annotation for the target of compartment: '+str(target_compartment.getId()))
                    continue
                if self.compareAnnotations(source_annotation, target_annotation):
                    sourceCompartmentID_targetCompartmentID[source_compartment.getId()] = target_compartment.getId() 
                    found = True
                    break
            if not found:
                toAddNum.append(i)
        for i in toAddNum:
            source_compartment = self.model.getCompartment(i)
            self._checklibSBML(source_compartment, 'Getting target compartment')
            target_compartment = target_model.createCompartment()
            self._checklibSBML(target_compartment, 'Creating target compartment')
            self._checklibSBML(target_compartment.setMetaId(source_compartment.getMetaId()),
                    'setting target metaId')
            self._checklibSBML(target_compartment.setId(source_compartment.getId()),
                    'setting target id')
            self._checklibSBML(target_compartment.setName(source_compartment.getName()),
                    'setting target name')
            self._checklibSBML(target_compartment.setConstant(source_compartment.getConstant()),
                    'setting target constant')
            self._checklibSBML(target_compartment.setAnnotation(source_compartment.getAnnotation()),
                    'setting target annotation')
            self._checklibSBML(target_compartment.setSBOTerm(source_compartment.getSBOTerm()),
                    'setting target annotation')
        ################ PARAMETERS ###########
        #WARNING: here we compare by ID
        targetParametersID = [i.getId() for i in target_model.getListOfParameters()]
        for source_parameter in self.model.getListOfParameters():
            if not source_parameter.getId() in targetParametersID:
                target_parameter = target_model.createParameter()
                self._checklibSBML(target_parameter, 'creating target parameter')
                self._checklibSBML(target_parameter.setId(source_parameter.getId()), 'setting target parameter ID')
                self._checklibSBML(target_parameter.setSBOTerm(source_parameter.getSBOTerm()), 
                    'setting target parameter SBO')
                self._checklibSBML(target_parameter.setUnits(source_parameter.getUnits()),
                    'setting target parameter Units')
                self._checklibSBML(target_parameter.setValue(source_parameter.getValue()),
                    'setting target parameter Value') 
                self._checklibSBML(target_parameter.setConstant(source_parameter.getConstant()),
                    'setting target parameter ID')
        ################ MODEL FBC ########################
        if not target_model.isPackageEnabled('fbc'):
            self._checklibSBML(target_model.enablePackage(
                'http://www.sbml.org/sbml/level3/version1/fbc/version2',
                'fbc',
                True),
                    'Enabling the FBC package')
        target_fbc = target_model.getPlugin('fbc')
        source_fbc = self.model.getPlugin('fbc')
        ################ FBC GENE PRODUCTS ########################
        #make a list of all the gene product
        #WARNING: here we compare by ID
        targetGenProductID = [i.getId() for i in target_fbc.getListOfGeneProducts()]
        for source_geneProduct in source_fbc.getListOfGeneProducts():
            if not source_geneProduct.getId() in targetGenProductID:
                target_geneProduct = target_fbc.createGeneProduct()
                self._checklibSBML(target_geneProduct, 'creating target gene product')
                self._checklibSBML(target_geneProduct.setId(source_geneProduct.getId()), 
                    'setting target gene product id')
                self._checklibSBML(target_geneProduct.setLabel(source_geneProduct.getLabel()), 
                    'setting target gene product label')
                self._checklibSBML(target_geneProduct.setName(source_geneProduct.getName()), 
                    'setting target gene product name')
                self._checklibSBML(target_geneProduct.setMetaId(source_geneProduct.getMetaId()),
                    'setting target gene product metaID')
        ############### FBC OBJECTIVES ############
        #WARNING: here we compare the Objective by ID, and we add the downstream fluxObjectives
        targetObjectiveID = [i.getId() for i in target_fbc.getListOfObjectives()]
        for source_objective in source_fbc.getListOfObjectives():
            if not source_objective.getId() in targetObjectiveID:
                target_objective = target_fbc.createObjective()
                self._checklibSBML(target_objective, 'creating target objective')
                self._checklibSBML(target_objective.setId(source_objective.getId()), 'setting target objective')
                self._checklibSBML(target_objective.setName(source_objective.getName()), 'setting target objective')
                self._checklibSBML(target_objective.setType(source_objective.getType()), 
                        'setting target objective type')
                for source_fluxObjective in source_objective.getListOfFluxObjectives():
                    target_fluxObjective = target_objective.createFluxObjective()
                    self._checklibSBML(target_fluxObjective, 'creating target flux objective')
                    self._checklibSBML(target_fluxObjective.setName(source_fluxObjective.getName()),
                        'setting target flux objective name')
                    self._checklibSBML(target_fluxObjective.setCoefficient(source_fluxObjective.getCoefficient()),
                        'setting target flux objective coefficient')
                    self._checklibSBML(target_fluxObjective.setReaction(source_fluxObjective.getReaction()),
                        'setting target flux objective reaction')
        ################ SPECIES ####################
        #TODO: modify the name to add rpPathway
        sourceSpeciesID_targetSpeciesID = {}
        toAddNum = []
        for i in range(self.model.getNumSpecies()):
            found = False
            source_species = self.model.getSpecies(i)
            self._checklibSBML(source_species, 'Getting source species')
            source_annotation = source_species.getAnnotation()
            self._checklibSBML(source_annotation, 'Getting source annotation')
            if not source_annotation:
                logging.warning('No annotation for the source of compartment '+str(source_compartment.getId()))
                continue
            for y in range(target_model.getNumSpecies()):
                target_species = target_model.getSpecies(y)
                self._checklibSBML(target_species, 'Getting target species')
                target_annotation = target_species.getAnnotation()
                #self._checklibSBML(target_annotation, 'Getting target annotation')
                if not target_annotation:    
                    logging.warning('Cannot find target number: '+str(y))
                    continue
                if self.compareAnnotations(source_annotation, target_annotation):
                    #save the speciesID as being the same
                    sourceSpeciesID_targetSpeciesID[self.model.species[i].getId()] = target_model.species[y].getId()
                    found = True
                    break
            #if it has not been found then add it to the target_model
            if not found:
                toAddNum.append(i)
        for i in toAddNum:
            source_species = self.model.getSpecies(i)
            self._checklibSBML(source_species, 'fetching source species')
            target_species = target_model.createSpecies()
            self._checklibSBML(target_species, 'creating species')
            self._checklibSBML(target_species.setMetaId(source_species.getMetaId()),
                    'setting target metaId')
            self._checklibSBML(target_species.setId(source_species.getId()),
                    'setting target id')
            self._checklibSBML(target_species.setCompartment(source_species.getCompartment()),
                    'setting target compartment')
            self._checklibSBML(target_species.setInitialConcentration(
                source_species.getInitialConcentration()),
                    'setting target initial concentration')
            self._checklibSBML(target_species.setBoundaryCondition(
                source_species.getBoundaryCondition()),
                    'setting target boundary concentration')
            self._checklibSBML(target_species.setHasOnlySubstanceUnits(
                source_species.getHasOnlySubstanceUnits()),
                    'setting target has only substance units')
            self._checklibSBML(target_species.setBoundaryCondition(
                source_species.getBoundaryCondition()),
                    'setting target boundary condition')
            self._checklibSBML(target_species.setConstant(source_species.getConstant()),
                'setting target constant')
            self._checklibSBML(target_species.setAnnotation(source_species.getAnnotation()),
                'setting target annotation')
        ################ REACTIONS ###################
        #Find the ID's of the similar target_model reactions
        #need to create a new instance of reactions to add to the model
        #test to see if the target model has the FBC package and if not add it
        sourceReactionsID_targetReactionsID = {}
        toAddNum = []
        if not target_model.isPackageEnabled('fbc'):
            self._checklibSBML(target_model.enablePackage(
                'http://www.sbml.org/sbml/level3/version1/fbc/version2',
                'fbc',
                True),
                    'Enabling the FBC package')
        #note sure why one needs to set this as False
        self._checklibSBML(self.document.setPackageRequired('fbc', False), 'enabling FBC package')
        for i in range(self.model.getNumReactions()):
            found = False
            source_reaction = self.model.getReaction(i)
            self._checklibSBML(source_reaction, 'fetching source reaction')
            source_annotation = source_reaction.getAnnotation()
            #self._checklibSBML(source_annotation, 'fetching source reaction annotation')
            if not source_annotation:
                logging.warning('No annotation for the source of reaction: '+str(source_reaction.getId()))
                continue
            for y in range(target_model.getNumReactions()):
                target_reaction = target_model.getReaction(y)
                self._checklibSBML(target_reaction, 'fetching target reaction annotation') 
                target_annotation = target_reaction.getAnnotation()
                #self._checklibSBML(target_annotation, 'fetching target reaction annotation') 
                if not target_annotation:    
                    logging.warning('No annotation for the target of reaction: '+str(target_reaction.getId()))
                    continue
                if self.compareAnnotations(source_annotation, target_annotation):
                    sourceReactionsID_targetReactionsID[self.model.reactions[i].getId()] = target_model.reactions[y].getId()
                    found = True
                    break
            if not found:
                toAddNum.append(i)
        for i in toAddNum:
            source_reaction = self.model.getReaction(i)
            self._checklibSBML(source_reaction, 'fetching source reaction')
            target_reaction = target_model.createReaction()
            self._checklibSBML(target_reaction, 'create reaction')
            target_fbc = target_reaction.getPlugin('fbc')
            self._checklibSBML(target_fbc, 'fetching target FBC package')
            source_fbc = source_reaction.getPlugin('fbc')
            self._checklibSBML(source_fbc, 'fetching source FBC package')
            source_upperFluxBound = source_fbc.getUpperFluxBound()
            self._checklibSBML(source_upperFluxBound, 'fetching upper flux bound')
            self._checklibSBML(target_fbc.setUpperFluxBound(source_upperFluxBound), 
                    'setting upper flux bound')
            source_lowerFluxBound = source_fbc.getLowerFluxBound()
            self._checklibSBML(source_lowerFluxBound, 'fetching lower flux bound')
            self._checklibSBML(target_fbc.setLowerFluxBound(source_lowerFluxBound), 
                    'setting lower flux bound')
            self._checklibSBML(target_reaction.setId(source_reaction.getId()), 'set reaction id')
            self._checklibSBML(target_reaction.setName(source_reaction.getName()), 'set name')
            self._checklibSBML(target_reaction.setSBOTerm(source_reaction.getSBOTerm()), 
                    'setting the reaction system biology ontology (SBO)') #set as process
            #TODO: consider having the two parameters as input to the function
            self._checklibSBML(target_reaction.setReversible(source_reaction.getReversible()), 
                    'set reaction reversibility flag')
            self._checklibSBML(target_reaction.setFast(source_reaction.getFast()), 
                    'set reaction "fast" attribute')
            self._checklibSBML(target_reaction.setMetaId(source_reaction.getMetaId()), 'setting species metaID')
            self._checklibSBML(target_reaction.setAnnotation(source_reaction.getAnnotation()), 
                    'setting annotation for source reaction')
            #reactants_dict
            for y in range(source_reaction.getNumReactants()):
                target_reactant = target_reaction.createReactant()
                self._checklibSBML(target_reactant, 'create target reactant')
                source_reactant = source_reaction.getReactant(y)
                self._checklibSBML(source_reactant, 'fetch source reactant')
                try:
                    #try to get the reactant from the dictionnary if annotations comparison 
                    #elects them to be the same
                    reactantID = sourceSpeciesID_targetSpeciesID[source_reactant.species]
                except KeyError:
                    #if not found in dictionnary then muct be part of the added ones
                    reactantID = source_reactant.species
                self._checklibSBML(target_reactant.setSpecies(reactantID), 'assign reactant species')
                #TODO: check to see the consequences of heterologous parameters not being constant
                self._checklibSBML(target_reactant.setConstant(source_reactant.getConstant()), 
                        'set "constant" on species '+str(source_reactant.getConstant()))
                self._checklibSBML(target_reactant.setStoichiometry(source_reactant.getStoichiometry()),
                        'set stoichiometry ('+str(source_reactant.getStoichiometry)+')')
            #products_dict
            for y in range(source_reaction.getNumProducts()):
                target_product = target_reaction.createProduct()
                self._checklibSBML(target_product, 'create target product')
                source_product = source_reaction.getProduct(y)
                self._checklibSBML(source_product, 'fetch source product')
                try:
                    #try to get the reactant from the dictionnary if annotations comparison 
                    #elects them to be the same
                    productID = sourceSpeciesID_targetSpeciesID[source_product.species]
                except KeyError:
                    #if not found in dictionnary then muct be part of the added ones
                    productID = source_product.species
                self._checklibSBML(target_product.setSpecies(productID), 'assign reactant species')
                #TODO: check to see the consequences of heterologous parameters not being constant
                self._checklibSBML(target_product.setConstant(source_product.getConstant()), 
                        'set "constant" on species '+str(source_product.getConstant()))
                self._checklibSBML(target_product.setStoichiometry(source_product.getStoichiometry()),
                        'set stoichiometry ('+str(source_product.getStoichiometry)+')')
        #### GROUPS #####
        if not target_model.isPackageEnabled('groups'):
            self._checklibSBML(target_model.enablePackage(
                'http://www.sbml.org/sbml/level3/version1/groups/version1',
                'groups',
                True),
                    'Enabling the GROUPS package')
        #!!!! must be set to false for no apparent reason
        self._checklibSBML(self.document.setPackageRequired('groups', False), 'enabling groups package')
        #self._checklibSBML(self.sbmlns.addPkgNamespace('groups',1), 'Add groups package')
        source_groups = self.model.getPlugin('groups')
        self._checklibSBML(source_groups, 'fetching the source model groups')
        target_groups = target_model.getPlugin('groups')
        self._checklibSBML(target_groups, 'fetching the target model groups')
        self._checklibSBML(target_groups.addGroup(source_groups.getGroup('rp_pathway')),
                'copying the source groups "rp_pathway" to the target groups')
        

    #########################################################################
    ############################# MODEL CREATION FUNCTIONS ##################
    #########################################################################


    ## Create libSBML model instance
    #
    # Function that creates a new libSBML model instance and initiates it with the appropriate packages. Creates a cytosol compartment
    #
    # @param name The name of the model
    # @param modelID The id of the mode
    # @param metaID metaID of the model. Default None means that we will generate a hash from the modelID
    def createModel(self, name, modelID, metaID=None):
        ## sbmldoc
        self.sbmlns = libsbml.SBMLNamespaces(3,1)
        self._checklibSBML(self.sbmlns, 'generating model namespace')
        self._checklibSBML(self.sbmlns.addPkgNamespace('groups',1), 'Add groups package')
        self._checklibSBML(self.sbmlns.addPkgNamespace('fbc',2), 'Add FBC package')
        #sbmlns = libsbml.SBMLNamespaces(3,1,'groups',1)
        self.document = libsbml.SBMLDocument(self.sbmlns)
        self._checklibSBML(self.document, 'generating model doc')
        #!!!! must be set to false for no apparent reason
        self._checklibSBML(self.document.setPackageRequired('fbc', False), 'enabling FBC package') 
        #!!!! must be set to false for no apparent reason
        self._checklibSBML(self.document.setPackageRequired('groups', False), 'enabling groups package')
        ## sbml model
        self.model = self.document.createModel()
        self._checklibSBML(self.model, 'generating the model')
        self._checklibSBML(self.model.setId(modelID), 'setting the model ID')
        model_fbc = self.model.getPlugin('fbc')
        model_fbc.setStrict(True)
        if metaID==None:
            metaID = self._genMetaID(modelID)
        self._checklibSBML(self.model.setMetaId(metaID), 'setting model metaID')
        self._checklibSBML(self.model.setName(name), 'setting model name')
        self._checklibSBML(self.model.setTimeUnits('second'), 'setting model time unit')
        self._checklibSBML(self.model.setExtentUnits('mole'), 'setting model compartment unit')
        self._checklibSBML(self.model.setSubstanceUnits('mole'), 'setting model substance unit')


    ## Create libSBML compartment 
    #
    # cytoplasm compartment TODO: consider seperating it in another function if another compartment is to be created
    #
    # @param model libSBML model object to add the compartment
    # @param size Set the compartement size
    # @return boolean Execution success
    #TODO: set the compName as None by default. To do that you need to regenerate the compXref to 
    #use MNX ids as keys instead of the string names
    def createCompartment(self, size, compId, compName, compXref, metaID=None):
        comp = self.model.createCompartment()
        self._checklibSBML(comp, 'create compartment')
        self._checklibSBML(comp.setId(compId), 'set compartment id')
        self.compartmentId = compId
        if compName:
            self._checklibSBML(comp.setName(compName), 'set the name for the cytoplam')
            self.compartmentName = compName
        self._checklibSBML(comp.setConstant(True), 'set compartment "constant"')
        self._checklibSBML(comp.setSize(size), 'set compartment "size"')
        self._checklibSBML(comp.setSBOTerm(290), 'set SBO term for the cytoplasm compartment')
        if metaID==None:
            metaID = self._genMetaID(compId)
        self._checklibSBML(comp.setMetaId(metaID), 'set the metaID for the compartment')
        #self.compartmentNames.append(compartmentName) #this assumes there is only one compartmentName
        annotation = '''<annotation>
  <rdf:RDF 
  xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" 
  xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" 
  xmlns:bqmodel="http://biomodels.net/model-qualifiers/">'''
        # if the name of the species is MNX then we annotate it using MIRIAM compliance
        #TODO: need to add all known xref from different databases (not just MetaNetX)
        annotation += '''
    <rdf:Description rdf:about="#'''+str(metaID or '')+'''">
      <bqbiol:is>
        <rdf:Bag>'''
        #TODO: for yout to complete
        id_ident = {'mnx': 'metanetx.compartment/', 'bigg': 'bigg.compartment/'}
        #WARNING: compartmentNameID as of now, needs to be a MNX ID
        if compName in compXref: 
            for databaseId in compXref[compName]:
                for compartmentId in compXref[compName][databaseId]:
                    try:
                        annotation += '''
          <rdf:li rdf:resource="http://identifiers.org/'''+str(id_ident[databaseId])+str(compartmentId)+'''"/>'''
                    except KeyError:
                        continue
        annotation += '''
        </rdf:Bag>
      </bqbiol:is>
    </rdf:Description> 
  </rdf:RDF>
</annotation>'''
        self._checklibSBML(comp.setAnnotation(annotation), 'setting annotation for reaction '+str(compName))


    ## Create libSBML unit definition
    #
    # Function that creates a unit definition (composed of one or more units)
    #
    # @param model libSBML model to add the unit definition
    # @param unit_id ID for the unit definition
    # @param metaID metaID for the unit definition. If None creates a hash from unit_id
    # @return Unit definition
    def createUnitDefinition(self, unit_id, metaID=None):
        unitDef = self.model.createUnitDefinition()
        self._checklibSBML(unitDef, 'creating unit definition')
        self._checklibSBML(unitDef.setId(unit_id), 'setting id')
        if metaID==None:
            metaID = self._genMetaID(unit_id)
        self._checklibSBML(unitDef.setMetaId(metaID), 'setting metaID')
        #self.unitDefinitions.append(unit_id)
        return unitDef


    ## Create libSBML unit
    #
    # Function that created a unit
    #
    # @param unitDef libSBML unit definition
    # @param libsmlunit libSBML unit parameter
    # @param exponent Value for the exponent (ex 10^5 mol/sec)
    # @param scale Value for the scale 
    # @param multiplier Value for the multiplie
    # @return Unit
    def createUnit(self, unitDef, libsbmlunit, exponent, scale, multiplier):
        unit = unitDef.createUnit()
        self._checklibSBML(unit, 'creating unit')
        self._checklibSBML(unit.setKind(libsbmlunit), 'setting the kind of unit')
        self._checklibSBML(unit.setExponent(exponent), 'setting the exponenent of the unit')
        self._checklibSBML(unit.setScale(scale), 'setting the scale of the unit')
        self._checklibSBML(unit.setMultiplier(multiplier), 'setting the multiplier of the unit')


    ## Create libSBML parameters
    #
    # Parameters, in our case, are used for the bounds for FBA analysis. Unit parameter must be an instance of unitDefinition
    #
    # @param parameter_id SBML id
    # @param value Float value for this parameter
    # @param unit libSBML unit parameter
    # @param metaID String Optional parameter for SBML metaID
    # @return libSBML parameter object
    def createParameter(self, parameter_id, value, unit, metaID=None):
        newParam = self.model.createParameter()
        self._checklibSBML(newParam, 'Creating a new parameter object')
        self._checklibSBML(newParam.setConstant(True), 'setting as constant')
        self._checklibSBML(newParam.setId(parameter_id), 'setting ID')
        self._checklibSBML(newParam.setValue(value), 'setting value')
        self._checklibSBML(newParam.setUnits(unit), 'setting units')
        self._checklibSBML(newParam.setSBOTerm(625), 'setting SBO term')
        if metaID==None:
            metaID = self._genMetaID(parameter_id)
        self._checklibSBML(newParam.setMetaId(metaID), 'setting meta ID')
        #self.parameters.append(parameter_id)
        return newParam


    ## Create libSBML reaction
    #
    # Create a reaction. fluxBounds is a list of libSBML.UnitDefinition, length of exactly 2 with the first position that is the upper bound and the second is the lower bound. reactants_dict and reactants_dict are dictionnaries that hold the following parameters: name, compartment, stoichiometry
    #
    # @param name Name for the reaction
    # @param reaction_id Reaction ID
    # @param fluxUpperBounds FBC id for the upper flux bound for this reaction
    # @param fluxLowerBounds FBC id for the lower flux bound for this reaction
    # BILAL check the lower
    # @param step 2D dictionnary with the following structure {'left': {'name': stoichiometry, ...}, 'right': {}}
    # @param reaction_smiles String smiles description of this reaction (added in IBISBA annotation)
    # @param compartmentId String Optinal parameter compartment ID
    # @param hetero_group Groups Optional parameter object that holds all the heterologous pathways
    # @param metaID String Optional parameter reaction metaID
    # @return metaID meta ID for this reaction
    def createReaction(self,
            reacId,
            fluxUpperBound,
            fluxLowerBound,
            step,
            reaction_smiles,
            reacXref,
            compartmentId=None,
            hetero_group=None,
            metaID=None):
        '''
        #### format of a step
        [{'rule_id': 'RR-01-f145fec21961-49-F',
          'right': {'CMPD_0000000003': 1},
          'left': {'MNXM10': 1, 'MNXM188': 1, 'MNXM4': 1},
          'step': 1,
          'path_id': 1,
          'transformation_id': 'TRS_0_1_11'},
         {'rule_id': 'RR-01-a0cc0be463ff-49-F',
          'right': {'TARGET_0000000001': 1},
          'left': {'CMPD_0000000003': 1, 'MNXM4': 1},
          'step': 0,
          'path_id': 1,
          'transformation_id': 'TRS_0_0_0'}]
        '''
        reac = self.model.createReaction()
        self._checklibSBML(reac, 'create reaction')
        ################ FBC ####################
        reac_fbc = reac.getPlugin('fbc')
        self._checklibSBML(reac_fbc, 'extending reaction for FBC')
        #bounds
        self._checklibSBML(reac_fbc.setUpperFluxBound(fluxUpperBound), 'setting '+str(reacId)+' upper flux bound')
        self._checklibSBML(reac_fbc.setLowerFluxBound(fluxLowerBound), 'setting '+str(reacId)+' lower flux bound')
        #########################################
        #reactions
        self._checklibSBML(reac.setId(reacId), 'set reaction id') #same convention as cobrapy
        #self._checklibSBML(reac.setName(str(reacId)+), 'set name') #same convention as cobrapy
        self._checklibSBML(reac.setSBOTerm(185), 'setting the system biology ontology (SBO)') #set as process
        #TODO: consider having the two parameters as input to the function
        self._checklibSBML(reac.setReversible(True), 'set reaction reversibility flag')
        self._checklibSBML(reac.setFast(False), 'set reaction "fast" attribute')
        if metaID==None:
            metaID = self._genMetaID(reacId)
        self._checklibSBML(reac.setMetaId(metaID), 'setting species metaID')
        #reactants_dict
        for reactant in step['left']:
            spe = reac.createReactant()
            self._checklibSBML(spe, 'create reactant')
            #use the same writing convention as CobraPy
            if compartmentId:
                self._checklibSBML(spe.setSpecies(str(reactant)+'__64__'+str(compartmentId)), 'assign reactant species')
            else:
                self._checklibSBML(spe.setSpecies(str(reactant)+'__64__'+str(self.compartmentId)), 'assign reactant species')
            #TODO: check to see the consequences of heterologous parameters not being constant
            self._checklibSBML(spe.setConstant(True), 'set "constant" on species '+str(reactant))
            self._checklibSBML(spe.setStoichiometry(float(step['left'][reactant])),
                'set stoichiometry ('+str(float(step['left'][reactant]))+')')
        #products_dict
        for product in step['right']:
            pro = reac.createProduct()
            self._checklibSBML(pro, 'create product')
            if compartmentId:
                self._checklibSBML(pro.setSpecies(str(product)+'__64__'+str(compartmentId)), 'assign product species')
            else:
                self._checklibSBML(pro.setSpecies(str(product)+'__64__'+str(self.compartmentId)), 'assign product species')
            #TODO: check to see the consequences of heterologous parameters not being constant
            self._checklibSBML(pro.setConstant(True), 'set "constant" on species '+str(product))
            self._checklibSBML(pro.setStoichiometry(float(step['right'][product])),
                'set the stoichiometry ('+str(float(step['right'][product]))+')')
        #annotation
        annotation = '''<annotation>
  <rdf:RDF
  xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
  xmlns:bqbiol="http://biomodels.net/biology-qualifiers/"
  xmlns:bqmodel="http://biomodels.net/model-qualifiers/">'''
        # if the name of the species is MNX then we annotate it using MIRIAM compliance
        #TODO: need to add all known xref from different databases (not just MetaNetX)
        annotation += '''
    <rdf:Description rdf:about="#'''+str(metaID or '')+'''">
      <bqbiol:is>
        <rdf:Bag>'''
        id_ident = {'mnx': 'metanetx.reaction/', 'rhea': 'rhea/', 'reactome': 'reactome/', 'bigg': 'bigg.reaction/', 'sabiork': 'sabiork.reaction/', 'ec-code': 'ec-code/', 'biocyc': 'biocyc/'}
        if reacId in reacXref:
            for dbId in reacXref[reacId]:
                for cid in reacXref[reacId][dbId]:
                    try:
                        annotation += '''
          <rdf:li rdf:resource="http://identifiers.org/'''+str(id_ident[dbId])+str(cid)+'''"/>'''
                    except KeyError:
                        continue
        annotation += '''
        </rdf:Bag>
      </bqbiol:is>
    </rdf:Description>'''
        annotation += '''
    <rdf:Ibisba rdf:about="#'''+str(metaID or '')+'''">
      <ibisba:ibisba xmlns:ibisba="http://ibisba.eu">
        <ibisba:smiles>'''+str(reaction_smiles or '')+'''</ibisba:smiles>
        <ibisba:rule_id>'''+str(step['rule_id'] or '')+'''</ibisba:rule_id>
        <ibisba:rule_score value="'''+str(step['rule_score'] or '')+'''" />
      </ibisba:ibisba>
    </rdf:Ibisba>
  </rdf:RDF>
</annotation>'''
        # Took this out of above
        #<ibisba:fba_biomass_score value="" />
        #<ibisba:fba_target_score value="" />
        #<ibisba:fba_splitObj_score value="" />
        self._checklibSBML(reac.setAnnotation(annotation), 'setting annotation for reaction '+str(reacId))
        #### GROUPS #####
        if not hetero_group==None:
            newM = hetero_group.createMember()
            self._checklibSBML(newM, 'Creating a new groups member')
            self._checklibSBML(newM.setIdRef(reacId), 'Setting name to the groups member')
        elif not self.hetero_group==None:
            newM = self.hetero_group.createMember()
            self._checklibSBML(newM, 'Creating a new groups member')
            self._checklibSBML(newM.setIdRef(reacId), 'Setting name to the groups member')
        else:
            logging.warning('This pathway is not added to a particular group')


    ## Create libSBML reaction
    #
    # Create a reaction. fluxBounds is a list of libSBML.UnitDefinition, length of exactly 2 with the first position that is the upper bound and the second is the lower bound. reactants_dict and reactants_dict are dictionnaries that hold the following parameters: name, compartmentId, stoichiometry
    #
    # @param chemId Species name (as of now, can only handle MNX ids)
    # @param metaID Name for the reaction
    # @param inchi String Inchi associated with this species
    # @param smiles String SMILES associated with this species
    # @param compartmentId String Set this species to belong to another compartmentId than the one globally set by self.compartmentId
    # @param charge Optional parameter describing the charge of the molecule of interest
    # @param chemForm Optional chemical formulae of the substrate (not SMILES or InChI)
    # @param dG Optinal Thermodynamics constant for this species
    # @param dG_uncert Optional Uncertainty associated with the thermodynamics of the reaction 
    def createSpecies(self, 
            chemId,
            chemXref, 
            metaID=None, 
            inchi=None,
            inchiKey=None,
            smiles=None,
            compartmentId=None):
            #charge=0,
            #chemForm=''):
        spe = self.model.createSpecies()
        self._checklibSBML(spe, 'create species')
        ##### FBC #####
        spe_fbc = spe.getPlugin('fbc')
        self._checklibSBML(spe_fbc, 'creating this species as an instance of FBC')
        #spe_fbc.setCharge(charge) #### These are not required for FBA 
        #spe_fbc.setChemicalFormula(chemForm) #### These are not required for FBA
        if compartmentId:
            self._checklibSBML(spe.setCompartment(compartmentId), 'set species spe compartment')
        else:
            self._checklibSBML(spe.setCompartment(self.compartmentId), 'set species spe compartment')
        #ID same structure as cobrapy
        #TODO: determine if this is always the case or it will change
        self._checklibSBML(spe.setHasOnlySubstanceUnits(False), 'set substance units')
        self._checklibSBML(spe.setBoundaryCondition(False), 'set boundary conditions')
        self._checklibSBML(spe.setConstant(False), 'set constant')
        #useless for FBA (usefull for ODE) but makes Copasi stop complaining
        self._checklibSBML(spe.setInitialConcentration(1.0), 'set an initial concentration')
        #same writting convention as COBRApy
        self._checklibSBML(spe.setId(str(chemId)+'__64__'+str(compartmentId)), 'set species id')
        if metaID==None:
            metaID = self._genMetaID(chemId)
        self._checklibSBML(spe.setMetaId(metaID), 'setting reaction metaID')
        self._checklibSBML(spe.setName(chemId), 'setting name for the namebolites')
        ###### annotation ###
        annotation = '''<annotation>
  <rdf:RDF 
  xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" 
  xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" 
  xmlns:bqmodel="http://biomodels.net/model-qualifiers/">'''
        # if the name of the species is MNX then we annotate it using MIRIAM compliance
        #TODO: need to add all known xref from different databases (not just MetaNetX)
        annotation += '''
    <rdf:Description rdf:about="#'''+str(metaID or '')+'''">
      <bqbiol:is>
        <rdf:Bag>'''
        id_ident = {'mnx': 'metanetx.chemical/', 'chebi': 'chebi/CHEBI:', 'bigg': 'bigg.metabolite/', 'hmdb': 'hmdb/', 'kegg_c': 'kegg.compound/', 'kegg_d': 'kegg.drug/', 'biocyc': 'biocyc/META:', 'seed': 'seed.compound/', 'metacyc': 'metacyc/', 'sabiork': 'seed.compound/', 'reactome': 'reactome.compound/'}
        if chemId in chemXref: 
            for dbId in chemXref[chemId]:
                for cid in chemXref[chemId][dbId]:
                    try:
                        if dbId == 'kegg' and cid[0] == 'C':
                            annotation += '''        
          <rdf:li rdf:resource="http://identifiers.org/'''+id_ident['kegg_c']+str(cid)+'''"/>'''
                        elif dbId == 'kegg' and cid[0] == 'D':
                            annotation += '''        
          <rdf:li rdf:resource="http://identifiers.org/'''+id_ident['kegg_d']+str(cid)+'''"/>'''
                        else:
                            annotation += '''        
          <rdf:li rdf:resource="http://identifiers.org/'''+str(id_ident[dbId])+str(cid)+'''"/>'''
                    except KeyError:
                        continue
        annotation += '''
        </rdf:Bag>
      </bqbiol:is>
    </rdf:Description>'''   
        ###### IBISBA additional information ########
        annotation += '''
    <rdf:Ibisba rdf:about="#'''+str(metaID or '')+'''">
      <ibisba:ibisba xmlns:ibisba="http://ibisba.eu/qualifiers">
        <ibisba:smiles>'''+str(smiles or '')+'''</ibisba:smiles>
        <ibisba:inchi>'''+str(inchi or '')+'''</ibisba:inchi>
        <ibisba:inchikey>'''+str(inchiKey or '')+'''</ibisba:inchikey>
      </ibisba:ibisba>
    </rdf:Ibisba>'''
        annotation += '''
  </rdf:RDF>
</annotation>'''
        self._checklibSBML(spe.setAnnotation(annotation), 'setting the annotation for new species')

############
#Removed from the above
#        <ibisba:dG_prime_o units="kj_per_mol" value=""/>
#        <ibisba:dG_prime_m units="kj_per_mol" value=""/>
#        <ibisba:dG_uncert units="kj_per_mol" value=""/>

    ## Create libSBML pathway
    #
    # Create the collection of reactions that constitute the pathway using the Groups package and create the custom IBIBSA annotations
    #
    # @param model libSBML model to add the unit definition
    # @param reaction_id Reaction ID
    # @param name Name for the reaction
    # @param fluxBounds list of size 2 that describe the FBC upper and lower bounds for this reactions flux
    # @param reactants list of species that are the reactants of this reaction
    # @param products list of species that are the products of this reaction
    # @param reaction_smiles String smiles description of this reaction (added in IBISBA annotation)
    # @return hetero_group The number libSBML groups object to pass to createReaction to categorise the new reactions
    def createPathway(self, path_id, metaID=None):
        groups_plugin = self.model.getPlugin('groups')
        self.hetero_group = groups_plugin.createGroup()
        self.hetero_group.setId('rp_pathway')
        if metaID==None:
            metaID = self._genMetaID('rp_pathway')
        self.hetero_group.setMetaId(metaID)
        self.hetero_group.setKind(libsbml.GROUP_KIND_COLLECTION)
        annotation = '''<annotation>
  <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" 
  xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
    <rdf:Ibisba rdf:about="#'''+str(metaID or '')+'''">
      <ibisba:ibisba xmlns:ibisba="http://ibisba.eu">
      </ibisba:ibisba>
    </rdf:Ibisba>
  </rdf:RDF>
</annotation>'''
        self.hetero_group.setAnnotation(annotation)
        ##Taken from above
        #<ibisba:dG_prime_o units="kj_per_mol" value=""/>
        #<ibisba:dG_prime_m units="kj_per_mol" value=""/>
        #<ibisba:dG_uncert units="kj_per_mol" value=""/>
        #<ibisba:fba_biomass_score value=""/>
        #<ibisba:fba_target_score value=""/>
        #<ibisba:fba_splitObj_score value=""/>


    ## Create libSBML gene
    #
    # Create the list of genes in the model including its custom IBISBA annotatons
    #
    # @param model libSBML model to add the unit definition
    # @param reac libSBML reaction object
    # @param step_id The step for the number of 
    # @return libSBML gene object
    def createGene(self, reac, step_id, metaID=None):
        #TODO: pass this function to Pablo for him to fill with parameters that are appropriate for his needs
        geneName = 'RP'+str(step_id)+'_gene'
        fbc_plugin = self.model.getPlugin('fbc')
        #fbc_plugin = reac.getPlugin("fbc")
        gp = fbc_plugin.createGeneProduct()
        gp.setId(geneName)
        if metaID==None:
            metaID = self._genMetaID(str(geneName))
        gp.setMetaId(metaID)
        gp.setLabel('gene_'+str(step_id))
        gp.setAssociatedSpecies('RP'+str(step_id))
        ##### NOTE: The parameters here require the input from Pablo to determine what he needs
        annotation = '''<annotation>
  <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" 
        xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
    <rdf:Ibisba rdf:about="#'''+str(metaID or '')+'''">
      <ibisba:ibisba xmlns:ibisba="http://ibisba.eu">
        <ibisba:fasta value="" />
      </ibisba:ibisba>
    </rdf:Ibisba>
  </rdf:RDF>
</annotation>'''
        gp.setAnnotation(annotation)
    

    ## Create libSBML flux objective 
    #
    # Using the FBC package one can add the FBA flux objective directly to the model. This function sets a particular reaction as objective with maximization or minimization objectives
    #
    # @param model libSBML model to add the unit definition
    # @param fluxObjID The id given to this particular objective
    # @param reactionName The name or id of the reaction that we are setting a flux objective
    # @param coefficient FBA coefficient 
    # @param isMax Boolean to determine if we are maximizing or minimizing the objective
    # @param metaID Set the metaID
    # @return Boolean exit code
    def createFluxObj(self, fluxObjID, reactionName, coefficient, isMax=True):#, metaID=None):
        #TODO: define more complex objectives and test the simulation using a libSBML object
        fbc_plugin = self.model.getPlugin('fbc')
        target_obj = fbc_plugin.createObjective()
        target_obj.setId(fluxObjID)
        #if metaID==None:
        #    metaID = self._genMetaID(fluxObjID)
        #target_obj.setMetaId(metaID)
        if isMax:
            target_obj.setType('maximize')
        else:
            target_obj.setType('minimize')
        fbc_plugin.setActiveObjectiveId(fluxObjID) # this ensures that we are using this objective when multiple
        target_flux_obj = target_obj.createFluxObjective()
        target_flux_obj.setReaction(reactionName)
        target_flux_obj.setCoefficient(coefficient)


    ## Generate a generic model 
    #
    # Since we will be using the same type of parameters for the RetroPath model, this function
    # generates a libSBML model with parameters that will be mostly used
    #
    #
    #
    def genericModel(self, modelName, modelID, compXref):
        self.createModel(modelName, modelID)
        # mmol_per_gDW_per_hr
        unitDef = self.createUnitDefinition('mmol_per_gDW_per_hr')
        self.createUnit(unitDef, libsbml.UNIT_KIND_MOLE, 1, -3, 1)
        self.createUnit(unitDef, libsbml.UNIT_KIND_GRAM, 1, 0, 1)
        self.createUnit(unitDef, libsbml.UNIT_KIND_SECOND, 1, 0, 3600)
        # kj_per_mol
        gibbsDef = self.createUnitDefinition('kj_per_mol')
        self.createUnit(gibbsDef, libsbml.UNIT_KIND_JOULE, 1, 3, 1)
        self.createUnit(gibbsDef, libsbml.UNIT_KIND_MOLE, 1, 1, 1)
        # infinity parameters (FBA)
        #upInfParam = self.createParameter('B_INF', float('inf'), 'kj_per_mol')
        #lowInfParam = self.createParameter('B__INF', float('-inf'), 'kj_per_mol')
        upNineParam = self.createParameter('B__999999', -999999, 'mmol_per_gDW_per_hr')
        lowNineParam = self.createParameter('B_999999', 999999, 'mmol_per_gDW_per_hr')
        lowZeroParam = self.createParameter('B_0', 0, 'mmol_per_gDW_per_hr')
        #compartment
        #TODO: create a new compartment 
        self.createCompartment(1, 'MNXC3', 'cytoplasm', compXref) 


    ##################################################################################################
    ########################################### TESTING ##############################################
    ##################################################################################################


    ## Test writing of an SBML model
    #
    # Private test function that will write a given pathway to file
    #
    # @param cofactors_rp_paths Dictionnary of the outputs from RP2paths
    # @param path_id Integer The cofactors_rp_paths pathway to write to SBML
    def _testWriter(self, cofactors_rp_paths=None, path_id=None):
        #### Define the if global or local parameter is to be used
        try:
            if cofactors_rp_paths==None and self.cofactors_rp_paths==None:
                raise TypeError 
            if cofactors_rp_paths==None and not self.cofactors_rp_paths==None:
                cofactors_rp_paths = self.cofactors_rp_paths
        except EmptyOutRPpaths:
            logging.error('Both class and input cofactors_rp_paths are empty')
        #### Check what is the best path to use
        if path_id==None:
            path_id = 1
        return None

