import cobra
import libsbml

#import rpSBML

## Class to simulate an rpsbml object using different FBA types and objective functions
#
# At this point we want to have the BIOMASS, target and shared objective
#TODO: add the pareto frontier optimisation as an automatic way to calculate the optimal fluxes
class rpFBA:
    def __init__(self, rpsbml):
        self.rpsbml = rpsbml
        self.cobraModel = None
  

    ##
    #
    #   
    def allObj(self):
        fbc_plugin = self.rpsbml.model.getPlugin('fbc')
        for objId in [i.getId() for i in fbc_plugin.getListOfObjectives()]:
            groups = self.rpsbml.model.getPlugin('groups')
            rp_pathway = groups.getGroup('rp_pathway')
            #find all the species of the reactions in the rp_pathway to return the shadow price
            ''' TO BE DETERMINED IF USED 
            mem = []
            for member in rp_pathway.getListOfMembers():
                reac = self.rpsbml.model.getReaction(member.getIdRef())
                for pro in reac.getListOfProducts():
                    mem.append(pro.getSpecies())
                for rea in reac.getListOfReactants():
                    mem.append(rea.getSpecies())
            '''
            #run the FBA
            #fbc_plugin = rpsbml.model.getPlugin('fbc')
            fbc_plugin.setActiveObjectiveId(objId)
            cobraModel = cobra.io.read_sbml_model(self.rpsbml.document.toXMLNode().toXMLString(), use_fbc_package=True)
            res = cobraModel.optimize()
            #update the annotations with the reaction flux results
            #reac.getChild('RDF').getChild('Ibisba').getChild('ibisba').getChild(child).addAttr('value', str(value))
            #update the annotations with the species shadow prices
            for member in rp_pathway.getListOfMembers():
                reac = self.rpsbml.model.getReaction(member.getIdRef())
                reac_annot = reac.getAnnotation()
                ibisba_annot = reac_annot.getChild('RDF').getChild('Ibisba').getChild('ibisba')
                #NOTE: for some reaseon one neesds to wrap a child XMLnode into a parent node for it to be valid
                #here we just pass the child node to be added
                tmpAnnot = libsbml.XMLNode.convertStringToXMLNode('<ibisba:ibisba xmlns:ibisba="http://ibisba.eu"> <ibisba:fba_'+str(objId)+' units="mmol_per_gDW_per_hr" value="'+str(res.fluxes.get(member.getIdRef()))+'" /> </ibisba:ibisba>')
                ibisba_annot.addChild(tmpAnnot.getChild('fba_'+str(objId)))
            ''' TO BE DETERMINED IF USED
            #update the shadow prices for species
            for speName in list(set(mem)): #remoce duplicates
                #only choose the heterologous species
                if len([x for x in speName.split('_') if x])==4:
                    spe = self.rpsbml.model.getSpecies(speName)
                    spe_annot = spe.getAnnotation()
                    ibisba_annot = spe_annot.getChild('RDF').getChild('Ibisba').getChild('ibisba')
                    tmpAnnot = libsbml.XMLNode.convertStringToXMLNode('<ibisba:ibisba xmlns:ibisba="http://ibisba.eu"> <ibisba:shadow_price_'+str(objId)+' units="mmol_per_gDW_per_hr" value="'+str(res.shadow_prices.get(speName))+'" /> </ibisba:ibisba>')
                    ibisba_annot.addChild(tmpAnnot.getChild('shadow_price_'+str(objId)))
            '''


    ## function to add the splitObj or any other objective that does not exist in the model
    #
    #
    def addObjective():
        pass       


    ## Given that there can be multiple objectives defined in the SBML, this function switches between different ones
    #
    # TODO: 
    def switchObjective(self, objId):
        fbc_plugin = self.rpsbml.model.getPlugin('fbc')
        listObj = fbc_plugin.getListOfObjectives()
        if objId in [i.getId() for i in listObj]:
            listObj.setActiveObjective(objId)
            return objId
        else:
            logger.warning('The objective Id '+str(objId)+' does not exist in rpsbml')
            return None
        self.libsbml_to_cobra()


    ##
    #
    #
    def libsbml_to_cobra(self):
        self.cobraModel = cobra.io.read_sbml_model(document.toXMLNode().toXMLString())
        #for an old version of cobrapy (0.4)
        #return cobra.io.sbml3.parse_xml_into_model(lxml.etree.fromstring(rpsbml.document.toXMLNode().toXMLString()))
        
    
    ##
    #
    #TODO: should update the 
    def simulate(self):
        res = self.cobraModel.optimize()
        return self.cobraModel.summary()
        
    
    ##
    #
    # child: rule_score, fba_biomass_score, fba_target_score, fba_splitObj_score, global_score
    def updateFBAscore(self, child, value):
        self.rpsbml()
        #reaction
        for reac in self.rpsbml.getListOfReactions():
            reac.getChild('RDF').getChild('Ibisba').getChild('ibisba').getChild(child).addAttr('value', str(value))
        #pathway


    ########################################################################
    ############################### FBA pathway ranking ####################
    ########################################################################

    #1) Number of interventions
    # need to calculate the number of steps that are not native to know the number of interventions

    #2) Maximal growth rate

    #3) Minimum product yeild at maximal growth rate

    #4) Minimum product yeild

    #5) Anaerobic condition

    #6) Number of potentially disruptive products

        #Toxicity?

    #7) Number of accessible metabolites (avoid intermediate accumulation)

    #8) Thermodynamics (MDF)

    #9) The overlap of the same changes --> might not be applicable in our case

    #10) Reduced model

    #11) ECM


