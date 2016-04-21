# -*- coding: utf-8 -*-
"""
Definition of general helper functions to create the various
objects in the SBML. 
These functions are called with the information dictionaries 
during the generation of cell and tissue model.
"""
# TODO: refactor code in the sbmlutils

import libsbml
from libsbml import UNIT_KIND_DIMENSIONLESS, UnitKind_toString

SBML_LEVEL = 3
SBML_VERSION = 1
  
def createUnitDefinition(model, uid, units):
    ''' Creates the defined unit definitions. '''
    unitdef = model.createUnitDefinition()
    unitdef.setId(uid)
    for data in units:
        kind = data[0]
        exponent = data[1]
        _createUnit(unitdef, kind, exponent)
             
def _createUnit(unitdef, kind, exponent, scale=0, multiplier=1.0):
    unit = unitdef.createUnit()
    unit.setKind(kind)
    unit.setExponent(exponent)
    unit.setScale(scale)
    unit.setMultiplier(multiplier)
            
            
def setMainUnits(model, main_units):
    ''' Sets the main units for the model. '''
    for key in ('time', 'extent', 'substance', 'length', 'area', 'volume'):
        unit = main_units[key]
        unit = getUnitString(unit)
        # set the values
        if key == 'time':
            model.setTimeUnits(unit)
        elif key == 'extent':
            model.setExtentUnits(unit)
        elif key == 'substance':
            model.setSubstanceUnits(unit)
        elif key == 'length':
            model.setLengthUnits(unit)
        elif key == 'area':
            model.setAreaUnits(unit)
        elif key == 'volume':
            model.setVolumeUnits(unit)
    
def getUnitString(unit):
    if type(unit) is int:
        unit = UnitKind_toString(unit)
    if unit == '-':
        unit = UnitKind_toString(UNIT_KIND_DIMENSIONLESS)
    return unit

def createParameters(model, pdict):
    for pdata in pdict.values():
        # id, name, value, unit, constant
        pid = pdata[0]
        name = pdata[1]
        value = pdata[2]
        unit = getUnitString(pdata[3])
        constant = pdata[4]
        createParameter(model, pid=pid, unit=unit, name=name, value=value, constant=constant)
        

def createParameter(model, pid, unit, name=None, value=None, constant=True):
    p = model.createParameter()
    p.setId(pid)
    p.setUnits(unit)
    if name != None:
        p.setName(name)
    if value != None:
        p.setValue(value)
    p.setConstant(constant)
    return p
    
def createCompartments(model, comps):
    for cid in sorted(comps):
        # comps['PV'] = ('[PV] perivenious', 3, 'm3', True, 'value)
        data = comps[cid]
        name = data[0]
        dims = data[1]
        units = data[2]
        constant = data[3]
        value = data[4]
        _createCompartment(model, cid=cid, name=name, dims=dims, 
                                  units=units, constant=constant, value=value)
           
def _createCompartment(model, cid, name, dims, units, constant, value):
    c = model.createCompartment()
    c.setId(cid)
    if name:
        c.setName(name)
    c.setSpatialDimensions(dims)            
    c.setUnits(units)
    c.setConstant(constant)
    if type(value) is str:
        # _createInitialAssignment(model, sid=cid, formula=value)
        _createAssignmentRule(model, sid=cid, formula=value)
    else:
        c.setValue(value)
    
def createSpecies(model, sdict):
    for sid in sorted(sdict):
        data = sdict[sid]
        _createSpecie(model, sid, name=data[0], init=data[1], units=data[2], 
                      compartment=data[3], boundaryCondition=data[4])
    
def _createSpecie(model, sid, name, init, units, compartment, boundaryCondition):
    s = model.createSpecies()
    s.setId(sid)
    if name:
        s.setName(name)
    s.setInitialConcentration(init)
    s.setUnits(units)
    s.setCompartment(compartment)
        
    s.setSubstanceUnits(model.getSubstanceUnits());
    s.setHasOnlySubstanceUnits(False);
    s.setConstant(False)
    s.setBoundaryCondition(boundaryCondition)

## InitialAssignments ##      
def createInitialAssignments(model, assignments, names):
    for data in assignments:
        # id, assignment, unit
        pid = data[0]
        unit = getUnitString(data[2])
        # Create parameter if not existing
        if (not model.getParameter(pid)) and (not model.getSpecies(pid)):
            createParameter(model, pid, unit, name=names.get(pid, None), value=None, constant=True)
        _createInitialAssignment(model, sid=pid, formula=data[1])
    
def _createInitialAssignment(model, sid, formula):
    assignment = model.createInitialAssignment();
    assignment.setSymbol(sid);
    astnode = libsbml.parseL3FormulaWithModel(formula, model)
    assignment.setMath(astnode);   
         
## AssignmentRules ##
def createAssignmentRules(model, rules, names):
    for data in rules:
        # id, rule, unit
        pid = data[0]
        unit = getUnitString(data[2])
        # Create parameter if not existing
        if (not model.getParameter(pid)) and (not model.getSpecies(pid)):
            createParameter(model, pid, unit, name=names.get(pid, None), value=None, constant=False)
        if (not model.getRule(pid)):
            _createAssignmentRule(model, sid=pid, formula=data[1])
            
def _createAssignmentRule(model, sid, formula):
    rule = model.createAssignmentRule();
    rule.setVariable(sid);
    astnode = libsbml.parseL3FormulaWithModel(formula, model)
    rule.setMath(astnode);   

######################################
# Deficiency Events (Galactosemias)
######################################
def getDeficiencyEventId(deficiency):
    return 'EDEF_{:0>2d}'.format(deficiency)
    
def createDeficiencyEvent(model, deficiency):
    eid = getDeficiencyEventId(deficiency)
    e = model.createEvent();
    e.setId(eid);
    e.setUseValuesFromTriggerTime(True);
    t = e.createTrigger()
    t.setInitialValue(False) # ! not supported by Copasi -> lame fix via time
    t.setPersistent(True)    # ! not supported by Copasi -> careful with usage
    formula = '(time>0) && (deficiency=={:d})'.format(deficiency);
    astnode = libsbml.parseL3FormulaWithModel(formula, model)
    t.setMath(astnode);
    return e;

##########################################
# Simulation Events (Peaks & Challenges)
##########################################
def createSimulationEvents(model, elist):
    for edata in elist:
        createEventFromEventData(model, edata)

def createEventFromEventData(model, edata):
    e = model.createEvent();
    e.setId(edata.eid)
    e.setName(edata.name);
    e.setUseValuesFromTriggerTime(True);
    t = e.createTrigger();
    t.setInitialValue(False)
    t.setPersistent(True)
    astnode = libsbml.parseL3FormulaWithModel(edata.trigger, model)
    t.setMath(astnode)
    # assignments
    for key, value in edata.assignments.iteritems(): 
        astnode = libsbml.parseL3FormulaWithModel(value, model)
        ea = e.createEventAssignment()
        ea.setVariable(key)
        ea.setMath(astnode)