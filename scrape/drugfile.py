#!/usr/bin/python

# import pprint


class DrugModel:
    def __init__(self, soup):

        self.source = soup

        self.drugId = soup.drugModel.drugId.string
        self.drugModelId = soup.drugModel.drugModelId.string
        self.covariates = []
        for cov in soup.drugModel.covariates.find_all('covariate'):
            self.covariates.append(Covariate(cov))
        self.formulationAndRoutes = []
        for cov in soup.drugModel.formulationAndRoutes.find_all('formulationAndRoute'):
            self.formulationAndRoutes.append(FormulationAndRoute(cov))

        for analyteGroup in soup.drugModel.analyteGroups.find_all('analyteGroup'):
            analytes = analyteGroup.analytes
            for analyte in analytes.find_all('analyte'):
                self.errorModel = analyte.errorModel.errorModelType.string

        # self.atc = soup.head.atc.string
        # self.drugid = soup.head.drugid.string
        # self.modelid = soup.head.modelid.string
        # self.drugnames = []
        # for dn in soup.head.drugnames.find_all('name'):
        # 	self.drugnames.append(dn.string)
        # self.intake = soup.adme.intake.string
        # self.distribution = soup.adme.distribution
        # self.elimination = soup.adme.elimination
        # self.doses = []
        # for d in soup.dosages.doses.find_all('dose'):
        # 	self.doses.append(d.string)
        # self.defaultdose = soup.doses['default']
        # self.defaultunit = soup.doses['unit']
        # self.intervals = []
        # for i in soup.intervals.find_all('interval'):
        # 	self.intervals.append(i.string)
        # self.defaultinterval = soup.intervals['default']
        # self.infusions = []
        # for inf in soup.infusions.find_all('infusion'):
        # 	self.infusions.append(inf.string)
        # self.defaultinfusion = soup.infusions['default']
        # self.targets = []
        # for tgt in soup.targets.find_all('target'):
        # 	self.targets.append(Target(tgt))
        # self.globaleps = soup.errormodel.additive.string
        # self.globaleta = soup.errormodel.proportional.string
        # self.parameters = []
        # for prm in soup.parameters.find_all('parameter'):
        # 	self.parameters.append(Parameter(prm))
        # self.operations = []
        # for op in soup.operations.find_all('operation'):
        # 	self.operations.append(Operation(op))


class FormulationAndRoute:
    def __init__(self, f):
        self.formulationAndRouteId = f.formulationAndRouteId.string
        self.formulation = f.formulation.string
        self.administrationName = f.administrationName.string
        self.administrationRoute = f.administrationRoute.string
        self.absorptionModel = f.absorptionModel.string
        self.defaultdose = float(f.dosages.availableDoses.findChild("default").standardValue.string)
        if self.absorptionModel == 'infusion':
            self.defaultinfusion = float(f.dosages.availableInfusions.findChild("default").standardValue.string)


class Parameter:
    def __init__(self, prm):
        self.id = prm.id  # check
        self.unit = prm.unit  # check
        self.value = prm.value  # check
        self.eps = prm.bsv.additive.string  # check
        self.eta = prm.bsv.proportional.string  # check
        self.min = prm.min.string
        self.max = prm.max.string


class Operation:
    def __init__(self, op):
        self.parameter = op.parameter.string
        self.formula = op.formula.string


class Covariate:
    def __init__(self, cov):
        self.covariateId = cov.covariateId.string
        self.unit = cov.unit.string
# check for valid unit
        self.covariateType = cov.covariateType.string
        self.dataType = cov.dataType.string
# do a check for valid type
        self.defaultvalue = cov.covariateValue.standardValue.string
# check for valid defaultvalue


class Target:
    def __init__(self, tgt):
        self.type = tgt.type.string
        self.concentrations = []
        for conc in tgt.find_all('concentrations'):
            self.concentrations.append(Concentration(conc))
        self.times = []
        for time in tgt.find_all('times'):
            self.times.append(Time(time))


class Concentration:
    def __init__(self, conc):
        self.min = conc.min.value.string
        self.max = conc.max.value.string
        self.best = conc.best.value.string
        self.unit = conc['unit']


class Time:
    def __init__(self, time):
        self.min = time.min.value.string
        self.max = time.max.value.string
        self.best = time.best.value.string
        self.unit = time['unit']
