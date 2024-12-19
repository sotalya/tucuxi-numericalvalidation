#!/usr/bin/python

# import datetime
# import random

from scrape.drugfile import *
from sotalya.data.query import Query
import os
from bs4 import BeautifulSoup
from colorama import init, Fore, Back, Style

# For colorama, to reset the color on every new print
init(autoreset=True)


class SystemTester:
    """
    * This test program is supposed to be greedy about which tests to run. There
      are overlaps in that the specifying arguments may specify the same conditions
      more than once, and we need to filter them out. One argument may affect the
      effect of another. Mathematically, its the outer product of the argument criteria
    * e.g. If only a drug is specified, test the drug with all patients, models greedily.
           But if only a patient is specified, test the patient with
           all applicable drugs and models etc.
    """
    def __init__(self):
        """
        * System tester is the builder/manager pattern for system tests.
        * This class needs access to datasets and drugs, so it uses dictionaries
          with ids as keys.
        """
        self.query_dict = {}
        self.drug_dict = {}

    def check_conflicts(self, patients, drugs, models):
        """
        * The command line is greedy so it wants the outer product of all the arguments
          to test under all conditions specified. This function determines which test cases
          to run.

            Args:

            * patients ([str]): The list of parsed patient ids.
            * drugs ([str]): The list of parsed drug ids.
            * models ([str]): The list of models.
        """
        if patients:
            for p in patients:
                if p in self.query_dict:
                    if self.query_dict[p].drugid not in self.drug_dict:
                        return 'No data for drug: {drug} of patient {patient}'.format(drug=self.query_dict[p].drugid,
                                                                                      patient=p)
                else:
                    return 'No data for patient: {patient}'.format(patient=p)

        if drugs:
            for d in drugs:
                if d not in self.drug_dict:
                    return 'No data for drug: {drug}'.format(drug=d)

        if models:
            for m in models:
                flag = 0
                for k, v in iter(self.drug_dict.items()):
                    if v.modelid.replace('ezechiel.model', '') in m:
                        flag = 1
                if flag == 0:
                    return 'No drugfiles loaded with model: {model}'.format(model=m)
        return 'none'

    def import_query(self, filenamepath: str, filename: str):
        """
        Imports query file by trying to parse the tqf file.

        :param str filenamepath: The file name containing the query
        :return: A query
        """

        base, ext = os.path.splitext(filenamepath)
        if ext == '.tqf':

            content = open(filenamepath).read()
            soup = BeautifulSoup(content, 'xml')

            if soup.query:
                try:
                    query = Query(soup)
                    query.sourceFile = os.path.abspath(filenamepath)
                    self.query_dict[filename] = query
                    print("Imported query file " + filenamepath)
                except Exception as e:
                    print(Fore.RED + 'Can not import the following query file : ' + filenamepath)
                    print(e)
            else:
                print(Fore.RED + "Can not import the following query file : " + filenamepath)

    def import_queries(self, folderpaths):
        """
        * Imports query files by trying to parse each xml file in the folderpaths.
            Args:
                folderpaths ([str]): The list of directories to search for xml files.
        """
        if folderpaths:
            for fpath in folderpaths:
                for dirname, dirnames, filenames in os.walk(fpath):
                    for filename in sorted(filenames):
                        self.import_query(os.path.join(dirname, filename), filename)

    def import_drugs(self, folderpaths):
        """
        * Imports drugfiles by trying to parse each xml file in the folderpaths.
            Args:
                folderpaths ([str]): The list of directories to search for xml files.
        """

        # print(folderpaths)
        # pdb.set_trace()
        for fpath in folderpaths:
            for dirname, dirnames, filenames in os.walk(fpath):
                for filename in filenames:
                    base, ext = os.path.splitext(filename)
                    if ext == '.tdd':
                        content = open(os.path.join(dirname, filename)).read()
                        soup = BeautifulSoup(content, 'xml')
                        if soup.drugModel:
                            # print(soup.drugModel.drugModelId.text)
                            soup = BeautifulSoup(content, 'xml')
                            drug = DrugModel(soup)
                            drug.sourceFile = os.path.abspath(os.path.join(dirname, filename))
                            self.drug_dict[drug.drugModelId] = drug
                        else:
                            print("Can not import the following drug file : " + filename)

    def get_drug_model_by_drug_id(self, drug_id):
        """
        * Iterate over the drug models to select the ones dealing with the required drug.
        """
        result = []
        for k, v in list(self.drug_dict.items()):
            if v.drugId == drug_id:
                result.append(v)
        return result

    # def generate_datasets(self, rp, rs, ss):
    #     """
    #     * DEPRECATED
    #     * Generates fake datasets for a certain number of random patients, random samples
    #       for whether its steady state or not.
    #
    #         Args:
    #             * rp (int): The number of random patients to generate
    #             * rs (int): The number of random samples to generate
    #             * ss (bool): If is steady state
    #
    #     """
    #     for i in range(0, rp):
    #         for k, v in iter(self.drug_dict.items()):
    #             ds = self.randomize_dataset_template(v, rs, ss)
    #             self.query_dict[ds.patient.id] = ds

    def filter_datasets(self, patients, drugs, models):
        """
        * After determining which cases to run with the check_conflicts
          method, this method filters the dataset dictionary according
          to which test cases to run by removing duplicates.

            Args:
                * patients ([str]): The list of parsed patient ids.
                * drugs ([str]): The list of parsed drug ids.
                * models ([str]): The list of models.

        """
        tobedel = []
        for k, v in iter(self.query_dict.items()):
            flag = 0
            if patients:
                for p in patients:
                    if p == k:
                        flag = 1
            if drugs:
                for d in drugs:
                    if d == v.drug.id:
                        flag = 1
            if models:
                for m in models:
                    if m == self.drug_dict[v.drug.id].modelid:
                        flag = 1
            if not patients and not drugs and not models:
                flag = 1

            if flag == 0:
                tobedel.append(k)

        for p in tobedel:
            del self.query_dict[p]

#     def randomize_dataset_template(self, drug, rs, ss):
#         """
#         * DEPRECATED
#         * Starting with an xml template ('template.xml') in the
#           test/system folder, this method generates random data
#           for a random patient and returns a dataset from it.
#
#             Args:
#                 * drug (DrugModel): The drug the random patient will take.
#                 * rs (int): The number of random samples to gen for this dataset.
#                 * ss (bool): If its at steadystate.
#
#             Returns:
#                 * ds (Query): The dataset filled with generated data.
#         """
#         bs = BeautifulSoup(open('template.xml'), 'xml')
#         ds = Query(bs)
#         ds.isgenerated = True
#         ds.drug.id = drug.drugid
#         ds.drug.atc = drug.atc
#
#         stime = datetime.datetime.now()
#         dosage = ds.dosages[0]
#         sample = ds.samples[0]
#         self.randomize_patient(ds, stime)
#         self.randomize_dosage(dosage, drug, stime)
#         self.randomize_sample(sample, stime, drug.defaultinterval, drug.targets, ss)
#         ds.dosages.append(dosage)
#
#         for i in range(1, rs):
#             # print i
#             sample = Sample(bs.sample)
#             self.randomize_sample(sample, stime, drug.defaultinterval, drug.targets, ss)
#             ds.samples.append(sample)
#
#         covariate = ds.covariates[0]
#         covariate.name = 'Weight'
#         covariate.date = stime.strftime('%Y-%m-%dT%H:%M:%S')
#         covariate.value = random.gauss(70, 5)
#
#         for cov in drug.covariates:
#             if cov.id == 'weight':
#                 ds.covariates.pop(0)
#             covariate = DatasetCovariate(bs.covariate)
#             covariate.name = cov.id
#             covariate.date = stime.strftime('%Y-%m-%dT%H:%M:%S')
#             if cov.type == 'int':
#                 covariate.value = str(int(round(random.gauss(float(cov.defaultvalue), float(cov.defaultvalue)
#                                                 * 0.1))))
#             elif cov.type == 'bool':
#                 covariate.value = random.choice(['0', '1'])
#             elif cov.type == 'double':
#                 covariate.value = str(random.gauss(float(cov.defaultvalue), float(cov.defaultvalue) * 0.1))
#                 if 'pna' in covariate.name:
#                     covariate.value = str(float(cov.defaultvalue))
#             ds.covariates.append(covariate)
#         return ds
#
#     @staticmethod
#     def randomize_dosage(dosage, drug, stime):
#         """
#         * DEPRECATED
#         * This function generates a dosage from the default values of the
#           dosage object.
#
#             Args:
#                 * dosage (Dosage): The dosage the random patient will take.
#                 * drug (DrugModel): The drug the random patient will take.
#                 * stime (Date): The time to use for the start of the dosage.
#         """
#         dosage.startdate = stime.strftime('%Y-%m-%dT%H:%M:%S')
#         dosage.lastdate = stime + datetime.timedelta(days=10)
#         dosage.dose = drug.defaultdose
#         dosage.doseunit = drug.defaultunit
#         dosage.interval = drug.defaultinterval
#         if 'infu' in drug.modelid or drug.drugid == 'ch.heig-vd.ezechiel.sunitinib':
#             dosage.infusion = 0 + int(round(random.gauss(30, 7)))
#             if drug.defaultinfusion != "0":
#                 dosage.infusion = int(round(random.gauss(float(drug.defaultinfusion),
#                                                          0.1 * (float(drug.defaultinfusion)))))
#         else:
#             dosage.infusion = '0'
#
#     @staticmethod
#     def randomize_patient(ds, starttime):
#         """
#         * DEPRECATED
#         * This function generates a patient from the default values of the
#           dosage object.
#
#             Args:
#                 * drug (DrugModel): The drug the random patient will take.
#                 * rs (int): The number of random samples to gen for this dataset.
#                 * starttime (Date): The time to use for the birthdate of the patient.
#         """
#         fstime = starttime.strftime('%d/%m/%Y-%H:%M')
#         ds.patient.id = 'GenPatient_ATC_{atc}_{stime}'.format(atc=ds.drug.atc, stime=starttime)
# #       ds.patient.id = 'GenPatient_ATC_{atc}_{stime}'.format(atc = ds.drug.atc, stime =
#                          int(round(random.gauss(30,7))))
#         ds.patient.firstname = 'GenPatient_ATC_{atc}_{stime}'.format(atc=ds.drug.atc, stime=starttime)
#         ds.patient.lastname = 'GenPatient_ATC_{atc}_{stime}'.format(atc=ds.drug.atc, stime=starttime)
#         ds.patient.birthdate = fstime
#
#     @staticmethod
#     def randomize_sample(sample, stime, interval, targets, ss):
#         """
#         * DEPRECATED
#         * This function generates a sample from the default values of the
#           sample object.
#
#             Args:
#                 * sample (Sample): The sample the random patient will take.
#                 * stime (Date): The time of the start of the dosage. A random
#                     quantity is added to this time to get the time of the sample
#                 * interval (double): The length of an interval in the dosage
#                 * targets ([Target]): The list of targets for this dataset
#                 * ss (bool): Whether its at steadystate
#         """
#         if ss == 2:
#             samtime = stime + datetime.timedelta(hours=random.gauss(round(int(interval) * 10.5), 2))
#         else:
#             samtime = stime + datetime.timedelta(hours=random.gauss(round(int(interval) * 1.5), 2))
#         sample.id = samtime.strftime('%Y-%m-%dT%H:%M:%S')
#         sample.sampledate = samtime.strftime('%Y-%m-%dT%H:%M:%S')
#         if targets:
#             if targets[0].type == 'mean':
#                 sample.concentration = round(random.gauss(float(targets[0].concentrations[0].best),
#                                                           float(targets[0].concentrations[0].best) * 0.1), 1)
#                 sample.unit = targets[0].concentrations[0].unit
#             if targets[0].type == 'peak':
#                 sample.concentration = round(random.gauss(float(targets[0].concentrations[0].best) * 0.8,
#                                                           float(targets[0].concentrations[0].best) * 0.1), 1)
#                 sample.unit = targets[0].concentrations[0].unit
#             if targets[0].type == 'residual':
#                 sample.concentration = round(random.gauss(float(targets[0].concentrations[0].best) * 1.2,
#                                                           float(targets[0].concentrations[0].best) * 0.1), 1)
#                 sample.unit = targets[0].concentrations[0].unit
#         else:
#             sample.concentration = round(random.gauss(1500, 200), 1)
#             sample.unit = 'ug/l'

    def filter_drugs(self):
        """
        * After determining which cases to run with the check_conflicts
          method, this method filters the drug dictionary according
          to which test cases to run by removing duplicates.
        """
        for k, v in list(self.drug_dict.items()):
            flag = 0
            for dk, dv in iter(self.query_dict.items()):
                for request in dv.requests:
                    if request.drugModelId == k:
                        flag = 1
            if flag == 0:
                del self.drug_dict[k]

    def print_drug(self):
        """
        * Prints all the drugs in the dictionary
        """
        print(Back.MAGENTA + "Printing Drugs")
        print(Back.MAGENTA + "**************")
        for k, v in iter(self.drug_dict.items()):
            # print('\n{name}\t\t:\t{value}'.format(name='ATC', value=v.atc))
            print(Fore.CYAN + "*********************************************************")
#            print('{name}\t\t:\t{value}'.format(name='drugId', value=v.drugId))
            print('{name}\t\t:\t{value}'.format(name='drugModelId', value=v.drugModelId))
#            print('{name}\t\t:\t{value}'.format(name='drugnames', value=v.drugnames))
#            print('{name}\t\t:\t{value}'.format(name='doses', value=v.doses))
#            print('{name}\t\t:\t{value}'.format(name='defaultdose', value=v.defaultdose))
#            print('{name}\t\t:\t{value}'.format(name='intervals', value=v.intervals))
#            print('{name}\t\t:\t{value}'.format(name='defaultinterval', value=v.defaultinterval))
#            print('{name}\t\t:\t{value}'.format(name='targets', value=len(v.targets)))
#            self.print_targets(v.targets)
            print("\n")

    def print_dataset(self):
        """
        * Prints all the datasets in the dictionary
        """
        print(Back.MAGENTA + "Printing Datasets")
        print(Back.MAGENTA + "*****************")
        for k, v in iter(self.query_dict.items()):
            print('\n{name}\t\t:\t{value}'.format(name='drugid', value=v.drug.id))
            print(Fore.CYAN + "*********************************************************")
            print('{name}\t\t:\t{value}'.format(name='ATC', value=v.drug.atc))
            print('{name}\t\t:\t{value}'.format(name='Dosages', value=len(v.dosages)))
            self.print_dosages(v.dosages)
            print('{name}\t\t:\t{value}'.format(name='Patient FirstName', value=v.patient.firstname))
            print('{name}\t\t:\t{value}'.format(name='Patient LastName', value=v.patient.lastname))
            print('{name}\t\t:\t{value}'.format(name='Patient Birthdate', value=v.patient.birthdate))
            print('{name}\t\t:\t{value}'.format(name='Covariates', value=len(v.covariates)))
            self.print_covariates(v.covariates)
            print('{name}\t\t:\t{value}'.format(name='Samples', value=len(v.samples)))
            self.print_samples(v.samples)
            print('{name}\t\t:\t{value}'.format(name='DrugResponseAnalyses', value=len(v.drugresponseanalyses)))
            self.print_responses(v.drugresponseanalyses)
            print('\n')

    @staticmethod
    def print_covariates(covs):
        """
        * Prints all the covariates in covs
            Args:
                covs ([DatasetCovariate]): The covariates from the dataset.
        """
        for c in covs:
            print('\n\t{name}\t\t:\t{value}'.format(name='id', value=c.covariateId))
            print('\t{name}\t\t:\t{value}'.format(name='date', value=c.date))
            print('\t{name}\t\t:\t{value}'.format(name='value', value=c.value))

    @classmethod
    def print_responses(cls, resps):
        """
        * Prints all the drugresponseanalayses in resps
            Args:
                resps ([DrugResponseAnalyses]): The drugresponseanalyses from the dataset.
        """
        for r in resps:
            print('\n\t{name}\t\t:\t{value}'.format(name='drugtreatmentname', value=r.drugtreatmentname))
            print('\t{name}\t\t:\t{value}'.format(name='drugmodelname', value=r.drugmodelname))
            print('\t{name}\t\t:\t{value}'.format(name='postenginename', value=r.postenginename))
            print('\t{name}\t\t:\t{value}'.format(name='reverseenginename', value=r.reverseenginename))
            print('\t{name}\t\t:\t{value}'.format(name='percentileenginename', value=r.percentileenginename))
            print('{name}\t\t:\t{value}'.format(name='Predictions', value=len(r.predictions)))
            cls.print_predictions(r.predictions)

    @staticmethod
    def print_predictions(preds):
        """
        * Prints all the predictions
            Args:
                preds ([Prediction]): The predictions from the dataset.
        """
        for p in preds:
            print('\t{name}\t\t:\t{value}'.format(name='name', value=p.name))
            print('\t{name}\t\t:\t{value}'.format(name='curvetype', value=p.curvetype))
            print('\t{name}\t\t:\t{value}'.format(name='paramstype', value=p.paramstype))
            print('\t{name}\t\t:\t{value}'.format(name='nbpoints', value=p.nbpoints))

    @staticmethod
    def print_samples(samples):
        """
        * Prints all the samples
            Args:
                samples ([Sample]): The samples from the dataset.
        """
        for s in samples:
            print('\n\t{name}\t\t:\t{value}'.format(name='id', value=s.id))
            print('\t{name}\t\t:\t{value}'.format(name='date', value=s.sampledate))
            print('\t{name}\t\t:\t{value}'.format(name='value', value=s.concentration))
            print('\t{name}\t\t:\t{value}'.format(name='unit', value=s.unit))

    @staticmethod
    def print_targets(targets):
        """
        * Prints all the targets
            Args:
                targets ([Target]): The targets from the dataset.
        """
        for t in targets:
            print('\n\t{name}\t\t:\t{value}'.format(name='targettype', value=t.type))
            for c in t.concentrations:
                print(Fore.CYAN + '\t{name}\t\t:'.format(name='max') + Style.RESET_ALL +
                      '\t{value}'.format(value=c.max))
                print('\t{name}\t\t:\t{value}'.format(name='min', value=c.min))
                print('\t{name}\t\t:\t{value}'.format(name='best', value=c.best))
                print('\t{name}\t\t:\t{value}'.format(name='unit', value=c.unit))

    @staticmethod
    def print_dosages(dosages):
        """
        * Prints all the dosages
            Args:
                dosages ([Dosage]): The dosages from the dataset.
        """
        for d in dosages:
            print('\n\t{name}\t\t:\t{value}'.format(name='startdate', value=d.startdate))
            print('\t{name}\t\t:\t{value}'.format(name='lastdate', value=d.lastdate))
            print('\t{name}\t\t:\t{value}'.format(name='dose', value=d.dose))
            print('\t{name}\t\t:\t{value}'.format(name='doseunit', value=d.doseunit))
            print('\t{name}\t\t:\t{value}'.format(name='interval', value=d.interval))
            print('\t{name}\t\t:\t{value}'.format(name='intervalunit', value=d.intervalunit))
            print('\t{name}\t\t:\t{value}'.format(name='infusion', value=d.infusion))
            print('\t{name}\t\t:\t{value}'.format(name='infusionunit', value=d.infusionunit))
            print('\t{name}\t\t:\t{value}'.format(name='intake', value=d.intake))
