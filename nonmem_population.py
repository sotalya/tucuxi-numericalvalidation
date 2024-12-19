import subprocess
from subprocess import check_output, CalledProcessError
import os
import sys
import csv
import shutil
import math
import numpy
import threading
from colorama import init

# Sotalya
from sotalya.data.query import *
from sotalya.data.requests import *


from sotalya.tucuxi.utils import get_platform, run_cmd, raise_exception
from scrape.drugfile import DrugModel
from dateutil.relativedelta import relativedelta


import queue

# For colorama, to reset the color on every new print
init(autoreset=True)


# Suffixes for NONMEM files
suffix_apriori_prediction = '.apriori'
suffix_apriori_percentiles = '.apriori.perc'
suffix_apriori_withVariability = '.apriori.withVariability'
suffix_aposteriori_prediction = '.aposteriori'
suffix_aposteriori_at_sample_time = '.aposteriori.sample'
suffix_aposteriori_percentiles = '.aposteriori.perc'  # not yet supported

if get_platform() == 'Windows':
    PATH_SEPARATOR = '\\'
    win = True
    linux = False
else:
    PATH_SEPARATOR = '/'
    linux = True
    win = False


class RequestNonmem(Request):
    def __init__(self, request:Request, source_file='', drug_model=None):
        super().__init__(request_id=request.requestId, drug_id=request.drugId, drug_model_id=request.drugModelId,
                         computing_traits=request.computingTraits)
        self.sourceFile = source_file
        self.drugModel = drug_model

class NonmemData:
    """
    * This class exposes methods to prepare the NONMEM dataset.
    * A CSV data input file is written from the datasets for NONMEM to execute the calculations. This file, called
      "data" represents the patient.
    * Familiarity with the NM-TRAN syntax and the NONMEM system can be attained by reading the manuals and the
      following site: `http://www.nonlin-model.org/ <http://www.nonlin-model.org/>`_.
    * It is really hard to figure out what is going on in the files without reading the site and the manuals.

    """
    def __init__(self, query: Query, request: Request, drug_model: DrugModel, times, ss,
                 foldername: str, with_variability: bool = False):
        """
        * An invocation of this constructor produces a dataset by NONMEM.

            1. As each query corresponds to a patient, this constructor creates the patient with its measures,
               dosages, and covariates in a data file.

        * All the NONMEM calls are made from the working directory corresponding to a particular sqlite database,
          passed in as "db".

            Args:
                * query(Query): The query from which the dataset will be generated
                * request (Request): The request this run of NONMEM calculates for.
                * drug_model (DrugModel): The drug of the Query.
                * times ([double]): the times provided from the analogous run of Tucuxi
                    in order to get the same point density.
                * ss (bool): if steadystate
                * with_variability : to indicate whether the apriori run should contain variability or not
        """
        print(Back.MAGENTA + '*************************************')
        print(Back.MAGENTA + 'RUNNING NONMEM')
        print(Back.MAGENTA + '*************************************')


        request_nonmem = RequestNonmem(request)
        request_nonmem.drugModel = drug_model

        drug_treatment = query.drugs[0]

        current_request = request_nonmem.computingTraits.computingTraits
        self.requestType = None
        self.parametersType = current_request.computingOption.parametersType
        self.query = query

        self.with_variability = with_variability

        # initiation of request dates according to the request type
        if current_request.requestType == RequestType.Prediction:
            self.requestType = RequestType.Prediction
            self.start = current_request.dateInterval.startDate
            self.endtime = current_request.dateInterval.endDate
            self.conc = int(current_request.nbPointPerHour * ((self.endtime - self.start).total_seconds() / 3600.0))

        elif current_request.requestType == RequestType.PredictionAtTimes:
            self.requestType = RequestType.PredictionAtTimes
            self.start = drug_treatment.dosageHistory.dosageTimeRanges[0].start
            self.endtime = datetime.strptime(current_request.dates.contents[-1].contents[0], '%Y-%m-%dT%H:%M:%S')


        elif current_request.requestType == RequestType.PredictionAtSampleTime:
            self.requestType = RequestType.PredictionAtSampleTime
        elif current_request.requestType == RequestType.Percentiles:
            self.requestType = RequestType.Percentiles
            self.start = current_request.dateInterval.startDate
            self.endtime = current_request.dateInterval.endDate
            self.conc = int(current_request.nbPointPerHour * ((self.endtime - self.start).total_seconds() / 3600.0))

        elif current_request.requestType == RequestType.Adjustment:
            self.requestType = RequestType.Adjustment

        # Dose information
        self.dend = drug_treatment.dosageHistory.dosageTimeRanges[0].end
        self.dstart = drug_treatment.dosageHistory.dosageTimeRanges[0].start
        dosageduration = self.dend - self.dstart
        inter = drug_treatment.dosageHistory.dosageTimeRanges[0].dosage.interval.total_seconds()
        self.cycles = int((dosageduration.total_seconds() / inter) / 3600.0)
        self.lastcycledur = (self.endtime - self.dend).total_seconds() / 3600.0
        self.times = times

        # create a folder for the nonmem results
        self.output_dir = foldername

        if not os.path.isdir(self.output_dir):
            try:
                subprocess.run('mkdir ' + self.output_dir, shell=True)
            except CalledProcessError as e:
                print(Fore.RED + e.output)
                raise_exception(e.output)

        self.output_dir = os.path.join(foldername, 'nonmem')

        if not os.path.isdir(self.output_dir):
            try:
                subprocess.run('mkdir ' + self.output_dir, shell=True)
            except CalledProcessError as e:
                print(Fore.RED + e.output)
                raise_exception(e.output)

        # This variable indicates if for a posteriori predictions, we shall get values for all times (True), or
        # only for times of samples (False). When True, the NONMEM file shall have an MDV field after DV, so the
        # NONMEM files have to be modified accordingly
        self.with_mdv = False
        if self.parametersType == ParametersTypeEnum.aposteriori and self.requestType == RequestType.Prediction:
            self.with_mdv = True
        if self.parametersType == ParametersTypeEnum.aposteriori and self.requestType == RequestType.Percentiles:
            self.with_mdv = True

        self.generate_tree()
        # generate the dataset
        self.generate_patient_file(query, request_nonmem, drug_model, ss, self.with_mdv)



    def generate_tree(self):
        """
        * This method sets the working directory for the calls to NONMEM from the "db" class member.
        * If the directory does not exist, it will create it.
        """
        # if not os.path.isdir(self.output_dir):
        #    clean_and_exit(-1)
        os.chdir(self.output_dir)

    def generate_patient_file(self, query: Query, request: RequestNonmem, drug: DrugModel, ss: bool, with_mdv: bool):
        """
        * This method is used to write the data.csv file with patient information.
        * The NONMEMdatasets and NONMEMrecords correspond to the format of the data files.
        * What happens here is really specific for the NM-TRAN format.

            Args:
                * query (Query): The query containing the current request
                * request (RequestNonmem): The request this run of NONMEM calculates for.
                * drug (DrugModel): The drug of the Query.
                * ss (bool): whether at steadystate
                * with_mdv (bool): whether we use MDV field in NONMEM files
        """

        [covariate_names, covariate_values] = self.load_covariate_names_and_defaults(drug)

        if with_mdv:
            hitems = ['ID', 'TIME', 'EVID', 'AMT', 'CMT', 'SS', 'II', 'DV', 'MDV'] + covariate_names
        else:
            hitems = ['ID', 'TIME', 'EVID', 'AMT', 'CMT', 'SS', 'II', 'DV'] + covariate_names
        nmds = NONMEMdataset(hitems)

        # Dose information
        self.load_dosages_as_records(query, request, drug, nmds, ss, covariate_values, with_mdv)
        intakes = sorted(nmds.records, key=self.get_time)

        # Concentration information
        rollingindex = 0
        if (self.parametersType != ParametersTypeEnum.aposteriori) or with_mdv:
            for intake_index in range(len(intakes) - 1):
                inf_covariate_values = intakes[intake_index].items[-len(covariate_values):]
                if drug.formulationAndRoutes[0].absorptionModel == 'infusion':
                    inf_covariate_values[-1] = "."

                rollingindex = self.load_concs_as_records(rollingindex,
                                                          intakes[intake_index + 1].time, nmds, inf_covariate_values,
                                                          with_mdv)
            intake_index = len(intakes) - 1
            inf_covariate_values = intakes[intake_index].items[-len(covariate_values):]
            if drug.formulationAndRoutes[0].absorptionModel == 'infusion':
                inf_covariate_values[-1] = "."

            self.load_concs_as_records(rollingindex,  self.times[-1] + 1, nmds,
                                                      inf_covariate_values, with_mdv)
            del rollingindex

        # inputting samples
        if self.parametersType == ParametersTypeEnum.aposteriori:
            inf_covariate_values = intakes[0].items[-len(covariate_values):]
            if drug.formulationAndRoutes[0].absorptionModel == 'infusion':
                inf_covariate_values[-1] = "."
            self.load_samples_as_records(query, request, nmds, inf_covariate_values, with_mdv)


        # If the data.csv file already exists, then we fill up the existing file
        # (useful if you want to execute a run on a population rather than a single patient)

        file_exist = os.path.isfile('data.csv')

        if not file_exist:
            f = open('data.csv', 'w')
            f.write(nmds.get_head())
        else:
            f = open('data.csv', 'a')

        sortedrecs = sorted(nmds.records, key=self.sort_record)
        for rec in sortedrecs:
            f.write('\n' + rec.get_rec())
        f.close()


    @staticmethod
    def sort_record(record):
        return record.items[1], record.items[2]

    @staticmethod
    def load_covariate_names_and_defaults(drug):
        """
        * This method is used to write the data.csv file with covariate information.
        * The NONMEMdatasets and NONMEMrecords correspond to the format of the data files.
        * What happens here is really specific for the NM-TRAN format.

            Args:
                * drug (DrugModel): The drug of the Query.
        """
        covariate_names = []
        covariate_defaults = []
        for cov in drug.covariates:
            base, ext = os.path.splitext(cov.covariateId)
            if ext != '':
                covariate_names.append(ext[1:])
            else:
                covariate_names.append(base)
            covariate_defaults.append(float(cov.defaultvalue))
        # check intake valid
        if drug.formulationAndRoutes[0].absorptionModel == 'infusion':
            covariate_names.append("RATE")
            covariate_defaults.append(float(drug.formulationAndRoutes[0].defaultdose) /
                                      float(drug.formulationAndRoutes[0].defaultinfusion)/60.0)
        return [covariate_names, covariate_defaults]

    @staticmethod
    def load_covariates_for_cycle(query, covariate_values, request, drug, dose_date):
        """
        * For a cycle, set the specific set of covariate values in order to adapt to changes

            Args:
                * query (Query): the query for this cycle
                * covariate_values ([double]): Set of covariate values.
                * request (Request): the request for this cycle
                * drug (DrugModel): the drug for this cycle
                * dose_date (date): the date of this cycle
        """

        ddate = dose_date
        if request.computingTraits.computingTraits.computingOption.parametersType != ParametersTypeEnum.population:
            for i, dvar in enumerate(drug.covariates):

                # select the first date of the covariate value
                date_first_cov = None
                for pvar in query.covariates:
                    if pvar.covariateId == dvar.covariateId and date_first_cov is None:
                        date_first_cov = datetime.strptime(pvar.date, "%Y-%m-%dT%H:%M:%S")
                    elif pvar.covariateId == dvar.covariateId and date_first_cov > datetime.strptime(pvar.date, "%Y-%m-%dT%H:%M:%S"):
                        date_first_cov = datetime.strptime(pvar.date, "%Y-%m-%dT%H:%M:%S")


                for j, pvar in enumerate(query.covariates):
                    date_cov = datetime.strptime(pvar.date, "%Y-%m-%dT%H:%M:%S")
                    # If the measured value is the only one, or if it occurs before or at the same time as the dose,
                    # then the value is retained. Otherwise, leave the default value.
                    if (pvar.covariateId == dvar.covariateId) and ((date_cov <= ddate) or date_cov == date_first_cov):
                        covariate_values[i] = float(pvar.value)

                    # Age covariate:
                    if (pvar.covariateId == 'birthdate') and (dvar.covariateType == 'ageInYears'):
                        startdate = request.computingTraits.computingTraits.dateInterval.startDate
                        try:
                            birthdate = datetime.strptime(pvar.value, '%Y-%m-%dT%H:%M:%S')
                        except ValueError:
                            birthdate = datetime.strptime(pvar.value, '%Y-%m-%d')
                        difference_in_years = relativedelta(startdate, birthdate).years
                        covariate_values[i] = difference_in_years

                    if (pvar.covariateId == 'birthdate') and (dvar.covariateType == 'ageInDays'):
                        startdate = request.computingTraits.computingTraits.dateInterval.startDate
                        try:
                            birthdate = datetime.strptime(pvar.value, '%Y-%m-%dT%H:%M:%S')
                        except ValueError:
                            birthdate = datetime.strptime(pvar.value, '%Y-%m-%d')
                        delta = startdate - birthdate
                        difference_in_days = delta.days
                        covariate_values[i] = difference_in_days

                    if (pvar.covariateId == 'birthdate') and (dvar.covariateType == 'ageInWeeks') :
                        try:
                            birthdate = datetime.strptime(pvar.value, '%Y-%m-%dT%H:%M:%S')
                        except ValueError:
                            birthdate = datetime.strptime(pvar.value, '%Y-%m-%d')
                        delta = dose_date - birthdate
                        difference_in_weeks = round(delta.days / 7, 2)
                        covariate_values[i] = difference_in_weeks

        return covariate_values

    def load_samples_as_records(self, query, request, nmds, covariate_values, with_mdv):
        """
        * This method translates a list of samples as records in the NONMEM data file

            Args:
                * query (Query): the query for this cycle
                * request (RequestNonmem): the request for this cycle
                * nmds ([NONMEMrecord]): the rolling list of records that make up a
                    NONMEM data file.
                * covariate_values ([double]): Set of covariate values for the samples
                    (WARNING! cannot be variable here)
                * with_mdv (bool): If yes, add an MDV field after DV
        """
        cmt = '.'

        times_for_samples_prediction = self.times.copy()

        for s in query.drugs[0].samples:
            sampletime = s.sampledate
            elapsedtime = sampletime - self.dstart
            if s.unit == 'mg/l':
                dv1 = float(s.concentration) * 1000
            else:
                dv1 = float(s.concentration)

            if request.drugModel.errorModel == 'exponential':
                dv1 = math.log(dv1)

            patient_id = self.query.queryId.split("_")[0]
            if with_mdv:
                items = [patient_id, elapsedtime.total_seconds()/3600.0, 0, '.', cmt, 0, '.', dv1, '0'] + covariate_values
            else:
                items = [patient_id, elapsedtime.total_seconds()/3600.0, 0, '.', cmt, 0, '.', dv1] + covariate_values
            nmds.records.append(NONMEMrecord(items))

            times_for_samples_prediction.remove(elapsedtime.total_seconds()/3600.0)

        for i in times_for_samples_prediction:
            patient_id = self.query.queryId.split("_")[0]
            if with_mdv:
                items = [patient_id, i, 2, '.', cmt, 0, '.', ".",
                         '0'] + covariate_values
            else:
                items = [patient_id, i, 2, '.', cmt, 0, '.', "."] + covariate_values
            nmds.records.append(NONMEMrecord(items))


    def load_dosages_as_records(self, query, request, drug, nmds, ss, covariate_values, with_mdv):
        """
        * This method translates a list of dosages as records in the NONMEM data file

            Args:
                * query (Query): Query containing the request
                * request (RequestNonmem): the request for this cycle
                * with_mdv: If yes, add an MDV field after DV
                * covariate_values ([double]): Set of covariate values for the samples
                    (WARNING! cannot be variable here)
                * ss: steady-state
                * nmds ([NONMEMrecord]): the rolling list of records that make up a
                    NONMEM data file.
                * drug: the drug for this cycle
        """
        for d in query.drugs[0].dosageHistory.dosageTimeRanges:
            # pdb.set_trace()
            dend = d.end
            dstart = d.start
            try:
                request_end = request.computingTraits.computingTraits.dateInterval.endDate
                if request_end < dend:
                    dend = request_end
            except AttributeError:
                dend = dend

            duration = (dend - dstart).total_seconds()/3600.0
            sincestart = (dstart - self.start).total_seconds()/3600.0
            cycle = 0
            # update covariate values
            interval = d.dosage.interval.total_seconds() / 3600.0
            while cycle * float(interval) < duration:
                cmt = '.'
                inkstart = sincestart + cycle * float(interval)
                if drug.formulationAndRoutes[0].absorptionModel == 'infusion':
                    covariate_values[-1] = float(d.dosage.dose.value) / \
                                           (float(d.dosage.dose.infusionTimeInMinutes.total_seconds())/3600.0)


                dose_date = dstart + timedelta(hours=cycle * float(interval))
                self.load_covariates_for_cycle(query, covariate_values, request, drug, dose_date)

                patient_id = self.query.queryId.split("_")[0]
                if ss:
                    if with_mdv:
                        items = [patient_id, inkstart, 4, d.dosage.dose.value, '.', 1, d.dosage.interval, '.', '1'] \
                                + covariate_values
                    else:
                        items = [patient_id, inkstart, 4, d.dosage.dose.value, '.', 1, d.dosage.interval, '.'] \
                                + covariate_values
                else:
                    if with_mdv:
                        items = [patient_id, inkstart, 1, d.dosage.dose.value, cmt, 0, '.', '.', '1'] + covariate_values
                    else:
                        items = [patient_id, inkstart, 1, d.dosage.dose.value, cmt, 0, '.', '.'] + covariate_values
                nmds.records.append(NONMEMrecord(items))
                cycle = cycle + 1

    def load_concs_as_records(self, rollingindex, end, nmds, covariate_values, with_mdv):
        """
        * This method translates a list of times as records in the NONMEM data file.
            When NONMEM runs with this job, it fills out concentrations at these times.

            Args:
                * nmds ([NONMEMrecord]): the rolling list of records that make up a
                    NONMEM data file.
                * with_mdv: If yes, add an MDV field after DV
                * covariate_values ([double]): Set of covariate values for the samples
                    (WARNING! cannot be variable here)
                * rollingindex: the number of records
                * end: the end of the cycle
        """


        cmt = '.'

        patient_id = self.query.queryId.split("_")[0]
        while rollingindex < len(self.times) and self.times[rollingindex] < end:
            if with_mdv:
                items = [patient_id, self.times[rollingindex], 0, '.', cmt, 0, '.', '.', '1'] + covariate_values
            else:
                items = [patient_id, self.times[rollingindex], 0, '.', cmt, 0, '.', '.'] + covariate_values
            rollingindex = rollingindex + 1
            nmds.records.append(NONMEMrecord(items))
        return rollingindex


    @staticmethod
    def get_time(item):
        return item.time




class NONMEMdataset:
    def __init__(self, hitems):
        self.records = []
        self.hitems = hitems

    def get_head(self):
        outstring = str(self.hitems[0])
        for hitem in self.hitems[1:]:
            outstring = outstring + "," + str(hitem)
        return outstring


class NONMEMrecord:
    def __init__(self, items):
        self.time = float(items[1])
        self.items = items

    def get_rec(self):
        outstring = str(self.items[0])
        for item in self.items[1:]:
            outstring = outstring + "," + str(item)
        return outstring



class NonmemExecute:
    """
        * This class exposes methods to prepare and run the NONMEM executable with arguments.
        * A CSV data input file is written from the datasets for NONMEM to execute the calculations. This file, called
          "data" represents the patient.
        * For each drug, a coresponding model file containing NM-TRAN code must exist with the same name in the folder
          "models". The model files for apriori, aposteriori, and percentiles are slightly different.
        * This class requires feedback from the NONMEM output files which are parsed from the working directory after the
          calculations.
        * The sequence of commands and the flow of data is managed in the constructor, which accepts arguments from
        * GlobalTester to configure what NONMEM will do.
        * Each command and the contents of the data file is echoed to the console such that to manually reproduce what this
          class automates, one could copy/paste the commands from the console output and recreate the data file
          (useful in cases of error).
        * Familiarity with the NM-TRAN syntax and the NONMEM system can be attained by reading the manuals and the
          following site: `http://www.nonlin-model.org/ <http://www.nonlin-model.org/>`_.
        * It is really hard to figure out what is going on in the files without reading the site and the manuals.
        * The PRED routines are used as they correspond to the ADME models Tucuxi uses for both steadystate and non.

        """
    def __init__(self,
                 drug_model: DrugModel, list_queries: [Query],
                 use_cache: bool, clean: bool,nonmem_models_path: str, foldername: str,
                 cmd: str, steady_state: bool, with_variability: bool =False,
                 nb_patient_per_thread: int = 1000, nb_threads: int = 8):
        """
                * An invocation of this constructor produces a complete test of one dataset by NONMEM.

                    1. this constructor creates the dataset, containing each patient with its measures,
                       dosages, and covariates in a data file.
                    2. Then it runs a calculation with the NM-TRAN model file corresponding to the drug of the patient.
                       It chooses the model file based on the type of calculation will be done: apriori, aposteriori,
                       percentiles.
                    3. The same return data types are used for this class and Tucuxi class. The data is parsed from the NONMEM
                       output files in the working directy after the run.

                * All the NONMEM calls are made from the working directory corresponding to a particular sqlite database,
                  passed in as "db".

                    Args:
                        * drug_model (DrugModel): The drug of the Query.
                        * list_queries (Query): a list containing one or several queries (one per patient)
                        * nonmem_models_path (str): path of the nonmem models
                        * foldername: name of the folder to save the results
                        * clean (bool): if to not clean up (remove all the working files after)
                        * steady-state (bool): if steadystate
                        * with_variability (bool): if we want to predict apriori concentrations with variability or not
                        * cmd (str): command to execute NONMEM
                        * use_cache (bool): if true, try to get data from existing cache files
                        * nb_patient_per_thread (int): number of patients to be calculated by a thread for percentiles
                        * nb_threads (int): number of threads to launch for percentiles calculation
                """

        self.drug_model = drug_model

        self.list_queries = list_queries
        self.list_requests = []

        self.nb_patient_per_thread = nb_patient_per_thread
        self.nb_threads = nb_threads
        self.nonmem_models_path = nonmem_models_path
        self.cmd = cmd

        self.status = True

        self.ss= steady_state
        self.with_variability = with_variability

        list_with_mdv = []
        list_request_type = []
        list_parameters_type = []

        self.nb_samples = 0

        # Dataset creation
        for i in range(len(list_queries)):
            for ri, response in enumerate(list_queries[i].requests):
                self.list_requests.append(response)

                times_samples = []
                for d in response.computingTraits.computingTraits.dates.contents:
                    date = datetime.strptime(d.contents[0], '%Y-%m-%dT%H:%M:%S')
                    times_samples.append((date - list_queries[i].drugs[0].dosageHistory.dosageTimeRanges[
                        0].start).total_seconds() / 3600)

                nm = NonmemData(list_queries[i], response,
                                self.drug_model,
                                times_samples, self.ss, with_variability=self.with_variability,
                                foldername= foldername)

                self.nb_samples = self.nb_samples + len(times_samples)
                list_with_mdv.append(nm.with_mdv)
                list_request_type.append(nm.requestType)
                list_parameters_type.append(nm.parametersType)

        # If several queries are used to run the same model, then we check that the query
        # type is the same for each patient.
        try:
            self.check_similar_list(list_with_mdv)
        except ValueError as e:
            print(str(e))
        self.with_mdv = list_with_mdv[0]

        try:
            self.check_similar_list(list_request_type)
        except ValueError as e:
            print(str(e))
        self.request_type_population = list_request_type[0]

        try:
            self.check_similar_list(list_parameters_type)
        except ValueError as e:
            print(str(e))
        self.parameters_type_population = list_parameters_type[0]


        if not os.path.isfile('{nmodelspath}{separator}{drug}{suffix}'.format(nmodelspath=self.nonmem_models_path,
                                                                              separator=PATH_SEPARATOR,
                                                                              drug=self.drug_model.drugModelId,
                                                                              suffix=suffix_apriori_prediction)):
            raise_exception('No model file with drug name \"{drug}{suffix}\"'
                            ' in model folder'.format(drug=self.drug_model.drugModelId, suffix=suffix_apriori_prediction))

        # requests are defined: if there is only one patient, then the request is used directly.
        # In the case of a population, we define a new request with the population values

        if len(list_queries)==1:
            request = list_queries[0].request[0]
            query_id = list_queries[0].get_id()
        else:
            request_population = Request()
            request_population.requestId = self.drug_model.drugModelId + self.request_type_population.value
            request_population.drugId = self.drug_model.drugId
            request_population.drugModelId = self.drug_model.drugModelId

            request = request_population
            query_id = "pop"

        is_filed = False

        # Nonmem execute
        if (self.request_type_population == RequestType.Prediction) or (self.request_type_population == RequestType.PredictionAtTimes):
            if ((self.parameters_type_population == ParametersTypeEnum.population) or
                    (self.parameters_type_population == ParametersTypeEnum.apriori)):
                if use_cache:
                    is_filed = self.read_apriori_from_file(drug_model, request, query_id)
                if not is_filed:
                    if self.with_variability:
                        self.apriori_with_variability(request, query_id)
                    else:
                        self.apriori(request, query_id)
                    self.results = self.extract_concentration()
                    is_filed = self.write_apriori_to_file(drug_model, query_id, self.results, request)
            elif self.parameters_type_population == ParametersTypeEnum.aposteriori:
                if use_cache:
                    is_filed = self.read_aposteriori_from_file(drug_model, request, query_id)
                if not is_filed:
                    self.aposteriori_prediction(request, query_id)
                    self.results = self.extract_aposteriori_concentration(request)
                    is_filed = self.write_aposteriori_to_file(drug_model, query_id, self.results, request)
            else:
                print("Very strange, that shouldn't happen")


        elif self.request_type_population == RequestType.PredictionAtSampleTime:
            if self.parameters_type_population == ParametersTypeEnum.aposteriori:
                if use_cache:
                    is_filed = self.read_aposteriori_from_file(drug_model, request, query_id)
                if not is_filed:
                    self.aposteriori_only_samples(request, query_id)
                    self.results = self.extract_aposteriori_only_samples(request)
                    is_filed = self.write_aposteriori_to_file(drug_model, query_id, self.results, request)
            else:
                print("Very strange, that shouldn't happen")

        elif self.request_type_population == RequestType.Percentiles:
            if use_cache:
                is_filed = self.read_perc_from_file(drug_model, request, query_id)
            if not is_filed:
                self.results = self.mt_percentiles(request, self.with_mdv)
                if self.results:
                    is_filed = self.write_perc_to_file(drug_model, query_id, self.results, request)
                else:
                    self.status = False
                    print(Fore.RED + "Error with percentiles calculation")

        if not is_filed:
            print(Fore.RED + "IO Error writing to csv file.")


        if clean:
            # This check allows to let the folder untouched if something bad happened
            if hasattr(self, 'results'):
                try:
                    subprocess.run('rm -r * ', shell=True)
                except CalledProcessError as e:
                    print(Fore.RED + e.output)
                    raise_exception(e.output)

        os.chdir(os.path.join("..", "..", ".."))

    @staticmethod
    def check_similar_list(list_to_check : []):
            """
              * This function is used to check whether the list contains the same elements.

                  Args:
                    * list_to_check (list): list to be checked
            """
            if len(list_to_check) == 0:
                # if the list is empty
                raise ValueError("This list is empty")

            # compare all the element in the list with the first element
            first_value = list_to_check[0]
            for element in list_to_check:
                if element != first_value:
                    raise ValueError("This list is not identical")


    def system_exit(self):
        pass

    @staticmethod
    def is_obsolete(filename: str, request: RequestNonmem):
        """

        :param filename:
        :param request:
        :return:
        """
        try:
            filetime = os.path.getmtime(filename)
            queryfiletime = os.path.getmtime(request.sourceFile)
            if filetime < queryfiletime:
                return True

        except os.error:
            return True

        return False

    def read_apriori_from_file(self, drug, request, query_id: str):
        """
        * This method reads apriori cached values from a file with a predefined
          name made of the drugid and patientid. This is done to avoid calculating
          with nonmem because it takes forever, or because it is not installed.

            Args:
                * request (RequestNonmem): The request this run of NONMEM calculates for.
                * drug (DrugModel): The drug of the Query.
                * query_id(str): the name of the query
            Returns:
                * anonymous (bool): whether the read was successful or not.
        """
        times = []
        dv = []

        filename = os.path.join('..', '..', '..', 'NONMEM_stored_apriori', '{d}_{id}.csv'.format(d=drug.drugModelId,
                                                                                                 id=query_id +
                                                                                                 request.get_id()))
        if self.is_obsolete(filename, request):
            return False

        try:
            f = open(filename, 'r')
            reader = csv.reader(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            for row in reader:
                # pdb.set_trace()
                times.append(float(row[0]))
                dv.append(float(row[1]))
            f.close()
            self.results = [times, dv]

        except IOError:
            return False
        return True

    @staticmethod
    def write_apriori_to_file(drug: DrugModel, query_id: str, results, request):
        """
        * This method writes apriori cached values from a file with a predefined
          name made of the drugid and patientid. This is done to avoid calculating
          with nonmem because it takes forever, or because it is not installed.

            Args:
                * request (Request): The request this run of NONMEM calculates for.
                * drug (DrugModel): The drug of the Query.
                * results ([[][]]): The data to cache
                * query_id (str): name of the query
        """
        try:
            f = open(os.path.join('..', '..', '..', 'NONMEM_stored_apriori', '{d}_{id}.csv'.
                                  format(d=drug.drugModelId, id=query_id + request.get_id())), 'w')
            writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            for i, v in enumerate(results[0]):
                # pdb.set_trace()
                writer.writerow([results[0][i]] + [results[1][i]])
            f.close()

        except IOError:
            return False
        return True

    def read_aposteriori_from_file(self, drug, request, query_id):
        """
        * This method reads apost cached values from a file with a predefined
          name made of the drugid and patientid. This is done to avoid calculating
          with nonmem because it takes forever, or because it is not installed.

            Args:
                * request (RequestNonmem): The request this run of NONMEM calculates for.
                * drug (DrugModel): The drug of the Query.
                * query_id (str): name of the query
        """
        times = []
        dv = []

        filename = os.path.join('..', '..', '..', 'NONMEM_stored_apost', '{d}_{id}.csv'.format(d=drug.drugModelId,
                                                                                               id= query_id
                                                                                               + request.get_id()))
        if self.is_obsolete(filename, request):
            return False

        try:
            f = open(filename, 'r')
            reader = csv.reader(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            for row in reader:
                # pdb.set_trace()
                times.append(float(row[0]))
                dv.append(float(row[1]))
            f.close()
            self.results = [times, dv]

        except IOError:
            return False
        return True

    @staticmethod
    def write_aposteriori_to_file(drug, query_id: str, results, request):
        """
        * This method writes apost cached values from a file with a predefined
          name made of the drugid and patientid. This is done to avoid calculating
          with nonmem because it takes forever, or because it is not installed.

            Args:
                * query_id (str): The ID of the query containing the request
                * request (Request): The request this run of NONMEM calculates for.
                * drug (DrugModel): The drug of the Query.
                * results ([[][]]): The data to cache
        """

        try:
            f = open(os.path.join('..', '..', '..', 'NONMEM_stored_apost', '{d}_{id}.csv'.
                                  format(d=drug.drugModelId, id=query_id + request.get_id())), 'w')
            writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            for i, v in enumerate(results[0]):
                # pdb.set_trace()
                writer.writerow([results[0][i]] + [results[1][i]])
            f.close()

        except IOError:
            return False
        return True

    def read_perc_from_file(self, drug, request, query_id):
        """
        * This method reads apost cached values from a file with a predefined
          name made of the drugid and patientid. This is done to avoid calculating
          with nonmem because it takes forever, or because it is not installed.

            Args:
                * request (RequestNonmem): The request this run of NONMEM calculates for.
                * drug (DrugModel): The drug of the Query.
                * query_id (str): name of the query
        """

        cperc = []
        percentiles = request.computingTraits.computingTraits.computingOption.ranks
        ctimes = []

        filename = os.path.join('..', '..', '..',
                                'NONMEM_stored_percentiles', '{d}_{id}.csv'.format(d=drug.drugModelId,
                                                                                   id= query_id +
                                                                                   request.get_id()))

        if self.is_obsolete(filename, request):
            return False

        try:

            f = open(filename, 'r')
            reader = csv.reader(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            for row in reader:
                # pdb.set_trace()
                ctimes.append(float(row[0]))
                cp = []
                for i in range(1, len(percentiles) + 1):
                    cp.append(float(row[i]))
                cperc.append(cp)
            f.close()
            if len(ctimes) == 0:
                return False
            self.results = [ctimes, cperc]

        except IOError:
            return False
        return True

    @staticmethod
    def write_perc_to_file(drug, query_id, results, request):
        """
        * This method writes percentiles cached values from a file with a predefined
          name made of the drugid and patientid. This is done to avoid calculating
          with nonmem because it takes forever, or because it is not installed.

            Args:
                * query_id (str): name of the query
                * request (Request): The request this run of NONMEM calculates for.
                * drug (DrugModel): The drug of the Query.
                * results ([[][]]): The data to cache
        """

        print(os.getcwd())
        try:
            f = open(os.path.join('..', '..', '..', 'NONMEM_stored_percentiles', '{d}_{id}.csv'
                     .format(d=drug.drugModelId, id=query_id + request.get_id())), 'w')
            writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            for i, v in enumerate(results[0]):
                # pdb.set_trace()
                writer.writerow([results[0][i]] + results[1][i])
            f.close()

        except IOError:
            return False
        return True

    @staticmethod
    def copy_files(drug_model_id, suffix, source_folder,
                   dest_folder='.'):

        folder = os.path.abspath(source_folder)

        shutil.copy2('{folder}{separator}{drug}{suffix}'.format(folder=folder, separator=PATH_SEPARATOR,
                                                                drug=drug_model_id, suffix=suffix),
                     '{dest}{separator}{drug}{suffix}'.format(dest=dest_folder, separator=PATH_SEPARATOR,
                                                              drug=drug_model_id, suffix=suffix))

        shutil.copy2('{folder}{separator}{drug}.pkmodel'.format(folder=folder, separator=PATH_SEPARATOR,
                                                                drug=drug_model_id, suffix=suffix),
                     '{dest}{separator}{drug}.pkmodel'.format(dest=dest_folder, separator=PATH_SEPARATOR,
                                                              drug=drug_model_id, suffix=suffix))

        shutil.copy2('{folder}{separator}general.header'.format(folder=folder, separator=PATH_SEPARATOR,
                                                                drug=drug_model_id, suffix=suffix),
                     '{dest}{separator}general.header'.format(dest=dest_folder, separator=PATH_SEPARATOR,
                                                              drug=drug_model_id, suffix=suffix))

        shutil.copy2('{folder}{separator}general.apriori.footer'.format(folder=folder, separator=PATH_SEPARATOR,
                                                                        drug=drug_model_id, suffix=suffix),
                     '{dest}{separator}general.apriori.footer'.format(dest=dest_folder, separator=PATH_SEPARATOR,
                                                                      drug=drug_model_id, suffix=suffix))

        shutil.copy2(
            '{folder}{separator}general.apriori.withVariability.footer'.format(folder=folder, separator=PATH_SEPARATOR,
                                                                               drug=drug_model_id, suffix=suffix),
            '{dest}{separator}general.apriori.withVariability.footer'.format(dest=dest_folder, separator=PATH_SEPARATOR,
                                                                             drug=drug_model_id, suffix=suffix))

        shutil.copy2('{folder}{separator}general.aposteriori.footer'.format(folder=folder, separator=PATH_SEPARATOR,
                                                                            drug=drug_model_id, suffix=suffix),
                     '{dest}{separator}general.aposteriori.footer'.format(dest=dest_folder, separator=PATH_SEPARATOR,
                                                                          drug=drug_model_id, suffix=suffix))

        shutil.copy2('{folder}{separator}general.aposteriori.sample.footer'.format(folder=folder,
                                                                                   separator=PATH_SEPARATOR,
                                                                                   drug=drug_model_id,
                                                                                   suffix=suffix),
                     '{dest}{separator}general.aposteriori.sample.footer'.format(dest=dest_folder,
                                                                                 separator=PATH_SEPARATOR,
                                                                                 drug=drug_model_id,
                                                                                 suffix=suffix))

        shutil.copy2('{folder}{separator}general.apriori.perc.footer'.format(folder=folder, separator=PATH_SEPARATOR,
                                                                             drug=drug_model_id, suffix=suffix),
                     '{dest}{separator}general.apriori.perc.footer'.format(dest=dest_folder, separator=PATH_SEPARATOR,
                                                                           drug=drug_model_id, suffix=suffix))

        shutil.copy2('{folder}{separator}general.dataimport'.format(folder=folder, separator=PATH_SEPARATOR,
                                                                    drug=drug_model_id, suffix=suffix),
                     '{dest}{separator}general.dataimport'.format(dest=dest_folder, separator=PATH_SEPARATOR,
                                                                  drug=drug_model_id, suffix=suffix))

    def aposteriori_prediction(self, request, query_id):
        """
        * This method runs NONMEM from the specified path (default is /opt/nm72/run/nmfe72 but can be overridden by
        * the option).
        * The model file corresponding to the drug for aposteriori calculations is selected. This model does posthoc
        * estimates.
        * I tested the default calculation engine in NONMEM (FOCE), and Laplacian, but in fact the implementation in the
        * foce plugin is laplacian. There is little difference (Aziz is using FOCE for all).

            Args:
                * request (Request): the request being run
                * query_id: name of the query
        """

        self.copy_files(request.drugModelId, suffix_aposteriori_prediction, self.nonmem_models_path)

        print(
            Fore.CYAN + '{nm} {drug}{suffix} {jobname}'.format(nm=self.cmd, drug=request.drugModelId,
                                                               suffix=suffix_aposteriori_prediction,
                                                               jobname= query_id + request.get_id() +
                                                                       '.aposteriori.log'))

        try:
            file_in = request.drugModelId + suffix_aposteriori_prediction
            file_out = query_id + request.get_id() + '.aposteriori.log'
            # check if platform is windows. if so, execute the subprocess as a pipeline.
            if sys.platform in ["win32", "cygwin"]:
                subprocess.run([self.cmd, file_in, file_out], capture_output=True)
            else:
                check_output([self.cmd, file_in, file_out])

            # check_output(['{nm} {drug}{suffix} {jobname}'.format(nm=self.cmd,
            #                                                      drug=request.drugModelId,
            #                                                      suffix=suffix_aposteriori_prediction,
            #                                                      jobname=request.get_full_id()+'.aposteriori.log'),
            #               ''], shell=True)
        except CalledProcessError as e:
            print(Fore.RED, e.output)
            self.system_exit()

    def aposteriori_only_samples(self, request, query_id):
        """
        * This method runs NONMEM from the specified path (default is /opt/nm72/run/nmfe72 but can be overridden by
        * the option).
        * The model file corresponding to the drug for aposteriori calculations is selected. This model does posthoc
        * estimates.
        * I tested the default calculation engine in NONMEM (FOCE), and Laplacian, but in fact the implementation in the
        * foce plugin is laplacian. There is little difference (Aziz is using FOCE for all).

            Args:
                * request (Request): the request being run
                * query_id (str): name of the query
        """

        self.copy_files(request.drugModelId, suffix_aposteriori_at_sample_time, self.nonmem_models_path)

        print(
            Fore.CYAN + '{nm} {drug}{suffix} {jobname}'.format(nm=self.cmd,
                                                               drug=request.drugModelId,
                                                               suffix=suffix_aposteriori_at_sample_time,
                                                               jobname=query_id + request.get_id() +
                                                                       '.aposteriori.sample.log'))

        try:
            file_in = request.drugModelId + suffix_aposteriori_at_sample_time
            file_out = query_id + request.get_id() + '.aposteriori.sample.log'
            # check if platform is windows. if so, execute the subprocess as a pipeline.
            if sys.platform in ["win32", "cygwin"]:
                subprocess.run([self.cmd, file_in, file_out], capture_output=True)
            else:
                check_output([self.cmd, file_in, file_out])

        except CalledProcessError as e:
            print(Fore.RED + e.output)
            raise_exception(e.output)

    def apriori(self, request, query_id):
        """
        * This method runs NONMEM with the model file corresponding to the apriori calculations. A simulation is done
        * with one single individual with zero variance.

            Args:
                * request (Request): the request being run
                * query_id : name of the query
        """

        self.copy_files(request.drugModelId, suffix_apriori_prediction, self.nonmem_models_path)

        print(
            Fore.CYAN + '{nm} {drug}{suffix} {jobname}'
            .format(nm=self.cmd, drug=self.drug_model.drugModelId, suffix=suffix_apriori_prediction,
                    jobname=query_id+ request.get_id() + '.apriori.log'))

        try:
            file_in = request.drugModelId + suffix_apriori_prediction
            file_out = query_id + request.get_id() + '.apriori.log'

            # check if platform is windows. if so, execute the subprocess as a pipeline.
            if sys.platform in ["win32", "cygwin"]:
                subprocess.run([self.cmd, file_in, file_out], capture_output=True)
            else:
                check_output([self.cmd, file_in, file_out])

        except CalledProcessError as e:
            print(Fore.RED, e.output)
            raise_exception(e.output)

    def apriori_with_variability(self, request, query_id):
        """
        * This method runs NONMEM with the model file corresponding to the apriori calculations. A simulation is done
        * with one single individual with variance.

            Args:
                * request (Request): the request being run
                * query_id (str): name of the query
        """

        self.copy_files(request.drugModelId, suffix_apriori_withVariability, self.nonmem_models_path)

        print(
            Fore.CYAN + '{nm} {drug}{suffix} {jobname}'
            .format(nm=self.cmd, drug=request.drugModelId, suffix=suffix_apriori_withVariability,
                    jobname=query_id + request.get_id() + '.apriori.withVariability.log'))

        try:
            file_in = request.drugModelId + suffix_apriori_withVariability
            file_out = query_id + request.get_id() + '.apriori.withVariability.log'

            # check if platform is windows. if so, execute the subprocess as a pipeline.
            if sys.platform in ["win32", "cygwin"]:
                subprocess.run([self.cmd, file_in, file_out], capture_output=True)
            else:
                check_output([self.cmd, file_in, file_out])

        except CalledProcessError as e:
            print(Fore.RED, e.output)
            raise_exception(e.output)

    def mt_percentiles(self, request, with_mdv):
        """
        * Calculating Percentiles (as far as I know) is not something NONMEM normally does for the user.
            The way it is done for this project takes a long time, so instead of running 10,000 patients
            all one after the other, we chop up the work into chunks of 500 patients, and then bring them all
            together afterward to calculate percentiles using numpy.

            Args:
                * request (Request): the request being run
                * with_mdv (bool): Whether to use or not the MDV field in NONMEM files
        """
        qu = queue.Queue()
        tresults = {}

        nb_jobs = self.nb_threads
        job_ids = []
        for i in range(nb_jobs):
            job_ids.append(i)
            qu.put(i)

        try:
            threads = [threading.Thread(target=self.nm_perc_worker, args=(request, qu, tresults, with_mdv))
                       for _i in job_ids]
            for thread in threads:
                try:
                    thread.daemon = True
                    thread.start()
                    # join here for test purpose
                    # thread.join()
                except threading.ThreadError:
                    print(Fore.RED + "Threading error")
                    # thread.exit()
                    raise_exception("Threading error")
                # pdb.set_trace()
                # qu.put(None)
            qu.join()
        except RuntimeError as err:
            raise_exception(err.args[0])

        # pdb.set_trace()
        print(Fore.CYAN + 'threads returned\n')
        if tresults == {}:
            return False

        # Get the number of values from the data file. -1 because of the header row
        # nb_values = len(open("data.csv").readlines()) - 1
        nb_values = len(self.nb_samples)

        ret = self.taggregate(tresults, request, nb_values)
        return ret

    def nm_perc_worker(self, request, qu, results, with_mdv, query_id):
        """
        * For the mt NONMEM monte carlo, this is what a worker thread does
                :param with_mdv: whether to use or not the MDV field in NONMEM files
                :param request: name of the request
                :param results: results of the request
                :param qu
                :param query_id: name of the query

            Args:
                * dataset (Query): the dataset being run
        """
        if self.parameters_type_population == ParametersTypeEnum.aposteriori:
            print(Fore.RED + "A posteriori percentiles not implemented in NONMEM")

        print(Fore.YELLOW + 'starting worker\n', )

        while True:
            # pdb.set_trace()
            jobno = qu.get()
            print(Fore.YELLOW + 'worker initialized with jobno: {jn}\n'.format(jn=jobno), )
            # pdb.set_trace())
            job = 'job' + str(jobno)
            if jobno is None:
                # pdb.set_trace()
                return
            try:
                if not os.path.isdir('job{jn}'.format(jn=jobno)):
                    # check_output(['mkdir job{jn}'.format(jn=jobno)], shell=True)

                    # Linux compatibility ?
                    os.makedirs(job)
            except IOError as e:
                print(e.error + '\n', )
                raise_exception(e.error)
            try:
                cmd1_info = "data.csv job{jn}{separator}data.csv".format(jn=jobno, separator=PATH_SEPARATOR)

                if with_mdv:
                    cmd2_info = "{nmpath}{separator}{drug}.sim.mdv job{jn}{separator}{drug}.sim.mdv" \
                        .format(nmpath=self.nonmem_models_path, jn=jobno,
                                separator=PATH_SEPARATOR, drug=request.drugModelId)
                    cmd3_info = "job{jn}".format(jn=jobno)
                    cmd4_info = "{nm} {drug}.sim.mdv {drug}_out.txt".format(nm=self.cmd, drug=request.drugModelId)

                    if win:
                        command_list = ["copy " + cmd1_info, "copy " + cmd2_info, "cd " + cmd3_info, cmd4_info]
                    else:
                        command_list = ["cp " + cmd1_info, "cp " + cmd2_info, "cd " + cmd3_info, cmd4_info]

                    run_cmd(command_list)
                else:
                    cmd2_info = "job{jn}".format(jn=jobno)
                    cmd3_info = " {nm} {drug}{suffix} {jobname}".format(nm=self.cmd, drug=request.drugModelId,
                                                                        suffix=suffix_apriori_percentiles,
                                                                        jobname= query_id +
                                                                                request.get_id() +
                                                                                '.apriori.perc.log')
                    if win:
                        command_list = ["copy " + cmd1_info, "cd " + cmd2_info, cmd3_info]
                    else:
                        command_list = ["cp " + cmd1_info, "cd " + cmd2_info, cmd3_info]

                    self.copy_files(request.drugModelId, suffix_apriori_percentiles, self.nonmem_models_path,
                                    'job{jn}'.format(jn=jobno))

                    # We modify the NONMEM footer to have a different seed for each job
                    data_file = "job{jn}{separator}general.apriori.perc.footer".format(jn=jobno,
                                                                                       separator=PATH_SEPARATOR)
                    fin = open(data_file, "rt")
                    data = fin.read()
                    # The SEED seems to be used for other types of computations
                    # #new_seed = "SEED=" + str(jobno + 10)
                    # data = data.replace("SEED=2", new_seed)
                    new_seed_sim = "$SIM(" + str(123456 + int(jobno + 10)) + ")"
                    data = data.replace("$SIM (123456)", new_seed_sim)

                    new_sub = "SUB=" + str(self.nb_patient_per_thread)
                    data = data.replace("SUB=1000", new_sub)

                    fin.close()
                    fin = open(data_file, "wt")
                    fin.write(data)
                    fin.close()

                    run_cmd(command_list)

            except CalledProcessError as e:
                print(e.output + '\n', )
                # raise_exception(e.output)
                return

            print(Fore.YELLOW + 'worker on job {jn} ran nonmem\n'.format(jn=jobno), )
            run_cmd(['cd ..'])
            # pdb.set_trace()
            try:
                self.extractp(jobno, results)
            except RuntimeError as err:
                print(Fore.RED + "Percentile thread error: " + err.args[0])
            qu.task_done()
            print(Fore.YELLOW + 'worker on job {jn} extracted values\n'.format(jn=jobno), )

    def extractp(self, jobno, results):
        """
        * Extracts the results for the job from file

                :param results: results of the request
                :param jobno: asdf

            Args:
                * file (str): the file containing the results for the job

            Returns:
                * vals ([]): the values from the files
        """
        try:
            f = open('job{jn}'.format(jn=jobno) + PATH_SEPARATOR + 'apriori.perc.out.dat', 'r')
        except IOError as e:
            print(Fore.RED + e.strerror + ' ' + e.filename + '\n', )
            raise_exception(e.strerror + ' ' + e.filename)
            return
        # pdb.set_trace()
        res = self.create_base_sim_generator(f)
        # res = self.get_perc_data(f)
        results['job{jn}'.format(jn=jobno)] = res

    @staticmethod
    def get_perc_data(file):
        """
        * Extracts the results for the job from file

            Args:
                * file (File): the file containing the results for the job

            Returns:
                * vals ([]): the values from the files
        """
        allpercs = []
        for line in file:
            # pdb.set_trace()
            vals = line.split()
            if vals[0] == "License":
                break
            # vals.pop(0)
            allpercs.append(vals)
        return allpercs

    @staticmethod
    def create_base_sim_generator(file):
        """
        * Extracts the results for the job from file

            Args:
                * file (File): the file containing the results for the job

            Returns:
                * vals ([]): the values from the files
        """
        for line in file:
            # pdb.set_trace()
            vals = line.split()
            if vals[0] == "License":
                break
            # vals.pop(0)
            yield vals

    @staticmethod
    def taggregate(results, request, nb_values):
        """
        * To aggregate, we add all results to the first job results

            Args:
                * results ({job, []}): the dataset being run
                * request (Request): The request of interest, to retrieve the percentiles ranks
                * nb_values (int): number of values request
        """
        cperc = []
        percentiles = request.computingTraits.computingTraits.ranks
        if len(results) == 0:
            return [[], []]

        numtimes = nb_values
        timesdict = []  # {k: [] for k in range(numtimes)}
        valsdict = {k: [] for k in range(numtimes)}

        time_index = 0
        val_index = 0
        for i, v in enumerate(results):
            for j, k in enumerate(results[v]):
                # pdb.set_trace()
                if (time_index < numtimes) and (float(k[2]) == 0):
                    timesdict.append(float(k[1]))
                    time_index = time_index + 1
                if float(k[2]) == 0:
                    valsdict[val_index % numtimes].append(float(k[4]))
                    val_index = val_index + 1

        for index in range(numtimes):
            valset = numpy.array(valsdict[index])
            res = numpy.percentile(valset, percentiles)
            cperc.append(list(res))
        ctimes = timesdict
        return [ctimes, cperc]


    @staticmethod
    def extract_aposteriori_only_samples(request):
        """
        * This method extracts the DV and the IPRED from the results file called "tab" for comparison. The residual is
        * just the difference between these two values.
        * The format of the "tab" file is defined in the table set at the end of each model file. The position of DV and
        * IPRED in this table has to be maintained to get the values out.
        * Returns: List of lists [[cobs],[cpred]] (measures, concentrations)
        """
        # open out.txt and get results
        try:
            f = open('aposteriori.sample.out.dat', 'r')
            cobs = []
            cpred = []
            f.readline()
            f.readline()
            for line in f:
                linetokens = line.split()
                if float(linetokens[1]) == 0:
                    # pdb.set_trace()
                    if request.drugModel.errorModel == 'exponential':
                        cobs.append(math.exp(float(linetokens[3])))
                        cpred.append(float(linetokens[4]))
                    else:
                        cobs.append(float(linetokens[3]))
                    cpred.append(float(linetokens[4]))

            f.close()

            return [cobs, cpred]

        except IOError as e:
            print(e.strerror + ' ' + e.filename)
            raise_exception(e.strerror + ' ' + e.filename)



    @staticmethod
    def extract_aposteriori_concentration(request):
        """
        * This method extracts the DV and the IPRED from the results file called "tab" for comparison. The residual is
        * just the difference between these two values.
        * The format of the "tab" file is defined in the table set at the end of each model file. The position of DV and
        * IPRED in this table has to be maintained to get the values out.
        * Returns: List of lists [[cobs],[cpred]] (measures, concentrations)
        """
        # open out.txt and get results
        try:
            f = open('aposteriori.out.dat', 'r')
            times = []
            cpred = []
            f.readline()
            f.readline()
            for line in f:
                linetokens = line.split()
                if float(linetokens[3]) == 0 and float(linetokens[1]) == 0:
                    times.append(float(linetokens[2]))
                    if request.drugModel.errorModel == 'exponential':
                        cpred.append(math.exp(float(linetokens[4])))
                    else:
                        cpred.append(float(linetokens[4]))
            f.close()

            return [times, cpred]

        except IOError as e:
            print(e.strerror + ' ' + e.filename)
            raise_exception(e.strerror + ' ' + e.filename)

    def extract_concentration(self):
        """
        * This method extracts the concentration values over many times throughout the treatment of one individual of a
        * subpopulation with zero inter/intra variance if self.with_variability = False
        * The NONMEM output file in this case is the "base_sim.fit".
        * Returns: List of lists [[times],[dataset]] (times, concentrations)
        """
        # pdb.set_trace()
        try:
            if self.with_variability:
                f = open('apriori.withVariability.out.dat', 'r')
            else:
                f = open('apriori.out.dat', 'r')
            patient_id = []
            dv = []
            times = []
            lcount = 1
            # pdb.set_trace()

            for line in f:
                linetokens = line.split()
                if float(linetokens[2]) == 0:
                    patient_id.append(str(int(float(linetokens[0]))))
                    times.append(float(linetokens[1]))
                    dv.append(float(linetokens[4]))
                    lcount = lcount + 1
            f.close()
            return [patient_id, times, dv]

        except IOError as e:
            print(e.strerror + ' ' + e.filename)
            raise_exception(e.strerror + ' ' + e.filename)

    @staticmethod
    def percentile(dataset):
        """
        * This method extracts the concentration values at deciles of monte carlo simulations run by NONMEM at each
        * measure time throughout the treatment.
        * The NONMEM output file in this case is the "base_sim.fit".
        * There is a line here where the number of patients must be manually maintained between ezechiel and NONMEM
        * Returns: List of lists [[ctimes],[cperc]] (times, concentrations at deciles at each measure)
        """
        cperc = []
        percentiles = [10, 20, 30, 40, 50, 60, 70, 80, 90]

        lines = []
        try:
            f = open('base_sim.fit', 'r')
            lines = f.readlines()
            f.close()

        except IOError as e:
            print(Fore.RED + e.strerror + ' ' + e.filename)
            raise_exception(e.strerror + ' ' + e.filename)
        dv = []
        times = []
        lcount = 0
        # here is the number of samples, must be manually maintained
        while lcount < (len(dataset.samples) + dataset.conc * dataset.cycles) * 1000:
            # line = f.readline()
            # pdb.set_trace()
            line = lines[lcount]
            linetokens = line.split()
            if float(linetokens[2]) == 0:
                times.append(float(linetokens[1]))
                val = float(linetokens[4])
                dv.append(val)
                lcount = lcount + 1
        # pdb.set_trace()
        uniquetimes = sorted(set(times))

        for utime in uniquetimes:
            indexes = [i for i, x in enumerate(times) if x == utime]
            dvsubset = []
            for ind in indexes:
                dvsubset.append(dv[ind])
            res = []
            dvsubsetarray = numpy.array(dvsubset)
            for perc in percentiles:
                p = numpy.percentile(dvsubsetarray, perc)
                # pdb.set_trace()
                res.append(float(p))
            #                  print p
            cperc.append(res)
        ctimes = list(uniquetimes)
        # pdb.set_trace()
        return [ctimes, cperc]
