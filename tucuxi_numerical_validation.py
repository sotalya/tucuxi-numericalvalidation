import os.path

import seaborn as sns

import configparser
import io
import glob


from argparse import ArgumentParser

import pandas as pd
import matplotlib.pyplot as plt

# With tucucli
#from tucuxi.processing.tucuxirun import TucuCliRun
#from tucuxi.data.query import Query
#from tucuxi.utils import *
#from tucuxi.importexport.exporttqf import  *

# With sotalya
from sotalya.tucuxi.utils import *
from sotalya.processing.tucuxirun import TucuPycliRun
from sotalya.importexport.exporttqf import *

from core.utils import set_output_folder_name
from core.systemtest import SystemTester

from core.report import Report
from population_numerical_validation import *
from nonmem_population import *


# reset all
init(autoreset=True)

sys.path.append(os.path.abspath(".."))

# globals
mgr = SystemTester()
SoftSelect = dict(both=0, onlytucuxi=1, onlynonmem=2)


dict_parameters = {"Ka": "KA", # tucuxi : nonmem
                   "CL": "CL",
                   "V": "VC",
                   "V1": "VC",
                   "V2": "VP",
                   "Q": "Q",
                   "Q2": "Q",
                   "Q3": "Q3",
                   "V3": "V3",
                   "Tlag": "ALAG1",
                   "F": "F1",
                   "Vmax":"VMAX",
                   "Km": "KM",
                   "Ktr":"KTR",

                   "AUC0_24": "AUC0_24",
                   "CMIN": "CMIN",
                   "CMAX": "CMAX",
                   "IPRED": "IPRED"}


dict_model = {"IV": ["210000", "220000", "230000",
                     "210100", "220100", "210200", "220200", "210300", "220300",
                     "210400", "220400", "210010", "210020", "210009",
                     "210001", "210002", "210003", "210004", "210007", "210008",
                     "220001", "220002", "220003", "220004", "220007","220009",
                     "230001", "230002", "230003", "230004", "230007", "230009"] ,
              "Oral": ["110000", "120000", "130000", "111000", "112000", "123000", "124000", "125000",
                       "126000", "127000", "128000", "110100", "120100", "110200", "120200",
                       "110300", "120300", "110400", "120400", "110010", "110020",
                       "110001", "110002", "110003", "110004", "110005", "110006", "110007", "110008",
                       "110009",
                       "120001", "120002", "120003", "120004", "120007", "120009",
                       "130001", "130002", "130003", "130004", "130007", "130009"],
              "Bolus": ["310000", "320000", "330000", "310100", "320100", "310200", "320200", "310300", "320300",
                        "310400", "320400", "310010", "310020", "310009",
                        "310001", "310002", "310003", "310004", "310007", "310008",
                        "320001", "320002", "320003", "320004", "320007", "320009",
                        "330001", "330002", "330003", "330004", "330007", "330009"]}


"""
*************************************** Args ***************************************************************************
"""

time_tucuxi=0

config = configparser.ConfigParser()
parser = ArgumentParser(description='vancomycin',
                        prog='experiment.py')


configValues = {"tucucli": "",
                "drugspath_validation": "",
                "nonmemmodelspath": "",
                "queriespath": "",
                "queryfile": "",
                "output": "",
                "inputfilename": "",
                "requesttemplate": "",
                "listtemplate": "",
                "runonserver": False,
                "exportprf": False,
                "computation": ""}



drugname = 'virtualdrug'

errorCounter = 0

errorfile: io.TextIOWrapper


def parse_the_args():
    """
    * This method parses the command line arguments. This project relies on windows having unix tools
      available in the console. We get this by installing git and allowing to use it via cmd.exe.
    """

    parser.add_argument('-q', '--quiet',
                        help='Dont show all the data, kind of like release mode.', action='store_true')

    parser.add_argument('-queriespath', type=str, default='',
                        help='Import folder paths for tucuxi queries (tqf)')

    parser.add_argument('-nonmemmodelspath', type=str,
                        help='Import folder paths for NonMem models')

    parser.add_argument('-noclean',
                        help='If true will not cleanup Tucuxi afterwards.',
                        action="store_true")

    parser.add_argument('-queryfile', type=str, default="", dest='queryfile',
                        help='Import file path for query (tqf). If this option is set,'
                             ' then only this specific query is executed')

    parser.add_argument('-drugspath_validation', type=str, dest='drugspath_validation',
                        help='Import folder paths for Tucuxi drugfiles (tdd)')

    parser.add_argument('-d', type=str, nargs='+',
                        help='Drug IDs of drugs to be tested. (e.g. ch.heig-vd.tucuxi.imatinib)')

    parser.add_argument('-t', '--tucucli', type=str, dest='tucucli',
                        help='Command to execute Tucuxi cli, default is tucucli')

    parser.add_argument('-inputfilename', type=str, dest='inputfilename',
                        help='Import data file (.csv)', default="")

    parser.add_argument('-inputdoses', type=str, dest='inputdoses',
                        help='Import doses file (.csv)', default="")

    parser.add_argument('-inputlevels', type=str, dest='inputlevels',
                        help='Import levels file (.csv)', default="")

    parser.add_argument('-inputscr', type=str, dest='inputscr',
                        help='Import sources file (.csv)', default="")

    parser.add_argument('-output', type=str, dest='output',
                        help='Export folder path for results', default="")

    parser.add_argument('-whichsoft', type=int,
                        help='Set to run both (0), only tucuxi (1), or only nonmem (2).',
                        choices=[0, 1, 2], default=0)

    parser.add_argument('-n', '--nonmem', type=str,
                        help='Command to execute nonmem. default = /opt/nm72/run/nmfe72',
                        default='/opt/nm72/run/nmfe72')

    parser.add_argument('-ss', '--steadystate', type=int,
                        help='Whether or not to calculate at steady state (0 or 1),'
                             ' value of 2 will approx steady state with 10.5 cycles',
                        default=0)

    parser.add_argument('-g', '--graph', help='If set, then graphs are generated.',
                        action="store_true")

    parser.add_argument('-outputlist', type=str, dest='outputlist',
                        help='Export folder path for list results', default="")

    parser.add_argument('-cache', dest='use_cache', help='Use cache instead of calculating with NONMEM.',
                        action="store_true", default=False)

    parser.add_argument('-requesttemplate', type=str, dest='requesttemplate',
                        help='Request template file (.xml)', default="")

    parser.add_argument('-listtemplate', type=str, dest='listtemplate',
                        help='List template file (.xml)', default="")

    parser.add_argument('-runonserver', type=bool, dest='runonserver',
                        help='Run the experiment on server ?', default="")

    parser.add_argument('-exportprf', type=bool, dest='exportprf',
                        help='Export pending request (.prf) ?', default="")

    parser.add_argument('-computation', type=str, dest='computation',
                        help='Computation API URL', default="")

    return parser.parse_args()


def args_to_dictionary(args):
    configValues["tucucli"] = args.tucucli
    configValues["drugspath_validation"] = args.drugspath_validation
    configValues["nonmemmodelspath"] = args.nonmemmodelspath
    configValues["queriespath"] = args.queriespath
    configValues["queryfile"] = args.queryfile
    configValues["output"] = args.output
    configValues["requesttemplate"] = args.requesttemplate
    configValues["inputfilename"] = args.inputfilename
    configValues["listtemplate"] = args.listtemplate
    configValues["runonserver"] = args.runonserver
    configValues["exportprf"] = args.exportprf
    configValues["computation"] = args.computation


def config_ini_file(args):
    args_to_dictionary(args)

    try:
        with open('config.ini') as f:
            config.read_file(f)
    except IOError:
        # Le fichier n'existe pas et chemins sont donnés en cmd
        if are_args_given(configValues):

            config['PATH'] = {'tucucli': configValues["tucucli"],
                              'drugspath_validation': configValues["drugspath_validation"],
                              'nonmemmodelspath': configValues["nonmemmodelspath"],
                              'queriespath': configValues["queriespath"],
                              'output': configValues["output"],
                              'queryfile': configValues["queryfile"],
                              'inputfilename': configValues["inputfilename"],
                              'listtemplate': configValues["listtemplate"],
                              'requesttemplate': configValues["requesttemplate"]}

            config['BOOLEAN'] = {'runonserver': str(configValues["runonserver"]),
                                 'exportprf': str(configValues["exportprf"])}

            config['URL'] = {'computation': configValues["computation"]}
            with open('config.ini', 'w') as iniConfigFile:
                config.write(iniConfigFile)
        else:
            print(Back.RED + '*************************************')
            print(Fore.RED + 'Error with the execution of experiment.py')
            print(Fore.RED + 'Not found : config.ini')
            print(Back.RED + '*************************************')
            sys.exit()

    config.read('config.ini')

    choose_data_from(config, configValues, PATH, "tucucli")
    choose_data_from(config, configValues, PATH, "drugspath_validation")
    choose_data_from(config, configValues, PATH, "nonmemmodelspath")
    choose_data_from(config, configValues, PATH, "queriespath")
    choose_data_from(config, configValues, PATH, "output")
    choose_data_from(config, configValues, PATH, "queryfile")
    choose_data_from(config, configValues, PATH, "requesttemplate")
    choose_data_from(config, configValues, PATH, "listtemplate")
    choose_data_from(config, configValues, PATH, "inputfilename")

    choose_data_from(config, configValues, BOOLEAN, "runonserver")
    choose_data_from(config, configValues, BOOLEAN, "exportprf")

    choose_data_from(config, configValues, URL, "computation")

    with open('config.ini', 'w') as iniConfigFile:
        config.write(iniConfigFile)

    if are_paths_wrong(configValues):
        print(Back.RED + '*************************************')
        print(Fore.RED + 'Error with the execution of experiment.py')
        print(Fore.RED + 'Wrong or missing path in config.ini')
        print(Back.RED + '*************************************')
        sys.exit()


"""
*************************************** Create Data and Compare ********************************************************
"""


def create_data(pop: Population, administration_route: str) -> List[ExpPatient]:
    """
            * This function generates a list of ExpPatient, wich contains information on dosage,
            measured concentrations and TDM results for each patient

                Args:
                * pop (Population): the population of patients

                return:
                * List of ExpPatient
    """

    patients = []
    for i in range(pop.nb_patients):
        current_patient = ExpPatient()
        current_patient.patient_id = str(pop.population[i].patient_id)

        # Information about TDM
        tdm_episode = EpisodeVirtualDrug(episode_id=str(pop.population[i].nb_tdm),
                                               birthdate=pop.population[i].birthdate,
                                            sampling_group = pop.population[i].sampling_group,
                                               etacl=str(pop.population[i].etacl),
                                               startonlydate=unis_str_to_datetime(
                                                   pop.population[i].dose.start_date.strftime("%m/%d/%Y") +
                                                   'T00:00:00'))

        # Information about dosage regimen
        for j in range(0, pop.population[i].dose.nb_doses):
            day = DayNumericalValidation()
            day.dose = str(pop.population[i].dose.amt)
            day.number_of_doses = 1

            start_time = timedelta(hours=pop.population[i].dose.start_date.time().hour,
                                   minutes=pop.population[i].dose.start_date.time().minute)

            day.start_datetime = tdm_episode.startonlydate + start_time + timedelta(
                hours=j * pop.population[i].dose.interval)

            day.interval = str(pop.population[i].dose.interval)

            if administration_route=="Bolus":
                day.infusion_time = "0"
            else:
                day.infusion_time = str(pop.population[i].dose.infusion_time)

            day.end_datetime = day.start_datetime + timedelta(hours=float(day.interval)) * float(day.number_of_doses)

            day.weight = str(pop.population[i].bw)
            day.is_male = str(pop.population[i].is_male)

            tdm_episode.days.append(day)

        tdm_episode.samples += pop.population[i].sample

        current_patient.tdm.append(tdm_episode)

        patients.append(current_patient)

    return patients



"""
*************************************** Create Query *******************************************************************
"""


def create_query(patient: ExpPatient, model_id: str, folder_name: str, software: str,
                     administration_route: str):
    """
                        * This function creates the .tqf file required by tucuxi and nonmem

                            Args:

                            * patient (ExpPatient): list of patients
                            * model_id (str): the name of the popPK model to be used
                            * folder_name (str): the name of folder in which to save file
                            * software (str): nonmem or tucuxi
                            * administration_route (str): IV, oral, bolus
    """



    tdm = len(patient.tdm) -1
    episode = patient.tdm[tdm]

    templatefilename = configValues['queryfile']

    # Initiation of the Query
    query = Query()
    query.queryId = patient.patient_id + '_' + model_id
    query.patientId = patient.patient_id
    query.date = datetime.now().strftime('%Y-%m-%dT%H:%M:%S')

    drug = Drug()
    drug.activePrinciple = drugname
    drug.drugId = drugname

    # The TDM first dose protocol consists of taking measurements directly after the first dose and adjusting
    # according to the prediction made at steady-state. The episode must therefore be recopied in order to add
    # doses for adjustment at steady-state.

    episodes = []
    for ep in range(tdm + 1):
        episodes.append(patient.tdm[ep])

    # Add covariates information
    create_covariates_virtualdrug(query, episodes, model_id)

    # Add dosage history
    if administration_route == "IV":
        _ = create_dosage_history(drug.dosageHistory, episodes)
    elif administration_route == "Bolus":
        _ = create_dosage_history_bolus(drug.dosageHistory, episodes)
    else:
        abs_model = "extra"
        if model_id == "ch.tucuxi.virtualdrug.mod111000":
            abs_model = "extra.lag"
        _ = create_dosage_history_peros(drug.dosageHistory, episodes, abs_model)

    # Add samples
    if episode.samples[len(episode.samples) - 1].value != '':

        for episode in episodes:

            for sample_index, sample in enumerate(episode.samples):
                the_date = sample.sample_date

                if sample.value != '' :
                    sample_value = float(sample.value)
                    if sample_value < 0.1:
                        sample_value = 0.1

                    drug.samples.append(Sample.create_sample('c1', the_date, 'virtualdrug',
                                                             sample_value, 'ug/l'))


    query.drugs.append(drug)

    if software == "nonmem":
        if episode.samples[len(episode.samples) - 1].value == '':  # query for nonmem
            # Prediction a priori of the samples
            suffixe_tqf = "_nonmem"
            sample_date = []
            for s in episode.samples:
                sample_date.append(s.sample_date)

            computing_option = ComputingOption(ParametersTypeEnum.apriori,
                                               CompartmentOptionEnum.allActiveMoieties, True, True, True)

            request_apriori = Request()
            request_apriori.requestId = model_id + '.apriori.point'
            request_apriori.drugId = drugname
            request_apriori.drugModelId = model_id
            request_apriori.computingTraits = PredictionAtTimesTraits.create_prediction_at_times_traits(sample_date,
                                                                                                       computing_option)

            query.requests.append(request_apriori)
        else:
            # Prediction a posteriori
            suffixe_tqf = "_nonmem"
            sample_date = []
            for s in episode.samples:
                sample_date.append(s.sample_date)

            computing_option = ComputingOption(ParametersTypeEnum.aposteriori,
                                               CompartmentOptionEnum.allActiveMoieties, True, True, True)

            request_apriori = Request()
            request_apriori.requestId = model_id + '.aposteriori.point'
            request_apriori.drugId = drugname
            request_apriori.drugModelId = model_id
            request_apriori.computingTraits = PredictionAtTimesTraits.create_prediction_at_times_traits(sample_date,
                                                                                                        computing_option)

            query.requests.append(request_apriori)

    else: # query for tucuxi
        suffixe_tqf = "_tucuxi"

        nb_points_per_hour = 20

        computing_option = ComputingOption(ParametersTypeEnum.aposteriori,
                                           CompartmentOptionEnum.allActiveMoieties, True, True, True)

        # A posteriori evaluation
        request_aposterori_stat = Request()
        request_aposterori_stat.requestId = model_id + '.stats'
        request_aposterori_stat.drugId = drugname
        request_aposterori_stat.drugModelId = model_id
        aposteriori_date = episodes[0].days[0].start_datetime
        aposteriori_end = episodes[-1].days[-1].end_datetime + timedelta(hours=24)
        request_aposterori_stat.computingTraits = PredictionTraits.create_prediction_traits(
            nb_points_per_hour, aposteriori_date, aposteriori_end, computing_option)

        query.requests.append(request_aposterori_stat)

        # A posteriori prediction of samples
        sample_date = []
        sample_times = []
        for s in episode.samples:
            sample_date.append(s.sample_date)
            sample_times.append((s.sample_date - aposteriori_date).seconds/3600)

        request_aposterori_pred = Request()
        request_aposterori_pred.requestId = model_id + '.pred'
        request_aposterori_pred.drugId = drugname
        request_aposterori_pred.drugModelId = model_id
        request_aposterori_pred.computingTraits = PredictionAtTimesTraits.create_prediction_at_times_traits(sample_date,
                                                                                                   computing_option)

        query.requests.append(request_aposterori_pred)

    print('Computation run locally')
    output_tqf_filename = os.path.join(folder_name, query.queryId + suffixe_tqf + '.tqf')
    print('.tqf ', output_tqf_filename)
    exporter_tqf = ExportTqf()
    exporter_tqf.export_to_file(query, output_tqf_filename, templatefilename)


"""
*************************************** NONMEM *************************************************************************
"""


class RunNonmem:
    """
                            * This class executes the nonmem run

                                Args:

                                * args
                                * nmresult []: list of predicted concentrations : [ID], [Time of samples], [
                                Predicted concentrations].
                                * report (Report): reports run problems
        """

    def __init__(self, args, report: Report):

        self.args = args
        self.nmresult = pd.DataFrame.empty
        self.report = report

    def run_query(self, list_queries: List[Query], folder_name: str, nm_cmd: str, model_name: str):
        """
                            * This function executes a nonmem run based on a list of queries

                                Args:

                                * list_queries (Query): list of queries for each patient
                                * folder_name (str): name of the folder to save results
                                * nm_cmd (str): command used to execute run on nonmem
                                * model_name (str): model to be executed
        """

        if self.args.whichsoft != SoftSelect['onlynonmem']:

            self.report.start_query(list_queries[0].sourceFile + " : " + list_queries[0].queryId)

            if True:
                try:
                    list_requests= []
                    list_time_cmin = pd.DataFrame(columns=[" ID"] + ["TIME_CMIN"])

                    for i in range(len(list_queries)):
                        for ri, response in enumerate(list_queries[i].requests):
                            list_requests.append(response)

                            times_samples = []
                            for d in response.computingTraits.computingTraits.dates.contents:
                                date_sample = datetime.strptime(d.contents[0], '%Y-%m-%dT%H:%M:%S')
                                times_samples.append((date_sample - list_queries[i].drugs[0].dosageHistory.dosageTimeRanges[
                                    0].start).total_seconds() / 3600)

                            # add Cmin:
                            times_samples.append((list_queries[i].drugs[0].dosageHistory.dosageTimeRanges[-1].end -
                                                  list_queries[i].drugs[0].dosageHistory.dosageTimeRanges[0].start).total_seconds() / 3600)

                            id_patient = float(list_queries[i].queryId.split("_")[0])
                            temp_cmin = pd.DataFrame([[id_patient,(list_queries[i].drugs[0].dosageHistory.dosageTimeRanges[-1].end -
                                                  list_queries[i].drugs[0].dosageHistory.dosageTimeRanges[0].start).total_seconds() / 3600]], columns=[' ID','TIME_CMIN'])

                            if list_time_cmin.empty:
                                list_time_cmin = temp_cmin
                            else:
                                list_time_cmin = pd.concat([list_time_cmin,temp_cmin], ignore_index=True)

                            times_samples.sort()

                            # nonmem data file
                            _ = NonmemData(list_queries[i], response,
                                            drug_model=mgr.drug_dict[list_queries[0].requests[0].drugModelId],
                                            times=times_samples, ss=self.args.steadystate,
                                            with_variability=True,
                                            foldername=folder_name)


                    name_file=''
                    if sys.platform == "linux" or sys.platform == "linux2" or sys.platform == "darwin":
                        name_file = folder_name+"/nonmem/"
                    elif sys.platform == "win32":
                        name_file = folder_name + "\\nonmem\\"

                    # execute nonmem
                    print(Fore.CYAN + nm_cmd + ' ' + name_file +model_name + ' output.txt')
                    os.chdir(name_file)
                    subprocess.run(nm_cmd + ' ' + name_file +model_name + ' output.txt', shell=True)
                    os.chdir(folder_name)

                    # Results
                    file_result = glob.glob(name_file+"*.tab")
                    pd_results = pd.read_table(file_result[0], skiprows=1, delimiter=",")

                    for c in pd_results.columns:
                        c_new = c.split("   ")[0]
                        pd_results.rename(columns={c:c_new}, inplace=True)


                    # Add TIME_CMIN and TIME CMAX to the results
                    pd_results = pd_results.merge(list_time_cmin, on=' ID')

                    # Format results
                    pd_results = results_format(pd_results)

                    self.nmresult = pd_results

                except RuntimeError as err:
                    self.report.jrep.add_testcase(False, 'nonmem_dataset', 'dataset', err.args[0])

        else:
            self.report.jrep.add_testcase(False, 'analysis', 'run', 'Computation not compared by Python')


    def select_soft_from_args(self):
        """
        Select which software to run: tucuxi, nonmem or both. Make sure the executables are present on the machine.
        :return:
        """
        if self.args.whichsoft != SoftSelect['onlytucuxi']:
            if not os.path.isfile(self.args.nonmem):
                if self.args.use_cache:
                    print(Fore.CYAN + 'Cache mode set. All NONMEM data will be exclusively loaded from cache.')
                else:
                    try:
                        Popen(['rm -rf ' + get_output_folder_name()], shell=True)
                    except CalledProcessError as err:
                        print(Fore.RED + err.output)
                    sys.exit('Nonmem executable not found at: ' + self.args.nonmem)
        if self.args.whichsoft != SoftSelect['onlynonmem']:
            if not os.path.isfile(configValues['tucucli']) and not os.path.isfile(configValues['tucucli'] + '.exe'):
                try:
                    Popen(['rm -rf ' + get_output_folder_name()], shell=True)
                except CalledProcessError as err:
                    print(Fore.RED + err.output)
                sys.exit('Tucuxi cli executable not found at: ' + configValues['tucucli'] + '.exe or ' + configValues['tucucli'])

def results_format(pd_results):

    # Calcul AUC on 24h
    max_time_rows = pd_results[pd_results['TIME'] == pd_results['TIME_CMIN']].copy()
    max_time_rows.drop_duplicates(subset=" ID", keep='last', inplace=True)
    max_time_minus_24 = pd_results.groupby(' ID')['TIME_CMIN'].max() - 24
    max_time_minus_24_df = pd.DataFrame(
        {' ID': max_time_minus_24.index, 'Max_TIME_minus_24': max_time_minus_24.values})
    pd_results = pd_results.merge(max_time_minus_24_df, on=' ID')
    max_time_minus_24_rows = pd_results[pd_results['TIME'] == pd_results['Max_TIME_minus_24']].copy()
    max_time_minus_24_rows.drop_duplicates(subset=" ID", keep='last', inplace=True)
    result = max_time_rows['AUC'].values - max_time_minus_24_rows['AUC'].values
    result = pd.DataFrame({' ID': max_time_rows[' ID'].values, 'AUC0_24': result})
    pd_results = pd_results.merge(result, on=' ID')
    pd_results.drop(columns="Max_TIME_minus_24", inplace=True)

    # Add column CMIN
    result = pd_results[pd_results['TIME'] == pd_results['TIME_CMIN']].copy()
    result.drop_duplicates(subset=" ID", keep='last', inplace=True)
    result = result[[" ID", "IPRED"]]
    result.rename(columns={"IPRED": "CMIN"}, inplace=True)
    pd_results = pd_results.merge(result, on=' ID')


    # Correct column Cmax
    result = pd_results[pd_results['CMAX'] == pd_results['CMAX'].groupby(pd_results[' ID']).transform('last')].copy()
    result.drop_duplicates(subset=" ID", keep='last', inplace=True)
    result = result[[" ID", "CMAX"]]
    pd_results.drop(columns="CMAX", inplace=True)
    pd_results = pd_results.merge(result, on=' ID')

    # Correct column Tmax
    result = pd_results[
        pd_results['TMAX'] == pd_results['TMAX'].groupby(pd_results[' ID']).transform('last')].copy()
    result.drop_duplicates(subset=" ID", keep='last', inplace=True)
    result = result[[" ID", "TMAX"]]
    pd_results.drop(columns="TMAX", inplace=True)
    pd_results = pd_results.merge(result, on=' ID')

    # drop row of Cmin
    pd_results['is_CMIN'] = np.where(pd_results['TIME'] == pd_results['TIME_CMIN'], True, False)
    pd_results['count'] = pd_results.groupby(['TIME', ' ID'])['TIME'].transform('count')
    pd_results = pd_results[(pd_results['is_CMIN'] == False) | (pd_results['count'] > 1)]

    # Delete columns is_CMIN, TIME_CMIN, count
    pd_results = pd_results.drop(columns='count')
    pd_results = pd_results.drop(columns='TIME_CMIN')
    pd_results = pd_results.drop(columns='is_CMIN')

    pd_results.drop_duplicates(subset=(" ID", "TIME", "EVID"), keep='last', inplace=True)

    # keep only sample rows
    pd_results.drop(pd_results[pd_results['EVID'] == 1].index, inplace=True)

    # delete columns EVID, PRED, ETA
    eta_to_delete = [col for col in pd_results.columns if 'ETA' in col]
    pd_results = pd_results.drop(eta_to_delete, axis=1)
    pd_results = pd_results.drop(columns='EVID')
    pd_results = pd_results.drop(columns='PRED')
    pd_results = pd_results.drop(columns='AUC')

    pd_results.reset_index(drop=True, inplace=True)
    pd_results.sort_values(by=[" ID", "TIME"], inplace=True)

    return pd_results




def transfer_models_to_local(foldername_nonmem, filepaths_model_origin,
                             model_name):
    """
    * Transfers the population file and models from the original place, to the local ouput folder, in a new sub-folder.
    * Returns nothing

        Args:
        * foldername_nonmem (str): path to create the new sub-folder
        * filepaths_model_origin (str): path giving the location of the models files
        * model_name (str): name of the model
    """
    # get working directory
    filepath = get_output_folder_name()

    # if a nonmem folder does not exist already, create one to run the nonmem command inside
    print(Fore.CYAN + filepath+foldername_nonmem)
    if not os.path.isdir(filepath+foldername_nonmem):
        try:
            subprocess.run('mkdir ' + filepath+foldername_nonmem, shell=True)
        except CalledProcessError as err:
            print(Fore.RED + err.output)
            sys.exit()

    if not os.path.isfile(model_name):
        try:
            print(Fore.CYAN + model_name)
            origin = ''
            dest = ''
            if sys.platform == "linux" or sys.platform == "linux2" or sys.platform == "darwin":
                origin = filepaths_model_origin+"/"+model_name
                dest = filepath+foldername_nonmem+"/"+model_name
            elif sys.platform == "win32":
                origin = filepaths_model_origin+"\\"+model_name
                dest = filepath+foldername_nonmem+"\\"+model_name

            shutil.copyfile(origin, dest)
        except CalledProcessError as err:
            print(Fore.RED + err.output)
            sys.exit()

"""
*************************************** Tucuxi *************************************************************************
"""


def run_tucuxi_query(query: Query):
    """
             * This function executes tucucli

                Args:

                * query (Query): query to be executed on tucucli
                * patient (ExpPatient): to update tucucli results
    """

    # If we use tucucli instead of sotalya :
    # output_dir = os.path.join(get_output_folder_name())
    # output_file_name = ''
    # if sys.platform == "linux" or sys.platform == "linux2" or sys.platform == "darwin":
    #     output_file_name = output_dir + '/' + query.queryId + '.xml'
    # elif sys.platform == "win32":
    #     output_file_name = output_dir + '\\' + query.queryId + '.xml'

    #tucuxi_run = TucuCliRun(configValues["tucucli"], configValues["drugspath_validation"], output_file_name)

    patient_id = query.queryId.split("_")[0]

    global time_tucuxi
    start_time = datetime.now()

    tucuxi_run = TucuPycliRun(configValues["drugspath_validation"])

    # For tucucli :
    #tucuxi_run = TucuCliRun(configValues["tucucli"], configValues["drugspath_validation"], output_file_name)

    # Sotalya:

    query_response = tucuxi_run.run_tucuxi_from_file(query.sourceFile)

    end_time = datetime.now()

    time_tucuxi = time_tucuxi + (end_time - start_time).seconds/60

    if query_response is None:
        print(Fore.RED + "No result from Tucuxi")
        return

    # Results
    result_aposteriori = {}

    for response in query_response.responses:

        # A posteriori results
        if response.requestType == "prediction":

            apostsinglepredictiondata = response.computingTrait.cycleDatas

            # PK parameters
            for p in apostsinglepredictiondata[-1].parameters:
                result_aposteriori[p.id] = p.value

            result_aposteriori["AUC0_24"] = apostsinglepredictiondata[-2].statistics.auc24
            result_aposteriori["mean"] = apostsinglepredictiondata[-2].statistics.mean
            result_aposteriori["CMAX"] = apostsinglepredictiondata[-2].statistics.peak
            result_aposteriori["CMIN"] = apostsinglepredictiondata[-2].statistics.residual

        # Prediction of concentrations at sampling time
        if response.requestType == "singlePoints":
            dosage_start = query.drugs[0].dosageHistory.dosageTimeRanges[0].start

            result_aposteriori["TIME"] = []
            result_aposteriori["IPRED"] = []


            for d in response.computingTrait.points:
                result_aposteriori["TIME"].append(((datetime.strptime(d.time, '%Y-%m-%dT%H:%M:%S') - dosage_start).days*24 +
                                              (datetime.strptime(d.time, '%Y-%m-%dT%H:%M:%S') - dosage_start).seconds/3600))
                result_aposteriori["IPRED"].append(d.value)

    pd_temp = pd.DataFrame(columns=["ID"] + list(result_aposteriori.keys()))

    pd_temp["ID"] = [patient_id] * len(result_aposteriori["TIME"])

    for key, value in result_aposteriori.items():

        if (key=="Vmax" and ("mod110200" in query_response.queryId or "mod210200" in query_response.queryId or
        "mod120200" in query_response.queryId or "mod220200" in query_response.queryId or "mod110400" in query_response.queryId or
             "mod210400" in query_response.queryId or "mod120400" in query_response.queryId or "mod220400" in query_response.queryId or
                "mod310200" in query_response.queryId or "mod320200" in query_response.queryId or
                "mod310400" in query_response.queryId or "mod320400" in query_response.queryId )):
            pd_temp[key] = [value/1000] * len(result_aposteriori["TIME"])
        else :
            pd_temp[key] = [value] * len(result_aposteriori["TIME"])

    pd_temp["TIME"] = result_aposteriori["TIME"]
    pd_temp["IPRED"] = result_aposteriori["IPRED"]

    pd_temp.rename(columns=dict_parameters, inplace=True)

    result_tucuxi = pd_temp

    return result_tucuxi



"""
*************************************** Tucuxi numerical validation process **************************************************************************
"""


def find_administration_route(model_id: str):
    key_adm = None
    for key, value in dict_model.items():
        if model_id in value:
            key_adm = key
    return key_adm  # Returns None if the value is not found in the dictionary



def run_numerical_validation(pop: Population, model_id: str, foldername_init: str, args):
    """
             * This function performs the virtual TDM process

                Args:

                * pop (Population): the patient population to be used for virtual TDM
                * model_id (str): the name of the popPK model to be used
                * foldername_init: the name of the initial folder to save the results
                * args
    """

    # Oral or IV model:
    administration_route = find_administration_route(model_id)


    # creation of the data
    patients = create_data(pop, administration_route)


    folder_name = os.path.join(foldername_init, "run_" + model_id)

    print(Fore.CYAN + folder_name)
    if not os.path.isdir(folder_name):
        try:
            subprocess.run('mkdir ' + folder_name, shell=True)
        except CalledProcessError as err:
            print(Fore.RED + err.output)
            sys.exit()

    if sys.platform == "linux" or sys.platform == "linux2" or sys.platform == "darwin":
        set_output_folder_name(folder_name)
    elif sys.platform == "win32":
        set_output_folder_name(folder_name)

    global mgr
    mgr = SystemTester()


    file_drug = list()
    file_drug.append(configValues["drugspath_validation"])

    mgr.import_drugs(file_drug)

    model_tucuxi = "ch.tucuxi.virtualdrug.mod" + model_id


    # For each patient, we create the .tqf associated
    for i in range(len(patients)):
        create_query(patients[i], model_id=model_tucuxi, folder_name=folder_name, software="nonmem",
                     administration_route=administration_route)


    # nonmem : simulation model
    model_nonmem = "run" + model_id + "_simulations.mod"

    folder_nonmem=''
    if sys.platform == "linux" or sys.platform == "linux2" or sys.platform == "darwin":
        folder_nonmem = "/nonmem"
    elif sys.platform == "win32":
        folder_nonmem = "\\nonmem"

    filepath_model_origin = configValues['nonmemmodelspath']

    transfer_models_to_local(foldername_nonmem=folder_nonmem,
                             filepaths_model_origin= filepath_model_origin,
                             model_name = model_nonmem)

    report = Report(args.graph, args.whichsoft)
    nonmemrun_apriori = RunNonmem(args, report)

    # check if necessary executables exist where we expect them
    nonmemrun_apriori.select_soft_from_args()

    # import the queries
    mgr.import_queries([folder_name])
    mgr.print_drug()

    list_queries = []
    for(k, v) in mgr.query_dict.items():
        list_queries.append(v)

    # Run nonmem for all the queries
    nonmemrun_apriori.run_query(list_queries, folder_name, args.nonmem, model_nonmem)


    # Delete nonmem files
    directory_path = ''
    if sys.platform == "linux" or sys.platform == "linux2" or sys.platform == "darwin":
        directory_path = folder_name + '/nonmem/'
    elif sys.platform == "win32":
        directory_path = folder_name + '\\nonmem\\'

    delete_files_and_subdirectories(directory_path)

    # remove nonmem tqf
    name_path= ''
    if sys.platform == "linux" or sys.platform == "linux2" or sys.platform == "darwin":
        name_path = "/*_nonmem.tqf"
    elif sys.platform == "win32":
        name_path = "\\*_nonmem.tqf"

    files = glob.glob(folder_name + name_path)
    for f in files:
        os.remove(f)

    result_tucuxi = pd.DataFrame.empty

    # Tucucli run
    for i in range(len(patients)):
        # Updating concentration values
        patients[i].update_sampling_patient([list(nonmemrun_apriori.nmresult[" ID"].astype(int).astype(str)),
                                            list(nonmemrun_apriori.nmresult["TIME"]), list(nonmemrun_apriori.nmresult["DV"])])

        # creation of the query for tucucli
        create_query(patients[i], model_id=model_tucuxi, folder_name=folder_name, software="tucuxi",
                     administration_route=administration_route)

    start_time_total = datetime.now()

    # import the queries
    mgr.import_queries([folder_name])
    mgr.print_drug()

    item_patient = 0
    for (k, v) in mgr.query_dict.items():
        print("Patient ", str(v.queryId.split("_")[0]))

        if k[-10:] != "nonmem.tqf":

            result = run_tucuxi_query(v)
            if item_patient == 0:
                result_tucuxi = result
            else:
                result_tucuxi = pd.concat([result_tucuxi, result], ignore_index=True)

            item_patient +=1
            end_time = datetime.now()
            time_total = (end_time - start_time_total).seconds / 60

            print("Total time: " + str(time_total))
            print("Tucuxi time: " + str(time_tucuxi))

    # remove tucuxi tqf
    namefile = ''
    if sys.platform == "linux" or sys.platform == "linux2" or sys.platform == "darwin":
        namefile = "/*_tucuxi.tqf"
    elif sys.platform == "win32":
        namefile = "\\*_tucuxi.tqf"

    files = glob.glob(folder_name + namefile)
    for f in files:
        os.remove(f)

    # Run nonmem a posteriori
    for i in range(len(patients)):
        # create query for nonmem
        create_query(patients[i], model_id=model_tucuxi, folder_name=folder_name, software="nonmem",
                     administration_route=administration_route)

    if sys.platform == "linux" or sys.platform == "linux2" or sys.platform == "darwin":
        set_output_folder_name(folder_name)
    elif sys.platform == "win32":
        set_output_folder_name(folder_name)

    model_nonmem = "run"+ model_id +".mod"

    transfer_models_to_local(foldername_nonmem=folder_nonmem,
                             filepaths_model_origin= filepath_model_origin,
                             model_name = model_nonmem)

    nonmemrun_aposteriori = RunNonmem(args, report)

    # check if necessary executables exist where we expect them
    nonmemrun_aposteriori.select_soft_from_args()

    # import the queries
    mgr.query_dict = {}
    mgr.import_queries([folder_name])
    mgr.print_drug()

    list_queries = []
    for (k, v) in mgr.query_dict.items():
        list_queries.append(v)

    # Run nonmem for all the queries
    nonmemrun_aposteriori.run_query(list_queries, folder_name, args.nonmem, model_nonmem)

    # Delete nonmem files
    directory_path = ''
    if sys.platform == "linux" or sys.platform == "linux2" or sys.platform == "darwin":
        directory_path = folder_name + '/nonmem/'
    elif sys.platform == "win32":
        directory_path = folder_name + '\\nonmem\\'

    delete_files_and_subdirectories(directory_path)

    # remove nonmem tqf
    name_path = ''
    if sys.platform == "linux" or sys.platform == "linux2" or sys.platform == "darwin":
        name_path = "/*_nonmem.tqf"
    elif sys.platform == "win32":
        name_path = "\\*_nonmem.tqf"

    files = glob.glob(folder_name + name_path)
    for f in files:
        os.remove(f)


    # rename columns for nonmem or tucuxi
    for c in nonmemrun_apriori.nmresult.columns:
        if c != " ID" and c!="TIME" and c!="S2_SAMPLING":
            c_new = c + "_true"
            nonmemrun_apriori.nmresult.rename(columns={c: c_new}, inplace=True)
        if c==" ID":
            nonmemrun_apriori.nmresult.rename(columns={c: "ID"}, inplace=True)

    for c in nonmemrun_aposteriori.nmresult.columns:
        if c != " ID" and c!="TIME" and c!="S2_SAMPLING":
            c_new = c + "_nonmem"
            nonmemrun_aposteriori.nmresult.rename(columns={c: c_new}, inplace=True)
        if c==" ID":
            nonmemrun_aposteriori.nmresult.rename(columns={c: "ID"}, inplace=True)

    for c in result_tucuxi.columns:
        if c != "ID" and c != "TIME":
            c_new = c + "_tucuxi"
            result_tucuxi.rename(columns={c: c_new}, inplace=True)

    # convert nonmem ID in str
    nonmemrun_apriori.nmresult["ID"] = nonmemrun_apriori.nmresult["ID"].astype(int).astype(str)
    nonmemrun_aposteriori.nmresult["ID"] = nonmemrun_aposteriori.nmresult["ID"].astype(int).astype(str)

    # correct times
    result_tucuxi["TIME"] = round(result_tucuxi["TIME"], ndigits=4)

    # merge results of nonmem a priori/a posteriori and tucuxi
    global_results = pd.merge(nonmemrun_apriori.nmresult.copy(), result_tucuxi, on=["ID", "TIME"])
    global_results = global_results.merge(nonmemrun_aposteriori.nmresult.copy(), on=["ID", "TIME","S2_SAMPLING"])

    # Add comparison columns
    global_results = compare_results(global_results, "_nonmem", "_tucuxi", list(dict_parameters.values()))
    global_results = compare_results(global_results, "_true", "_nonmem", list(dict_parameters.values()))
    global_results = compare_results(global_results, "_true", "_tucuxi", list(dict_parameters.values()))


    global_results.to_csv("results_"+model_nonmem+".csv", sep=",", index=False)
    create_graph_resultats(global_results, folder_name, list(dict_parameters.values()))





def compare_results(data: pd.DataFrame, soft1:str, soft2:str, list_param:[]):

    for i in range(len(list_param)):
        param_soft1 = list_param[i] + soft1
        if soft1 == "_true" and list_param[i]=="IPRED":
            param_soft1 = "DV" + soft1

        param_soft2 = list_param[i] + soft2
        if soft2 == "_true" and list_param[i] == "IPRED":
            param_soft2 = "DV" + soft2

        if param_soft1 in data.columns and param_soft2 in data.columns:
            column_name = list_param[i] + soft1 + soft2

            data[column_name] = (data[param_soft2] - data[param_soft1]) / data[param_soft1] * 100

            grouped = data.groupby('S2_SAMPLING')

            # MPE_comparison
            data["MPE_comparison" + column_name] = grouped[param_soft2].transform(
                lambda group: np.log(group) - np.log(data.loc[group.index, param_soft1])
            )

            # Calcul MPE_moyenne_IC
            mpe_moyenne_ic = grouped['MPE_comparison' + column_name].apply(lambda x:
                f"{(np.exp(x.mean()) - 1) * 100} [ {((np.exp(x.mean() - 1.96 * x.std()/np.sqrt(len(x))) - 1) * 100)}, "
                f"{((np.exp(x.mean() + 1.96 * x.std()/np.sqrt(len(x)))- 1) * 100)}]"
            )

            # Calcul MSE_comparison
            data['MSE_comparison' + column_name] = data['MPE_comparison' + column_name] ** 2

            # Calcul RMSE_moyenne_IC
            rmse_moyenne_ic = grouped['MSE_comparison' + column_name].apply(lambda x:
                f"{(np.exp(np.sqrt(x.mean())) - 1) * 100} [ {(((1 / np.exp(np.sqrt(x.mean()))) - 1) * 100)}, "
                f"{((np.exp(np.sqrt(x.mean())) - 1) * 100)}]"
            )


            new_columns = pd.DataFrame({
                'MPE_moyenne_IC' + column_name: data['S2_SAMPLING'].map(mpe_moyenne_ic),
                'RMSE_moyenne_IC_' + column_name: data['S2_SAMPLING'].map(rmse_moyenne_ic)
            })

            data = pd.concat([data, new_columns], axis=1)
    return data




def delete_files_and_subdirectories(directory_path):
   try:
     with os.scandir(directory_path) as entries:
       for entry in entries:
         if entry.is_file():
            os.unlink(entry.path)
         else:
            shutil.rmtree(entry.path)
     print("All files and subdirectories deleted successfully.")
   except OSError:
     print("Error occurred while deleting files and subdirectories.")

   os.rmdir(directory_path)



"""
*************************************** Graph **************************************************************************
"""


def create_graph_resultats(pd_results: pandas.DataFrame, folder_name: str, list_param:[]):

    df = pd_results.drop_duplicates(subset=["ID"], keep='first')


    name_folder = ''
    if sys.platform == "linux" or sys.platform == "linux2" or sys.platform == "darwin":
        name_folder = folder_name + "/"
    elif sys.platform == "win32":
        name_folder = folder_name + '\\'

    for i in range(len(list_param)):

        if (list_param[i]+ "_nonmem") in pd_results.columns and (list_param[i]+ "_tucuxi") in pd_results.columns :
            sns.lmplot(x=list_param[i]+ "_nonmem", y=list_param[i]+"_tucuxi", data = df, hue= "S2_SAMPLING", fit_reg=True,
                       ci=None)
            x = [df[list_param[i]+ "_nonmem"].min(), df[list_param[i]+ "_nonmem"].max()]
            y = [df[list_param[i]+ "_nonmem"].min(), df[list_param[i]+ "_nonmem"].max()]
            plt.plot(x, y, linestyle='--', color='gray')

            plt.grid(True)
            plt.savefig(name_folder + list_param[i] + "_nonmem_tucuxi")
            plt.close()

            sns.lmplot(x=list_param[i] + "_true", y=list_param[i] + "_nonmem", data=df, hue="S2_SAMPLING", fit_reg=True,
                       ci=None)
            x = [df[list_param[i] + "_nonmem"].min(), df[list_param[i] + "_nonmem"].max()]
            y = [df[list_param[i] + "_nonmem"].min(), df[list_param[i] + "_nonmem"].max()]
            plt.plot(x, y, linestyle='--', color='gray')

            plt.grid(True)
            plt.savefig(name_folder + list_param[i] + "_true_nonmem")
            plt.close()

            sns.lmplot(x=list_param[i] + "_true", y=list_param[i] + "_tucuxi", data=df, hue="S2_SAMPLING", fit_reg=True,
                       ci = None)
            x = [df[list_param[i] + "_tucuxi"].min(), df[list_param[i] + "_tucuxi"].max()]
            y = [df[list_param[i] + "_tucuxi"].min(), df[list_param[i] + "_tucuxi"].max()]
            plt.plot(x, y, linestyle='--', color='gray')

            plt.grid(True)
            plt.savefig(name_folder + list_param[i] + "_true_tucuxi")
            plt.close()





"""
*************************************** Main ***************************************************************************
"""


def main():
    args = parse_the_args()

    config_ini_file(args)

    global mgr

    global dict_parameters

    global dict_model

    global time_tucuxi



    """
    pop = Population.virtualdrug_from_nb_patients(nb_patients=400, nb_samples=4,
                                                  study_design= {"10": {"amount": 10,
                                                                "nb_doses": 1,
                                                                "sample_times": [[0.5,1.5], [3,5], [6,10], [18, 30]]},
                                                                 "11": {"amount": 30,
                                                                       "nb_doses": 1,
                                                                       "sample_times": [[0.5, 1.5], [3, 5], [6, 10],
                                                                                        [18, 30]]},
                                                                 "12": {"amount": 60,
                                                                       "nb_doses": 1,
                                                                       "sample_times": [[0.5, 1.5], [3, 5], [6, 10],
                                                                                        [18, 30]]},
                                                                 "13": {"amount": 80,
                                                                       "nb_doses": 1,
                                                                       "sample_times": [[0.5, 1.5], [3, 5], [6, 10],
                                                                                        [18, 30]]},
                                                                 "14": {"amount": 120,
                                                                       "nb_doses": 1,
                                                                       "sample_times": [[0.5, 1.5], [3, 5], [6, 10],
                                                                                        [18, 30]]},
                                                                 "20": {"amount": 10,
                                                                        "nb_doses": 1,
                                                                        "sample_times": [[14, 34]]},
                                                                 "21": {"amount": 30,
                                                                       "nb_doses": 1,
                                                                       "sample_times": [[14, 34]]},
                                                                 "22": {"amount": 60,
                                                                       "nb_doses": 1,
                                                                       "sample_times": [[14, 34]]},
                                                                 "23": {"amount": 80,
                                                                       "nb_doses": 1,
                                                                       "sample_times": [[14, 34]]},
                                                                 "24": {"amount": 120,
                                                                       "nb_doses": 1,
                                                                       "sample_times": [[14, 34]]},

                                                                 "30": {"amount": 10,
                                                                        "nb_doses": 10,
                                                                        "sample_times": [[214.5, 215.5], [216, 218],
                                                                                         [219, 221], [222,226]]},
                                                                 "31": {"amount": 30,
                                                                       "nb_doses": 10,
                                                                       "sample_times": [[214.5, 215.5], [216, 218],
                                                                                        [219, 221], [222, 226]]},
                                                                 "32": {"amount": 60,
                                                                       "nb_doses": 10,
                                                                       "sample_times": [[214.5, 215.5], [216, 218],
                                                                                        [219, 221], [222, 226]]},
                                                                 "33": {"amount": 80,
                                                                       "nb_doses": 10,
                                                                       "sample_times": [[214.5, 215.5], [216, 218],
                                                                                        [219, 221], [222, 226]]},
                                                                 "34": {"amount": 120,
                                                                       "nb_doses": 10,
                                                                       "sample_times": [[214.5, 215.5], [216, 218],
                                                                                        [219, 221], [222, 226]]},
                                                                 "40": {"amount": 10,
                                                                        "nb_doses": 10,
                                                                        "sample_times":[[230, 250]]},
                                                                 "41": {"amount": 30,
                                                                       "nb_doses": 10,
                                                                       "sample_times": [[230, 250]]},
                                                                 "42": {"amount": 60,
                                                                       "nb_doses": 10,
                                                                       "sample_times": [[230, 250]]},
                                                                 "43": {"amount": 80,
                                                                       "nb_doses": 10,
                                                                       "sample_times": [[230, 250]]},
                                                                 "44": {"amount": 120,
                                                                       "nb_doses": 10,
                                                                       "sample_times": [[230, 250]]}})
    """
    name_file = ''
    if sys.platform == "linux" or sys.platform == "linux2" or sys.platform == "darwin":
        name_file = 'run_{time}/'.format(time=datetime.now().strftime('%Y.%m.%d_%H.%M.%S'))
    elif sys.platform == "win32":
        name_file = 'run_{time}\\'.format(time=datetime.now().strftime('%Y.%m.%d_%H.%M.%S'))


    foldername_init = os.path.join(configValues['output'],
                                   name_file)


    print(Fore.CYAN + foldername_init)
    if not os.path.isdir(foldername_init):
        try:
            subprocess.run('mkdir ' + foldername_init, shell=True)
        except CalledProcessError as err:
            print(Fore.RED + err.output)
            sys.exit()


    list_model_id =  ["210000", "220000", "230000",
                     "210100", "220100", "210200", "220200", "210300", "220300",
                     "210400", "220400", "210010", "210020", "210009",
                     "210001", "210002", "210003", "210004", "210007", "210008",
                     "220001", "220002", "220003", "220004", "220007","220009",
                     "230001", "230002", "230003", "230004", "230007", "230009",

                      "110000", "120000", "130000", "111000", "112000", "123000", "124000", "125000",
                       "126000", "127000", "128000", "110100", "120100", "110200", "120200",
                       "110300", "120300", "110400", "120400", "110010", "110020",
                       "110001", "110002", "110003", "110004", "110005", "110006", "110007", "110008",
                       "110009",
                       "120001", "120002", "120003", "120004", "120007", "120009",
                       "130001", "130002", "130003", "130004", "130007", "130009",

                      "310000", "320000", "330000", "310100", "320100", "310200", "320200", "310300", "320300",
                        "310400", "320400", "310010", "310020", "310009",
                        "310001", "310002", "310003", "310004", "310007", "310008",
                        "320001", "320002", "320003", "320004", "320007", "320009",
                        "330001", "330002", "330003", "330004", "330007", "330009"]


    for model_id in list_model_id:
        pop = Population.virtualdrug_from_csv_file(configValues['inputfilename'],
                                                   model_id=model_id)
        run_numerical_validation(pop, model_id, foldername_init, args)






if __name__ == "__main__":
    main()