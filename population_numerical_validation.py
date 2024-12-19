import pandas

from colorama import init
from dataclasses import dataclass, field

# with sotalya
from sotalya.data.computingqueryresponse import *
from sotalya.processing.querytopendingrequest import *
from sotalya.data.query import *
from sotalya.data.requests import *

# with tucucli
#from tucuxi.data.computingqueryresponse import *
#from tucuxi.processing.querytopendingrequest import *
#from tucuxi.data.query import *
#from tucuxi.data.requests import *


from math import *
from core.extractxml import *
import re

# reset all
init(autoreset=True)


"""
*************************************** Useful functions ***************************************************************
"""


def fact():
    """
    * function returning an empty list

        Args:

        * no Args
    """
    return []


def unis_str_to_datetime(string: str):
    """
    * This function returns a datetime in the right format.

        Args:

        * string (str): a datetime as a string.
    """
    return datetime.strptime(string, '%m/%d/%YT%H:%M:%S')


def create_excel_results(dict_data: list, folder_name: str, file_name: str, nb_samples: int, list_sample_times: list):
    """
    * Once a population is created, this will create an Excel to be able to use those results further on.
    * It helps to be able to use the old results again in the code, to try new doses for example.

        Args:

        * dict_data (list): list containing the population data
        * folder_name (str): the name of the folder where we want the population to be stored
        * file_name (str): the name of the file where we want the population to be stored
        * nb_samples (int): number of sample
        * list_sample_times (list): list of sampling times
    """
    if len(dict_data) == 0:
        return False

    data = pandas.DataFrame(columns=dict_data[0].keys())

    # saving each patient in 1 line, the columns being the co-variates, dates, and other relevant information
    # from the population class
    i = 0
    for d in dict_data:
        list_data = []
        for key, value in d.items():
            list_data.append(value)

        data.loc[i] = list_data
        i += 1
    data["nb_samples"] = nb_samples
    data["list_sample_times"] = [list_sample_times] * len(data.index)

    data.to_csv(folder_name + file_name, index=False)
    return data


"""
*************************************** Day and Date *******************************************************************
"""


@dataclass
class Day:
    """
    * The class contains all information useful for 1 patient, in 1 Episode / cycle, and in 1 day of this cycle:
        The dose, number of doses given that day, the weight, the start date and end date, the dose interval
        and the infusion's time.
    * This class has 1 child class for each type of Day we want to have, depending on the study and covariates used:
        ** DayNumericalValidation
    """
    dose: str = ''
    number_of_doses: int = 0
    weight: str = ''
    start_datetime: datetime = datetime(1900, 1, 1, 0, 0)
    end_datetime: datetime = datetime(1900, 1, 1, 0, 0)
    interval: str = ''
    infusion_time: str = ''


@dataclass(init=False)
class DayNumericalValidation(Day):
    """
    * The class contains all information useful for 1 patient, in 1 Episode / cycle, and in 1 day of this cycle:
        The dose, number of doses given that day, the weight and creatinine measured, as well as the post-natal age
        (pna), the start date and end date, the dose interval and the infusion's time.
    * This Day creation is related to the numerical validation between tucuxi and nonmem using ch.tucuxi.virtualdrug.modX.
    * This class has 1 parent class: Day.
    """
    def __init__(self, dose: str = '', number_of_doses: int = 0, weight: str = '', is_male:str ='',
                 start_datetime: datetime = datetime(1900, 1, 1, 0, 0),
                 end_datetime: datetime = datetime(1900, 1, 1, 0, 0),
                 interval: str = '', infusion_time: str = ''):
        """
        * Initialisation of the DayNumericalValidation class.
        * This class needs access to the dose, number of doses given that day, the weight measured,
            as well as the sex, the start date and end date, the dose interval and the infusion's time.
        """
        super().__init__(dose,  number_of_doses, weight, start_datetime, end_datetime, interval, infusion_time)
        self.is_male = is_male


    def create_dict(self):
        """
        * Creating a dictionary based on the results in this class.
        * It is used when want to print the output while the code is running.

            Args:

            * no Args
        """
        dict_dose = {"bodyweight": self.weight,
                     "isMale": self.is_male,
                     "amount": self.dose,
                     "infusion_time": self.infusion_time,
                     "interval": self.interval
                     }
        return dict_dose



"""
*************************************** Dose & Concentrations **********************************************************
"""


@dataclass(init=False)
class DosePatient:
    """
    * This class contains the information of when the treatment started, the amount of drug a patient was given,
        the number of times the patient was given the dose and at which interval, including the infusion time of the
        given drug.
    * Of course, if there is no infusion time, it can be set to 0.
    """
    start_date: datetime = datetime(1900, 1, 1, 0, 0, 0)
    amt: float = 0.0
    nb_doses: int = 0
    infusion_time: float = 0.0
    interval: float = 0.0

    def __init__(self, start_date, amt, nb_doses, infusion_time, interval):
        """
        * Initialisation of the DosePatient class.
        * This class needs access to the start date, amount of drug per dose and number of doses to give,
            the infusion time, and the interval between each dose.
        """
        self.start_date = start_date
        self.amt = amt
        self.nb_doses = nb_doses
        self.infusion_time = infusion_time
        self.rate = self.amt / self.infusion_time
        self.interval = interval


    @classmethod
    def numerical_validation_initial_dose(cls, nb_doses: int, amt:float):
        """
        * This is a Constructor of DosePatient.
        * This returns the initial dose given to a patient
        * This dose calculation is related to the ch.tucuxi.virtualdrug study.

            Args:

            * nb_doses: int, the number of doses intake
        """
        start_date = datetime(2023, 1, 1, 8, 00, 00)
        infusion_time = 1

        interval = 24

        return cls(start_date=start_date, amt=amt, nb_doses=nb_doses, infusion_time=infusion_time, interval=interval)


@dataclass(init=False)
class SamplePatient:
    """
    * This class represents a sampling time for 1 patient at 1 specific time.
    * It will therefore contain a datetime at which this sampling was done, as well as the concentration / value
        measured at that time.
    * e.g. 40mg/L (value) of VirtualDrug were measured the 14.02.2022 (sample_date)
    * The information needed are initialized to 01.01.1900 and an empty string, if none is given.
    """
    sample_date: datetime = datetime(1900, 1, 1)
    value: str = ''

    def __init__(self, sample_date=datetime(1900, 1, 1), value=""):
        """
        * Initialisation of the SamplePatient class.
        * This class needs access to a date at which a sampling was done, as well as the value measured while doing
            the sampling.
        * The information needed are initialized to 01.01.1900 and an empty string, if none is given.
        """
        self.sample_date = sample_date
        self.value = value

    @classmethod
    def random_sample(cls, day_sample: datetime, min_hour: float, max_hour: float):
        """
        * initialization of a sample taken at a random time, between 'min_hour' and 'max_hour'
        """
        random_sample_hour = np.round(np.random.uniform(min_hour, max_hour, 1), 2)
        time_sample = timedelta(hours=random_sample_hour[0])
        sample_date = day_sample + time_sample
        return cls(sample_date=sample_date, value='')

    def create_dict(self):
        """
        * Creating a dictionary based on the results in this class.
        * It is used when want to print the output while the code is running.

            Args:

            * no Args
        """
        dict_sample = {"sample_date": self.sample_date,
                       "sample_value": self.value}
        return dict_sample


@dataclass(init=False)
class DosageTime:
    """
    * This class represents the time for which a special drug dosage is given.
    * It contains only the start and the end date for the specific dosage, and the dosage itself.
    * The dosage can be of the class DosageRepeat, DosageSequence, LastingDosage, DailyDosage or WeeklyDosage.
    """
    start: datetime
    end: datetime

    def __init__(self, start: datetime = str_to_datetime("1111-11-11T11:11:11"), dosage=None,
                 end: datetime = str_to_datetime("1111-11-11T11:11:11")):
        """
        * Initialisation of the DosageTime class.
        * This class needs access to the start and the end date for the specific dosage, and the dosage itself
        * The information needed are initialized, if nothing is given.
        """
        self.start = start
        self.end = end
        self.dosage = dosage

    @classmethod
    def copy_dosage_time(cls, d=None):
        """
        * This is a Constructor of DosageTime.
        * This returns a new instance of the class DosageTime based on a variable given that contains all the
            information needed to create one.
        * The variable given can be also of the class DosageTime. Then, it will simply do a copy of this one.
        * If no variable is given, it will act as the __init__() of the class.

            Args:

            * d (unknown): a class containing a variable "start", "end" and "dosage". Initialized as None.
        """
        start = str_to_datetime("1111-11-11T11:11:11")
        end = str_to_datetime("1111-11-11T11:11:11")
        dosage = None

        if d is not None:
            start = str_to_datetime(d.start.string)
            end = str_to_datetime(d.end.string)

            if d.dosage.dosageLoop:
                dosage = choose_dosage(d.dosage.dosageLoop)
            else:
                dosage = choose_dosage(d.dosage)

        return cls(start=start, dosage=dosage, end=end)

    def is_valid(self):
        """
        * This function checks if the start and end date are considered valid or not.
        * If they are not of the type datetime, it returns an error.

            Args:

            * no Args
        """
        if not (type(self.start) is datetime):
            print(Fore.RED + "Dosage start date invalid")
            return False
        if not (type(self.end) is datetime):
            print(Fore.RED + "Dosage start date invalid")
            return False
        return True


"""
*************************************** Results and Comparison *********************************************************
"""


@dataclass
class ResultsAposteriori:
    """
    * These results are the assessment of the patient's current dosing regimen.
    * This class registers tdm_id, the start and end date of the cycle, AUC24, AUC, trough, if the target was achieved
        or not, as well as the graph values.
    """
    tdm_id: str
    date_start: datetime
    date_end: datetime
    auc24: float
    auc: float
    trough: float
    target_achieved: bool
    list_graph_values: []

    def create_dict(self):
        """
        * Creating a dictionary based on the results in this class.
        * It is used when want to print the output while the code is running.

            Args:

            * no Args
        """
        dict_results = {"TDM_ID": self.tdm_id,
                        "AUC0-24": self.auc24,
                        "Cmin": self.trough,
                        "target_achieved": self.target_achieved}
        return dict_results


@dataclass
class ResultsAdj:
    """
    * These are the results that tucuxi offers in adjustment, with the new dosage to be taken.
    * This class registers the tdm_id, the date of the adjustment of the dose, the new dose, its interval, infusion
        time, as well as the score and the target for the TDM in question and the list of the values for the graph.
    """
    tdm_id: str
    date_first_dose_adj: datetime
    dose: float
    interval: float
    infusion: float
    score: float
    target: float
    list_graph_values: []


@dataclass
class CompareDose:
    """
    * This class helps to find out if the proposed dose is the same as the dose already received.
    * It contains the tdm id, the dose difference in mg, the dose difference relative, the interval difference in hours,
        the interval difference relative to one another, and the infusion difference, relative and in hours.
    """
    tdm_id: str = ''
    dose_diff_inmg: float = 0.0
    dose_diff_relative: float = 0.0
    interval_diff_inhours: float = 0.0
    interval_diff_relative: float = 0.0
    infusion_diff_inhours: float = 0.0
    infusion_diff_relative: float = 0.0


@dataclass(init=False)
class PKParameters:
    """
    * This class generates a new set of PK parameters for a patient.
    * PK parameters: Clearance and Volume.
    * How those parameters are measured: Nonmem and Tucuxi.
    """
    CL_dict = {}
    V_dict = {}

    def __init__(self, cl_nonmem, cl_tucuxi, v_nonmem, v_tucuxi):
        """
        * Initialisation of the PKParameters class.
        * This class needs access to CL and V parameters, both for nonmem and tucuxi.
        """
        self.CL_dict = {"Nonmem_CL": cl_nonmem, "Tucuxi_CL": cl_tucuxi}
        self.V_dict = {"Nonmem_V": v_nonmem, "Tucuxi_V": v_tucuxi}


"""
*************************************** Cycle **************************************************************************
"""


@dataclass(init=False)
class Episode:
    """
    * This class generates a cycle for a patient.
    * A cycle is a period of time for which the patient is treated a special dose and has a special weight
    * This class has 1 child class for each type of Day we want to have, depending on the study and covariates used:
        ** EpisodeVirtualDrug
    """
    episode_id: int = 0
    birthdate: datetime = datetime(1900, 1, 1, 0, 0, 0)
    etacl: float = 0.0
    samples: List[SamplePatient] = field(default_factory=fact)
    days: [] = field(default_factory=fact)
    startonlydate: datetime = datetime(1900, 1, 1)

    def __init__(self, episode_id=0, birthdate=datetime(1900, 1, 1, 0, 0, 0), etacl=0.0, samples=None, days=None,
                 startonlydate=datetime(1900, 1, 1)):
        """
        * Initialisation of the Episode class.
        * This class needs access to the episode id, the birthdate, the inter variability of the clearance (etacl),
            the list of the sampling time (samples), the list of the information for the patient for 1 cycle (day),
            and the start only date.
        """
        if samples is None:
            samples = fact()
        if days is None:
            days = fact()
        self.episode_id = episode_id
        self.birthdate = birthdate
        self.etacl = etacl
        self.samples = samples
        self.days = days
        self.startonlydate = startonlydate


@dataclass(init=False)
class EpisodeVirtualDrug(Episode):
    """
    * This class generates a cycle for a patient.
    * A cycle is a period of time for which the patient is treated a special dose and has a special weight
    * This Episode creation is related to the numerical validation between tucuxi and nonmem using ch.tucuxi.virtualdrug.modX
    * This class has 1 parent class: Episode.
    """
    def __init__(self, episode_id=0, birthdate=datetime(1900, 1, 1, 0, 0, 0),
                 is_male=0, sampling_group= 1, etacl=0.0, samples=None, days=None, startonlydate=datetime(1900, 1, 1)):
        """
        * This class needs access to the episode id, the sex of the patient, the sampling group, the inter
            variability of the clearance (etacl), the list of the sampling time (samples), the list of the information
            for the patient for 1 cycle (day), and the start only date.
        """
        if samples is None:
            samples = fact()
        if days is None:
            days = fact()
        super().__init__(episode_id, birthdate, etacl, samples, days, startonlydate)
        self.is_male = is_male
        self.sampling_group = sampling_group

    def create_dict(self):
        """
        * Creating a dictionary based on the results in this class.
        * It is used when want to print the output while the code is running.

            Args:

            * no Args
        """
        dict_episode = {"isMale": self.is_male,
                        "etaCL": self.etacl}
        dict_episode.update(self.days[-1].create_dict())
        dict_episode.update(self.samples[-1].create_dict())
        return dict_episode


"""
*************************************** Patient ************************************************************************
"""


def create_dosage_history(history: DosageHistory, episodes: []) -> datetime:
    """
    * Function checking what is the last dosage date and then manage the tqf requests.

        Args:

        * history (DosageHistory): dosage history
        * episodes (List[Episode]): list of the episodes for 1 patient
    """
    very_end_date = episodes[0].days[0].start_datetime

    for episode in episodes:

        for day in episode.days:
            lasting_dosage = LastingDosage()
            lasting_dosage.dose.value = day.dose
            lasting_dosage.dose.infusionTimeInMinutes =  timedelta(minutes=(float(day.infusion_time) * 60.0))
            lasting_dosage.dose.unit = 'mg'
            lasting_dosage.formulationAndRoute.absorptionModel = 'infusion'
            lasting_dosage.formulationAndRoute.formulation = 'parenteralSolution'
            lasting_dosage.formulationAndRoute.administrationRoute = 'intravenousDrip'
            lasting_dosage.formulationAndRoute.administrationName = 'foo bar'
            lasting_dosage.interval = timedelta(hours=float(day.interval))

            start_date = day.start_datetime
            end_date = day.end_datetime

            if end_date > very_end_date:
                very_end_date = end_date

            if start_date + lasting_dosage.interval * int(day.number_of_doses) > very_end_date:
                very_end_date = start_date + lasting_dosage.interval * int(day.number_of_doses)

            dosage_time_range = DosageTime(start_date, lasting_dosage, end_date)

            history.dosageTimeRanges.append(dosage_time_range)

    return very_end_date



def create_dosage_history_peros(history: DosageHistory, episodes: [], abs_model: str) -> datetime:
    """
    * Function checking what is the last dosage date and then manage the tqf requests.

        Args:

        * history (DosageHistory): dosage history
        * episodes (List[Episode]): list of the episodes for 1 patient
        * abs_model: to indicate if absorptionModel = extra or extra lag
    """
    very_end_date = episodes[0].days[0].start_datetime

    for episode in episodes:

        for day in episode.days:
            lasting_dosage = LastingDosage()
            lasting_dosage.dose.value = day.dose
            lasting_dosage.dose.infusionTimeInMinutes = timedelta(minutes=(float(day.infusion_time) * 60.0))
            lasting_dosage.dose.unit = 'mg'
            lasting_dosage.formulationAndRoute.absorptionModel = abs_model
            lasting_dosage.formulationAndRoute.formulation = 'oralSolution'
            lasting_dosage.formulationAndRoute.administrationRoute = 'oral'
            lasting_dosage.formulationAndRoute.administrationName = 'foo bar'
            lasting_dosage.interval = timedelta(hours=float(day.interval))

            start_date = day.start_datetime
            end_date = day.end_datetime

            if end_date > very_end_date:
                very_end_date = end_date

            if start_date + lasting_dosage.interval * int(day.number_of_doses) > very_end_date:
                very_end_date = start_date + lasting_dosage.interval * int(day.number_of_doses)

            dosage_time_range = DosageTime(start_date, lasting_dosage, end_date)

            history.dosageTimeRanges.append(dosage_time_range)

    return very_end_date

def create_dosage_history_bolus(history: DosageHistory, episodes: []) -> datetime:
    """
    * Function checking what is the last dosage date and then manage the tqf requests.

        Args:

        * history (DosageHistory): dosage history
        * episodes (List[Episode]): list of the episodes for 1 patient
    """
    very_end_date = episodes[0].days[0].start_datetime

    for episode in episodes:

        for day in episode.days:
            lasting_dosage = LastingDosage()
            lasting_dosage.dose.value = day.dose
            lasting_dosage.dose.infusionTimeInMinutes =  timedelta(minutes=(0.0 * 60.0))
            lasting_dosage.dose.unit = 'mg'
            lasting_dosage.formulationAndRoute.absorptionModel = 'bolus'
            lasting_dosage.formulationAndRoute.formulation = 'parenteralSolution'
            lasting_dosage.formulationAndRoute.administrationRoute = 'intravenousBolus'
            lasting_dosage.formulationAndRoute.administrationName = 'foo bar'
            lasting_dosage.interval = timedelta(hours=float(day.interval))

            start_date = day.start_datetime
            end_date = day.end_datetime

            if end_date > very_end_date:
                very_end_date = end_date

            if start_date + lasting_dosage.interval * int(day.number_of_doses) > very_end_date:
                very_end_date = start_date + lasting_dosage.interval * int(day.number_of_doses)

            dosage_time_range = DosageTime(start_date, lasting_dosage, end_date)

            history.dosageTimeRanges.append(dosage_time_range)

    return very_end_date


def create_covariates_virtualdrug(query: Query, episodes: List[EpisodeVirtualDrug], model_name:str):
    """
    * This function adds the episodes and their information to the given query, thus writing the tqf file.
    * The function is related to the ch.tucuxi.virtualdrug.modX study.

        Args:

        * query (Query): enter a query.
        * episodes (List[EpisodeVirtualDrug]): list of episodes of only 1 person.
            This argument is like a patient.
    """
    first_episode = episodes[0]
    episode_start_date = episodes[0].startonlydate

    query.covariates.append(
        Covariate.create_covariate('birthdate', episode_start_date,
                                   first_episode.birthdate.strftime('%Y-%m-%dT%H:%M:%S'), '', 'date', 'discrete'))

    query.covariates.append(
        Covariate.create_covariate('sampling_group', episode_start_date,
                                   first_episode.sampling_group, '', 'int', 'discrete'))

    if model_name=="ch.tucuxi.virtualdrug.mod110009" or model_name=="ch.tucuxi.virtualdrug.mod210009" :
        for episode in episodes:

            for d in range(len(episode.days)):
                dose = episode.days[d]

                if d == 0 or (d != 0 and dose.start_datetime.date() != episode.days[d - 1].start_datetime.date()):
                    the_date = dose.start_datetime

                    if dose.weight != '':
                        bodyweight = str(float(dose.weight))
                        query.covariates.append(Covariate.create_covariate('bodyweight', the_date, bodyweight, 'kg',
                                                                           'double', 'continuous'))

                    if dose.is_male != '':
                        is_male = "0"
                        if dose.is_male == "True" :
                            is_male = "1"
                        query.covariates.append(Covariate.create_covariate('sex', the_date, is_male,
                                                                           '', 'double', 'continuous'))




@dataclass(init=False)
class Patient:
    """
    * This class generates a patient with an ID, a bodyweight, a birthdate, a parameter of inter variability on the
        clearance, a number of TDM already done, a dose and sample.
    * This class has 1 child class for each type of Patient we want to have, depending on the study and covariates used:
        ** PatientVirtualDrug
    """
    patient_id: int
    bw: float
    birthdate: datetime
    etacl: float
    nb_tdm: int
    dose: DosePatient
    sample: List[SamplePatient]

    def __init__(self, patient_id, bw, birthdate, etacl, dose, sample):
        """
        * Initialisation of the Patient class.
        * This class needs access to an ID, a bodyweight, a birthdate, a parameter of inter variability on the
            clearance, a number of TDM already done, a dose and sample.
        """
        self.patient_id = patient_id
        self.bw = bw
        self.birthdate = birthdate
        self.etacl = etacl
        self.nb_tdm = 0
        self.dose = dose
        self.sample = []

        if isinstance(sample, list):
            self.sample += sample
        else:
            self.sample.append(sample)


@dataclass(init=False)
class PatientVirtualDrug(Patient):
    """
    * This class generates a patient with an ID, a bodyweight, a birthdate, a parameter of inter variability on the
        clearance, a number of TDM already done, a dose and sample, as well as sex.
    * This patient creation is related to the ch.tucuxi.virtualdrug study.
    * This class has 1 parent class: Patient.
    """
    def __init__(self, patient_id, birthdate, bw, is_male,  sampling_group, dose, sample):
        """
        * Initialisation of the PatientVirtualDrug class.
        * This class needs access to ID, a bodyweight, a birthdate, a parameter of inter variability on the
            clearance, a number of TDM already done, a dose and sample, as well as birthweight.
        """

        super().__init__(patient_id, bw, birthdate, 0, dose, sample)
        self.is_male = is_male
        self.sampling_group = sampling_group

    @classmethod
    def from_random(cls, patient_id, is_male, bw, sampling_group, list_sample_times, amount, nb_doses):
        """
        * This is a Constructor of PatientVirtualDrug
        * This constructor creates a new Patient based on the desired input.
        * It automaticaly creates the bodyweight, the sex, the dose and the sample.

            Args:

            * id (int): patient's id.
            * is_male (bool): if patient is a male -> 1, if not -> 0.
            * bw : bodyweight of the patient
            * etacl (float): the parameter of inter variability on the clearance.
            * list_sample_times (list) : time of the sample ("random" or "residual" or float)
        """
        nb_samples = len(list_sample_times)

        # birthdate
        birthdate = datetime(2023, 1, 1, 0, 0, 0)


        # measure initial dose that the patient will receive
        dose_patient = DosePatient.numerical_validation_initial_dose(nb_doses = nb_doses, amt=amount)

        # samples
        list_samples = []

        for s in range(nb_samples):
            if list_sample_times[s] == "residual":
                time_sample = timedelta(hours=dose_patient.interval) - timedelta(minutes=5)
                sample_date = dose_patient.start_date + time_sample
                sample_patient = SamplePatient(sample_date=sample_date, value='')

            elif (list_sample_times[s] == "random") | isinstance(list_sample_times[s], list):
                if isinstance(list_sample_times[s], list):
                    sample_patient = SamplePatient.random_sample(day_sample=dose_patient.start_date,
                                                                 min_hour=list_sample_times[s][0],
                                                                 max_hour=list_sample_times[s][1])
                else:
                    sample_patient = SamplePatient.random_sample(day_sample=dose_patient.start_date,
                                                                 min_hour=0,
                                                                 max_hour=dose_patient.interval)

            else:
                time_sample = timedelta(hours=list_sample_times[s])
                sample_date = dose_patient.start_date + time_sample
                sample_patient = SamplePatient(sample_date=sample_date, value='')

            list_samples.append(sample_patient)


        return cls(patient_id=patient_id, birthdate=birthdate, bw=bw, is_male=is_male,
                   sampling_group= sampling_group,
                    dose=dose_patient, sample=list_samples)

    def create_dict_patient(self, nb_samples_max=int) -> dict:
        dict_patient = {"patient_id": self.patient_id,
                        "sampling_group": self.sampling_group,
                        "birthdate": self.birthdate,
                        "bodyweight": self.bw,
                        "isMale": self.is_male,
                        "etacl": self.etacl,
                        "start_date": self.dose.start_date,
                        "amt": self.dose.amt,
                        "infusion_time": self.dose.infusion_time,
                        "interval": self.dose.interval,
                        "nb_doses": self.dose.nb_doses,
                        "rate": self.dose.rate}

        for i in range(0, len(self.sample)):
            dict_patient["sample_date_{key}".format(key=i)] = self.sample[i].sample_date

        for i in range(len(self.sample), nb_samples_max):
            dict_patient["sample_date_{key}".format(key=i)] = datetime(1900,1,1,0,0)


        return dict_patient


@dataclass(init=False)
class ExpPatient:
    """
    * This class generates a patient with an ID, the list of pk_parameters for the patient,
        as well as the list of Cycles he/she goes through (Cycle = Episode), the list of
        TDM results it has generated previously, the list of compare_dosage and
        the list of results a posteriori.
    """
    tdm: []  # List[Episode]
    results: List[ResultsAdj]
    compare_dosage: List[CompareDose]
    results_Aposteriori: List[ResultsAposteriori]
    PK_Parameters: List[PKParameters]
    patient_id: str = ''

    def __init__(self):
        """
        * Initialisation of the ExpPatient class.
        * This class does not need any special access.
        """
        self.results = []
        self.tdm = []
        self.compare_dosage = []
        self.results_Aposteriori = []
        self.PK_Parameters = []

    def update_sampling_patient(self, nmresult: [[str], [float], [float]]):
        """
        * Updating the ExpPatient class with the new values from nonmem.

            Args:

            * nmresult ([float]): The nonmem results ([ID], [Time of samples], [samples]).
        """
        index_last_tdm = len(self.tdm) - 1

        episodes = self.tdm[index_last_tdm].samples

        samples_id = [index for index, value in enumerate(nmresult[0]) if value == self.patient_id]
        if len(samples_id)==len(episodes):
            for i in range(len(episodes)):
                if nmresult[2][samples_id[i]] >= 0:
                    episodes[i].value = nmresult[2][samples_id[i]]
                else:
                    episodes[i].value = 0

        self.tdm[index_last_tdm].samples = episodes


    def create_dict(self):
        """
        * Creating a dictionary based on the results in this class.
        * It is used when want to print the output while the code is running.

            Args:

            * no Args
        """
        dict_patient = {"ID": self.patient_id}
        dict_patient.update(self.tdm[-1].create_dict())
        dict_patient.update(self.results_Aposteriori[-1].create_dict())
        dict_patient.update(self.PK_Parameters[-1].CL_dict)
        dict_patient.update(self.PK_Parameters[-1].V_dict)
        return dict_patient


"""
*************************************** Population *********************************************************************
"""


@dataclass(init=False)
class Population:
    """
    * This class generates a list of patients based on a certain number of patient given upfront.
    """
    nb_patients: int = 0
    population: [] = field(default_factory=fact)  # List[Patient]
    nb_samples: int = 1
    list_sample_times: [] = field(default_factory=fact)

    def __init__(self, nb_patients, population, nb_samples, list_sample_times):
        """
        * Initialisation of the Population class.
        * This class needs access to a number of patient and a list of patients.
        """
        self.nb_patients = nb_patients
        self.population = population
        self.nb_samples = nb_samples
        self.list_sample_times = list_sample_times

    @classmethod
    def virtualdrug_from_nb_patients(cls, nb_patients, nb_samples, study_design:dict):
        """
        * This is a Constructor of Population.
        * This returns a population containing a certain number of patient, as well as a list of patients.
        * This population creation is related to the ch.tucuxi.virtualdrug study.

            Args:

            * nb_patients (int): Number of patients to add. If the number is odd, the number will be
                incremented of 1 to be divisible by 2 so that the population contains as many male
                as female
        """
        # If the number is odd, the number will be incremented of 1 to be divisible by 2 so that the population
        # contains as many male as female
        if (nb_patients % len(study_design)) != 0:
            nb_patients += nb_patients % len(study_design)


        # CO-VARIATES

        # bodyweight
        bw_pop = np.floor(np.random.uniform(20, 100, nb_patients))

        # study group
        group_pop=[]
        list_sample_times_pop = []
        amount_pop = []
        nb_doses_pop = []
        is_male_pop = []

        for key, value in study_design.items():
            group_pop = group_pop + [key] * round(nb_patients/len(study_design))
            is_male_pop = (is_male_pop + [True] * ceil((nb_patients/len(study_design))/2) +
                           [False] * floor((nb_patients/len(study_design))/2))
            list_sample_times_pop = list_sample_times_pop + [value["sample_times"]] * round(nb_patients/len(study_design))
            amount_pop = amount_pop + [value["amount"]] * round(nb_patients/len(study_design))
            nb_doses_pop = nb_doses_pop + [value["nb_doses"]] * round(nb_patients/len(study_design))

        # future list of patients
        population = []
        # for each patient
        for i in range(nb_patients):

            # create and save the patient
            population.append(PatientVirtualDrug.from_random(patient_id=str(i + 1), is_male=is_male_pop[i],
                                                                   bw=bw_pop[i], sampling_group=group_pop[i],
                                                                   list_sample_times=list_sample_times_pop[i],
                                                             amount=amount_pop[i],
                                                             nb_doses=nb_doses_pop[i]))

        # create the population
        return cls(nb_patients=nb_patients, population=population, nb_samples=nb_samples,
                   list_sample_times=[])

    @classmethod
    def virtualdrug_from_csv_file(cls, filename: str, model_id: str):
        """
                * This is a Constructor of Population.
                * This returns a population containing a certain number of patient, as well as a list of patients.
                * This population creation is related to the ch.tucuxi.virtualdrug study.
                * This population is based on a csv file

                    Args:

                    * filename (str): name of the file containing the population
        """

        df = pandas.read_csv(filename)
        if model_id == "111000":
            df.loc[df['sampling_group']<=15, 'sample_date_0'] = "1900-01-01 00:00:00"
            df.loc[(df['sampling_group'] >= 30) & (df['sampling_group'] <= 35), 'sample_date_1'] = "1900-01-01 00:00:00"
        list_col = df.columns.tolist()

        population = []
        list_sample_times = df["list_sample_times"][0]

        for i in df.index:
            dose_patient = DosePatient(start_date=datetime.strptime(df["start_date"][i], "%Y-%m-%d %H:%M:%S"),
                                       amt=df["amt"][i], nb_doses=int(df["nb_doses"][i]),
                                       interval=int(df["interval"][i]), infusion_time=int(df["infusion_time"][i]))

            sample = []

            list_col_sample = [m.string if m is not None else 'No Match' for m in [re.search("sample_date",
                                                                                             col) for col in list_col]]

            for j in range(0, len(list_col_sample)):
                if list_col_sample[j] != 'No Match':
                    date = datetime.strptime(df[list_col_sample[j]][i], "%Y-%m-%d %H:%M:%S")

                    if date.year != 1900:
                        sample.append(SamplePatient(sample_date=date, value=''))

            population.append(PatientVirtualDrug(patient_id=str(df["patient_id"][i]),
                                                 birthdate=datetime.strptime(df["birthdate"][i],
                                                                                   "%Y-%m-%d"),
                                                 bw=df["bodyweight"][i],
                                                 is_male=df["isMale"][i],
                                                  dose=dose_patient, sample=sample,
                                                 sampling_group=str(df["sampling_group"][i])
                                                ))

        return cls(nb_patients=len(df), population=population, nb_samples=len(list_sample_times),
                   list_sample_times=list_sample_times)

    def save_population_in_file(self, folder_name: str, file_name: str):
        """
                * This function saves the created population in a csv file

                    Args:

                    * folder_name (str): the name of the folder to save the file
                    * file_name (str): the name of the file containing the population
        """
        list_dict_patients = []
        for patient in self.population:
            list_dict_patients.append(patient.create_dict_patient())

        create_excel_results(list_dict_patients, folder_name, file_name, self.nb_samples, self.list_sample_times)

    def virtualdrug_save_population_in_file(self, folder_name: str, file_name: str):
        """
                * This function saves the created population in a csv file

                    Args:

                    * folder_name (str): the name of the folder to save the file
                    * file_name (str): the name of the file containing the population
        """
        list_dict_patients = []
        for patient in self.population:
            list_dict_patients.append(patient.create_dict_patient(self.nb_samples))

        create_excel_results(list_dict_patients, folder_name, file_name, self.nb_samples, self.list_sample_times)


