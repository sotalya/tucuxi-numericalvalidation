

from core.utils import get_output_folder_name
from sotalya.data.plotter import Plotter
from sotalya.data.requests import Request, RequestType, ParametersTypeEnum
from sotalya.data.query import Query
from core.junit import JUnitReport
from colorama import Fore


plotargs = ['bo', 'go', 'po', 'oo', 'co']
postcount = 0
perccount = 0


SoftSelect = dict(both=0, onlytucuxi=1, onlynonmem=2)

generalReport = ''


class Report:

    def __init__(self, generate_graphs, whichsoft):
        global generalReport
        generalReport = self
        self.generateGraphs = generate_graphs
        self.whichsoft = whichsoft
        self.jrep = JUnitReport("testsuite")

    def write_report(self):
        print('Report file : {}'.format(get_output_folder_name() + '/junit_report'))
        self.jrep.write_junit_report(get_output_folder_name() + '/junit_report')

    def start_query(self, query_id):
        self.jrep.add_testsuite(query_id)

    def report_singlecurve(self, tuc_results: [[float], [float]], nm_results: [[float], [float]],
                           request: Request, query: Query):
        """
        * This method reports the results from calculations of Tucuxi and NONMEM with the same data.
        * Its for single curves for apriori and aposteriori. It calls the plotter to make graphs and print
          them to images. It exports xUnit compatible xml for Jenkins
            Args:
                tuc_results ([[],[]]): results from tucucli in the format [[times],[values]]
                nm_results ([[],[]]): results from nonmem in the format [[times],[values]]
                request (Request): the dataset for which to report
                query (Query): The query containing the request
        """

        if self.whichsoft == SoftSelect['both']:
            if request.computingTraits.computingTraits.requestType == RequestType.Prediction:

                title = 'Undefined title'
                if request.computingTraits.computingTraits.computingOption.parametersType \
                        == ParametersTypeEnum.population:
                    title = 'Population'
                elif request.computingTraits.computingTraits.computingOption.parametersType \
                        == ParametersTypeEnum.apriori:
                    title = 'A priori'
                elif request.computingTraits.computingTraits.computingOption.parametersType \
                        == ParametersTypeEnum.aposteriori:
                    title = 'A posteriori'

                if self.generateGraphs:
                    pltr = Plotter(get_output_folder_name(), query.get_id() + request.get_id(), plotargs[postcount])
                    pltr.plot_single_prediction(tuc_results, query.get_id() + request.get_id(), title)
                    pltr.plot_single_prediction_nonmem(nm_results, query.get_id() + request.get_id() + "_nm", title)
                    pltr.plot_c(nm_results, tuc_results, query.get_id() + request.get_id() + "_comp", title)
                list_results = tuc_results[0][1]
                has_passed, failtext = self.compute_error(list_results, nm_results[1], title)

                self.jrep.add_testcase(has_passed, request.requestId, 'prediction', failtext)

            # elif args.whichsoft == SoftSelect['onlynonmem']:
            #     pltr.plot_single_obs_vs_pred(tuc_results, request.get_full_id())
            elif request.computingTraits.computingTraits.requestType == RequestType.PredictionAtSampleTime:
                tuc_results[0][0] = nm_results[0]

                if self.generateGraphs:

                    pltr = Plotter(get_output_folder_name(), query.get_id() + request.get_id(), plotargs[postcount])
                    pltr.plot_obs_vs_pred(nm_results, tuc_results, query.get_id() + request.get_id())

                    samples = [[], []]
                    for i in range(0, len(query.drugs[0].samples)):
                        samples[0].append(nm_results[0][i])
                        samples[1].append(float(query.drugs[0].samples[i].concentration) * 1000.0)
                    # pltr.plot_obs_vs_pred_with_samples(nm_results, tuc_results, samples, request.get_full_id())
                    pltr.plot_pred_vs_pred(nm_results, tuc_results, request.get_id())

    def report_adjustments(self, tuc_results: [[float], [float]], request: Request, query: Query):
        """
        * This method reports the results from adjustments of Tucuxi.
        * It is for single curves for apriori and aposteriori. It calls the plotter to make graphs and print
          them to images. It exports xUnit compatible xml for Jenkins
            Args:
                tuc_results ([[],[]]): results from tucucli in the format [[times],[values]]
                request (Request): the dataset for which to report
                query (Query): The query containing the request
        """
        if self.generateGraphs:
            pltr = Plotter(get_output_folder_name(), query.get_id() + request.get_id(), plotargs[postcount])

            title = 'Adjustments'
            pltr.plot_adjustments(tuc_results, query.get_id() + request.get_id(), title)

    def report_single(self, tuc_results: [[float], [float]], request: Request, query: Query):
        """
        * This method reports the results from calculations of Tucuxi.
        * Its for single curves for apriori and aposteriori. It calls the plotter to make graphs and print
          them to images. It exports xUnit compatible xml for Jenkins
            Args:
                tuc_results ([[],[]]): results from tucucli in the format [[times],[values]]
                request (Request): the dataset for which to report
                query (Query): The query containing the request
        """

        if self.generateGraphs:
            pltr = Plotter(get_output_folder_name(), query.get_id() + request.get_id(), plotargs[postcount])

            if request.computingTraits.computingTraits.requestType == RequestType.Prediction:

                title = 'Undefined title'
                if request.computingTraits.computingTraits.computingOption.parametersType == \
                        ParametersTypeEnum.population:
                    title = 'Population'
                elif request.computingTraits.computingTraits.computingOption.parametersType == \
                        ParametersTypeEnum.apriori:
                    title = 'A priori'
                elif request.computingTraits.computingTraits.computingOption.parametersType == \
                        ParametersTypeEnum.aposteriori:
                    title = 'A posteriori'

                pltr.plot_single_prediction(tuc_results, query.get_id() + request.get_id(), title)

    @staticmethod
    def compute_error(tuc_results, nm_results, title, abs_threshold=100.0, rel_threshold=0.01):
        sum_squares = 0.0
        relative = 0.0

        alldiff = []

        has_passed = True
        failtext = ''

        if len(tuc_results) != len(nm_results):
            has_passed = False
            failtext += title + \
                ' Test Fail : Tucuxi and Nonmem number of points differ: {t} in Tucuxi vs {n} in Nonmem.' \
                .format(t=len(tuc_results), n=len(nm_results))
            print(Fore.RED + failtext)
            return has_passed, failtext

        for i in range(0, len(tuc_results)):
            diff = tuc_results[i] - nm_results[i]
            alldiff.append(diff)
            mean_result = (tuc_results[i] + nm_results[i]) / 2
            if mean_result != 0:
                relative += abs(diff) / mean_result
            sum_squares += (diff * diff)
        err = sum_squares / len(tuc_results)
        relative_error = relative / len(tuc_results)
        # slope, intercept, r_value, p_value, std_err = stats.linregress(tuc_results, nm_results)

        # if r_value ** 2 < thresh:
        #     has_passed = False
        #     failtext = title + ' Test Fail. Tucuxi vs Nonmem r-squared < {thr} % : {rval} %' \
        #         .format(thr=thresh * 100, rval=(r_value ** 2) * 100)
        #     print(Fore.RED + failtext)
        # else:
        #     print(Fore.GREEN + title + ' r_squared: {rval} %'.format(rval=(r_value ** 2) * 100))

        if abs_threshold > 0:
            if err > abs_threshold:
                if has_passed:
                    failtext += title + ' Test Fail : '
                failtext += 'Tucuxi vs Nonmem mean squared error > {thr} : {err}.' \
                    .format(thr=abs_threshold, err=err)
                print(Fore.RED + failtext)
                has_passed = False
            else:
                print(Fore.GREEN + title + ' mean squared error: {err}'.format(err=err))

        if rel_threshold > 0:
            if relative_error > rel_threshold:
                if has_passed:
                    failtext += title + ' Test Fail : '
                failtext += 'Tucuxi vs Nonmem relative mean error > {thr} : {err}.' \
                    .format(thr=rel_threshold, err=relative_error)
                print(Fore.RED + failtext)
                has_passed = False
            else:
                print(Fore.GREEN + title + ' mean relative error: {err}'.format(err=relative_error))

        return has_passed, failtext

    def report_single_percentiles(self, query, tuc_results, request):

        if self.generateGraphs:
            # Plot only Tucuxi results
            pltr = Plotter(get_output_folder_name(), query.get_id() + request.get_id() + '_tucuxi', plotargs[perccount])
            pltr.new_plot_percentiles_conc_time(tuc_results)

    def report_percentiles(self, query: Query, request: Request, tuc_results: [[float], [float]],
                           nm_results: [[float], [float]]):
        """
        * This method reports the results from percentiles calculations of Tucuxi and NONMEM with the same data.
        * Its for apriori percentiles. It calls the plotter to make graphs and print
          them to images. It exports xUnit compatible xml for Jenkins
            Args:
                tuc_results ([[],[]]): results from Tucuxi percentiles in the format [[times],[values]]
                nm_results ([[],[]]): results from nonmem percenitles in the format [[times],[values]]
                request (Request): the dataset for which to report
                query (Query): The query containing the request
        """

        if self.generateGraphs:
            # Plot only Tucuxi results
            pltr = Plotter(get_output_folder_name(), query.get_id() + '_tucuxi', plotargs[perccount])
            pltr.new_plot_percentiles_conc_time(tuc_results)

            pltr1 = Plotter(get_output_folder_name(), query.get_id() + '_nonmem',
                            plotargs[postcount])
            pltr1.new_plot_percentiles_conc_time(nm_results)

            pltr.plot_all_percentiles_compare(tuc_results, nm_results,
                                              query.get_id() + request.get_id() + ".full_percentiles_compare")

        for i in range(0, len(nm_results[1])):
            has_passed, failtext = self.compute_error(tuc_results[1][i], nm_results[1][i], 'Percentile {p}'
                                                      .format(
                p=request.computingTraits.computingTraits.ranks[i]), 0, 0.05)

            self.jrep.add_testcase(has_passed, request.requestId + '_p{p}'.format(
                p=request.computingTraits.computingTraits.ranks[i]), 'percentiles', failtext)
