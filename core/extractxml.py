
from scrape.response import *
from datetime import datetime
from core.utils import get_output_folder_name
from sotalya.data.requests import RequestType
import numpy as np


import os
import subprocess


class QueryResponseXmlExtractor:

    def __init__(self, queryresponse):
        self.queryResponse = queryresponse
        self.output_dir = os.path.join(get_output_folder_name(), 'tucuxi')

    def run_extractor(self, query):
        responses = Responses(self.queryResponse.queryId)
        index = 0
        for r in self.queryResponse.responses:

            response = Response(r.requestId, r.requestType)
            response_request_type = response.requestType

            if response_request_type == RequestType.Adjustment.value:
                response.results = []
                adjustment = r.computingTrait.adjustments.adjustment
                score = adjustment.score.string
                times = []
                values = []
                feed_times_values(times, values, adjustment)
                response.results.append([score, times, values])
                response.request = query.requests[index]
                responses.add_response(response)
            elif response_request_type == RequestType.Prediction.value:
                response.results = []
                prediction = r.computingTrait
                times = []
                values = []
                feed_times_values(times, values, prediction)
                response.results.append([times, values])
                response.request = query.requests[index]
                responses.add_response(response)
            elif response_request_type == RequestType.Percentiles.value:

                times = []
                values = []
                temp_list = []
                response.results = [[], []]

                first_cycle_data = r.computingTrait.percentileList[0].cycleDatas[0]
                first_time = first_cycle_data.start

                for percentile in r.computingTrait.percentileList:
                    for cycleData in percentile.cycleDatas:
                        for j in range(len(cycleData.values)):
                            temp_list.append(float(cycleData.values[j]))

                    values.append(temp_list)
                    temp_list = []
                response.results[1] = values

                for cycleData in r.computingTrait.percentileList[0].cycleDatas:
                    for i in range(len(cycleData.times)):
                        per = (((cycleData.start - first_time).total_seconds()) / 3600.0) \
                              + cycleData.times[i]
                              # + first_cycle_data.times[i]
                        times.append(float(per))

                response.results[0] = times

                response.request = query.requests[index]
                responses.add_response(response)

            elif response_request_type == RequestType.PredictionAtTimes.value:
                start_dose = query.drugs[0].dosageHistory.dosageTimeRanges[0].start
                response.results = []
                prediction = r.computingTrait
                times = []
                values = []
                feed_times_values_singlepoints(times, values, prediction, start_dose)
                response.results.append([times, values])
                response.request = query.requests[index]
                responses.add_response(response)

            index = index + 1

        return responses


def feed_times_values_singlepoints(realtime, realvalue, xmlroot, start_dose):
    for points in xmlroot.points:
        date = datetime.strptime(points.time, '%Y-%m-%dT%H:%M:%S')
        absolutetime = date - start_dose
        realtime.append(absolutetime.total_seconds() / 3600.0)
        realvalue.append(float(points.value))

def feed_times_values(realtime, realvalue, xmlroot):
    for cycleData in xmlroot.cycleDatas:
        for time in cycleData.times:
            absolutetime = cycleData.start - xmlroot.cycleDatas[0].start
            realtime.append(absolutetime.total_seconds() / 3600.0 + time)
        for value in cycleData.values:
            realvalue.append(float(value))


def str_to_datetime(string):
    return datetime.strptime(string, '%Y-%m-%dT%H:%M:%S')