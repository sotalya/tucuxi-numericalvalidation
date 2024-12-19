#!/usr/bin/python


class Responses:
    def __init__(self, query_id):
        self.responses = []
        self.queryId = query_id

    def add_response(self, response):
        self.responses.append(response)


class Response:
    def __init__(self, request_id, request_type):
        self.requestId = request_id
        self.requestType = request_type
