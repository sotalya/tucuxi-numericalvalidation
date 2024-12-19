#!/usr/bin/env python
# coding=utf-8

import xml.etree.cElementTree as elementTree

# from bs4 import BeautifulSoup


class JUnitReport:

    def __init__(self, testsuites_name='testsuites'):
        self.suites = elementTree.Element(testsuites_name)
        self.suite = None
        self.currentSuiteName = ''

    def write_junit_report(self, filename):
        tree = elementTree.ElementTree(self.suites)
        # with open(filename + '.xml', "w", encoding='utf-8') as file:
        #     file.write(str(BeautifulSoup(elementTree.tostring(tree), "xml").prettify()))
        tree.write(filename + '.xml')

    def add_testsuite(self, numtests):
        self.currentSuiteName = '{nt}'.format(nt=numtests)
        self.suite = elementTree.SubElement(self.suites, "testsuite", tests=self.currentSuiteName)

    def add_testcase(self, has_passed, clname, testname, failtext):
        classname = clname
#        classname = '{} : {}'.format(self.currentSuiteName, clname)
        if has_passed:
            elementTree.SubElement(self.suite, "testcase", classname=classname, name=testname)
        else:
            fail = elementTree.SubElement(self.suite, "testcase", classname=classname, name=testname)
            elementTree.SubElement(fail, "failure", type="").text = failtext
