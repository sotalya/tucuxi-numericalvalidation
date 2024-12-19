#!/usr/bin/python


foldername = ''


def get_output_folder_name():
    global foldername
    return foldername


def set_output_folder_name(name):
    global foldername
    foldername = name

# I removed this functionality as rm -rf can be a bit dangerous...
# def clean_and_exit(message):
    #     global foldername
    #     Popen('rm -rf ' + foldername, shell=True)
#     exit(message)


# I removed this functionality as rm -rf can be a bit dangerous...
# def clean():
    #     global foldername
    # Popen('rm -rf ' + foldername, shell=True)


