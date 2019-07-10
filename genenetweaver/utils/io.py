from time import strftime, localtime


def get_time_string():
    return strftime("%Y-%m-%d_%H:%M:%S", localtime())


def get_time_stamp():
    return strftime("%Y%m%d_%H%M%S", localtime())

