import logging


def get_default_logger(name=...):
    logger_args = []
    if name is not ...: logger_args.append(name)
    logger = logging.getLogger(*logger_args)
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)
    return logger
