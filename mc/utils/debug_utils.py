import inspect
import logging

def debug_locals(frames_back=1, logger=logging):
    current_frame = inspect.currentframe()
    frame = current_frame
    try:
        for i in range(frames_back): frame = frame.f_back
        fn = inspect.stack()[1 + frames_back][3]
        msg = "fn: {fn}, locals: {locals_}".format(
            fn=fn,
            locals_=frame.f_back.f_locals
        )
        logger.debug(msg)
    finally: del current_frame
