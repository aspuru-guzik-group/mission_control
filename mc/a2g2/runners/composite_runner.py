import time


class CompositeRunner(object):
    def __init__(self, *args, runners=None, tick_interval=120):
        self.runners = runners or {}
        self.tick_counter = 0
        self.tick_interval = tick_interval

    def run(self, ntimes=None, tick_interval=None):
        if ntimes:
            for i in range(ntimes):
                self._tick_and_sleep(tick_interval=tick_interval)
        else:
            while self._ticking:
                self._tick_and_sleep(tick_interval=tick_interval)

    def _tick_and_sleep(self, tick_interval=None):
        if tick_interval is None: tick_interval = self.tick_interval
        self.tick()
        time.sleep(tick_interval)

    def tick(self):
        self.tick_counter += 1
        for runner in self.runners.values(): runner.tick()
