def run_until_completed(self, args=None, kwargs=None, unparsed_args=None):
    max_ticks = 1e3
    tick_interval = .5
    mc_dao = self._get_mc_dao()
    mc_clients = self._get_mc_clients(mc_dao=mc_dao)
    flow_runner = self._get_flow_runner(mc_clients=mc_clients)
    job_runner = self._get_job_runner(mc_clients=mc_clients)
    tick_counter = 0
    while self._has_unfinished_items(mc_dao=mc_dao):
        tick_counter += 1
        self.logger.info('Tick #%s' % (tick_counter))
        if tick_counter > max_ticks: raise Exception("exceed max ticks")
        flow_runner.tick()
        job_runner.tick()
        time.sleep(tick_interval)
    failed_flows = self._get_failed_flows(mc_dao=mc_dao)
    if failed_flows:
        raise Exception("Failed flows:\n{failed_flows_dump}".format(
            failed_flows_dump=pprint.pformat(failed_flows)))
    self.logger.info('Completed.')
