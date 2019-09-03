"""Module in charge of parallelizing the execution of the GSO algorithm in different clusters."""

from multiprocessing import Process, cpu_count
import cProfile
from lightdock.util.logger import LoggingManager


class Tentacle(Process):
    """A Kraken without tentacles would be a sea serpent, right?"""
    def __init__(self, tasks, profiling=False):
        super(Tentacle, self).__init__()
        self.tasks = tasks
        self.profiling = profiling
        self.log = LoggingManager.get_logger('kraken')
        self.log.info("Tentacle ready with %d tasks" % len(self.tasks))
        
    def run(self):
        for task in self.tasks:
            if not self.profiling:
                task.run()
            else:
                cProfile.runctx('task.run()', globals(), locals(), 'process_%s.out' % self.name)
        self.log.info("folding tentacle %s" % self.name)
        

class Kraken(object):
    """Below the thunders of the upper deep;
    Far, far beneath in the abysmal sea, 
    His ancient, dreamless, uninvaded sleep 
    The Kraken sleepeth: faintest sunlights flee
    
    The Kraken 1830, Alfred Tennyson 
    """    
    def __init__(self, tasks, num_cpus=0, profiling=False):
        self.log = LoggingManager.get_logger('kraken')
        try:
            self.num_processes = int(num_cpus)
            if self.num_processes < 1:
                raise ValueError()
            if self.num_processes > cpu_count():
                self.log.warning("Number of cores (%d) larger than available." % self.num_processes)
                raise ValueError()
        except (ValueError, TypeError):
            self.log.warning("Number of cores has not been specified or is incorrect. Using available cores.")
            self.num_processes = cpu_count()
        
        self.log.info("Kraken has %d tentacles (cpu cores)" % self.num_processes)
        
        self.tasks = tasks
        self.num_tasks = len(tasks)
        self.tentacles = []
        tentacle_tasks = [tasks[i::self.num_processes] for i in range(self.num_processes)]
        
        for i in range(self.num_processes):
            tentacle = Tentacle(tentacle_tasks[i], profiling)
            self.tentacles.append(tentacle)
        
        self.log.info("%d ships ready to be smashed" % self.num_tasks)
    
    def release(self):
        """Unleash the wrath of this monster"""
        self.log.info("Release the Kraken!")
        for tentacle in self.tentacles:
            tentacle.start()
            
        for tentacle in self.tentacles:
            tentacle.join()
        
        self.log.info("%d ships destroyed" % self.num_tasks)
        
        reports = [task.gso.report() for task in self.tasks]
        
        return reports

    def sink(self):
        """Sink this monster"""
        for tentacle in self.tentacles:
            tentacle.terminate()
        self.log.warning("Kraken sunk to the bottom of the ocean")
